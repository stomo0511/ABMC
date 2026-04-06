#include <queue>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cmath>
#include <type_traits>
#include <vector>
#include <string>

#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/Coloring.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"
#include "abmc.hpp"

// =============================================================================
// COMPATIBILITY SHIM FOR RABBIT ORDER
// =============================================================================
// 1. 名前空間と構造体の前方宣言
namespace rabbit_order { namespace aux { struct atom; } }

// 2. 標準ライブラリのトレイトを std 名前空間内で特殊化
namespace std {
    template<> 
    struct is_trivially_copyable<rabbit_order::aux::atom> : std::true_type {};
}

// 3. Boost 独自のトレイトを boost 名前空間内で特殊化
#include <boost/type_traits/is_trivially_copyable.hpp>
namespace boost {
    template<> 
    struct is_trivially_copyable<rabbit_order::aux::atom> : true_type {};
}

// 4. オリジナルヘッダの読み込み
#include "rabbit_order.hpp"
// =============================================================================

using namespace rabbit_order::aux;

/**
 * Rabbit Order と ABMC の統合処理を管理するクラス
 */
class RabbitABMC {
public:
    // コンストラクタ
    RabbitABMC(const Graph& G, int nb) 
        : G_(G), N_(static_cast<int>(boost::num_vertices(G))) {
        block_size_ = (N_ + nb - 1) / nb; // ターゲットとするブロックサイズ s
        block_of_.assign(N_, -1);
    }

    /**
     * ハイブリッド順序付けのメイン実行フロー
     */
    BlockPartition run() {
        // Step 1: Rabbit 用の隣接リスト作成 (重みの絶対値化により Assertion 回避)
        auto adj = convert_to_rabbit_adj();

        // Step 2: Rabbit Order の実行
        rg_ = std::make_unique<graph>(aggregate(std::move(adj)));

        // Step 3: 事前情報の計算 (絶縁判定・サイズ・階層情報)
        is_insular_ = identify_insular_nodes();
        subtree_sizes_.assign(N_, 0);
        is_pure_insular_subtree_.assign(N_, false);

        if (rg_->tops) {
            for (vint top : *rg_->tops) {
                preprocess_subtree_info(top);
            }
        }

        // Step 4: パティショニング（不均一性解消ロジック付き）
        int current_block_id = 0;
        residual_nodes_.clear();

        if (rg_->tops) {
            for (vint top : *rg_->tops) {
                partition_with_packing(top, current_block_id);
            }
        }

        // 最後にバッファに残ったノードを最終ブロックとして放出
        flush_residual_nodes(current_block_id, true);

        // 結果のパッケージング
        BlockPartition part;
        part.n = N_;
        part.s = block_size_;
        part.nb = current_block_id;
        part.block_of = block_of_;
        part.blocks.resize(part.nb);
        for (int i = 0; i < N_; ++i) {
            if (part.block_of[i] != -1) {
                part.blocks[part.block_of[i]].push_back(i);
            }
        }
        return part;
    }

private:
    const Graph& G_;
    int N_;
    int block_size_;
    std::unique_ptr<graph> rg_;
    std::vector<int> block_of_;
    std::vector<bool> is_insular_;
    std::vector<int> subtree_sizes_;
    std::vector<bool> is_pure_insular_subtree_;
    std::vector<int> residual_nodes_; // サイズ s に満たないノードの集約バッファ

    // --- 内部ロジック関数 ---

    // 重みを絶対値にして Rabbit 用の隣接リストを生成
    std::vector<std::vector<edge>> convert_to_rabbit_adj() {
        auto weightmap = boost::get(boost::edge_weight, G_);
        std::vector<std::vector<edge>> adj(N_);
        for (int u = 0; u < N_; ++u) {
            boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
            for (boost::tie(ei, ei_end) = boost::out_edges(u, G_); ei != ei_end; ++ei) {
                int v = (int)boost::target(*ei, G_);
                float weight = std::abs(static_cast<float>(weightmap[*ei]));
                adj[u].push_back({static_cast<vint>(v), weight});
            }
        }
        return adj;
    }

    // 絶縁ノード (自分のコミュニティ内にエッジが閉じているノード) の判定
    std::vector<bool> identify_insular_nodes() {
        std::vector<bool> res(N_, true);
        #pragma omp parallel for
        for (int u = 0; u < N_; ++u) {
            vint u_com = rg_->coms[u];
            for (const auto& e : rg_->es[u]) {
                if (rg_->coms[e.first] != u_com) { res[u] = false; break; }
            }
        }
        return res;
    }

    // サブツリーのサイズと絶縁性をボトムアップで計算
    void preprocess_subtree_info(int v) {
        if (v == vmax) return;
        int size = 1;
        bool all_insular = is_insular_[v];

        vint child = rg_->vs[v].a->child;
        while (child != vmax) {
            preprocess_subtree_info(child);
            size += subtree_sizes_[child];
            if (!is_pure_insular_subtree_[child]) all_insular = false;
            child = rg_->vs[child].sibling;
        }
        subtree_sizes_[v] = size;
        is_pure_insular_subtree_[v] = all_insular;
    }

    // 均一化（パッキング）と絶縁優先を考慮した分割ロジック
    void partition_with_packing(int v, int& current_block_id) {
        if (v == vmax || block_of_[v] != -1) return;

        // 絶縁ノード優先: 十分なサイズ (0.8s以上) かつ純粋絶縁なら独立ブロック化
        if (subtree_sizes_[v] <= block_size_ && subtree_sizes_[v] >= block_size_ * 0.8 && is_pure_insular_subtree_[v]) {
            assign_subtree_to_block(v, current_block_id++);
            return;
        }

        if (subtree_sizes_[v] > block_size_) {
            vint child = rg_->vs[v].a->child;
            while (child != vmax) {
                partition_with_packing(child, current_block_id);
                child = rg_->vs[child].sibling;
            }
            // 親自身が未割当ならバッファへ
            if (block_of_[v] == -1) {
                residual_nodes_.push_back(v);
                flush_residual_nodes(current_block_id, false);
            }
        } else {
            // サイズ s 未満のサブツリーは集約バッファに入れて結合を図る
            collect_subtree_nodes(v);
            flush_residual_nodes(current_block_id, false);
        }
    }

    // サブツリー内の未割当ノードをバッファに収集
    void collect_subtree_nodes(int v) {
        if (v == vmax || block_of_[v] != -1) return;
        residual_nodes_.push_back(v);
        vint child = rg_->vs[v].a->child;
        while (child != vmax) {
            collect_subtree_nodes(child);
            child = rg_->vs[child].sibling;
        }
    }

    // バッファに溜まったノードをブロックサイズ s 単位でブロック化
    void flush_residual_nodes(int& current_block_id, bool force) {
        while (residual_nodes_.size() >= (size_t)block_size_ || (force && !residual_nodes_.empty())) {
            int count = 0;
            while (!residual_nodes_.empty() && count < block_size_) {
                int node = residual_nodes_.back();
                residual_nodes_.pop_back();
                if (block_of_[node] == -1) {
                    block_of_[node] = current_block_id;
                    count++;
                }
            }
            if (count > 0) current_block_id++;
            if (force && residual_nodes_.empty()) break;
            if (!force && residual_nodes_.size() < (size_t)block_size_) break;
        }
    }

    // 指定したサブツリー以下を全て現在のブロックに割り当て
    void assign_subtree_to_block(int v, int block_id) {
        if (v == vmax || block_of_[v] != -1) return;
        block_of_[v] = block_id;
        vint child = rg_->vs[v].a->child;
        while (child != vmax) {
            assign_subtree_to_block(child, block_id);
            child = rg_->vs[child].sibling;
        }
    }
};

// =============================================================================
// MAIN FUNCTION
// =============================================================================

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <matrix.mtx> <# blocks>\n", argv[0]);
        return 1;
    }

    // 1. データの読み込み
    Graph G = Read_MM_UD(argv[1]);
    int nb = std::atoi(argv[2]);

    // 2. Rabbit-ABMC ハイブリッド順序付けの実行
    RabbitABMC solver(G, nb);
    BlockPartition part = solver.run();

    // 3. ブロックグラフの構築と彩色
    Graph T = BuildBlockGraph(G, part.block_of, BlockEdgeWeight::Binary);
    std::vector<int> block_color;
    int nc = Greedy_Coloring(T, block_color);
    RelabelColorsByClassSize(block_color);

    // 4. 結果の出力 (ブロック数をファイル名に付与)
    std::string stem = file_stem(argv[1]);
    std::string blk_path = stem + "_rabbit_B" + std::to_string(nb) + ".blk";
    std::string bcol_path = stem + "_rabbit_B" + std::to_string(nb) + ".bcol";

    WriteBlockInfo_1Based(part.block_of, blk_path);
    WriteBlockColor_1Based(block_color, nc, bcol_path);

    // 5. 評価値の表示
    double Q = Modularity_Unweighted(G, part.block_of);
    std::printf("Rabbit-ABMC Completed.\n#Blocks: %d, #Colors: %d, Modularity: %.6f\n", part.nb, nc, Q);

    return 0;
}