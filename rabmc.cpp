#include <queue>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cmath>
#include <type_traits>

#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/Coloring.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"

// =============================================================================
// COMPATIBILITY SHIM FOR RABBIT ORDER
// =============================================================================
// オリジナルの rabbit_order.hpp を修正せずに現代の Boost でビルドするための処置。
// boost::atomic<atom> が要求する TriviallyCopyable 特性を明示的に付与します。

// 1. 名前空間と構造体の前方宣言
namespace rabbit_order { namespace aux { struct atom; } }

// 2. 標準ライブラリのトレイトを std 名前空間内で特殊化
namespace std {
    template<> 
    struct is_trivially_copyable<rabbit_order::aux::atom> : std::true_type {};
}

// 3. Boost 独自のトレイトを boost 名前空間内で特殊化 (Boost の内部チェック用)
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
 * Rabbit Order と ABMC の統合処理をカプセル化したクラス
 */
class RabbitABMC {
public:
    RabbitABMC(const Graph& G, int nb) 
        : G_(G), N_(static_cast<int>(boost::num_vertices(G))) {
        block_size_ = (N_ + nb - 1) / nb;
        block_of_.assign(N_, -1);
    }

    /**
     * 再順序付けアルゴリズムのメイン実行フロー
     */
    BlockPartition run() {
        // Step 1: Rabbit Order 用の隣接リスト作成 (重みの絶対値化)
        auto adj = convert_to_rabbit_adj();

        // Step 2: Rabbit Order (Parallel Incremental Aggregation) の実行
        rg_ = std::make_unique<graph>(aggregate(std::move(adj)));

        // Step 3: 絶縁ノードの同定とサブツリーサイズの計算
        is_insular_ = identify_insular_nodes();
        subtree_sizes_.assign(N_, 0);
        if (rg_->tops) {
            for (vint top : *rg_->tops) {
                calculate_subtree_sizes(top);
            }
        }

        // Step 4: デンドログラムの再帰的切り出し
        int current_block_id = 0;
        if (rg_->tops) {
            for (vint top : *rg_->tops) {
                partition_dendrogram(top, current_block_id);
            }
        }

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

    // Rabbit Order 入力形式への変換 (符号反転による Assertion 回避)
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

    // 各ノードの絶縁性 (所属コミュニティ内に全エッジが閉じているか) の判定
    std::vector<bool> identify_insular_nodes() {
        std::vector<bool> res(N_, true);
        #pragma omp parallel for
        for (int u = 0; u < N_; ++u) {
            int u_com = rg_->coms[u];
            for (const auto& e : rg_->es[u]) {
                if (rg_->coms[e.first] != u_com) { res[u] = false; break; }
            }
        }
        return res;
    }

    // 各サブツリーに含まれる原子ノード数の合計を計算
    int calculate_subtree_sizes(int v) {
        if (v == vmax) return 0;
        int size = 1;
        vint child = rg_->vs[v].a->child;
        while (child != vmax) {
            size += calculate_subtree_sizes(child);
            child = rg_->vs[child].sibling;
        }
        return subtree_sizes_[v] = size;
    }

    // デンドログラムを走査し、サイズ s 以下のサブツリーをブロックとして抽出
    void partition_dendrogram(int v, int& current_block_id) {
        if (v == vmax) return;

        if (subtree_sizes_[v] <= block_size_) {
            assign_subtree_to_block(v, current_block_id++);
            return;
        }

        vint child = rg_->vs[v].a->child;
        if (child == vmax) {
            if (block_of_[v] == -1) block_of_[v] = current_block_id++;
        } else {
            while (child != vmax) {
                partition_dendrogram(child, current_block_id);
                child = rg_->vs[child].sibling;
            }
            if (block_of_[v] == -1) block_of_[v] = current_block_id++;
        }
    }

    // 特定のサブツリー以下の全ノードを指定ブロックIDに割り当てる
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

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <matrix.mtx> <# blocks>\n", argv[0]);
        return 1;
    }

    // 1. Matrix Market データの読み込み
    Graph G = Read_MM_UD(argv[1]);
    int nb = std::atoi(argv[2]);

    // 2. ハイブリッド順序付けの実行
    RabbitABMC solver(G, nb);
    BlockPartition part = solver.run();

    // 3. ブロックグラフの構築と彩色
    Graph T = BuildBlockGraph(G, part.block_of, BlockEdgeWeight::Binary);
    std::vector<int> block_color;
    int nc = Greedy_Coloring(T, block_color);
    RelabelColorsByClassSize(block_color);

    // 4. 結果の出力
    std::string stem = file_stem(argv[1]);
    std::string blk_path = stem + "_rabbit_B" + std::to_string(nb) + ".blk";
    std::string bcol_path = stem + "_rabbit_B" + std::to_string(nb) + ".bcol";

    WriteBlockInfo_1Based(part.block_of, blk_path);
    WriteBlockColor_1Based(block_color, nc, bcol_path);

    // 5. 評価値の表示
    double Q = Modularity_Unweighted(G, part.block_of);
    std::printf("Rabbit-ABMC Execution Completed.\n#Blocks: %d, #Colors: %d, Modularity: %.6f\n", part.nb, nc, Q);

    return 0;
}