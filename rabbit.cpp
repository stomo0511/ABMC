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
class RabbitReorder {
public:
    explicit RabbitReorder(const Graph& G)
        : G_(G), N_(static_cast<int>(boost::num_vertices(G))) {}

    std::vector<int> run() {
        // Step 1: Rabbit Order 用の隣接リスト作成 (重みの絶対値化)
        auto adj = convert_to_rabbit_adj();

        // Step 2: Rabbit Order (Parallel Incremental Aggregation) の実行
        rg_ = std::make_unique<graph>(aggregate(std::move(adj)));

        return build_permutation_by_dfs();
    }

private:
    const Graph& G_;
    int N_;
    std::unique_ptr<graph> rg_;

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

    // top-level vertices を集める
    std::vector<int> collect_toplevel_vertices() const {
        std::vector<int> tops;
        tops.reserve(N_);

        for (int v = 0; v < N_; ++v) {
            // dest[v] == v なら top-level とみなす
            // rabbit_order.hpp の graph 構造が論文どおりならこれでよい
            if (rg_->coms[v] == v) {
                tops.push_back(v);
            }
        }

        // 上の判定が実装依存で外れる場合に備えたフォールバック:
        // coms の代表元を集める
        if (tops.empty()) {
            std::unordered_set<int> reps;
            for (int v = 0; v < N_; ++v) reps.insert(rg_->coms[v]);
            tops.assign(reps.begin(), reps.end());
            std::sort(tops.begin(), tops.end());
        }

        return tops;
    }

    void dfs_dendrogram(int v, std::vector<int>& order) const {
        if (v == static_cast<int>(vmax)) return;

        vint child = rg_->vs[v].a->child;
        if (child == vmax) {
            order.push_back(v);
            return;
        }

        vint cur = child;
        while (cur != vmax) {
            dfs_dendrogram(static_cast<int>(cur), order);
            cur = rg_->vs[cur].sibling;
        }

        order.push_back(v);
    }

    std::vector<int> build_permutation_by_dfs() const {
        std::vector<int> order;
        order.reserve(N_);

        auto tops = collect_toplevel_vertices();
        for (int top : tops) {
            dfs_dendrogram(top, order);
        }

        std::vector<char> used(N_, 0);
        std::vector<int> unique_order;
        unique_order.reserve(N_);

        for (int v : order) {
            if (0 <= v && v < N_ && !used[v]) {
                used[v] = 1;
                unique_order.push_back(v);
            }
        }
        for (int v = 0; v < N_; ++v) {
            if (!used[v]) unique_order.push_back(v);
        }

        std::vector<int> perm(N_, -1);
        for (int new_id = 0; new_id < N_; ++new_id) {
            perm[unique_order[new_id]] = new_id;
        }
        return perm;
    }
};

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <matrix.mtx> [output.mtx]\n", argv[0]);
        return 1;
    }

    const std::string in_path  = argv[1];
    const std::string stem     = file_stem(in_path);
    const std::string out_mtx  = (argc >= 3) ? argv[2] : (stem + "_rabbit.mtx");
    const std::string out_perm = stem + "_rabbit.perm";

    // Matrix Market の無向重み付き読み込み
    Graph G = Read_MM_UD(in_path.c_str());

    // 2. ハイブリッド順序付けの実行（現状はRabbit ordering）
    RabbitReorder reorder(G);
    std::vector<int> perm = reorder.run();

    // 3. 並べ替えたグラフを Matrix Market 形式で出力
    WriteReorderedMatrixMarket_UD(G, perm, out_mtx);
    // WritePermutation_1Based(perm, out_perm);

    std::cout << "Wrote reordered matrix : " << out_mtx  << "\n";
    // std::cout << "Wrote permutation      : " << out_perm << "\n";

    return 0;
}