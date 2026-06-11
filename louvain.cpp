#include <queue>
#include <algorithm>
#include <chrono>
#include <numeric>
#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/Coloring.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"

extern "C" {
    #include <igraph/igraph.h>
}

using Clock = std::chrono::steady_clock;

static double elapsed_seconds(Clock::time_point begin, Clock::time_point end)
{
    return std::chrono::duration<double>(end - begin).count();
}

// Louvain法によるクラスタリング
std::vector<int> Louvain_from_Boost(const Graph& G) {
    igraph_rng_seed(igraph_rng_default(), 42u);

    igraph_t ig;
    if (igraph_empty(&ig, static_cast<igraph_integer_t>(num_vertices(G)), IGRAPH_UNDIRECTED)) {
        throw std::runtime_error("igraph_empty failed");
    }

    std::vector<igraph_integer_t> edge_list;
    edge_list.reserve(2 * num_edges(G));

    for (auto e : boost::make_iterator_range(edges(G))) {
        edge_list.push_back(static_cast<igraph_integer_t>(source(e, G)));
        edge_list.push_back(static_cast<igraph_integer_t>(target(e, G)));
    }

    int rc = IGRAPH_SUCCESS;

    igraph_vector_int_t edge_vec =
        igraph_vector_int_view(
            edge_list.data(),
            static_cast<igraph_integer_t>(edge_list.size())
        );
    rc = igraph_add_edges(&ig, &edge_vec, /*attr=*/nullptr);

    if (rc != IGRAPH_SUCCESS) {
            igraph_destroy(&ig);
            throw std::runtime_error("igraph_add_edges failed: " + std::to_string(rc));
    }

    std::vector<int> result(num_vertices(G), 0);

    igraph_vector_int_t membership;
    igraph_matrix_int_t memberships;
    igraph_vector_t modularity;

    igraph_vector_int_init(&membership, 0);
    igraph_matrix_int_init(&memberships, 0, 0);
    igraph_vector_init(&modularity, 0);

    rc = igraph_community_multilevel(
        &ig, /*weights*/ nullptr, /*resolution*/ 1.0,
        &membership, &memberships, &modularity);

    if (rc != IGRAPH_SUCCESS) { throw std::runtime_error("community_multilevel failed"); }

    const igraph_integer_t levels = igraph_vector_size(&modularity);
    const igraph_integer_t comm_count = igraph_vector_int_max(&membership) + 1;

    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&membership); ++i)
        result[static_cast<size_t>(i)] = static_cast<int>(VECTOR(membership)[i]);

    igraph_vector_int_destroy(&membership);
    igraph_matrix_int_destroy(&memberships);
    igraph_vector_destroy(&modularity);
    igraph_destroy(&ig);

    return result;
}

// 返り値: 各ノードの新コミュニティID（0..C-1、サイズの大きい順）
// オプション: sizes_out に新ID順のサイズ（長さC）を返す
std::vector<int> relabel_dense(
    const std::vector<int>& comm,
    std::vector<int>* sizes_out = nullptr)
{
    const size_t N = comm.size();

    // 旧ラベルごとのサイズを集計
    struct Stat { int label; int count; };
    std::unordered_map<int, int> cnt;
    cnt.reserve(N * 2);
    for (int c : comm) {
        if (c >= 0) ++cnt[c];
    }
    if (cnt.empty()) { // すべて未所属など
        if (sizes_out) sizes_out->clear();
        return std::vector<int>(N, -1);
    }

    // 並べ替え: 大きい順、同数は旧ラベルの昇順
    std::vector<Stat> order;
    order.reserve(cnt.size());
    for (auto& kv : cnt) order.push_back({kv.first, kv.second});
    std::sort(order.begin(), order.end(), [](const Stat& a, const Stat& b){
        if (a.count != b.count) return a.count > b.count;
        return a.label < b.label;
    });

    // 旧→新ラベル写像を作成（大きい順に 0,1,2,... を付与）
    std::unordered_map<int,int> old2new;
    old2new.reserve(order.size());
    if (sizes_out) sizes_out->assign(order.size(), 0);
    for (int new_id = 0; new_id < static_cast<int>(order.size()); ++new_id) {
        old2new[order[new_id].label] = new_id;
        if (sizes_out) (*sizes_out)[new_id] = order[new_id].count;
    }

    // 各ノードのラベルを詰め直し
    std::vector<int> dense(N, -1);
    for (size_t i = 0; i < N; ++i) {
        int c = comm[i];
        if (c >= 0) dense[i] = old2new[c];
    }
    return dense;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <matrix.mtx>\n", argv[0]);
        return 1;
    }

    Graph G = Read_MM_UD(argv[1]);    // 疎行列の隣接グラフ（無向グラフ）

    auto after_read_to_write_begin = Clock::now();

    std::vector<int> block_of = Louvain_from_Boost(G);
    std::vector<int> sizes;
    auto dense = relabel_dense(block_of, &sizes);
    int nb = sizes.size();   // コミュニティの数（= ブロック数）
    
    //////////////////////////////////////////////
    // ブロックグラフの作成
    Graph T = BuildBlockGraph(G, dense, BlockEdgeWeight::Binary);

    //////////////////////////////////////////////
    // ブロック内結合度の評価
    auto internal = CountInternalEdges(G, dense, nb);
    double total_avg = 0.0;

    for (int b = 0; b < nb; ++b) {
        double avg_deg = (sizes[b] > 0) ? 2.0 * internal[b] / sizes[b] : 0.0;
        total_avg += avg_deg;
    }

    //////////////////////////////////////////////
    // ブロックグラフの彩色
    std::vector<int> block_color;
    int nc = Greedy_Coloring(T, block_color);

    // 色ラベルを頻度順に付け替え
    RelabelColorsByClassSize(block_color);

    auto after_read_to_write_end = Clock::now();

    // 出力ファイル名は <入力行列のstem>.blk, <stem>.bcol
    std::string stem = file_stem(argv[1]);
    stem += "_louvain";
    std::string blk_path  = stem + ".blk";
    std::string bcol_path = stem + ".bcol";

    // ブロック情報データの出力
    WriteBlockInfo_1Based(dense, blk_path);

    // ブロック色情報データの出力
    WriteBlockColor_1Based(block_color, nc, bcol_path);

    // --- モジュラリティ（未加重）
    double Q = Modularity_Unweighted(G, dense);

    std::cout << "time = " << elapsed_seconds(after_read_to_write_begin, after_read_to_write_end)
              << " NB = " << nb
              << " NC = " << nc
              << " modularity = " << Q << "\n";
    return 0;
}
