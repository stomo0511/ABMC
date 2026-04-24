#include <queue>
#include <algorithm>
#include <numeric>
#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/Coloring.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"

extern "C" {
    #include <igraph/igraph.h>
}

// 0.10判定マクロ（0.9では未定義のことがあるので保険）
#ifndef IGRAPH_VERSION_MAJOR
#  define IGRAPH_VERSION_MAJOR 0
#  define IGRAPH_VERSION_MINOR 9
#endif
#if (IGRAPH_VERSION_MAJOR > 0) || (IGRAPH_VERSION_MAJOR==0 && IGRAPH_VERSION_MINOR>=10)
#  define IGRAPH_AT_LEAST_010 1
#else
#  define IGRAPH_AT_LEAST_010 0
#endif

// Louvain法によるクラスタリング
std::vector<int> Louvain_from_Boost(const Graph& G) {
    igraph_rng_seed(igraph_rng_default(), 42u);  // 任意の固定シードで実行時の再現性を確保

    // 1) igraph を頂点数指定で初期化（無向）
    igraph_t ig;
    if (igraph_empty(&ig, static_cast<igraph_integer_t>(num_vertices(G)), IGRAPH_UNDIRECTED)) {
        throw std::runtime_error("igraph_empty failed");
    }

    // 2) Boost.Graph の辺を igraph の edges ベクタ（u0, v0, u1, v1, ...）へ
    std::vector<igraph_integer_t> edge_list;
    edge_list.reserve(2 * num_edges(G));
    for (auto e : boost::make_iterator_range(edges(G))) {
        edge_list.push_back(static_cast<igraph_integer_t>(source(e, G)));
        edge_list.push_back(static_cast<igraph_integer_t>(target(e, G)));
    }

    // 3) view（所有権なし）を作って add_edges
    //    ※ 0.10系は「戻り値で view を返す」2引数API
    std::vector<int> result(num_vertices(G), 0);

    int rc = IGRAPH_SUCCESS;

#if IGRAPH_AT_LEAST_010
    igraph_vector_int_t edge_vec;
    igraph_vector_int_view(&edge_vec, edge_list.data(),
                            static_cast<igraph_integer_t>(edge_list.size()));
    rc = igraph_add_edges(&ig, &edge_vec, /*attr=*/nullptr);

#else
    igraph_vector_t evec;
    igraph_vector_init(&evec, static_cast<long>(edge_list.size()));
    for (long i = 0; i < static_cast<long>(edge_list.size()); ++i) {
        VECTOR(evec)[i] = static_cast<igraph_real_t>(edge_list[i]);
    }
    rc = igraph_add_edges(&ig, &evec, /*attr=*/nullptr);
    igraph_vector_destroy(&evec);
#endif

    if (rc != IGRAPH_SUCCESS) {
            igraph_destroy(&ig);
            throw std::runtime_error("igraph_add_edges failed: " + std::to_string(rc));
        }

    // 4) Louvain (multilevel) 実行
#if IGRAPH_AT_LEAST_010
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

    // std::cout << "vertices: " << igraph_vcount(&ig) << "\n";
    // std::cout << "communities: " << comm_count << "\n";
    // std::cout << "final modularity: " << VECTOR(modularity)[levels - 1] << "\n";

    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&membership); ++i)
        result[static_cast<size_t>(i)] = static_cast<int>(VECTOR(membership)[i]);

    igraph_vector_int_destroy(&membership);
    igraph_matrix_int_destroy(&memberships);
    igraph_vector_destroy(&modularity);

#else
    igraph_vector_t membership;   // 実数ベクタ（中身は整数ID）
    igraph_matrix_t memberships;  // 実数行列（階層ごとの membership）
    igraph_vector_t modularity;   // 実数ベクタ（階層ごとの Q）

    igraph_vector_init(&membership, 0);
    igraph_matrix_init(&memberships, 0, 0);
    igraph_vector_init(&modularity, 0);

    rc = igraph_community_multilevel(
        &ig, /*weights*/ nullptr, /*resolution*/ 1.0,
        &membership, &memberships, &modularity);
    if (rc != IGRAPH_SUCCESS) { throw std::runtime_error("community_multilevel failed"); }

    const long levels = igraph_vector_size(&modularity);
    const long comm_count = static_cast<long>(igraph_vector_max(&membership)) + 1;

    // std::cout << "vertices: " << igraph_vcount(&ig) << "\n";
    // std::cout << "communities: " << comm_count << "\n";
    std::cout << "NB = " << comm_count;
    // std::cout << "final modularity: " << VECTOR(modularity)[levels - 1] << "\n";

    for (long i = 0; i < igraph_vector_size(&membership); ++i)
        result[static_cast<size_t>(i)] = static_cast<int>(VECTOR(membership)[i]);

    igraph_vector_destroy(&membership);
    igraph_matrix_destroy(&memberships);
    igraph_vector_destroy(&modularity);
#endif

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

    std::vector<int> block_of = Louvain_from_Boost(G);

    std::vector<int> sizes;
    auto dense = relabel_dense(block_of, &sizes);

    int nb = sizes.size();   // コミュニティの数（= ブロック数）
    
    //////////////////////////////////////////////
    // ブロックグラフの作成
    // Graph T = BuildBlockGraph(G, block_of, BlockEdgeWeight::Binary);
    Graph T = BuildBlockGraph(G, dense, BlockEdgeWeight::Binary);

    //////////////////////////////////////////////
    // ブロック内結合度の評価
    // auto internal = CountInternalEdges(G, block_of, nb);
    auto internal = CountInternalEdges(G, dense, nb);
    double total_avg = 0.0;
    for (int b = 0; b < nb; ++b) {
        double avg_deg = (sizes[b] > 0) ? 2.0 * internal[b] / sizes[b] : 0.0;
        total_avg += avg_deg;
        // std::cout << "Block " << b << ": nodes=" << sizes[b]
        //           << ", internal_edges=" << internal[b]
        //           << ", avg_deg=" << avg_deg << "\n";
    }
    // std::cout << "Total average degree: " << (nb > 0 ? total_avg / nb : 0.0) << "\n";

    //////////////////////////////////////////////
    // ブロック間結合度の評価
    // ブロック間の複数エッジのエッジをカウントする場合
    // Graph T_count = BuildBlockGraph(G, part.block_of, BlockEdgeWeight::Count);
    // auto Tw = get(boost::edge_weight, T_count);
    // for (auto eIt = edges(T_count); eIt.first != eIt.second; ++eIt.first) {
    //     auto e = *eIt.first;
    //     int bu = (int)source(e, T_count);
    //     int bv = (int)target(e, T_count);
    //     double w = Tw[e];
    //     std::cout << "Block " << bu << " - Block " << bv
    //               << ": inter_edges=" << w << "\n";
    // }
    // ブロック間の複数エッジを1本にカウントする場合
    // Graph T_bin = BuildBlockGraph(G, part.block_of, BlockEdgeWeight::Binary);
    // 各ブロックの次数を計算
    total_avg = 0.0;
    for (int b = 0; b < nb; ++b) {
        int deg = boost::degree(b, T);
        // std::cout << "Block " << b << ": degree=" << deg << "\n";
        total_avg += deg;
    }
    // std::cout << "Block graph average degree: "
    //           << (nb > 0 ? total_avg / nb : 0.0) << "\n";

    //////////////////////////////////////////////
    // ブロックグラフの彩色
    std::vector<int> block_color;
    int nc = Greedy_Coloring(T, block_color);

    // 色ラベルを頻度順に付け替え
    RelabelColorsByClassSize(block_color);

    // 出力ファイル名は <入力行列のstem>.blk, <stem>.bcol
    std::string stem = file_stem(argv[1]);
    stem += "_louvain";
    std::string blk_path  = stem + ".blk";
    std::string bcol_path = stem + ".bcol";

    // ブロック情報データの出力
    // WriteBlockInfo_1Based(block_of, blk_path);
    WriteBlockInfo_1Based(dense, blk_path);

    // ブロック色情報データの出力
    WriteBlockColor_1Based(block_color, nc, bcol_path);
    std::cout << " NC = " << nc << "\n";

    // --- モジュラリティ（未加重）
    // double Q = Modularity_Unweighted(G, block_of);
    double Q = Modularity_Unweighted(G, dense);
    // std::printf("Modularity (unweighted)   = %.6f\n", Q);
    return 0;
}
