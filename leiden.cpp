#include <queue>
#include <algorithm>
#include <cstdlib>
#include <cerrno>
#include <numeric>
#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/Coloring.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"
#include "common/Timer.hpp"

extern "C" {
    #include <igraph/igraph.h>
}

// Leiden法によるクラスタリング
std::vector<int> Leiden_from_Graph(const Graph& G, double resolution, unsigned long seed, double& quality_out)
{
    igraph_rng_seed(igraph_rng_default(), static_cast<igraph_uint_t>(seed));

    igraph_t ig;
    if (igraph_empty(&ig,
                     static_cast<igraph_integer_t>(num_vertices(G)),
                     IGRAPH_UNDIRECTED)) {
        throw std::runtime_error("igraph_empty failed");
    }

    std::vector<igraph_integer_t> edge_list;
    edge_list.reserve(2 * num_edges(G));

    for_each_undirected_edge(G, [&](int u, int v, double) {
        edge_list.push_back(static_cast<igraph_integer_t>(u));
        edge_list.push_back(static_cast<igraph_integer_t>(v));
    });

    igraph_vector_int_t edge_vec = igraph_vector_int_view(
        edge_list.data(),
        static_cast<igraph_integer_t>(edge_list.size())
    );
    int rc = igraph_add_edges(&ig, &edge_vec, nullptr);

    if (rc != IGRAPH_SUCCESS) {
        igraph_destroy(&ig);
        throw std::runtime_error("igraph_add_edges failed: " + std::to_string(rc));
    }

    igraph_vector_int_t membership;
    igraph_vector_int_init(&membership, 0);

    igraph_int_t nb_clusters = 0;
    igraph_real_t quality = 0.0;

#ifdef CPM
    ////////////////////////////////////////////////
    // CPM
    rc = igraph_community_leiden_simple(
        &ig,
        /* weights */ nullptr,
        /* objective */ IGRAPH_LEIDEN_OBJECTIVE_CPM,
        /* resolution */ resolution,
        /* beta */ 0.01,
        /* start */ false,
        /* n_iterations */ 2,
        &membership,
        &nb_clusters,
        &quality
    );
#else
    ////////////////////////////////////////////////
    // modularity
    rc = igraph_community_leiden_simple(
        &ig,
        /* weights */ nullptr,
        /* objective */ IGRAPH_LEIDEN_OBJECTIVE_MODULARITY,
        /* resolution */ resolution,
        /* beta */ 0.01,
        /* start */ false,
        /* n_iterations */ 2,
        &membership,
        &nb_clusters,
        &quality
    );
#endif

    if (rc != IGRAPH_SUCCESS) {
        igraph_vector_int_destroy(&membership);
        igraph_destroy(&ig);
        throw std::runtime_error("community_leiden_simple failed: " + std::to_string(rc));
    }

    std::vector<int> result(num_vertices(G), 0);
    for (igraph_integer_t i = 0; i < igraph_vector_int_size(&membership); ++i) {
        result[static_cast<size_t>(i)] =
            static_cast<int>(VECTOR(membership)[i]);
    }

    quality_out = static_cast<double>(quality);

    igraph_vector_int_destroy(&membership);
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

std::string sanitize_resolution_for_filename(const std::string& s)
{
    std::string out;
    out.reserve(s.size());

    for (char ch : s) {
        if (ch == '.') {
            out += 'p';
        } else if (ch == '+') {
            // 通常は出ないが、+1.0 のような入力への対応
            continue;
        } else if (ch == '-') {
            // 念のため。解像度は0以上なので通常は不要
            out += 'm';
        } else {
            out += ch;
        }
    }

    return out;
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 5) {
        std::fprintf(stderr,
                    "usage: %s <matrix.mtx> <resolution> [coloring(1|2|3)] [seed]\n",
                    argv[0]);
        return 1;
    }

    char* end = nullptr;
    errno = 0;
    double resolution = std::strtod(argv[2], &end);
    if (errno != 0 || end == argv[1] || *end != '\0') {
        std::fprintf(stderr, "invalid resolution: %s\n", argv[2]);
        return 1;
    }

    int c_in = (argc >= 4) ? std::atoi(argv[3]) : 1;  // 彩色ルーチン選択
    if (c_in < 1 || c_in > 3) c_in = 1;

    Graph G = Read_MM_UD(argv[1]);    // 疎行列の隣接グラフ（無向グラフ）

    unsigned long seed = 42u;  // デフォルトのシード値
    if (argc >= 5) {
        char* seed_end = nullptr;
        errno = 0;
        seed = std::strtoul(argv[4], &seed_end, 10);
        if (errno != 0 || seed_end == argv[4] || *seed_end != '\0') {
            std::fprintf(stderr, "invalid seed: %s\n", argv[4]);
            return 1;
        }
    }

    auto after_read_to_write_begin = Clock::now();

    double leiden_quality = 0.0;

    std::vector<int> block_of = Leiden_from_Graph(G, resolution, seed, leiden_quality);
    std::vector<int> sizes;
    auto dense = relabel_dense(block_of, &sizes);
    int nb = sizes.size();   // コミュニティの数（= ブロック数）
    
    //////////////////////////////////////////////
    // ブロックグラフの作成
    Graph T = BuildBlockGraph(G, dense, BlockEdgeWeight::Binary);

    //////////////////////////////////////////////
    // ブロックグラフの彩色
    std::vector<int> block_color;
    int nc = 0;

    switch (c_in) {
    case 1:
        nc = Greedy_Coloring(T, block_color);
        break;

    case 2:
        nc = Greedy_Coloring_Balanced(T, block_color);
        break;

    case 3:
        nc = Greedy_Coloring_FixedNc_Balanced(T, block_color);
        break;

    default:
        nc = Greedy_Coloring(T, block_color);
        break;
    }

    // 色ラベルを頻度順に付け替え
    RelabelColorsByClassSize(block_color);

    auto after_read_to_write_end = Clock::now();

    // 出力ファイル名は <入力行列のstem>.blk, <stem>.bcol
    std::string stem = file_stem(argv[1]);
#ifdef CPM
    stem += "_leiden_cpm";
#else
    stem += "_leiden";
#endif

    stem += "_r" + sanitize_resolution_for_filename(argv[2]);
    stem += "_c" + std::to_string(c_in);

    std::string blk_path  = stem + ".blk";
    std::string bcol_path = stem + ".bcol";

    // ブロック情報データの出力
    WriteBlockInfo_1Based(dense, blk_path);

    // ブロック色情報データの出力
    WriteBlockColor_1Based(block_color, nc, bcol_path);

    // --- モジュラリティ（未加重）
    double Q = Modularity_Unweighted(G, dense);

    std::cout << "time = " << elapsed_seconds(after_read_to_write_begin, after_read_to_write_end)
              << ", NB = " << nb
              << ", NC = " << nc
              << ", gamma = " << resolution
              << ", seed = " << seed
              << ", quality = " << leiden_quality << "\n";
            //   << ", modularity = " << Q << "\n";
    return 0;
}
