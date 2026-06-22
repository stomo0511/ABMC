#include <algorithm>
#include <iomanip>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

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
std::vector<int> Leiden_from_Graph(const Graph& G, double resolution) {
    igraph_rng_seed(igraph_rng_default(), 42u);

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

    struct Stat { int label; int count; };
    std::unordered_map<int, int> cnt;
    cnt.reserve(N * 2);
    for (int c : comm) {
        if (c >= 0) ++cnt[c];
    }
    if (cnt.empty()) {
        if (sizes_out) sizes_out->clear();
        return std::vector<int>(N, -1);
    }

    std::vector<Stat> order;
    order.reserve(cnt.size());
    for (auto& kv : cnt) order.push_back({kv.first, kv.second});
    std::sort(order.begin(), order.end(), [](const Stat& a, const Stat& b) {
        if (a.count != b.count) return a.count > b.count;
        return a.label < b.label;
    });

    std::unordered_map<int, int> old2new;
    old2new.reserve(order.size());
    if (sizes_out) sizes_out->assign(order.size(), 0);
    for (int new_id = 0; new_id < static_cast<int>(order.size()); ++new_id) {
        old2new[order[new_id].label] = new_id;
        if (sizes_out) (*sizes_out)[new_id] = order[new_id].count;
    }

    std::vector<int> dense(N, -1);
    for (size_t i = 0; i < N; ++i) {
        int c = comm[i];
        if (c >= 0) dense[i] = old2new[c];
    }
    return dense;
}

struct LeidenResult {
    double gamma;
    std::vector<int> dense;
    std::vector<int> block_color;
    std::vector<int> sizes;
    int nb;
    int nc;
    double modularity;
    double R;
};

LeidenResult RunLeidenAndColoring(const Graph& G, double gamma) {
    LeidenResult result;
    result.gamma = gamma;

    const std::vector<int> block_of = Leiden_from_Graph(G, gamma);
    result.dense = relabel_dense(block_of, &result.sizes);
    result.nb = static_cast<int>(result.sizes.size());

    Graph T = BuildBlockGraph(G, result.dense, BlockEdgeWeight::Binary);
    result.nc = Greedy_Coloring(T, result.block_color);
    if (result.nc == 0) {
        throw std::runtime_error("Greedy_Coloring returned zero colors");
    }
    RelabelColorsByClassSize(result.block_color);

    result.modularity = Modularity_Unweighted(G, result.dense);
    result.R = static_cast<double>(result.nb) /
               static_cast<double>(result.nc);
    return result;
}

int main(int argc, char** argv) {
    const int p = 6;              // # threads
    const double mu = 2.0;         // safety factor
    const int t_max = 4;           // max # of trials
    const double gamma0 = 1.0e-3;  // initial guess for gamma

    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <matrix.mtx>\n", argv[0]);
        return 1;
    }

    Graph G = Read_MM_UD(argv[1]);
    const auto search_begin = Clock::now();
    const double threshold = mu * static_cast<double>(p);

    double gamma = gamma0;
    int trial_count = 0;
    std::optional<LeidenResult> best;
    std::optional<LeidenResult> last;

    auto run_trial = [&]() -> bool {
        LeidenResult current = RunLeidenAndColoring(G, gamma);
        ++trial_count;
        const bool satisfied = current.R >= threshold;

        std::cout << "trial = " << trial_count << "\n"
                  << "gamma = " << std::setprecision(17) << current.gamma << "\n"
                  << "NB = " << current.nb << "\n"
                  << "NC = " << current.nc << "\n"
                  << "R = " << current.R << "\n"
                  << "condition = "
                  << (satisfied ? "satisfied" : "not_satisfied") << "\n";

        last = std::move(current);
        return satisfied;
    };

    if (run_trial()) {
        best = last;
        while (trial_count < t_max) {
            gamma /= 2.0;
            if (run_trial()) {
                best = last;
            } else {
                break;
            }
        }
    } else {
        while (trial_count < t_max) {
            gamma *= 2.0;
            if (run_trial()) {
                best = last;
                break;
            }
        }
    }

    const LeidenResult& selected = best ? *best : *last;
    const auto search_end = Clock::now();

    std::string stem = file_stem(argv[1]);
#ifdef CPM
    stem += "_ite_leiden_cpm";
#else
    stem += "_ite_leiden";
#endif
    WriteBlockInfo_1Based(selected.dense, stem + ".blk");
    WriteBlockColor_1Based(selected.block_color, selected.nc, stem + ".bcol");

    std::cout << "time = " << elapsed_seconds(search_begin, search_end) << "\n"
              << "selected_gamma = " << selected.gamma << "\n"
              << "trials = " << trial_count << "\n"
              << "NB = " << selected.nb << "\n"
              << "NC = " << selected.nc << "\n"
              << "R = " << selected.R << "\n"
              << "threshold = " << threshold << "\n"
              << "modularity = " << selected.modularity << "\n";
    return 0;
}
