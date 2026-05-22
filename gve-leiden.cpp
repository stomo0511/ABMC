#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "common/Types.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"
#include "common/Coloring.hpp"
#include "common/mm_io.hpp"
#include "gve-leiden-inc/main.hxx"

std::vector<int> relabel_dense(
    const std::vector<int>& comm,
    std::vector<int>* sizes_out = nullptr)
{
    const size_t N = comm.size();

    struct Stat {
        int label;
        int count;
    };

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
    for (const auto& entry : cnt) {
        order.push_back({entry.first, entry.second});
    }
    std::sort(order.begin(), order.end(), [](const Stat& a, const Stat& b) {
        if (a.count != b.count) return a.count > b.count;
        return a.label < b.label;
    });

    std::unordered_map<int, int> old_to_new;
    old_to_new.reserve(order.size());
    if (sizes_out) sizes_out->assign(order.size(), 0);

    for (int new_id = 0; new_id < static_cast<int>(order.size()); ++new_id) {
        old_to_new[order[new_id].label] = new_id;
        if (sizes_out) (*sizes_out)[new_id] = order[new_id].count;
    }

    std::vector<int> dense(N, -1);
    for (size_t i = 0; i < N; ++i) {
        const int c = comm[i];
        if (c >= 0) dense[i] = old_to_new[c];
    }
    return dense;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <matrix.mtx> [symmetric(0|1)] [weighted(0|1)]\n", argv[0]);
        return 1;
    }

    const char* file = argv[1];
    const bool already_symmetric = argc > 2 ? std::atoi(argv[2]) != 0 : false;
    const bool weighted = argc > 3 ? std::atoi(argv[3]) != 0 : false;

    std::ifstream in(file);
    if (!in) {
        std::cerr << "Error: Cannot open file " << file << "\n";
        return 1;
    }

    bool matrix_symmetric = false;
    size_t rows = 0;
    size_t cols = 0;
    size_t nnz = 0;
    readMtxHeader(in, matrix_symmetric, rows, cols, nnz);

    if (rows == 0 || cols == 0) {
        std::cerr << "Error: Invalid Matrix Market coordinate file: " << file << "\n";
        return 1;
    }

    using K = uint32_t;
    using V = float;

    DiGraph<K, None, V> graph;
    readMtxOmpW(graph, file, weighted);
    if (!already_symmetric) {
        graph = symmetricizeOmp(graph);
    }

    std::cout << "rows = " << rows
              << " cols = " << cols
              << " nnz = " << nnz
              << " matrix_symmetric = " << (matrix_symmetric ? 1 : 0) << "\n";
    std::cout << "graph_order = " << graph.order()
              << " graph_size = " << graph.size()
              << " graph_span = " << graph.span() << "\n";

    auto result = leidenStaticOmp(graph);

    Graph boost_graph = Read_MM_UD(file);
    const size_t boost_order = num_vertices(boost_graph);

    std::vector<int> block_of(boost_order, -1);
    for (size_t u = 0; u < boost_order; ++u) {
        const size_t gve_vertex = u + 1;
        if (gve_vertex < result.membership.size()) {
            block_of[u] = static_cast<int>(result.membership[gve_vertex]);
        }
    }

    std::vector<int> block_sizes;
    std::vector<int> dense = relabel_dense(block_of, &block_sizes);

    const size_t community_count = block_sizes.size();
    size_t min_community_size = community_count > 0
        ? static_cast<size_t>(block_sizes[0])
        : 0;
    size_t max_community_size = 0;
    for (int block_size : block_sizes) {
        const size_t size = static_cast<size_t>(block_size);
        if (size < min_community_size) min_community_size = size;
        if (size > max_community_size) max_community_size = size;
    }

    std::cout << "community_count = " << community_count
              << " max_community_size = " << max_community_size
              << " min_community_size = " << min_community_size << "\n";

    const double total_edge_weight = edgeWeightOmp(graph) / 2.0;
    const double modularity = total_edge_weight > 0.0
        ? modularityByOmp(graph, [&](K u) { return result.membership[u]; }, total_edge_weight)
        : 0.0;

    std::cout << "modularity = " << modularity << "\n";

    Graph block_graph = BuildBlockGraph(boost_graph, dense, BlockEdgeWeight::Binary);

    std::vector<int> block_color;
    const int color_count = Greedy_Coloring(block_graph, block_color);
    RelabelColorsByClassSize(block_color);

    std::string stem = file_stem(file);
    stem += "_gveleiden";
    const std::string blk_path = stem + ".blk";
    const std::string bcol_path = stem + ".bcol";

    WriteBlockInfo_1Based(dense, blk_path);
    WriteBlockColor_1Based(block_color, color_count, bcol_path);

    std::cout << "color_count = " << color_count << "\n"
              << "blk_path = " << blk_path << "\n"
              << "bcol_path = " << bcol_path << "\n";

    return 0;
}
