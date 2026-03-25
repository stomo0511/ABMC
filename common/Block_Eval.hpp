#pragma once
#include <vector>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include "Types.hpp"
#include "Coloring.hpp"
#include "BlockIO.hpp"

// 内部エッジ数
inline std::vector<int>
CountInternalEdges(const Graph& G, const std::vector<int>& block_of, int nb)
{
    std::vector<int> internal(nb, 0);
    for (auto it = edges(G); it.first != it.second; ++it.first) {
        auto e = *it.first;
        int u = (int)source(e, G);
        int v = (int)target(e, G);
        int bu = block_of[u], bv = block_of[v];
        if (bu >= 0 && bu == bv)
            internal[bu] += 1;
    }
    return internal;
}

// ノード数
inline std::vector<int>
CountNodesPerBlock(const std::vector<int>& block_of, int nb)
{
    std::vector<int> cnt(nb, 0);
    for (int b : block_of) if (b >= 0) cnt[b] += 1;
    return cnt;
}

// ブロックグラフ上の次数 (Binary モデル)
inline std::vector<int> BlockDegreesBinary(const Graph& T)
{
    int nb = (int)num_vertices(T);
    std::vector<int> deg(nb, 0);
    for (int b = 0; b < nb; ++b)
        deg[b] = (int)boost::degree(b, T);
    return deg;
}

// モジュラリティ (未加重)
inline double
Modularity_Unweighted(const Graph& G, const std::vector<int>& block_of)
{
    const int nb = 1 + *std::max_element(block_of.begin(), block_of.end());
    const double m = static_cast<double>(num_edges(G));
    if (m == 0.0) return 0.0;

    std::vector<int> Kb(nb, 0);
    std::vector<int> Lb(nb, 0);

    for (auto v = 0u; v < num_vertices(G); ++v) {
        int b = block_of[(int)v];
        Kb[b] += (int)boost::degree(v, G);
    }
    for (auto it = edges(G); it.first != it.second; ++it.first) {
        auto e = *it.first;
        int u = (int)source(e, G);
        int v = (int)target(e, G);
        int bu = block_of[u], bv = block_of[v];
        if (bu == bv) Lb[bu] += 1;
    }

    double Q = 0.0;
    for (int b = 0; b < nb; ++b) {
        double L = (double)Lb[b];
        double K = (double)Kb[b];
        Q += (L / m) - (K / (2.0*m)) * (K / (2.0*m));
    }
    return Q;
}