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
inline double Modularity_Unweighted(const Graph& G, const std::vector<int>& block_of)
{
    // 1. ブロック数の確定
    const int nb = 1 + *std::max_element(block_of.begin(), block_of.end());

    // 2. 総エッジ数 m
    const double m = static_cast<double>(num_edges(G));
    if (m == 0.0) return 0.0;

    // 3. 各ブロックの次数合計 Kb と内部エッジ数 Lbのカウント
    std::vector<int> Kb(nb, 0);
    std::vector<int> Lb(nb, 0);

    // ノード毎の次数を積算
    for (auto v = 0u; v < num_vertices(G); ++v) {
        int b = block_of[(int)v];
        if (b >= 0)
            Kb[b] += (int)boost::degree(v, G);
    }
    // 内部エッジをカウント
    for (auto it = edges(G); it.first != it.second; ++it.first) {
        auto e = *it.first;
        int u = (int)source(e, G);
        int v = (int)target(e, G);
        int bu = block_of[u], bv = block_of[v];
        if (bu >= 0 && bu == bv)
            Lb[bu] += 1;
    }

    // 4. モジュラリティの計算
    double Q = 0.0;
    for (int b = 0; b < nb; ++b) {
        double L = (double)Lb[b];
        double K = (double)Kb[b];
        Q += (L / m) - (K / (2.0*m)) * (K / (2.0*m));
    }
    return Q;
}

// モジュラリティ (加重)
inline double Modularity_Weighted(const Graph& G, const std::vector<int>& block_of) {
    const int nb = 1+*std::max_element(block_of.begin(),block_of.end());

    std::vector<double> Kb(nb,0.0);
    std::vector<double> Lb(nb,0.0);

    auto wmap = boost::get(boost::edge_weight,G);
    double W = 0.0;

    boost::graph_traits<Graph>::edge_iterator
        ei,ei_end;

    for(boost::tie(ei,ei_end)=boost::edges(G); ei!=ei_end;++ei){
        int u=boost::source(*ei,G);
        int v=boost::target(*ei,G);

        if (u==v) continue; // 自己ループは無視

        double w = std::abs((double)wmap[*ei]);
        W+=w;
    }

    if(W == 0.0)
        return 0.0;

    // weighted degree
    for(int u=0; u<(int)num_vertices(G); ++u){
        int b = block_of[u];
        double deg = 0.0;
        auto nbrs = boost::out_edges(u,G);

        for(auto ei=nbrs.first; ei!=nbrs.second; ++ei){
            deg += std::abs((double)wmap[*ei]);
        }
        Kb[b]+=deg;
    }

    // internal weight
    for(boost::tie(ei,ei_end)=boost::edges(G); ei!=ei_end;++ei){
        int u=boost::source(*ei,G);
        int v=boost::target(*ei,G);

        if (u==v) continue; // 自己ループは無視

        if(block_of[u]==block_of[v]){
            double w=std::abs((double)wmap[*ei]);
            Lb[block_of[u]]+=w;
        }
    }

    double Q=0.0;

    for(int b=0;b<nb;++b) {
        double L=Lb[b];
        double K=Kb[b];

        Q += (L/W) - (K/(2.0*W))*(K/(2.0*W));
    }

    return Q;
}
