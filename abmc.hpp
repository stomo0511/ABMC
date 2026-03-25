#pragma once
#include <queue>
#include <algorithm>

static inline int pick_next_seed(const std::vector<int>& block_of) {
    for (int v = 0; v < (int)block_of.size(); ++v)
        if (block_of[v] == -1) return v;
    return -1;
}

//////////////////////////////////////////
// デバッグ用出力ルーチン
// ブロック内容（含まれる元ノード）をすべて出力
void DumpBlocks(const BlockPartition& part)
{
    std::printf("==== Blocks ====\n");
    for (int k = 0; k < part.nb; ++k) {
        std::printf("Block %d (size %zu):", k, part.blocks[k].size());
        for (int v : part.blocks[k]) std::printf(" %d", v);
        std::printf("\n");
    }
}

// ブロックグラフの全エッジ（bk, bl, weight）を列挙
void DumpBlockEdges(const Graph& T)
{
    std::printf("==== Block edges (undirected) ====\n");
    auto Tw = get(boost::edge_weight, T);
    // 収まりが良いように (min,max) で整列
    std::vector<std::tuple<int,int,double>> E;
    for (auto it = edges(T); it.first != it.second; ++it.first) {
        auto e = *it.first;
        int a = (int)source(e, T), b = (int)target(e, T);
        if (a > b) std::swap(a, b);
        E.emplace_back(a, b, get(Tw, e));
    }
    std::sort(E.begin(), E.end());
    for (auto& [a,b,w] : E)
        std::printf("(%d, %d)  w=%.12g\n", a, b, w);
}

// 各ブロックの隣接リスト（相手ブロックと重み）を出力
void DumpBlockAdjacency(const Graph& T)
{
    std::printf("==== Block adjacency lists ====\n");
    auto Tw = get(boost::edge_weight, T);
    int nb = (int)num_vertices(T);
    for (int k = 0; k < nb; ++k) {
        std::printf("Bk %d:", k);
        auto adj = boost::adjacent_vertices(k, T);
        for (auto it = adj.first; it != adj.second; ++it) {
            int l = (int)*it;
            auto ep = edge(k, l, T);
            double w = ep.second ? get(Tw, ep.first) : 0.0;
            std::printf(" %d(w=%.12g)", l, w);
        }
        std::printf("\n");
    }
}

