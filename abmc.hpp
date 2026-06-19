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
    // 収まりが良いように (min,max) で整列
    std::vector<std::tuple<int,int,double>> E;
    for_each_undirected_edge(T, [&](int a, int b, double w) {
        if (a > b) std::swap(a, b);
        E.emplace_back(a, b, w);
    });
    std::sort(E.begin(), E.end());
    for (auto& [a,b,w] : E)
        std::printf("(%d, %d)  w=%.12g\n", a, b, w);
}

// 各ブロックの隣接リスト（相手ブロックと重み）を出力
void DumpBlockAdjacency(const Graph& T)
{
    std::printf("==== Block adjacency lists ====\n");
    int nb = (int)num_vertices(T);
    for (int k = 0; k < nb; ++k) {
        std::printf("Bk %d:", k);
        for (const auto& e : T.adj[k]) {
            std::printf(" %d(w=%.12g)", e.to, e.weight);
        }
        std::printf("\n");
    }
}
