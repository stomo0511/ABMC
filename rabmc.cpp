#include <queue>
#include <algorithm>
#include <numeric>
#include <iostream>
#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/Coloring.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"
#include "abmc.hpp"
#include "rabbit_order.hpp"

using namespace rabbit_order::aux; // 型記述を簡略化するため aux 名前空間を使用

// --- ヘルパー関数: 指定したサブツリーの全ノードを特定のブロックに割り当てる ---
void AssignSubtreeToBlock(int v, int block_id, const graph& rg, std::vector<int>& block_of) {
    if (v == vmax || block_of[v] != -1) return;

    block_of[v] = block_id;

    // 子ノードを再帰的に処理
    vint child = rg.vs[v].a->child;
    while (child != vmax) {
        AssignSubtreeToBlock(child, block_id, rg, block_of);
        child = rg.vs[child].sibling;
    }
}

// --- ステップ 4: Rabbit Order の階層構造から ABMC ブロックを生成する ---
void PartitionDendrogram(
    int v, int& current_block_id, int block_size, 
    const graph& rg, const std::vector<int>& sizes, 
    std::vector<int>& block_of
) {
    if (v == vmax) return;

    // もしこのサブツリーが目標サイズ s 以下なら、一つのブロックとして確定
    if (sizes[v] <= block_size) {
        AssignSubtreeToBlock(v, current_block_id++, rg, block_of);
        return;
    }

    // サイズが s より大きい場合は、さらに深く探索
    vint child = rg.vs[v].a->child;
    if (child == vmax) {
        // 葉ノード（原子ノード）だが何らかの理由で単独で block_size を超えることはないはずだが、
        // 安全のためにここで割り当てる
        if (block_of[v] == -1) block_of[v] = current_block_id++;
    } else {
        // 子ノードを走査
        while (child != vmax) {
            PartitionDendrogram(child, current_block_id, block_size, rg, sizes, block_of);
            child = rg.vs[child].sibling;
        }
        // 子の処理が終わった後、親(v)自身が未割当なら現在のブロックに入れる（余り物処理）
        if (block_of[v] == -1) block_of[v] = current_block_id++;
    }
}

// --- Rabbit 用の変換関数（提供済み） ---
std::vector<std::vector<edge>> ConvertToRabbitAdj(const Graph& G) {
    const int n = (int)boost::num_vertices(G);
    auto weightmap = boost::get(boost::edge_weight, G);
    std::vector<std::vector<edge>> adj(n);
    for (int u = 0; u < n; ++u) {
        boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::out_edges(u, G); ei != ei_end; ++ei) {
            int v = (int)boost::target(*ei, G);
            float weight = std::abs(static_cast<float>(weightmap[*ei]));
            adj[u].push_back({static_cast<vint>(v), weight});
            // adj[u].push_back({static_cast<vint>(v), (float)weightmap[*ei]});
        }
    }
    return adj;
}

// --- 絶縁ノード判定（提供済み） ---
std::vector<bool> IdentifyInsularNodes(const graph& rg) {
    const int n = rg.n();
    std::vector<bool> is_insular(n, true);
    #pragma omp parallel for
    for (int u = 0; u < n; ++u) {
        int u_com = rg.coms[u];
        for (const auto& e : rg.es[u]) {
            if (rg.coms[e.first] != u_com) { is_insular[u] = false; break; }
        }
    }
    return is_insular;
}

// --- サブツリーサイズ計算（提供済み） ---
int CalculateSubtreeSizes(int v, const graph& rg, std::vector<int>& sizes) {
    if (v == vmax) return 0;
    int size = 1;
    vint child = rg.vs[v].a->child;
    while (child != vmax) {
        size += CalculateSubtreeSizes(child, rg, sizes);
        child = rg.vs[child].sibling;
    }
    return sizes[v] = size;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <matrix.mtx> <# blocks>\n", argv[0]);
        return 1;
    }

    // 1. 行列の読み込みとグラフ構築
    Graph G = Read_MM_UD(argv[1]);
    int N = (int)num_vertices(G);
    int B_target = std::atoi(argv[2]);
    int bs = (N + B_target - 1) / B_target;

    // int isolated_count = 0;
    // for (int i = 0; i < N; ++i) {
    //     if (boost::degree(i, G) == 0) isolated_count++;
    // }
    // std::cout << "Debug: Isolated nodes count = " << isolated_count << std::endl;

    // 2. Rabbit Order 実行フェーズ
    auto adj = ConvertToRabbitAdj(G);

    // for (int u = 0; u < N; ++u) {
    //     float row_sum = 0;
    //     for (const auto& e : adj[u]) row_sum += e.second;
    //     if (row_sum < 0) {
    //         std::cout << "Debug: Node " << u << " has negative strength: " << row_sum << std::endl;
    //     }
    // }

    graph rg = aggregate(std::move(adj));

    // 3. 絶縁ノード同定とサイズ計算
    auto is_insular = IdentifyInsularNodes(rg);
    std::vector<int> subtree_sizes(N, 0);
    if (rg.tops) {
        for (vint top : *rg.tops) {
            CalculateSubtreeSizes(top, rg, subtree_sizes);
        }
    }

    // 4. デンドログラム・パティショニング（ブロック化）
    BlockPartition part;
    part.n = N;
    part.s = bs;
    part.block_of.assign(N, -1);
    int current_block_id = 0;

    if (rg.tops) {
        for (vint top : *rg.tops) {
            PartitionDendrogram(top, current_block_id, bs, rg, subtree_sizes, part.block_of);
        }
    }
    part.nb = current_block_id;

    // BlockPartition の blocks 配列（vector<vector<int>>）を再構成
    part.blocks.resize(part.nb);
    for (int i = 0; i < N; ++i) {
        if (part.block_of[i] != -1) {
            part.blocks[part.block_of[i]].push_back(i);
        }
    }

    // 5. ブロックグラフの構築と彩色
    Graph T = BuildBlockGraph(G, part.block_of, BlockEdgeWeight::Binary);
    
    std::vector<int> block_color;
    int nc = Greedy_Coloring(T, block_color);
    RelabelColorsByClassSize(block_color);

    // 6. 結果出力
    std::string stem = file_stem(argv[1]);
    std::string blk_path = stem + "_rabbit.blk";
    std::string bcol_path = stem + "_rabbit.bcol";

    WriteBlockInfo_1Based(part.block_of, blk_path);
    WriteBlockColor_1Based(block_color, nc, bcol_path);

    double Q = Modularity_Unweighted(G, part.block_of);
    std::printf("Rabbit-ABMC Completed.\n#Blocks: %d, #Colors: %d, Modularity: %.6f\n", part.nb, nc, Q);

    return 0;
}
