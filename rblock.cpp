#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <numeric>

#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/Coloring.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"
#include "abmc.hpp"

// ------------------------------------------------------------
// Louvain ブロック情報を読む
// 想定形式:
//   - N 個の整数（各ノードの block id）
//   - または先頭にヘッダ1個 + N 個の整数
// block id は 0-based / 1-based のどちらでも可
// ------------------------------------------------------------
std::vector<int> ReadBlockInfoFlexible(const std::string& path, int N)
{
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error("cannot open block file: " + path);
    }

    int nb;
    if (!(ifs >> nb)) {
        throw std::runtime_error("cannot read number of blocks: " + path);
    }

    std::vector<int> block_of(N, -1);

    for (int i = 0; i < N; ++i) {
        int node_id, block_id;
        if (!(ifs >> node_id >> block_id)) {
            throw std::runtime_error("block file format error at line " + std::to_string(i + 2));
        }

        if (node_id != i + 1) {
            throw std::runtime_error("unexpected node id at line " + std::to_string(i + 2));
        }

        block_of[i] = block_id - 1;  // 1-based -> 0-based
    }

    return block_of;
}
// std::vector<int> ReadBlockInfoFlexible(const std::string& path, int N)
// {
//     std::ifstream ifs(path);
//     if (!ifs) {
//         throw std::runtime_error("cannot open block file: " + path);
//     }

//     std::vector<long long> raw;
//     long long x;
//     while (ifs >> x) raw.push_back(x);

//     if (raw.empty()) {
//         throw std::runtime_error("empty block file: " + path);
//     }

//     std::vector<int> block_of;

//     if ((int)raw.size() == N) {
//         block_of.resize(N);
//         for (int i = 0; i < N; ++i) block_of[i] = static_cast<int>(raw[i]);
//     } else if ((int)raw.size() == N + 1) {
//         // 先頭がヘッダとみなせる場合
//         block_of.resize(N);
//         for (int i = 0; i < N; ++i) block_of[i] = static_cast<int>(raw[i + 1]);
//     } else {
//         throw std::runtime_error(
//             "invalid block file length: expected N or N+1 integers, got " +
//             std::to_string(raw.size()));
//     }

//     int mn = *std::min_element(block_of.begin(), block_of.end());
//     if (mn >= 1) {
//         for (int& b : block_of) --b;   // 1-based -> 0-based
//     } else if (mn < 0) {
//         throw std::runtime_error("negative block id found in: " + path);
//     }

//     return block_of;
// }

// ------------------------------------------------------------
// block id を 0..nb-1 に詰め直す
// sizes_out[new_id] = そのブロックサイズ
// ------------------------------------------------------------
std::vector<int> relabel_dense(
    const std::vector<int>& block_of,
    std::vector<int>* sizes_out = nullptr)
{
    const int N = static_cast<int>(block_of.size());

    std::unordered_map<int, int> cnt;
    cnt.reserve(N * 2);
    for (int b : block_of) {
        if (b >= 0) ++cnt[b];
    }

    if (cnt.empty()) {
        if (sizes_out) sizes_out->clear();
        return std::vector<int>(N, -1);
    }

    struct Stat { int old_id; int count; };
    std::vector<Stat> order;
    order.reserve(cnt.size());
    for (const auto& kv : cnt) {
        order.push_back({kv.first, kv.second});
    }

    // 大きい順、同数なら旧 id 昇順
    std::sort(order.begin(), order.end(),
        [](const Stat& a, const Stat& b) {
            if (a.count != b.count) return a.count > b.count;
            return a.old_id < b.old_id;
        });

    std::unordered_map<int, int> old2new;
    old2new.reserve(order.size());

    if (sizes_out) sizes_out->assign(order.size(), 0);

    for (int new_id = 0; new_id < (int)order.size(); ++new_id) {
        old2new[order[new_id].old_id] = new_id;
        if (sizes_out) (*sizes_out)[new_id] = order[new_id].count;
    }

    std::vector<int> dense(N, -1);
    for (int i = 0; i < N; ++i) {
        auto it = old2new.find(block_of[i]);
        if (it != old2new.end()) dense[i] = it->second;
    }
    return dense;
}

// ------------------------------------------------------------
// block_of から各ブロックの頂点集合を作る
// ------------------------------------------------------------
std::vector<std::vector<int>> BuildMembers(const std::vector<int>& block_of, int nb)
{
    std::vector<std::vector<int>> members(nb);
    for (int v = 0; v < (int)block_of.size(); ++v) {
        int b = block_of[v];
        if (b >= 0) members[b].push_back(v);
    }
    return members;
}

// ------------------------------------------------------------
// 頂点集合 verts の誘導部分グラフを作る
// local vertex 0..k-1 <-> global vertex verts[local]
// ------------------------------------------------------------
Graph BuildInducedSubgraph(const Graph& G, const std::vector<int>& verts)
{
    const int k = static_cast<int>(verts.size());
    Graph SG(k);

    std::vector<int> g2l(num_vertices(G), -1);
    for (int i = 0; i < k; ++i) g2l[verts[i]] = i;

    auto wmapG = get(boost::edge_weight, G);
    auto edges_pair = edges(G);

    for (auto it = edges_pair.first; it != edges_pair.second; ++it) {
        auto e = *it;
        int u = static_cast<int>(source(e, G));
        int v = static_cast<int>(target(e, G));

        int lu = g2l[u];
        int lv = g2l[v];
        if (lu >= 0 && lv >= 0) {
            double w = get(wmapG, e);
            boost::add_edge(lu, lv, w, SG);
        }
    }

    return SG;
}

// ------------------------------------------------------------
// 再ブロック化
// ルール:
//   |B| <= 2s : そのまま 1 ブロック
//   |B| >  2s : 誘導部分グラフ上で ABMC(block_size=s) を実行
// ------------------------------------------------------------
std::vector<int> ReblockLargeLouvainBlocks(
    const Graph& G,
    const std::vector<int>& louvain_block_of,
    int s,
    BlockPolicy policy = BlockPolicy::FIFO)
{
    if (s < 1) {
        throw std::runtime_error("block size s must be >= 1");
    }

    std::vector<int> sizes;
    std::vector<int> dense_block_of = relabel_dense(louvain_block_of, &sizes);
    const int nb0 = static_cast<int>(sizes.size());

    auto members = BuildMembers(dense_block_of, nb0);

    std::vector<int> final_block_of(num_vertices(G), -1);
    int next_bid = 0;

    for (int b = 0; b < nb0; ++b) {
        const auto& verts = members[b];
        const int sz = static_cast<int>(verts.size());

        if (sz == 0) continue;

        if (sz <= 2 * s) {
            // しきい値以下はそのまま 1 ブロック
            for (int gv : verts) final_block_of[gv] = next_bid;
            ++next_bid;
        } else {
            // 大ブロックのみ局所 ABMC
            Graph SG = BuildInducedSubgraph(G, verts);
            BlockPartition local = ABMC_Blocking(SG, s, policy);

            for (int lv = 0; lv < (int)verts.size(); ++lv) {
                int gv = verts[lv];
                final_block_of[gv] = next_bid + local.block_of[lv];
            }
            next_bid += local.nb;
        }
    }

    return final_block_of;
}

// ------------------------------------------------------------
// main
// usage:
//   rblock <matrix.mtx> <louvain.blk> <s>
// ------------------------------------------------------------
int main(int argc, char** argv)
{
    if (argc < 4) {
        std::fprintf(stderr, "usage: %s <matrix.mtx> <louvain.blk> <block_size_s>\n", argv[0]);
        return 1;
    }

    const std::string matrix_path = argv[1];
    const std::string louvain_blk_path = argv[2];
    const int s = std::atoi(argv[3]);

    if (s < 1) {
        std::fprintf(stderr, "error: block_size_s must be >= 1\n");
        return 1;
    }

    // 入力グラフ
    Graph G = Read_MM_UD(matrix_path);
    const int N = static_cast<int>(num_vertices(G));

    // Louvain ブロック情報
    std::vector<int> louvain_block_of = ReadBlockInfoFlexible(louvain_blk_path, N);

    // 再ブロック化
    // ここでは ABMC の policy は FIFO に固定
    std::vector<int> final_block_of =
        ReblockLargeLouvainBlocks(G, louvain_block_of, s, BlockPolicy::FIFO);

    // 念のため dense に詰め直す
    std::vector<int> sizes;
    final_block_of = relabel_dense(final_block_of, &sizes);
    const int nb = static_cast<int>(sizes.size());

    // --------------------------------------------------------
    // 評価・彩色
    // --------------------------------------------------------
    Graph T = BuildBlockGraph(G, final_block_of, BlockEdgeWeight::Binary);

    auto internal = CountInternalEdges(G, final_block_of, nb);
    double total_avg = 0.0;
    for (int b = 0; b < nb; ++b) {
        double avg_deg = (sizes[b] > 0) ? 2.0 * internal[b] / sizes[b] : 0.0;
        total_avg += avg_deg;
    }
    std::cout << "Total average degree: "
              << (nb > 0 ? total_avg / nb : 0.0) << "\n";

    total_avg = 0.0;
    for (int b = 0; b < nb; ++b) {
        int deg = boost::degree(b, T);
        total_avg += deg;
    }
    std::cout << "Block graph average degree: "
              << (nb > 0 ? total_avg / nb : 0.0) << "\n";

    std::vector<int> block_color;
    int nc = Greedy_Coloring(T, block_color);
    RelabelColorsByClassSize(block_color);

    // --------------------------------------------------------
    // 出力
    // --------------------------------------------------------
    std::string stem = file_stem(matrix_path);
    stem += "_rblock";
    stem += "_s" + std::to_string(s);

    std::string blk_path  = stem + ".blk";
    std::string bcol_path = stem + ".bcol";

    WriteBlockInfo_1Based(final_block_of, blk_path);
    WriteBlockColor_1Based(block_color, nc, bcol_path);

    // --------------------------------------------------------
    // 統計出力
    // --------------------------------------------------------
    int louvain_nb = 0;
    {
        std::vector<int> tmp_sizes;
        auto dense0 = relabel_dense(louvain_block_of, &tmp_sizes);
        (void)dense0;
        louvain_nb = static_cast<int>(tmp_sizes.size());

        int n_large = 0;
        for (int sz : tmp_sizes) {
            if (sz > 2 * s) ++n_large;
        }

        std::cout << "Input Louvain blocks: " << louvain_nb << "\n";
        std::cout << "Reblocked large blocks (> 2s): " << n_large << "\n";
    }

    std::cout << "Output blocks: " << nb << "\n";
    std::cout << "Block colors: " << nc << "\n";

    double Q = Modularity_Unweighted(G, final_block_of);
    std::printf("Modularity (unweighted) = %.6f\n", Q);

    // Q = Modularity_Unweighted(G, final_block_of);
    // std::printf("Modularity (weighted)   = %.6f\n", Q);

    return 0;
}