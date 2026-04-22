#include <queue>
#include <algorithm>
#include <numeric>
#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/Coloring.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"
#include "abmc.hpp"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <matrix.mtx> <# blocks> [policy(1|2|3)]\n", argv[0]);
        return 1;
    }

    Graph G = Read_MM_UD(argv[1]);    // 疎行列の隣接グラフ（無向グラフ）
    int N  = num_vertices(G);         // ノード数
    int B = std::atoi(argv[2]);       // ブロック数
    int bs = N % B == 0 ? N / B : N / B + 1;  // ブロックサイズ

    int p_in = (argc >= 4) ? std::atoi(argv[3]) : 1;  // ブロック化ポリシー
    if (p_in < 1 || p_in > 3) p_in = 1;

    // ブロック化ポリシー設定
    BlockPolicy policy = static_cast<BlockPolicy>(p_in);

    // ブロック化
    BlockPartition part = ABMC_Blocking(G, bs, policy);

    //////////////////////////////////////////////
    // デバッグ用出力
    // std::printf("ABMC blocking with %s\n", policy_name(policy));
    // std::printf("#blocks=%d, block size=%d\n", B, part.s);
    std::cout << "NB = " << B;
    // DumpBlocks(part);

    //////////////////////////////////////////////
    // ブロックグラフの作成
    Graph T = BuildBlockGraph(G, part.block_of, BlockEdgeWeight::Binary);

    //////////////////////////////////////////////
    // ブロック内結合度の評価
    auto internal = CountInternalEdges(G, part.block_of, part.nb);
    double total_avg = 0.0;
    for (int b = 0; b < part.nb; ++b) {
        int n = (int)part.blocks[b].size();
        double avg_deg = (n > 0) ? 2.0 * internal[b] / n : 0.0;
        total_avg += avg_deg;
        // std::cout << "Block " << b << ": nodes=" << n
        //           << ", internal_edges=" << internal[b]
        //           << ", avg_deg=" << avg_deg << "\n";
    }
    // std::cout << "Total average degree: " << (part.nb > 0 ? total_avg / part.nb : 0.0) << "\n";

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
    Graph T_bin = BuildBlockGraph(G, part.block_of, BlockEdgeWeight::Binary);
    // 各ブロックの次数を計算
    total_avg = 0.0;
    for (int b = 0; b < part.nb; ++b) {
        int deg = boost::degree(b, T_bin);
        // std::cout << "Block " << b << ": degree=" << deg << "\n";
        total_avg += deg;
    }
    // std::cout << "Block graph average degree: "
    //           << (part.nb > 0 ? total_avg / part.nb : 0.0) << "\n";

    //////////////////////////////////////////////
    // ブロックグラフの彩色
    std::vector<int> block_color;
    int nc = Greedy_Coloring(T, block_color);

    // 色ラベルを頻度順に付け替え
    RelabelColorsByClassSize(block_color);

    //////////////////////////////////////////////
    // デバッグ用出力
    // DumpBlockEdges(T);
    // DumpBlockAdjacency(T);

    // 出力ファイル名は <入力行列のstem>.blk, <stem>.bcol
    std::string stem = file_stem(argv[1]);
    stem += "_abmc";
    stem += "_B" + std::to_string(B);  // # blocks
    stem += "_p" + std::to_string(p_in);    // policy
    std::string blk_path  = stem + ".blk";
    std::string bcol_path = stem + ".bcol";

    // ブロック情報データの出力
    WriteBlockInfo_1Based(part.block_of, blk_path);

    // ブロック色情報データの出力
    WriteBlockColor_1Based(block_color, nc, bcol_path);
    std::cout << " NC = " << nc << "\n";

    // --- モジュラリティ（未加重）
    double Q = Modularity_Unweighted(G, part.block_of);
    // std::printf("Modularity (unweighted)   = %.6f\n", Q);
    return 0;
}
