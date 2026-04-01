#include <queue>
#include <algorithm>
#include <numeric>
#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/Coloring.hpp"


int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <matrix.mtx> \n", argv[0]);
        return 1;
    }

    Graph G = Read_MM_UD(argv[1]);    // 疎行列の隣接グラフ（無向グラフ）
    int N = boost::num_vertices(G);   // ノード数
    int E = boost::num_edges(G);      // エッジ数

    // std::cout << "N = " << N << ", E = " << E << "\n";

    std::vector<int> coloring;
    int nc = Greedy_Coloring(G, coloring);

    // 色ラベルを頻度順に付け替え
    RelabelColorsByClassSize(coloring);

    std::cout << "nc = " << nc << "\n";

    // // 出力ファイル名は <入力行列のstem>.col
    std::string stem = file_stem(argv[1]);
    stem += "_gmc";
    std::string col_path  = stem + ".col";

    FILE* fp = std::fopen(col_path.c_str(), "w");
    std::fprintf(fp, "%d\n", nc); // 色の数を先頭に出力
    for (int i = 0; i < N; ++i) {
        std::fprintf(fp, "%d %d\n", i + 1, coloring[i] + 1);
    }
    std::fclose(fp);

    return 0;
}
