#include <queue>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cmath>
#include <type_traits>
#include <boost/graph/cuthill_mckee_ordering.hpp>

#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/Coloring.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <matrix.mtx> [output.mtx]\n", argv[0]);
        return 1;
    }

    const std::string in_path  = argv[1];
    const std::string stem     = file_stem(in_path);
    const std::string out_mtx  = (argc >= 3) ? argv[2] : (stem + "_rcm.mtx");
    const std::string out_perm = stem + "_rcm.perm";

    // Matrix Market の無向重み付き読み込み
    Graph G = Read_MM_UD(in_path.c_str());

    // Cuthill-McKee 順序の計算 <- RCMではない
    std::vector<int> perm(boost::num_vertices(G));
    boost::cuthill_mckee_ordering(G, perm.rbegin());

    // 3. 並べ替えたグラフを Matrix Market 形式で出力
    WriteReorderedMatrixMarket_UD(G, perm, out_mtx);
    // WritePermutation_1Based(perm, out_perm);

    std::cout << "Wrote reordered matrix : " << out_mtx  << "\n";
    // std::cout << "Wrote permutation      : " << out_perm << "\n";

    return 0;
}