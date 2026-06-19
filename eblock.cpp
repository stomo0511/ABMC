#include <iostream>
#include <vector>

#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"

int main(int argc, char** argv)
{
    if(argc < 4){
        std::cerr
            << "Usage: "
            << argv[0]
            << " <matrix.mtx> <block_file.blk> <color_file.bcol>"
            << std::endl;

        return 1;
    }

    Graph G = Read_MM_UD(argv[1]);
    int N = num_vertices(G);

    std::vector<int> block_of;
    int nb = ReadBlockInfo_1Based(argv[2], N, block_of);

    std::vector<int> block_color;
    int nc = ReadBlockColor_1Based(argv[3], nb, block_color);

    EvaluatePartitioning(G, block_of, block_color, nb, nc);

    return 0;
}
