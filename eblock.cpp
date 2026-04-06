#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iomanip>

#include "common/Types.hpp"
#include "common/mm_io.hpp"
// #include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"
// #include "abmc.hpp"

/**
 * ChatGPTの指摘による評価指標の追加
 */
void EvaluateEdgeLocality(
    const Graph& G,
    const std::vector<int>& block_of,
    int nb)
{
    int internal = 0;
    int external = 0;

    int N = boost::num_vertices(G);

    for(int u=0;u<N;++u){

        boost::graph_traits<Graph>::out_edge_iterator ei,ei_end;

        for(boost::tie(ei,ei_end)=boost::out_edges(u,G);
            ei!=ei_end;++ei){

            int v = boost::target(*ei,G);

            if(u<v){

                if(block_of[u]==block_of[v])
                    internal++;
                else
                    external++;
            }
        }
    }

    double ratio =
        (double)internal /
        (internal + external);

    std::cout<<"Internal edges : "<<internal<<"\n";
    std::cout<<"External edges : "<<external<<"\n";
    std::cout<<"Internal ratio : "<<ratio<<"\n";
}

void EvaluateBlockDensity(
    const Graph& G,
    const std::vector<int>& block_of,
    int nb)
{
    auto nodes_per_block =
        CountNodesPerBlock(block_of,nb);

    auto internal_edges =
        CountInternalEdges(G,block_of,nb);

    double avg_density=0;

    for(int b=0;b<nb;++b){

        int n = nodes_per_block[b];

        if(n<=1)
            continue;

        double max_edges =
            n*(n-1)/2.0;

        double density =
            internal_edges[b] /
            max_edges;

        avg_density+=density;
    }

    avg_density/=nb;

    std::cout
        <<"Avg block density : "
        <<avg_density<<"\n";
}

void EvaluateBlockGraph(const Graph& T)
{
    int nb = boost::num_vertices(T);

    int edges = boost::num_edges(T);

    double density = edges / (nb*(nb-1)/2.0);

    std::cout
        <<"Block graph edges : "
        <<edges<<"\n";

    std::cout
        <<"Block graph density : "
        <<density<<"\n";
}

void EvaluateLocalityScore(
    const Graph& G,
    const std::vector<int>& block_of)
{
    int internal=0;
    int external=0;

    int N =
        boost::num_vertices(G);

    for(int u=0;u<N;++u){

        boost::graph_traits
        <Graph>::out_edge_iterator ei,ei_end;

        for(boost::tie(ei,ei_end)=
            boost::out_edges(u,G);

            ei!=ei_end;++ei){

            int v =
                boost::target(*ei,G);

            if(block_of[u]==block_of[v])
                internal++;
            else
                external++;
        }
    }

    double score =
        internal /
        (double)(external+1);

    std::cout
        <<"Locality score : "
        <<score<<"\n";
}

void EvaluateLoadBalance(
    const std::vector<int>& nodes_per_block)
{
    double mean=0;

    for(int n:nodes_per_block)
        mean+=n;

    mean/=nodes_per_block.size();

    double var=0;

    for(int n:nodes_per_block)
        var+=(n-mean)*(n-mean);

    var/=nodes_per_block.size();

    double stddev=sqrt(var);

    std::cout
        <<"Block size stddev : "
        <<stddev<<"\n";
}

/**
 * ブロックの品質を多角的に評価する
 */
void EvaluatePartitioning(const Graph& G, const std::vector<int>& block_of, int nb, int nc) {
    int N = static_cast<int>(boost::num_vertices(G));

    // 1. ブロック数 (nb)
    // 引数として渡された nb を使用

    // 2. ブロック内のノード数統計 (CountNodesPerBlock を使用)
    std::vector<int> nodes_per_block = CountNodesPerBlock(block_of, nb);
    
    double sum_nodes = std::accumulate(nodes_per_block.begin(), nodes_per_block.end(), 0.0);
    double avg_nodes = sum_nodes / nb;
    int max_nodes = *std::max_element(nodes_per_block.begin(), nodes_per_block.end());
    int min_nodes = *std::min_element(nodes_per_block.begin(), nodes_per_block.end());

    // 3. ブロック内の平均次数の平均
    // 「ブロック内のノードの平均次数」 = (2 * 内部エッジ数) / ブロック内ノード数
    std::vector<int> internal_edges = CountInternalEdges(G, block_of, nb);
    double sum_internal_avg_degree = 0.0;
    for (int b = 0; b < nb; ++b) {
        if (nodes_per_block[b] > 0) {
            // 各ブロック b における内部平均次数
            double local_avg_deg = (2.0 * internal_edges[b]) / nodes_per_block[b];
            sum_internal_avg_degree += local_avg_deg;
        }
    }
    double mean_of_internal_avg_degrees = sum_internal_avg_degree / nb;

    // 4. ブロックグラフの平均次数 (BlockDegreesBinary を使用)
    Graph T = BuildBlockGraph(G, block_of, BlockEdgeWeight::Binary);
    std::vector<int> block_degrees = BlockDegreesBinary(T);
    double sum_block_degrees = std::accumulate(block_degrees.begin(), block_degrees.end(), 0.0);
    double avg_block_degree = sum_block_degrees / nb;

    // 5. モジュラリティ (Modularity_Unweighted を使用)
    double Q = Modularity_Unweighted(G, block_of);

    // 結果表示
    std::cout << "\n==============================================" << std::endl;
    std::cout << "  Partitioning Evaluation Report" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "1. Total Blocks          : " << nb << std::endl;
    std::cout << "   Total Colors          : " << nc << std::endl;
    
    std::cout << "2. Nodes per Block" << std::endl;
    std::cout << "   - Average             : " << std::fixed << std::setprecision(2) << avg_nodes << std::endl;
    std::cout << "   - Maximum             : " << max_nodes << std::endl;
    std::cout << "   - Minimum             : " << min_nodes << std::endl;
    
    std::cout << "3. Internal Avg Degree   : " << std::setprecision(4) << mean_of_internal_avg_degrees 
              << " (Mean of per-block avg)" << std::endl;
              
    std::cout << "4. Block Graph Avg Deg   : " << avg_block_degree << std::endl;
    
    std::cout << "5. Modularity (Q)        : " << std::setprecision(6) << Q << std::endl;
    std::cout << "==============================================" << std::endl;

    // 評価の実行
    EvaluateEdgeLocality(G,block_of,nb);
    EvaluateBlockDensity(G,block_of,nb);
    EvaluateBlockGraph(T);
    EvaluateLocalityScore(G,block_of);
    EvaluateLoadBalance(nodes_per_block);
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <matrix.mtx> <block_file.blk> <color_file.bcol>" << std::endl;
        return 1;
    }

    // 行列の読み込み
    Graph G = Read_MM_UD(argv[1]); 
    int N = static_cast<int>(boost::num_vertices(G));

    // ブロック分けファイルの読み込み (1-based -> 0-based 変換込み)
    std::vector<int> block_of;
    int nb = ReadBlockInfo_1Based(argv[2], N, block_of); 

    // 彩色情報の読み込み
    std::vector<int> block_color;
    int nc = ReadBlockColor_1Based(argv[3], nb, block_color);

    // 評価の実行
    EvaluatePartitioning(G, block_of, nb, nc);

    return 0;
}