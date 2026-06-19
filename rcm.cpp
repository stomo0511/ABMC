#include <queue>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cmath>
#include <type_traits>

#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/Coloring.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"

static std::vector<int> reverse_cuthill_mckee_ordering(const Graph& G)
{
    const int n = G.n;
    std::vector<int> order;
    order.reserve(n);
    std::vector<char> visited(n, 0);

    auto less_degree_then_id = [&](int a, int b) {
        const std::size_t da = degree(a, G);
        const std::size_t db = degree(b, G);
        if (da != db) return da < db;
        return a < b;
    };

    for (;;) {
        int start = -1;
        for (int v = 0; v < n; ++v) {
            if (!visited[v] && (start == -1 || less_degree_then_id(v, start)))
                start = v;
        }
        if (start == -1) break;

        std::queue<int> q;
        visited[start] = 1;
        q.push(start);

        while (!q.empty()) {
            const int u = q.front();
            q.pop();
            order.push_back(u);

            std::vector<int> nbrs;
            nbrs.reserve(G.adj[u].size());
            for (const auto& e : G.adj[u]) {
                if (!visited[e.to])
                    nbrs.push_back(e.to);
            }
            std::sort(nbrs.begin(), nbrs.end(), less_degree_then_id);

            for (int v : nbrs) {
                if (visited[v]) continue;
                visited[v] = 1;
                q.push(v);
            }
        }
    }

    std::reverse(order.begin(), order.end());

    std::vector<int> perm(n, -1);
    for (int new_id = 0; new_id < n; ++new_id)
        perm[order[new_id]] = new_id;
    return perm;
}

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

    std::vector<int> perm = reverse_cuthill_mckee_ordering(G);

    // 3. 並べ替えたグラフを Matrix Market 形式で出力
    WriteReorderedMatrixMarket_UD(G, perm, out_mtx);
    // WritePermutation_1Based(perm, out_perm);

    std::cout << "Wrote reordered matrix : " << out_mtx  << "\n";
    // std::cout << "Wrote permutation      : " << out_perm << "\n";

    return 0;
}
