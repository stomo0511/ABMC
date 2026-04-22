#include <queue>
#include <algorithm>
#include <numeric>
#include "common/Types.hpp"
#include "common/mm_io.hpp"
#include "common/Coloring.hpp"
#include "common/BlockIO.hpp"
#include "common/Block_Eval.hpp"
#include "abmc.hpp"

// ブロック化（ポリシー選択あり）
BlockPartition ABMC_Blocking(const Graph& G, int block_size, BlockPolicy policy)
{
    const int N = (int)boost::num_vertices(G);
    if (N == 0) return {};
    if (block_size < 1) block_size = 1;

    BlockPartition out;
    out.n = N;
    out.s = block_size;
    out.block_of.assign(N, -1);

    // WeightMap wmap = get(boost::edge_weight, G);
    auto wmap = get(boost::edge_weight, G);

    int assigned = 0;
    while (assigned < N) {
        const int seed0 = pick_next_seed(out.block_of); // 最小番号の未割当ノード（論文のseed）:contentReference[oaicite:2]{index=2}
        if (seed0 < 0) break;

        // このブロックの目標サイズ（末尾ブロックは小さくてもOK）:contentReference[oaicite:3]{index=3}
        const int target = std::min(block_size, N - assigned);

        const int k = (int)out.blocks.size();
        out.blocks.emplace_back();
        std::vector<char> in_block(N, 0);
        std::vector<char> in_cand(N, 0);
        std::vector<int>  cand;  // 候補集合 C（重複なし）:contentReference[oaicite:4]{index=4}

        auto add_to_block = [&](int v){
            out.block_of[v] = k;
            out.blocks.back().push_back(v);
            in_block[v] = 1;
            ++assigned;

            // C ← C ∪ adj(v) （未割当のみ追加）:contentReference[oaicite:5]{index=5}
            auto nbrs = boost::adjacent_vertices(v, G);
            for (AdjacencyIter it = nbrs.first; it != nbrs.second; ++it) {
                int u = (int)*it;
                if (out.block_of[u] == -1 && !in_cand[u]) {
                    in_cand[u] = 1;
                    cand.push_back(u);
                }
            }
        };

        // Step 2: Cが空ならseed、そうでなければCからpolicyに従って選ぶ:contentReference[oaicite:6]{index=6}
        add_to_block(seed0);

        while ((int)out.blocks.back().size() < target) {
            int chosen = -1;

            if (cand.empty()) {
                // 近傍が尽きたら未割当から新しいseedで継続（G'からの取り出しに相当）:contentReference[oaicite:7]{index=7}
                int seed = pick_next_seed(out.block_of);
                if (seed < 0) break;
                chosen = seed;
            } else if (policy == BlockPolicy::FIFO) {
                // Policy 1: FIFO:contentReference[oaicite:8]{index=8}
                // candは集合だが、FIFOっぽく先頭を使う
                // 既に割当済みのものはスキップ
                while (!cand.empty() && out.block_of[cand.front()] != -1) {
                    in_cand[cand.front()] = 0;
                    cand.erase(cand.begin());
                }
                if (!cand.empty()) {
                    chosen = cand.front();
                    in_cand[chosen] = 0;
                    cand.erase(cand.begin());
                }
            } else {
                // Policy 2/3 はスコア最大の候補を取る
                // 近傍数または重み合計を高速に得るため、adj(v)をナメつつ in_block[] を参照
                double best = -1.0;
                int best_idx = -1;

                for (int idx = 0; idx < (int)cand.size(); ++idx) {
                    int v = cand[idx];
                    if (out.block_of[v] != -1) continue;

                    double score = 0.0;
                    auto nbrs = boost::adjacent_vertices(v, G);
                    for (AdjacencyIter it = nbrs.first; it != nbrs.second; ++it) {
                        int u = (int)*it;
                        if (!in_block[u]) continue;

                        if (policy == BlockPolicy::MaxEdges) {
                            // 本数を数える（整数だがdoubleに載せる）:contentReference[oaicite:9]{index=9}
                            score += 1.0;
                        } else {
                            // Weighted: Σ_{u∈Vk} |a_{vu}| を近似（無向1本をそのまま利用）:contentReference[oaicite:10]{index=10}
                            auto ep = boost::edge(v, u, G);
                            if (ep.second) {
                                double w = std::abs(get(wmap, ep.first));
                                score += w;
                            } else {
                                // (グラフが無向であれば通常到達する)
                            }
                        }
                    }
                    if (score > best || (score == best && v < best_idx)) {
                        best = score;
                        best_idx = v;
                    }
                }

                if (best_idx != -1) {
                    chosen = best_idx;
                    // cand から削除
                    auto it = std::find(cand.begin(), cand.end(), chosen);
                    if (it != cand.end()) cand.erase(it);
                    in_cand[chosen] = 0;
                } else {
                    // すべて割当済み/スコア不可→seed fallback
                    int seed = pick_next_seed(out.block_of);
                    if (seed < 0) break;
                    chosen = seed;
                }
            }

            if (chosen == -1) break;
            add_to_block(chosen);
        }
    }

    out.nb = (int)out.blocks.size();
    return out;
}
