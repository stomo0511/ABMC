#include <numeric>
#include "Types.hpp"
#include "Coloring.hpp"

// グラフ G に対して貪欲法で彩色を行う
// 戻り値: 使用した色の数
// 色情報を格納した配列 color も返す
// 彩色方針は「次数降順＋ID昇順タイブレーク」
int Greedy_Coloring(const Graph& G, std::vector<int>& color)
{
    const std::size_t N = G.n;
    color.assign(N, -1); // -1: uncolored

    struct VertexDegree {
        int v;
        std::size_t deg;
    };

    std::vector<VertexDegree> order;
    order.reserve(N);
    for (int v = 0; v < static_cast<int>(N); ++v) {
        order.push_back({v, G.adj[v].size()});
    }

    // 次数降順、同次数ならID昇順
    std::sort(order.begin(), order.end(),
              [](const VertexDegree& a, const VertexDegree& b) {
                  if (a.deg != b.deg) return a.deg > b.deg; // 降順
                  return a.v < b.v;                         // タイブレーク：ID昇順
              });

    // 近傍色のマーキング（タイムスタンプ法）
    std::vector<int> mark(N + 1, -1); // 色番号は最大でも N-1
    int stamp = 0;
    int max_color = -1;

    for (const auto& vd : order) {
        int u = vd.v;
        ++stamp;

        for (const auto& e : G.adj[u]) {
            int c = color[e.to];
            if (c >= 0) mark[c] = stamp;
        }

        // 最小許容色を割り当て
        int c = 0;
        while (c <= max_color && mark[c] == stamp) ++c;
        if (c == max_color + 1) ++max_color;
        color[u] = c;
    }

    return max_color + 1; // 使用色数
}

// 頻度順に色ラベルを付け替える
// 同数なら旧ラベルの昇順を優先（安定なタイブレーク）
void RelabelColorsByClassSize(std::vector<int>& color) {
    if (color.empty()) return;

    int nc = 1 + *std::max_element(color.begin(), color.end());
    std::vector<int> cnt(nc, 0);
    for (int c : color) if (c >= 0) ++cnt[c];

    std::vector<int> order(nc);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b){
                  if (cnt[a] != cnt[b]) return cnt[a] > cnt[b]; // 大きい順
                  return a < b;                                  // タイブレーク
              });

    // old -> new の写像を作る
    std::vector<int> new_id(nc, -1);
    for (int i = 0; i < nc; ++i) new_id[ order[i] ] = i;

    // 配列を書き換え
    for (int& c : color) if (c >= 0) c = new_id[c];
}
