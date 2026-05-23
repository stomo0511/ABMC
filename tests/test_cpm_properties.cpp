#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "inc/cpm_properties.hxx"

namespace {

void require_close(double actual, double expected, double eps, const char* label)
{
  if (std::abs(actual - expected) > eps) {
    std::cerr << label << ": expected " << expected << ", got " << actual << "\n";
    throw std::runtime_error(label);
  }
}

struct TinyGraph {
  using key_type = unsigned;
  std::vector<std::vector<std::pair<unsigned, double>>> adj;

  std::size_t span() const { return adj.size(); }
  bool hasVertex(unsigned u) const { return u < adj.size(); }

  template <class F>
  void forEachEdge(unsigned u, F f) const {
    for (const auto& e : adj[u]) f(e.first, e.second);
  }
};

void add_undirected(TinyGraph& g, unsigned u, unsigned v, double w = 1.0)
{
  if (std::max(u, v) >= g.adj.size()) g.adj.resize(std::max(u, v) + 1);
  g.adj[u].push_back({v, w});
  g.adj[v].push_back({u, w});
}

void test_cpm_community()
{
  require_close(cpmCommunity(5.0, 3.0, 0.2), 5.0 - 0.2 * 9.0, 1e-12,
                "cpmCommunity");
}

void test_delta_cpm_matches_quality_difference()
{
  TinyGraph g;
  g.adj.resize(4);
  add_undirected(g, 0, 1);
  add_undirected(g, 1, 2);
  add_undirected(g, 2, 3);

  std::vector<unsigned> before = {0, 0, 1, 1};
  std::vector<unsigned> after  = {0, 1, 1, 1};
  const double gamma = 0.1;

  const double q_before = cpmByOmp(g, [&](unsigned u) { return before[u]; }, gamma);
  const double q_after = cpmByOmp(g, [&](unsigned u) { return after[u]; }, gamma);

  // Move vertex 1 from community 0 to community 1.
  // u_to_c = edge(1,2) = 1, u_to_d = edge(1,0) = 1,
  // usize = 1, csize = 2, dsize = 2.
  const double delta = deltaCPM(1.0, 1.0, 1.0, 2.0, 2.0, gamma);
  require_close(delta, q_after - q_before, 1e-12, "deltaCPM");
}

}  // namespace

int main()
{
  test_cpm_community();
  test_delta_cpm_matches_quality_difference();
  std::cout << "test_cpm_properties: ok\n";
  return 0;
}
