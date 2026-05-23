#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "inc/cpm_leiden.hxx"

namespace {

struct TinyGraph {
  using key_type = unsigned;
  std::vector<std::vector<std::pair<unsigned, double>>> adj;

  std::size_t span() const { return adj.size(); }
  std::size_t order() const { return adj.size(); }
  bool hasVertex(unsigned u) const { return u < adj.size(); }

  template <class F>
  void forEachVertexKey(F f) const {
    for (unsigned u = 0; u < adj.size(); ++u) f(u);
  }

  template <class F>
  void forEachEdge(unsigned u, F f) const {
    for (const auto& e : adj[u]) f(e.first, e.second);
  }

  template <class F>
  void forEachEdgeKey(unsigned u, F f) const {
    for (const auto& e : adj[u]) f(e.first);
  }
};

void add_undirected(TinyGraph& g, unsigned u, unsigned v, double w = 1.0)
{
  const unsigned n = std::max(u, v) + 1;
  if (n > g.adj.size()) g.adj.resize(n);
  g.adj[u].push_back({v, w});
  g.adj[v].push_back({u, w});
}

TinyGraph complete_graph(unsigned n)
{
  TinyGraph g;
  g.adj.resize(n);
  for (unsigned u = 0; u < n; ++u) {
    for (unsigned v = u + 1; v < n; ++v) {
      add_undirected(g, u, v);
    }
  }
  return g;
}

TinyGraph two_triangles_with_bridge()
{
  TinyGraph g;
  g.adj.resize(6);
  add_undirected(g, 0, 1);
  add_undirected(g, 1, 2);
  add_undirected(g, 2, 0);
  add_undirected(g, 3, 4);
  add_undirected(g, 4, 5);
  add_undirected(g, 5, 3);
  add_undirected(g, 2, 3);
  return g;
}

TinyGraph ring_of_cliques(unsigned cliques, unsigned clique_size)
{
  TinyGraph g;
  g.adj.resize(cliques * clique_size);

  for (unsigned c = 0; c < cliques; ++c) {
    const unsigned base = c * clique_size;
    for (unsigned i = 0; i < clique_size; ++i) {
      for (unsigned j = i + 1; j < clique_size; ++j) {
        add_undirected(g, base + i, base + j);
      }
    }
  }

  for (unsigned c = 0; c < cliques; ++c) {
    const unsigned next = (c + 1) % cliques;
    add_undirected(g, c * clique_size, next * clique_size);
  }

  return g;
}

template <class K>
bool same_community(const std::vector<K>& membership, unsigned a, unsigned b)
{
  return membership[a] == membership[b];
}

void require(bool condition, const char* label)
{
  if (!condition) throw std::runtime_error(label);
}

CPMLeidenOptions options(double gamma)
{
  CPMLeidenOptions o;
  o.gamma = gamma;
  o.tolerance = 1e-9;
  o.maxIterations = 50;
  o.maxPasses = 10;
  return o;
}

void test_complete_graph_small_gamma()
{
  const TinyGraph g = complete_graph(5);
  const auto r = cpmLeiden(g, options(0.01));
  require(r.communities == 1, "complete graph with small gamma should be one community");
}

void test_two_triangles_gamma_split_and_merge()
{
  const TinyGraph g = two_triangles_with_bridge();

  const auto merged = cpmLeiden(g, options(0.01));
  require(merged.communities == 1, "two triangles should merge for small gamma");

  const auto split = cpmLeiden(g, options(0.20));
  require(split.communities >= 2, "two triangles should split for larger gamma");
  require(same_community(split.membership, 0, 1) &&
          same_community(split.membership, 1, 2),
          "left triangle should stay together");
  require(same_community(split.membership, 3, 4) &&
          same_community(split.membership, 4, 5),
          "right triangle should stay together");
  require(!same_community(split.membership, 0, 3),
          "two triangles should be separated");
}

void test_ring_of_cliques()
{
  const unsigned clique_count = 4;
  const unsigned clique_size = 4;
  const TinyGraph g = ring_of_cliques(clique_count, clique_size);
  const auto r = cpmLeiden(g, options(0.10));

  require(r.communities == clique_count, "ring of cliques should split by clique");
  for (unsigned c = 0; c < clique_count; ++c) {
    const unsigned base = c * clique_size;
    for (unsigned i = 1; i < clique_size; ++i) {
      require(same_community(r.membership, base, base + i),
              "vertices in a clique should share a community");
    }
  }
}

}  // namespace

int main()
{
  test_complete_graph_small_gamma();
  test_two_triangles_gamma_split_and_merge();
  test_ring_of_cliques();
  std::cout << "test_cpm_leiden: ok\n";
  return 0;
}
