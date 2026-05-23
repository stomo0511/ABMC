#pragma once

#include <cstddef>
#include <vector>

// Constant Potts Model objective:
//   Q_CPM = sum_c (e_c - gamma * n_c^2)
//
// e_c is the internal edge weight of community c.
// n_c is the node-size of community c.
// gamma is the CPM resolution parameter.

inline double cpmCommunity(double cin, double csize, double gamma)
{
  return cin - gamma * csize * csize;
}

// Delta Q_CPM for moving vertex u from its current community d to candidate c.
//
// u_to_c: total edge weight from u to community c
// u_to_d: total edge weight from u to community d
// usize:  node-size of u; normally 1 for original vertices, and the super-vertex
//         size after aggregation if CPM is run on an aggregated graph.
// csize:  current node-size of candidate community c before moving u
// dsize:  current node-size of source community d before moving u
//
// This function assumes u_to_c and u_to_d are counted as undirected internal
// edge-weight changes, i.e. moving u adds u_to_c and removes u_to_d.
// In a graph implementation that stores each undirected edge in both directions,
// scanning only u's outgoing edges still sees each neighbor once, so no factor 2
// is needed for the edge term. If the caller instead supplies directed internal
// weights counted in both directions, pass edge_count_factor=2.0.
//
// The penalty delta for -gamma*n_c^2 is:
//   -gamma*((csize+usize)^2 + (dsize-usize)^2 - csize^2 - dsize^2)
// =  2*gamma*usize*(dsize - csize - usize)
// Therefore the penalty part does require the factor 2 shown above for the
// n^2 CPM objective.
inline double deltaCPM(
  double u_to_c,
  double u_to_d,
  double usize,
  double csize,
  double dsize,
  double gamma,
  double edge_count_factor = 1.0)
{
  const double edge_delta = edge_count_factor * (u_to_c - u_to_d);
  const double penalty_delta = 2.0 * gamma * usize * (dsize - csize - usize);
  return edge_delta + penalty_delta;
}

// Compute Q_CPM for a whole partition in parallel with OpenMP pragmas.
//
// G must provide:
//   key_type
//   span()
//   hasVertex(u)
//   forEachEdge(u, fn(v, w))
//
// FC must map a vertex id u to its community id.
//
// The default undirected_edges_stored_twice=true matches the common GVE setup
// where an undirected graph is represented by two directed arcs. In that case
// the internal edge accumulator sees each undirected internal edge twice, so it
// is divided by 2 before cpmCommunity() is applied.
template <class G, class FC>
inline double cpmByOmp(
  const G& x,
  FC community_of,
  double gamma,
  bool undirected_edges_stored_twice = true)
{
  using K = typename G::key_type;

  const std::size_t S = x.span();
  std::vector<double> vertex_internal(S, 0.0);
  std::vector<double> community_internal(S, 0.0);
  std::vector<double> community_size(S, 0.0);

  #pragma omp parallel for schedule(dynamic, 2048)
  for (K u = 0; u < static_cast<K>(S); ++u) {
    if (!x.hasVertex(u)) continue;

    const K cu = community_of(u);
    x.forEachEdge(u, [&](auto v, auto w) {
      const K cv = community_of(v);
      if (cu == cv) vertex_internal[u] += static_cast<double>(w);
    });
  }

  #pragma omp parallel for schedule(static, 2048)
  for (K u = 0; u < static_cast<K>(S); ++u) {
    if (!x.hasVertex(u)) continue;

    const K cu = community_of(u);

    #pragma omp atomic
    community_internal[cu] += vertex_internal[u];

    #pragma omp atomic
    community_size[cu] += 1.0;
  }

  const double internal_scale = undirected_edges_stored_twice ? 0.5 : 1.0;
  double q = 0.0;

  #pragma omp parallel for schedule(static) reduction(+:q)
  for (K c = 0; c < static_cast<K>(S); ++c) {
    const double n = community_size[c];
    if (n <= 0.0) continue;

    const double cin = internal_scale * community_internal[c];
    q += cpmCommunity(cin, n, gamma);
  }

  return q;
}
