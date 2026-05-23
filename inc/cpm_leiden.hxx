#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

#include "cpm_properties.hxx"

// Options for CPM-maximizing Leiden.
struct CPMLeidenOptions {
  double gamma = 1.0;
  double tolerance = 1.0e-2;
  double aggregationTolerance = 0.8;
  double toleranceDrop = 10.0;
  int maxIterations = 20;
  int maxPasses = 10;
  int repeat = 1;
  bool verbose = false;
};

// Result of CPM-maximizing Leiden.
template <class K = std::uint32_t>
struct CPMLeidenResult {
  std::vector<K> membership;
  std::size_t communities = 0;
  double quality = 0.0;
  int iterations = 0;
  int passes = 0;
  double elapsedSeconds = 0.0;
};

template <class K = std::uint32_t, class W = double>
class CPMAggregatedGraph {
public:
  using key_type = K;
  using vertex_value_type = double;
  using edge_value_type = W;

  std::vector<std::size_t> offsets;
  std::vector<K> degrees;
  std::vector<K> edgeKeys;
  std::vector<W> edgeValues;

  CPMAggregatedGraph() = default;

  CPMAggregatedGraph(
    std::size_t n,
    std::vector<std::unordered_map<K, W>>& adjacency)
  {
    offsets.resize(n + 1, 0);
    degrees.resize(n, K());

    std::size_t m = 0;
    for (std::size_t u = 0; u < n; ++u) {
      offsets[u] = m;
      degrees[u] = static_cast<K>(adjacency[u].size());
      m += adjacency[u].size();
    }
    offsets[n] = m;

    edgeKeys.reserve(m);
    edgeValues.reserve(m);
    for (std::size_t u = 0; u < n; ++u) {
      std::vector<std::pair<K, W>> edges(adjacency[u].begin(), adjacency[u].end());
      std::sort(edges.begin(), edges.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
      });
      for (const auto& edge : edges) {
        edgeKeys.push_back(edge.first);
        edgeValues.push_back(edge.second);
      }
    }
  }

  inline std::size_t span() const noexcept { return degrees.size(); }
  inline std::size_t order() const noexcept { return degrees.size(); }
  inline std::size_t size() const noexcept { return edgeKeys.size(); }
  inline bool empty() const noexcept { return degrees.empty(); }
  inline bool directed() const noexcept { return true; }
  inline bool hasVertex(K u) const noexcept {
    return static_cast<std::size_t>(u) < degrees.size();
  }
  inline std::size_t degree(K u) const noexcept {
    return hasVertex(u) ? degrees[u] : 0;
  }

  template <class FP>
  inline void forEachVertexKey(FP fp) const noexcept {
    for (K u = 0; static_cast<std::size_t>(u) < degrees.size(); ++u) {
      fp(u);
    }
  }

  template <class FP>
  inline void forEachEdge(K u, FP fp) const noexcept {
    if (!hasVertex(u)) return;
    const std::size_t begin = offsets[u];
    const std::size_t end = begin + degrees[u];
    for (std::size_t i = begin; i < end; ++i) {
      fp(edgeKeys[i], edgeValues[i]);
    }
  }

  template <class FP>
  inline void forEachEdgeKey(K u, FP fp) const noexcept {
    if (!hasVertex(u)) return;
    const std::size_t begin = offsets[u];
    const std::size_t end = begin + degrees[u];
    for (std::size_t i = begin; i < end; ++i) {
      fp(edgeKeys[i]);
    }
  }
};

template <class K = std::uint32_t, class W = double>
struct CPMAggregationResult {
  CPMAggregatedGraph<K, W> graph;
  std::vector<double> vertexSize;
  std::vector<K> denseMembership;
  std::size_t communities = 0;
  bool aggregated = false;

  explicit operator bool() const noexcept { return aggregated; }
};

template <class G, class K>
inline std::size_t cpmRenumberCommunities(
  const G& graph,
  std::vector<K>& membership,
  std::vector<K>& oldToNew)
{
  oldToNew.assign(membership.size(), K());
  std::vector<char> seen(membership.size(), 0);
  std::size_t count = 0;

  graph.forEachVertexKey([&](auto raw_u) {
    const K u = static_cast<K>(raw_u);
    const K c = membership[u];
    const std::size_t ci = static_cast<std::size_t>(c);
    if (ci >= seen.size()) return;
    if (!seen[ci]) {
      seen[ci] = 1;
      oldToNew[ci] = static_cast<K>(count++);
    }
  });

  graph.forEachVertexKey([&](auto raw_u) {
    const K u = static_cast<K>(raw_u);
    const K c = membership[u];
    membership[u] = oldToNew[static_cast<std::size_t>(c)];
  });

  return count;
}

template <class G, class K>
inline void initializeCPMSizes(
  const G& graph,
  std::vector<K>& membership,
  std::vector<double>& vsize,
  std::vector<double>& csize)
{
  const std::size_t S = graph.span();

  membership.assign(S, K());

  // vsize is initialized to 1.0 only for the original graph. After aggregation,
  // callers should pass the super-vertex sizes in vsize; this function preserves
  // those values and only reinitializes csize for the current pass.
  if (vsize.size() != S) {
    vsize.assign(S, 0.0);
    graph.forEachVertexKey([&](auto u) {
      vsize[u] = 1.0;
    });
  }

  csize.assign(S, 0.0);
  graph.forEachVertexKey([&](auto u) {
    membership[u] = static_cast<K>(u);
    csize[u] = vsize[u];
  });
}

template <class G>
inline bool validateCPMSizes(
  const G& graph,
  const std::vector<double>& vsize,
  const std::vector<double>& csize,
  double tolerance = 1.0e-9)
{
  double vertexTotal = 0.0;
  graph.forEachVertexKey([&](auto u) {
    if (static_cast<std::size_t>(u) < vsize.size()) {
      vertexTotal += vsize[u];
    }
  });

  double communityTotal = 0.0;
  for (double n : csize) {
    communityTotal += n;
  }

  return std::abs(vertexTotal - communityTotal) <= tolerance;
}

template <class G, class K>
inline std::size_t cpmCountCommunities(const G& graph, const std::vector<K>& membership)
{
  std::vector<char> seen(membership.size(), 0);
  std::size_t count = 0;

  graph.forEachVertexKey([&](auto u) {
    const K c = membership[u];
    if (static_cast<std::size_t>(c) >= seen.size()) return;
    if (seen[c]) return;
    seen[c] = 1;
    ++count;
  });

  return count;
}

// CPM local-moving phase.
//
// This first implementation is sequential by design. It keeps membership and
// communitySize updates race-free while the CPM objective path is brought up.
// A parallel version needs per-thread scan buffers plus atomic or otherwise
// coordinated updates for membership and communitySize.
template <class G, class K>
inline int cpmLeidenMove(
  const G& graph,
  std::vector<K>& membership,
  const std::vector<double>& vertexSize,
  std::vector<double>& communitySize,
  const CPMLeidenOptions& options)
{
  const std::size_t S = graph.span();
  if (membership.size() < S || vertexSize.size() < S || communitySize.size() < S) {
    return 0;
  }

  std::vector<K> candidates;
  std::vector<double> edgeToCommunity(S, 0.0);
  std::vector<int> mark(S, 0);
  int stamp = 0;
  int iterations = 0;

  for (; iterations < options.maxIterations; ++iterations) {
    double totalGain = 0.0;
    int moves = 0;

    graph.forEachVertexKey([&](auto raw_u) {
      const K u = static_cast<K>(raw_u);
      const std::size_t ui = static_cast<std::size_t>(u);
      if (ui >= S) return;

      const K oldCommunity = membership[u];
      const std::size_t oldIndex = static_cast<std::size_t>(oldCommunity);
      if (oldIndex >= S) return;

      const double usize = vertexSize[u];
      if (usize <= 0.0) return;

      candidates.clear();
      ++stamp;
      if (stamp == 0) {
        std::fill(mark.begin(), mark.end(), 0);
        stamp = 1;
      }

      graph.forEachEdge(u, [&](auto raw_v, auto raw_w) {
        const K v = static_cast<K>(raw_v);
        if (v == u) return;

        const std::size_t vi = static_cast<std::size_t>(v);
        if (vi >= membership.size()) return;

        const K c = membership[v];
        const std::size_t ci = static_cast<std::size_t>(c);
        if (ci >= S) return;

        if (mark[ci] != stamp) {
          mark[ci] = stamp;
          edgeToCommunity[ci] = 0.0;
          candidates.push_back(c);
        }
        edgeToCommunity[ci] += static_cast<double>(raw_w);
      });

      const double uToOld = mark[oldIndex] == stamp
        ? edgeToCommunity[oldIndex]
        : 0.0;
      const double oldSize = communitySize[oldCommunity];

      K bestCommunity = oldCommunity;
      double bestGain = 0.0;

      for (K c : candidates) {
        if (c == oldCommunity) continue;

        const std::size_t ci = static_cast<std::size_t>(c);
        const double gain = deltaCPM(
          edgeToCommunity[ci],
          uToOld,
          usize,
          communitySize[c],
          oldSize,
          options.gamma);

        if (gain > bestGain) {
          bestGain = gain;
          bestCommunity = c;
        }
      }

      if (bestCommunity == oldCommunity || bestGain <= 0.0) return;

      communitySize[oldCommunity] -= usize;
      communitySize[bestCommunity] += usize;
      membership[u] = bestCommunity;
      totalGain += bestGain;
      ++moves;
    });

    // CPM is not normalized like modularity. Its absolute gain scales with
    // edge weights, gamma, and community sizes, so options.tolerance should be
    // interpreted in the raw CPM objective scale.
    if (moves == 0 || totalGain <= options.tolerance) {
      break;
    }
  }

  return iterations;
}

// CPM refinement phase.
//
// The current membership after local-moving is captured as the community bound.
// Refinement then restarts from singleton communities and only allows moves
// between vertices whose bound is identical. This mirrors the Leiden refinement
// idea while keeping the CPM objective state independent from modularity.
template <class G, class K>
inline int cpmLeidenRefine(
  const G& graph,
  std::vector<K>& membership,
  const std::vector<double>& vertexSize,
  std::vector<double>& communitySize,
  const CPMLeidenOptions& options)
{
  const std::size_t S = graph.span();
  if (membership.size() < S || vertexSize.size() < S || communitySize.size() < S) {
    return 0;
  }

  const std::vector<K> communityBound = membership;
  communitySize.assign(S, 0.0);
  graph.forEachVertexKey([&](auto raw_u) {
    const K u = static_cast<K>(raw_u);
    membership[u] = u;
    communitySize[u] = vertexSize[u];
  });

  std::vector<K> candidates;
  std::vector<double> edgeToCommunity(S, 0.0);
  std::vector<int> mark(S, 0);
  int stamp = 0;
  int iterations = 0;

  const double sizeTolerance = std::max(1.0e-12, std::abs(options.tolerance) * 1.0e-9);

  for (; iterations < options.maxIterations; ++iterations) {
    double totalGain = 0.0;
    int moves = 0;

    graph.forEachVertexKey([&](auto raw_u) {
      const K u = static_cast<K>(raw_u);
      const std::size_t ui = static_cast<std::size_t>(u);
      if (ui >= S) return;

      const K oldCommunity = membership[u];
      const std::size_t oldIndex = static_cast<std::size_t>(oldCommunity);
      if (oldIndex >= S) return;

      const double usize = vertexSize[u];
      if (usize <= 0.0) return;

      // In Leiden refinement, only singleton or isolated communities are moved.
      // With CPM state this is a size test rather than a total-degree test.
      if (std::abs(communitySize[oldCommunity] - usize) > sizeTolerance) return;

      candidates.clear();
      ++stamp;
      if (stamp == 0) {
        std::fill(mark.begin(), mark.end(), 0);
        stamp = 1;
      }

      graph.forEachEdge(u, [&](auto raw_v, auto raw_w) {
        const K v = static_cast<K>(raw_v);
        if (v == u) return;

        const std::size_t vi = static_cast<std::size_t>(v);
        if (vi >= membership.size() || vi >= communityBound.size()) return;
        if (communityBound[v] != communityBound[u]) return;

        const K c = membership[v];
        const std::size_t ci = static_cast<std::size_t>(c);
        if (ci >= S) return;

        if (mark[ci] != stamp) {
          mark[ci] = stamp;
          edgeToCommunity[ci] = 0.0;
          candidates.push_back(c);
        }
        edgeToCommunity[ci] += static_cast<double>(raw_w);
      });

      const double uToOld = mark[oldIndex] == stamp
        ? edgeToCommunity[oldIndex]
        : 0.0;
      const double oldSize = communitySize[oldCommunity];

      K bestCommunity = oldCommunity;
      double bestGain = 0.0;

      for (K c : candidates) {
        if (c == oldCommunity) continue;

        const std::size_t ci = static_cast<std::size_t>(c);
        const double gain = deltaCPM(
          edgeToCommunity[ci],
          uToOld,
          usize,
          communitySize[c],
          oldSize,
          options.gamma);

        if (gain > bestGain) {
          bestGain = gain;
          bestCommunity = c;
        }
      }

      if (bestCommunity == oldCommunity || bestGain <= 0.0) return;

      communitySize[oldCommunity] -= usize;
      communitySize[bestCommunity] += usize;
      membership[u] = bestCommunity;
      totalGain += bestGain;
      ++moves;
    });

    // CPM gain is in raw objective units, not modularity's normalized scale.
    if (moves == 0 || totalGain <= options.tolerance) {
      break;
    }
  }

  return iterations;
}

// CPM aggregation phase.
//
// Communities are renumbered to dense ids and collapsed into super-vertices.
// new_vsize[c] is the sum of vsize[u] for vertices assigned to community c.
// Self-loops are intentionally kept, matching the existing Leiden aggregation
// style that scans with SELF=true. For undirected graphs stored as two directed
// arcs, a community-internal undirected edge becomes a self-loop weight counted
// twice; cpmByOmp(..., undirected_edges_stored_twice=true) divides it by 2.
template <class G, class K>
inline CPMAggregationResult<K> cpmLeidenAggregate(
  const G& graph,
  std::vector<K>& membership,
  std::vector<double>& vertexSize,
  std::vector<double>& communitySize,
  const CPMLeidenOptions& options)
{
  (void)options;

  CPMAggregationResult<K> result;
  const std::size_t S = graph.span();
  if (membership.size() < S || vertexSize.size() < S) {
    return result;
  }

  std::vector<K> oldToNew;
  const std::size_t communities = cpmRenumberCommunities(graph, membership, oldToNew);
  result.communities = communities;
  result.denseMembership = membership;

  if (communities == 0) {
    communitySize.clear();
    vertexSize.clear();
    return result;
  }

  result.vertexSize.assign(communities, 0.0);
  graph.forEachVertexKey([&](auto raw_u) {
    const K u = static_cast<K>(raw_u);
    const K c = membership[u];
    result.vertexSize[static_cast<std::size_t>(c)] += vertexSize[u];
  });

  std::vector<std::unordered_map<K, double>> adjacency(communities);
  graph.forEachVertexKey([&](auto raw_u) {
    const K u = static_cast<K>(raw_u);
    const K cu = membership[u];
    const std::size_t cui = static_cast<std::size_t>(cu);
    if (cui >= communities) return;

    graph.forEachEdge(u, [&](auto raw_v, auto raw_w) {
      const K v = static_cast<K>(raw_v);
      const std::size_t vi = static_cast<std::size_t>(v);
      if (vi >= membership.size()) return;

      const K cv = membership[v];
      const std::size_t cvi = static_cast<std::size_t>(cv);
      if (cvi >= communities) return;

      adjacency[cui][cv] += static_cast<double>(raw_w);
    });
  });

  result.graph = CPMAggregatedGraph<K, double>(communities, adjacency);
  result.aggregated = communities < graph.order();

  // Keep the CPM size state consistent for callers that inspect this phase.
  // The top-level skeleton does not yet swap to result.graph, but later passes
  // must use result.vertexSize rather than resetting super-vertex size to 1.0.
  communitySize = result.vertexSize;

  return result;
}

template <class G>
inline auto cpmLeidenRunLevel(
  const G& graph,
  const std::vector<double>& inputVertexSize,
  const CPMLeidenOptions& options,
  int remainingPasses,
  int depth = 0)
  -> CPMLeidenResult<typename G::key_type>
{
  using K = typename G::key_type;

  CPMLeidenResult<K> result;
  std::vector<double> vertexSize = inputVertexSize;
  std::vector<double> communitySize(graph.span(), 0.0);

  initializeCPMSizes(graph, result.membership, vertexSize, communitySize);

  const int moved = cpmLeidenMove(
    graph,
    result.membership,
    vertexSize,
    communitySize,
    options);

  // cpmLeidenRefine stores the local-moving partition as the community bound,
  // restarts from singleton communities, and restricts refinement moves to that
  // bound. This keeps the top-level flow aligned with Leiden's dendrogram step.
  const int refined = cpmLeidenRefine(
    graph,
    result.membership,
    vertexSize,
    communitySize,
    options);

  auto aggregation = cpmLeidenAggregate(
    graph,
    result.membership,
    vertexSize,
    communitySize,
    options);

  result.iterations = moved + refined;
  result.passes = 1;

  if (options.verbose) {
    std::cerr << "CPM Leiden level " << depth
              << ": moved=" << moved
              << " refined=" << refined
              << " communities=" << aggregation.communities
              << " aggregated=" << (aggregation ? 1 : 0)
              << "\n";
  }

  if (!aggregation || remainingPasses <= 1 || aggregation.communities == 0) {
    std::vector<K> oldToNew;
    result.communities = cpmRenumberCommunities(graph, result.membership, oldToNew);
    return result;
  }

  auto upper = cpmLeidenRunLevel(
    aggregation.graph,
    aggregation.vertexSize,
    options,
    remainingPasses - 1,
    depth + 1);

  graph.forEachVertexKey([&](auto raw_u) {
    const K u = static_cast<K>(raw_u);
    const K superVertex = aggregation.denseMembership[u];
    const std::size_t si = static_cast<std::size_t>(superVertex);
    if (si < upper.membership.size()) {
      result.membership[u] = upper.membership[si];
    }
  });

  std::vector<K> oldToNew;
  result.communities = cpmRenumberCommunities(graph, result.membership, oldToNew);
  result.iterations += upper.iterations;
  result.passes += upper.passes;

  return result;
}

// Top-level CPM Leiden entry point.
//
// The implementation is sequential and recursive over aggregation levels. Each
// level runs local-moving, refinement, dense renumbering, aggregation, and then
// composes the upper-level membership back onto the original vertices for that
// level. Existing modularity Leiden code is not modified or called.
template <class G>
inline auto cpmLeiden(const G& graph, const CPMLeidenOptions& options = {})
  -> CPMLeidenResult<typename G::key_type>
{
  using K = typename G::key_type;
  using Clock = std::chrono::steady_clock;

  const auto begin = Clock::now();
  std::vector<double> vertexSize(graph.span(), 0.0);
  graph.forEachVertexKey([&](auto u) {
    vertexSize[u] = 1.0;
  });

  CPMLeidenResult<K> result = cpmLeidenRunLevel(
    graph,
    vertexSize,
    options,
    options.maxPasses);

  result.quality = cpmByOmp(
    graph,
    [&](K u) { return result.membership[u]; },
    options.gamma,
    true);

  const auto end = Clock::now();
  result.elapsedSeconds = std::chrono::duration<double>(end - begin).count();

  return result;
}
