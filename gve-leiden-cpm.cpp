#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "common/Types.hpp"
#include "common/Timer.hpp"
#include "gve-leiden-inc/main.hxx"
#include "inc/cpm_properties.hxx"
#include "inc/cpm_leiden.hxx"

static std::string file_stem_local(const std::string& path)
{
  const std::size_t slash = path.find_last_of("/\\");
  const std::string name = slash == std::string::npos ? path : path.substr(slash + 1);
  const std::size_t dot = name.find_last_of('.');
  return dot == std::string::npos ? name : name.substr(0, dot);
}

static std::vector<int> relabel_dense(
  const std::vector<int>& comm,
  std::vector<int>* sizes_out = nullptr)
{
  const std::size_t n = comm.size();

  struct Stat {
    int label;
    int count;
  };

  std::unordered_map<int, int> counts;
  counts.reserve(n * 2);
  for (const int c : comm) {
    if (c >= 0) ++counts[c];
  }

  if (counts.empty()) {
    if (sizes_out) sizes_out->clear();
    return std::vector<int>(n, -1);
  }

  std::vector<Stat> order;
  order.reserve(counts.size());
  for (const auto& entry : counts) {
    order.push_back({entry.first, entry.second});
  }
  std::sort(order.begin(), order.end(), [](const Stat& a, const Stat& b) {
    if (a.count != b.count) return a.count > b.count;
    return a.label < b.label;
  });

  std::unordered_map<int, int> old_to_new;
  old_to_new.reserve(order.size());
  if (sizes_out) sizes_out->assign(order.size(), 0);

  for (int new_id = 0; new_id < static_cast<int>(order.size()); ++new_id) {
    old_to_new[order[new_id].label] = new_id;
    if (sizes_out) (*sizes_out)[new_id] = order[new_id].count;
  }

  std::vector<int> dense(n, -1);
  for (std::size_t i = 0; i < n; ++i) {
    const int c = comm[i];
    if (c >= 0) dense[i] = old_to_new[c];
  }
  return dense;
}

static void write_partition_1based(
  const std::vector<int>& block_of,
  const std::string& out_path)
{
  std::ofstream out(out_path);
  if (!out) {
    throw std::runtime_error("cannot open output: " + out_path);
  }

  const int nb = block_of.empty()
    ? 0
    : (*std::max_element(block_of.begin(), block_of.end()) + 1);

  out << nb << "\n";
  for (std::size_t i = 0; i < block_of.size(); ++i) {
    const int block = block_of[i];
    if (block < 0) {
      throw std::runtime_error("unassigned node in partition");
    }
    out << (i + 1) << " " << (block + 1) << "\n";
  }
}

static void write_block_color_1based(
  const std::vector<int>& block_color,
  int color_count,
  const std::string& out_path)
{
  std::ofstream out(out_path);
  if (!out) {
    throw std::runtime_error("cannot open output: " + out_path);
  }

  out << color_count << "\n";
  for (std::size_t block = 0; block < block_color.size(); ++block) {
    const int color = block_color[block];
    if (color < 0 || color >= color_count) {
      throw std::runtime_error("block color out of range");
    }
    out << (block + 1) << " " << (color + 1) << "\n";
  }
}

static std::string default_bcol_path(const std::string& blk_path)
{
  const std::size_t slash = blk_path.find_last_of("/\\");
  const std::size_t dot = blk_path.find_last_of('.');
  if (dot != std::string::npos && (slash == std::string::npos || dot > slash)) {
    return blk_path.substr(0, dot) + ".bcol";
  }
  return blk_path + ".bcol";
}

template <class GveGraph>
static Graph build_block_graph_from_gve(
  const GveGraph& graph,
  const std::vector<int>& block_of,
  int block_count)
{
  struct EdgeHash {
    std::size_t operator()(std::uint64_t key) const {
      return static_cast<std::size_t>(key ^ (key >> 32));
    }
  };

  std::unordered_set<std::uint64_t, EdgeHash> block_edges;
  block_edges.reserve(graph.size());

  graph.forEachVertexKey([&](auto u) {
    if (u == 0) return;
    const std::size_t node_u = static_cast<std::size_t>(u - 1);
    if (node_u >= block_of.size()) return;

    const int bu = block_of[node_u];
    if (bu < 0) return;

    graph.forEachEdgeKey(u, [&](auto v) {
      if (v == 0) return;
      const std::size_t node_v = static_cast<std::size_t>(v - 1);
      if (node_v >= block_of.size()) return;

      const int bv = block_of[node_v];
      if (bv < 0 || bu == bv) return;

      const std::uint32_t a = static_cast<std::uint32_t>(std::min(bu, bv));
      const std::uint32_t b = static_cast<std::uint32_t>(std::max(bu, bv));
      block_edges.insert((static_cast<std::uint64_t>(a) << 32) | b);
    });
  });

  Graph block_graph(block_count);
  auto weights = get(boost::edge_weight, block_graph);
  for (const std::uint64_t key : block_edges) {
    const int bu = static_cast<int>(key >> 32);
    const int bv = static_cast<int>(key & 0xffffffffu);
    auto edge = add_edge(bu, bv, block_graph).first;
    weights[edge] = 1.0;
  }

  return block_graph;
}

static int greedy_coloring(const Graph& graph, std::vector<int>& color)
{
  using Vertex = boost::graph_traits<Graph>::vertex_descriptor;

  struct VertexDegree {
    Vertex vertex;
    std::size_t degree;
  };

  const std::size_t n = boost::num_vertices(graph);
  color.assign(n, -1);

  std::vector<VertexDegree> order;
  order.reserve(n);
  for (Vertex v = 0; v < n; ++v) {
    order.push_back({v, boost::degree(v, graph)});
  }

  std::sort(order.begin(), order.end(), [](const VertexDegree& a, const VertexDegree& b) {
    if (a.degree != b.degree) return a.degree > b.degree;
    return a.vertex < b.vertex;
  });

  std::vector<int> mark(n + 1, -1);
  int stamp = 0;
  int max_color = -1;

  for (const auto& item : order) {
    const Vertex u = item.vertex;
    ++stamp;

    auto adjacent = boost::adjacent_vertices(u, graph);
    for (auto it = adjacent.first; it != adjacent.second; ++it) {
      const int c = color[*it];
      if (c >= 0) mark[c] = stamp;
    }

    int c = 0;
    while (c <= max_color && mark[c] == stamp) ++c;
    if (c == max_color + 1) ++max_color;
    color[u] = c;
  }

  return max_color + 1;
}

static void relabel_colors_by_class_size(std::vector<int>& color)
{
  if (color.empty()) return;

  const int color_count = 1 + *std::max_element(color.begin(), color.end());
  std::vector<int> counts(color_count, 0);
  for (const int c : color) {
    if (c >= 0) ++counts[c];
  }

  std::vector<int> order(color_count);
  for (int i = 0; i < color_count; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    if (counts[a] != counts[b]) return counts[a] > counts[b];
    return a < b;
  });

  std::vector<int> new_id(color_count, -1);
  for (int i = 0; i < color_count; ++i) {
    new_id[order[i]] = i;
  }
  for (int& c : color) {
    if (c >= 0) c = new_id[c];
  }
}

static bool read_bool_arg(const char* value)
{
  return std::atoi(value) != 0;
}

static void usage(const char* program)
{
  std::fprintf(
    stderr,
    "usage: %s <matrix.mtx> [--gamma <value>] [--tolerance <value>] "
    "[--max-iterations <n>] [--max-passes <n>] [--repeat <n>] "
    "[--output <path>] [--symmetric 0|1(default=1)] [--weighted 0|1] [--verbose]\n",
    program);
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    usage(argv[0]);
    return 1;
  }

  const char* file = argv[1];
  CPMLeidenOptions options;
  std::string output_path;
  bool already_symmetric = true;
  bool weighted = false;

  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];

    auto require_value = [&](const char* name) -> const char* {
      if (i + 1 >= argc) {
        throw std::runtime_error(std::string("missing value for ") + name);
      }
      return argv[++i];
    };

    if (arg == "--gamma") {
      options.gamma = std::atof(require_value("--gamma"));
    } else if (arg == "--tolerance") {
      options.tolerance = std::atof(require_value("--tolerance"));
    } else if (arg == "--max-iterations") {
      options.maxIterations = std::atoi(require_value("--max-iterations"));
    } else if (arg == "--max-passes") {
      options.maxPasses = std::atoi(require_value("--max-passes"));
    } else if (arg == "--repeat") {
      options.repeat = std::atoi(require_value("--repeat"));
    } else if (arg == "--output") {
      output_path = require_value("--output");
    } else if (arg == "--symmetric") {
      already_symmetric = read_bool_arg(require_value("--symmetric"));
    } else if (arg == "--weighted") {
      weighted = read_bool_arg(require_value("--weighted"));
    } else if (arg == "--verbose") {
      options.verbose = true;
    } else if (arg == "--help" || arg == "-h") {
      usage(argv[0]);
      return 0;
    } else {
      throw std::runtime_error("unknown option: " + arg);
    }
  }

  if (options.gamma <= 0.0) {
    throw std::runtime_error("--gamma must be positive");
  }
  if (options.maxIterations < 1) {
    throw std::runtime_error("--max-iterations must be positive");
  }
  if (options.maxPasses < 1) {
    throw std::runtime_error("--max-passes must be positive");
  }
  if (options.repeat < 1) {
    throw std::runtime_error("--repeat must be positive");
  }

  std::ifstream in(file);
  if (!in) {
    std::cerr << "Error: Cannot open file " << file << "\n";
    return 1;
  }

  bool matrix_symmetric = false;
  std::size_t rows = 0;
  std::size_t cols = 0;
  std::size_t nnz = 0;
  readMtxHeader(in, matrix_symmetric, rows, cols, nnz);

  if (rows == 0 || cols == 0) {
    std::cerr << "Error: Invalid Matrix Market coordinate file: " << file << "\n";
    return 1;
  }
  if (rows != cols) {
    std::cerr << "Error: CPM Leiden expects a square matrix: " << file << "\n";
    return 1;
  }

  using K = uint32_t;
  using V = float;

  DiGraph<K, None, V> graph;
  readMtxOmpW(graph, file, weighted);
  if (!already_symmetric && !matrix_symmetric) {
    graph = symmetricizeOmp(graph);
  }

  const auto begin = Clock::now();
  auto result = cpmLeiden(graph, options);
  const auto end = Clock::now();

  std::vector<int> block_of(rows, -1);
  for (std::size_t u = 0; u < rows; ++u) {
    const std::size_t gve_vertex = u + 1;
    if (gve_vertex < result.membership.size()) {
      block_of[u] = static_cast<int>(result.membership[gve_vertex]);
    }
  }

  std::vector<int> block_sizes;
  std::vector<int> dense = relabel_dense(block_of, &block_sizes);
  const std::size_t community_count = block_sizes.size();

  if (output_path.empty()) {
    output_path = file_stem_local(file) + "_gveleiden_cpm.blk";
  }

  Graph block_graph = build_block_graph_from_gve(graph, dense, static_cast<int>(community_count));
  std::vector<int> block_color;
  const int color_count = greedy_coloring(block_graph, block_color);
  relabel_colors_by_class_size(block_color);

  const std::string bcol_path = default_bcol_path(output_path);
  write_partition_1based(dense, output_path);
  write_block_color_1based(block_color, color_count, bcol_path);

  std::cout << "objective = CPM"
            << " gamma = " << options.gamma
            << " NB = " << community_count
            << " NC = " << color_count
            << " quality = " << result.quality
            // << " time = " << elapsed_seconds(begin, end)
            << "\n";

  return 0;
}
