#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gve-leiden-inc/main.hxx"
#include "inc/cpm_properties.hxx"
#include "inc/cpm_leiden.hxx"

using Clock = std::chrono::steady_clock;

static double elapsed_seconds(Clock::time_point begin, Clock::time_point end)
{
  return std::chrono::duration<double>(end - begin).count();
}

static std::string file_stem_local(const std::string& path)
{
  const std::size_t slash = path.find_last_of("/\\");
  const std::string name = slash == std::string::npos ? path : path.substr(slash + 1);
  const std::size_t dot = name.find_last_of('.');
  return dot == std::string::npos ? name : name.substr(0, dot);
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
    "[--output <path>] [--symmetric 0|1] [--weighted 0|1] [--verbose]\n",
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
  bool already_symmetric = false;
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

  if (output_path.empty()) {
    output_path = file_stem_local(file) + "_gveleiden_cpm.blk";
  }
  write_partition_1based(block_of, output_path);

  std::cout << "objective = CPM"
            << " gamma = " << options.gamma
            << " NB = " << result.communities
            << " quality = " << result.quality
            << " time = " << elapsed_seconds(begin, end)
            << " output = " << output_path
            << "\n";

  return 0;
}
