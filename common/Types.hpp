#pragma once
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

// 無向グラフ + エッジに weight プロパティを持つ
// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
//                   boost::no_property,
//                   boost::property<boost::edge_weight_t, double>> Graph;
// typedef boost::graph_traits<Graph>::edge_descriptor EdgeIterator;
// typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
// typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIter;

// typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;

// Block graph T を作る（binary/count/abs-sum を選べる）
enum class BlockEdgeWeight { Binary, Count, AbsSum };

// エッジのプロパティ：重み
using EdgeWeightProperty = boost::property<boost::edge_weight_t, double>;

// 無向グラフ（隣接リスト表現）
using Graph = boost::adjacency_list<
    boost::vecS,          // 頂点コンテナ
    boost::vecS,          // エッジコンテナ
    boost::undirectedS,   // 無向グラフ
    boost::no_property,   // 頂点プロパティなし
    EdgeWeightProperty    // エッジプロパティ
>;

// エッジ重みアクセス用マップ
using WeightMap = boost::property_map<Graph, boost::edge_weight_t>::type;

// グラフ関連の便利型
using Vertex       = boost::graph_traits<Graph>::vertex_descriptor;
using Edge         = boost::graph_traits<Graph>::edge_descriptor;
using AdjacencyIter= boost::graph_traits<Graph>::adjacency_iterator;
using EdgeIter     = boost::graph_traits<Graph>::edge_iterator;
