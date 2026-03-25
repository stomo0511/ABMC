#pragma once
#include <iostream>
#include <fstream>

void WriteBlockInfo_1Based(const std::vector<int>& block_of, const std::string& out_path);
void WriteBlockColor_1Based(const std::vector<int>& block_color, int nc, const std::string& out_path);
Graph BuildBlockGraph(const Graph& G,
                      const std::vector<int>& block_of,
                      BlockEdgeWeight mode);