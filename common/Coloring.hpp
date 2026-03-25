#pragma once
#include "Types.hpp"

int Greedy_Coloring(const Graph& G, std::vector<int>& color);
void RelabelColorsByClassSize(std::vector<int>& color);