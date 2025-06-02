#ifndef HELD_KARP_H
#define HELD_KARP_H

#include "tsp_parser.h" 

#include <vector>
#include <utility> 

namespace HeldKarpAlgorithm {

    std::pair<std::vector<int>, double> solve(const TSPData& problem_data);

    vector<int> reconstruct_tour(const vector<vector<int>>& predecessor, 
                                    int num_cities, 
                                    int all_cities_mask, 
                                    int optimal_end_city);
} 

#endif 