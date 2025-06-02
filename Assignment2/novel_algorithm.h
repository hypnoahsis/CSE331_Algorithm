#ifndef NOVEL_ALGORITHM_H
#define NOVEL_ALGORITHM_H

#include "tsp_parser.h"
#include <vector>
#include <utility>

namespace NovelAlgorithm {

    // V0: Base Centroid-Hull Erosion
    std::pair<std::vector<int>, double> solve_v0_base_erosion(const TSPData& problem_data);

    // V1: Base Centroid-Hull Erosion + 2-Opt Refinement
    std::pair<std::vector<int>, double> solve_v1_erosion_plus_2opt(const TSPData& problem_data);

    // V2: Hull-Centroid for Angular Sort + 2-Opt Refinement
    std::pair<std::vector<int>, double> solve_v2_hull_centroid_plus_2opt(const TSPData& problem_data); 

} 

#endif 