#ifndef NOVEL_ALGORITHM_UTILS_H
#define NOVEL_ALGORITHM_UTILS_H

#include "tsp_parser.h" 
#include <vector>
#include <utility> 

namespace NovelTSPUtils {

    struct PointWithAngle { 
        Point p;
        double angle;
    };

    double calculate_euc_2d_distance_util(const Point& p1, const Point& p2); 

    std::vector<Point> compute_convex_hull(std::vector<Point> points_input);

    int find_point_index_util(const Point& p_to_find, const std::vector<Point>& all_coords); 

    std::pair<std::vector<int>, double> apply_2_opt(
        std::vector<int> tour_indices, 
        const std::vector<std::vector<double>>& dist_matrix
    );

} 

#endif 