#include "novel_algorithm_utils.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>

namespace NovelTSPUtils {

    double calculate_euc_2d_distance_util(const Point& p1, const Point& p2) {
        double dx = p1.x - p2.x;
        double dy = p1.y - p2.y;
        return std::round(std::sqrt(dx * dx + dy * dy));
    }

    double cross_product_util(Point O, Point A, Point B) {
        return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
    }

    std::vector<Point> compute_convex_hull(std::vector<Point> points_input) {
        std::vector<Point> P = points_input;
        int n = P.size();
        if (n <= 2) {
            return P;
        }

        std::sort(P.begin(), P.end(), [](const Point& a, const Point& b) {
            return a.x < b.x || (a.x == b.x && a.y < b.y);
        });

        P.erase(std::unique(P.begin(), P.end(), [](const Point& a, const Point& b){
            return std::abs(a.x - b.x) < 1e-9 && std::abs(a.y - b.y) < 1e-9;
        }), P.end());
        n = P.size();
        if (n <= 2) return P;

        std::vector<Point> hull;
        hull.reserve(2 * n);

        for (int i = 0; i < n; ++i) {
            while (hull.size() >= 2 && cross_product_util(hull[hull.size()-2], hull.back(), P[i]) <= 0) {
                hull.pop_back();
            }
            hull.push_back(P[i]);
        }
        for (int i = n - 2, t = hull.size() + 1; i >= 0; --i) {
            while (hull.size() >= t && cross_product_util(hull[hull.size()-2], hull.back(), P[i]) <= 0) {
                hull.pop_back();
            }
            hull.push_back(P[i]);
        }
        if (hull.size() > 1) {
             hull.pop_back(); 
        }
        return hull;
    }

    int find_point_index_util(const Point& p_to_find, const std::vector<Point>&) {
        return p_to_find.original_idx;
    }

    std::pair<std::vector<int>, double> apply_2_opt(
        std::vector<int> tour_indices,
        const std::vector<std::vector<double>>& dist_matrix
    ) {
        if (tour_indices.size() < 4) {
            double current_cost = 0;
            if(tour_indices.size() >= 2){
                for (size_t i = 0; i < tour_indices.size(); ++i) {
                    current_cost += dist_matrix[tour_indices[i]][tour_indices[(i + 1) % tour_indices.size()]];
                }
            } else if (tour_indices.size() == 1) current_cost = 0;
            return {tour_indices, current_cost};
        }

        bool improvement_found = true;
        int n = tour_indices.size();
        double current_best_cost = 0;
        for (int i = 0; i < n; ++i) {
            current_best_cost += dist_matrix[tour_indices[i]][tour_indices[(i + 1) % n]];
        }

        while (improvement_found) {
            improvement_found = false;

            for (int i = 0; i < n - 1; ++i) {
                for (int j = i + 2; j < n; ++j) { 

                    int city_A = tour_indices[i];
                    int city_B = tour_indices[(i + 1)]; 
                    int city_C = tour_indices[j];
                    int city_D = tour_indices[(j + 1) % n]; 

   
                    if (city_D == city_A) continue;


                    double cost_delta = (dist_matrix[city_A][city_C] + dist_matrix[city_B][city_D]) -
                                        (dist_matrix[city_A][city_B] + dist_matrix[city_C][city_D]);

                    if (cost_delta < -1e-9) { 
                        std::reverse(tour_indices.begin() + (i + 1), tour_indices.begin() + j + 1);
                        current_best_cost += cost_delta;
                        improvement_found = true;
                    }
                }
            }
        }
        return {tour_indices, current_best_cost};
    }

} 