#include "novel_algorithm.h"
#include "novel_algorithm_utils.h" 
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace std; 

namespace NovelAlgorithm {

pair<vector<int>, double> construct_centroid_hull_erosion_tour(
    const TSPData& problem_data,
    bool use_hull_centroid
) {
    int num_cities = problem_data.dimension;
    const auto& coords = problem_data.coords;
    const auto& dist_matrix = problem_data.dist_matrix;

    if (num_cities == 0) return {{}, 0.0};
    if (num_cities == 1) { 
        if (coords.empty()) {
            cerr << "Error (Novel Algo Core): N=1 but coords vector is empty." << endl;
            return {{}, -1.0}; 
        }
        return {{ NovelTSPUtils::find_point_index_util(coords[0], coords) }, 0.0};
    }


    // 1. Compute Convex Hull
    vector<Point> all_points_copy = coords;
    vector<Point> hull_point_structs = NovelTSPUtils::compute_convex_hull(all_points_copy);

    vector<int> main_tour_indices;
    main_tour_indices.reserve(num_cities);

    if (hull_point_structs.empty() && num_cities > 0) {
         cerr << "Error (Novel Algo Core): Convex hull is empty for N=" << num_cities << ". Fallback to sequential." << endl;
         for(int i=0; i<num_cities; ++i) main_tour_indices.push_back(i);
    } else if (hull_point_structs.size() == 1) {
        main_tour_indices.push_back(NovelTSPUtils::find_point_index_util(hull_point_structs[0], coords));
    } else { 
        for(const auto& p_struct : hull_point_structs) {
            int original_idx = NovelTSPUtils::find_point_index_util(p_struct, coords);
            if (original_idx != -1) {
                main_tour_indices.push_back(original_idx);
            } else {
                return {{}, -2.0}; 
            }
        }
    }

    vector<bool> is_in_tour(num_cities, false);
    for (int city_idx : main_tour_indices) {
        if(city_idx >=0 && city_idx < num_cities) is_in_tour[city_idx] = true;
        else {return {{}, -3.0}; } 
    }

    // 2. Compute Centroid
    Point centroid = {0, 0.0, 0.0, -1}; 
    if (use_hull_centroid) { 
        if (!hull_point_structs.empty()) {
            for (const auto& p_struct : hull_point_structs) {
                centroid.x += p_struct.x;
                centroid.y += p_struct.y;
            }
            centroid.x /= hull_point_structs.size();
            centroid.y /= hull_point_structs.size();
        } else { 
             if (num_cities > 0) {
                for (const auto& p : coords) { centroid.x += p.x; centroid.y += p.y; }
                centroid.x /= num_cities; centroid.y /= num_cities;
                cout << "  (Hull empty, using ALL points centroid for angular sort)" << endl;
             }
        }
    } else { 
        if (num_cities > 0) {
            for (const auto& p : coords) {
                centroid.x += p.x;
                centroid.y += p.y;
            }
            centroid.x /= num_cities;
            centroid.y /= num_cities;
        }
    }


    // 3. Identify and Sort Interior Points by Angle
    vector<NovelTSPUtils::PointWithAngle> interior_points_to_process;
    for (int i = 0; i < num_cities; ++i) {
        if (!is_in_tour[i]) {
            double angle = atan2(coords[i].y - centroid.y, coords[i].x - centroid.x);
            interior_points_to_process.push_back({coords[i], angle});
        }
    }
    sort(interior_points_to_process.begin(), interior_points_to_process.end(),
            [](const NovelTSPUtils::PointWithAngle& a, const NovelTSPUtils::PointWithAngle& b) {
                return a.angle < b.angle;
            });

    // 4. Iteratively Insert Interior Points
    for (const auto& pwa : interior_points_to_process) {
        int point_to_insert_idx = NovelTSPUtils::find_point_index_util(pwa.p, coords);

        if (point_to_insert_idx == -1) { continue; }
        if (is_in_tour[point_to_insert_idx]) continue; 

        double min_insertion_cost_increase = numeric_limits<double>::infinity();
        size_t best_insertion_pos = 0;

        if (main_tour_indices.empty()) { 
             main_tour_indices.push_back(point_to_insert_idx);
        } else if (main_tour_indices.size() == 1) {
             main_tour_indices.push_back(point_to_insert_idx);
        } else {
            for (size_t i = 0; i < main_tour_indices.size(); ++i) {
                int city_u_idx = main_tour_indices[i];
                int city_v_idx = main_tour_indices[(i + 1) % main_tour_indices.size()];

                double cost_increase = dist_matrix[city_u_idx][point_to_insert_idx] +
                                       dist_matrix[point_to_insert_idx][city_v_idx] -
                                       dist_matrix[city_u_idx][city_v_idx];

                if (cost_increase < min_insertion_cost_increase) {
                    min_insertion_cost_increase = cost_increase;
                    best_insertion_pos = i;
                }
            }
            main_tour_indices.insert(main_tour_indices.begin() + best_insertion_pos + 1, point_to_insert_idx);
        }
        if(point_to_insert_idx >=0 && point_to_insert_idx < num_cities) is_in_tour[point_to_insert_idx] = true;
        else {return{{},-4.0};}
    }

    // Calculate final tour cost
    double final_tour_cost = 0;
    if (main_tour_indices.size() >= 2) {
        for (size_t i = 0; i < main_tour_indices.size(); ++i) {
            final_tour_cost += dist_matrix[main_tour_indices[i]][main_tour_indices[(i + 1) % main_tour_indices.size()]];
        }
    } else if (main_tour_indices.size() == 1 && num_cities == 1) { 
        final_tour_cost = 0;
    } else if (main_tour_indices.empty() && num_cities == 0) {
        final_tour_cost = 0;
    } else if (main_tour_indices.size() < 2 && num_cities > 1) {
        cerr << "Warning (Novel Algo Core): Tour has < 2 cities for N > 1. Cost might be invalid." << endl;
    }

    if (num_cities > 0 && main_tour_indices.size() != static_cast<size_t>(num_cities)) {
        cerr << "    Warning (Novel Algo Core): Final constructed tour size (" << main_tour_indices.size()
             << ") does not match num_cities (" << num_cities << ")." << endl;
    }

    return {main_tour_indices, final_tour_cost};
}


// V0: Base Centroid-Hull Erosion (uses centroid of ALL points)
pair<vector<int>, double> solve_v0_base_erosion(const TSPData& problem_data) {
    return construct_centroid_hull_erosion_tour(problem_data, false); 
}

// V1: Base Centroid-Hull Erosion + 2-Opt Refinement (uses centroid of ALL points)
pair<vector<int>, double> solve_v1_erosion_plus_2opt(const TSPData& problem_data) {
    pair<vector<int>, double> initial_solution = construct_centroid_hull_erosion_tour(problem_data, false); 

    if (initial_solution.first.empty() && problem_data.dimension > 0) {
        return initial_solution;
    }

    return NovelTSPUtils::apply_2_opt(initial_solution.first, problem_data.dist_matrix);
}

// V2: Hull-Centroid for Angular Sort + 2-Opt Refinement
pair<vector<int>, double> solve_v2_hull_centroid_plus_2opt(const TSPData& problem_data) { 
    pair<vector<int>, double> initial_solution = construct_centroid_hull_erosion_tour(problem_data, true); 

    if (initial_solution.first.empty() && problem_data.dimension > 0) {
        return initial_solution;
    }

    return NovelTSPUtils::apply_2_opt(initial_solution.first, problem_data.dist_matrix);
}

} 