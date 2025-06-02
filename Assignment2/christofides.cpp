#include "christofides.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>       
#include <set>      
#include <limits>    

using namespace std; 

namespace ChristofidesAlgorithm {

// 2. MST algorithm (Prim's - Array-based for dense graphs)
vector<Edge> find_minimum_spanning_tree(
    int num_cities,
    const vector<vector<double>>& distance_matrix
) {
    if (num_cities == 0) return {};

    vector<double> min_edge_cost_to_connect(num_cities, numeric_limits<double>::infinity()); // Cost to connect this city to the MST
    vector<int> parent_in_mst(num_cities, -1); // Parent of city[i] in the MST
    vector<bool> is_in_mst(num_cities, false);   // True if city[i] is already included in MST

    min_edge_cost_to_connect[0] = 0;

    vector<Edge> mst_edges;
    mst_edges.reserve(num_cities - 1);

    for (int count = 0; count < num_cities; ++count) {
        int u = -1;
        double min_cost = numeric_limits<double>::infinity();
        for (int v_candidate = 0; v_candidate < num_cities; ++v_candidate) {
            if (!is_in_mst[v_candidate] && min_edge_cost_to_connect[v_candidate] < min_cost) {
                min_cost = min_edge_cost_to_connect[v_candidate];
                u = v_candidate;
            }
        }

        if (u == -1) { 
            if (count < num_cities -1 && num_cities > 1) { 
                 cerr << "    Warning: MST could not connect all cities (Prim's). Graph might be disconnected." << endl;
            }
            break; 
        }

        is_in_mst[u] = true;

        if (parent_in_mst[u] != -1) {
            mst_edges.push_back({parent_in_mst[u], u, distance_matrix[u][parent_in_mst[u]]});
        }

        for (int v_neighbor = 0; v_neighbor < num_cities; ++v_neighbor) {
            if (!is_in_mst[v_neighbor] && distance_matrix[u][v_neighbor] < min_edge_cost_to_connect[v_neighbor]) {
                parent_in_mst[v_neighbor] = u;
                min_edge_cost_to_connect[v_neighbor] = distance_matrix[u][v_neighbor];
            }
        }
    }
    
    if (mst_edges.size() != static_cast<size_t>(num_cities - 1) && num_cities > 1) {
         cerr << "    Warning: Prim's MST found " << mst_edges.size() << " edges, expected " 
              << num_cities - 1 << "." << endl;
    }

    return mst_edges;
}

// 3. Odd-degree vertex identification (This function remains the same)
vector<int> find_odd_degree_vertices(
    int num_cities,
    const vector<Edge>& mst_edges
) {

    vector<int> city_degrees(num_cities, 0);
    for (const auto& edge : mst_edges) {
        city_degrees[edge.city1_idx]++;
        city_degrees[edge.city2_idx]++;
    }

    vector<int> odd_degree_city_indices;
    for (int city_idx = 0; city_idx < num_cities; ++city_idx) {
        if (city_degrees[city_idx] % 2 != 0) {
            odd_degree_city_indices.push_back(city_idx);
        }
    }

    if (odd_degree_city_indices.size() % 2 != 0) {
         cerr << "    Error: Number of odd-degree vertices is not even! This should not happen." << endl;
    }
    return odd_degree_city_indices;
}

// 4. MWPM (Greedy heuristic) 
vector<Edge> find_greedy_min_weight_perfect_matching(
    const vector<int>& odd_degree_city_indices, 
    const vector<vector<double>>& distance_matrix
) {
    if (odd_degree_city_indices.empty()) {
        cout << "    No odd-degree vertices, no matching needed." << endl;
        return {};
    }
    if (odd_degree_city_indices.size() % 2 != 0) {
        cerr << "    Error: Cannot find perfect matching for an odd number of vertices." << endl;
        return {}; 
    }

    vector<Edge> candidate_matching_edges;
    for (size_t i = 0; i < odd_degree_city_indices.size(); ++i) {
        for (size_t j = i + 1; j < odd_degree_city_indices.size(); ++j) {
            int city_u_idx = odd_degree_city_indices[i];
            int city_v_idx = odd_degree_city_indices[j];
            candidate_matching_edges.push_back({city_u_idx, city_v_idx, distance_matrix[city_u_idx][city_v_idx]});
        }
    }

    sort(candidate_matching_edges.begin(), candidate_matching_edges.end());

    vector<Edge> mwpm_edges;
    mwpm_edges.reserve(odd_degree_city_indices.size() / 2);
    set<int> matched_cities_set; 

    for (const auto& edge : candidate_matching_edges) {
        if (matched_cities_set.find(edge.city1_idx) == matched_cities_set.end() &&
            matched_cities_set.find(edge.city2_idx) == matched_cities_set.end()) {
            mwpm_edges.push_back(edge);
            matched_cities_set.insert(edge.city1_idx);
            matched_cities_set.insert(edge.city2_idx);
        }
        if (matched_cities_set.size() == odd_degree_city_indices.size()) {
            break; 
        }
    }

    if (mwpm_edges.size() * 2 != odd_degree_city_indices.size()) {
        cerr << "    Warning: Greedy matching heuristic did not match all odd vertices. Expected "
             << odd_degree_city_indices.size() / 2 << " edges, got " << mwpm_edges.size() << endl;
    }

    return mwpm_edges;
}

// 5. Eulerian circuit (Hierholzer's algorithm) (This function remains the same)
vector<int> find_eulerian_circuit(
    int num_cities,
    int start_city_idx, 
    const vector<Edge>& multigraph_edges_const 
) {
    if (multigraph_edges_const.empty() && num_cities > 1) {
        cerr << "    Error: Multigraph has no edges. Cannot find Eulerian Circuit." << endl;
        return {};
    }
    if (num_cities == 0) return {};
    if (num_cities == 1) return {0}; 

    map<int, vector<pair<int, size_t>>> adjacency_list;
    vector<Edge> current_multigraph_edges = multigraph_edges_const; 
    vector<bool> edge_is_used(current_multigraph_edges.size(), false);
    vector<int> city_degrees(num_cities, 0);


    for (size_t i = 0; i < current_multigraph_edges.size(); ++i) {
        const auto& edge = current_multigraph_edges[i];
        adjacency_list[edge.city1_idx].push_back({edge.city2_idx, i});
        adjacency_list[edge.city2_idx].push_back({edge.city1_idx, i});
        city_degrees[edge.city1_idx]++;
        city_degrees[edge.city2_idx]++;
    }

    for (int i = 0; i < num_cities; ++i) {
        if (city_degrees[i] % 2 != 0) {
            cerr << "    Error: Graph for Eulerian circuit is not Eulerian. City " << i 
                 << " has odd degree " << city_degrees[i] << "." << endl;
            return {}; 
        }
    }

    vector<int> eulerian_path_nodes;
    vector<int> current_path_stack;

    int actual_start_city_idx = start_city_idx;
    if (city_degrees[actual_start_city_idx] == 0 && !current_multigraph_edges.empty()) {
        bool found_valid_start = false;
        for(int i = 0; i < num_cities; ++i) {
            if(city_degrees[i] > 0) {
                actual_start_city_idx = i;
                found_valid_start = true;
                break;
            }
        }
        if (!found_valid_start) {
             cerr << "    Error: No node with positive degree found, but edges exist." << endl; return {};
        }
    }
    
    current_path_stack.push_back(actual_start_city_idx);
    int current_city_in_path = actual_start_city_idx;

    while (!current_path_stack.empty()) {
        current_city_in_path = current_path_stack.back();
        bool found_unused_edge = false;

        for (size_t i = 0; i < adjacency_list[current_city_in_path].size(); ++i) {
            pair<int, size_t> neighbor_info = adjacency_list[current_city_in_path][i];
            int neighbor_city_idx = neighbor_info.first;
            size_t edge_original_idx = neighbor_info.second;

            if (!edge_is_used[edge_original_idx]) {
                edge_is_used[edge_original_idx] = true; 

                adjacency_list[current_city_in_path].erase(adjacency_list[current_city_in_path].begin() + i);

                for (size_t j = 0; j < adjacency_list[neighbor_city_idx].size(); ++j) {
                    if (adjacency_list[neighbor_city_idx][j].second == edge_original_idx) {
                        adjacency_list[neighbor_city_idx].erase(adjacency_list[neighbor_city_idx].begin() + j);
                        break;
                    }
                }
                current_path_stack.push_back(neighbor_city_idx); 
                found_unused_edge = true;
                break; 
            }
        }
        
        if (!found_unused_edge) { 
            eulerian_path_nodes.push_back(current_city_in_path);
            current_path_stack.pop_back();
        }
    }

    reverse(eulerian_path_nodes.begin(), eulerian_path_nodes.end()); 

    if (!multigraph_edges_const.empty() && eulerian_path_nodes.size() != multigraph_edges_const.size() + 1) {
        cout << "    Note: Eulerian path node count (" << eulerian_path_nodes.size() 
             << ") vs expected (" << multigraph_edges_const.size() + 1 << ")." << endl;
    }

    return eulerian_path_nodes;
}

// 6. Hamiltonian circuit conversion (shortcutting) (This function remains the same)
pair<vector<int>, double> convert_eulerian_to_hamiltonian(
    const vector<int>& eulerian_path_nodes,
    const vector<vector<double>>& distance_matrix
) {
    if (eulerian_path_nodes.empty()) {
        if (distance_matrix.size() == 1) return {{0}, 0.0}; 
        cerr << "    Error: Eulerian path is empty, cannot convert." << endl;
        return {{}, 0.0};
    }

    vector<int> hamiltonian_tour_nodes;
    hamiltonian_tour_nodes.reserve(distance_matrix.size()); 
    vector<bool> city_is_visited(distance_matrix.size(), false);
    double total_tour_cost = 0.0;

    int first_city_in_tour = eulerian_path_nodes[0];
    hamiltonian_tour_nodes.push_back(first_city_in_tour);
    city_is_visited[first_city_in_tour] = true;
    int last_city_added_to_hamiltonian = first_city_in_tour;

    for (size_t i = 1; i < eulerian_path_nodes.size(); ++i) {
        int next_city_in_eulerian = eulerian_path_nodes[i];
        if (!city_is_visited[next_city_in_eulerian]) {
            hamiltonian_tour_nodes.push_back(next_city_in_eulerian);
            city_is_visited[next_city_in_eulerian] = true;
            total_tour_cost += distance_matrix[last_city_added_to_hamiltonian][next_city_in_eulerian];
            last_city_added_to_hamiltonian = next_city_in_eulerian;
        }
    }

    if (hamiltonian_tour_nodes.size() > 1) { 
        total_tour_cost += distance_matrix[last_city_added_to_hamiltonian][hamiltonian_tour_nodes[0]];
    }
    
    size_t num_cities_in_problem = distance_matrix.size();
    if (hamiltonian_tour_nodes.size() != num_cities_in_problem && num_cities_in_problem > 0) {
        cerr << "    Warning: Hamiltonian tour does not visit all cities! Visited: " 
             << hamiltonian_tour_nodes.size() << "/" << num_cities_in_problem << endl;
    }

    return {hamiltonian_tour_nodes, total_tour_cost};
}


// Main solve function for Christofides 
pair<vector<int>, double> solve(const TSPData& problem_data) {

    int num_cities = problem_data.dimension;

    if (num_cities == 0) {
        cerr << "  Error: Number of cities is 0. Cannot solve." << endl;
        return {{}, 0.0};
    }
    if (num_cities == 1) {
        cout << "  Solved: TSP for 1 city. Tour: {0}, Cost: 0.0" << endl;
        return {{0}, 0.0}; 
    }
    if (problem_data.dist_matrix.empty() || problem_data.dist_matrix.size() != static_cast<size_t>(num_cities)) {
        cerr << "  Error: Distance matrix is invalid or mismatched with dimension." << endl;
        return {{}, 0.0};
    }

    vector<Edge> mst = find_minimum_spanning_tree(num_cities, problem_data.dist_matrix);
    if (mst.empty() && num_cities > 1) { 
        cerr << "  Error: MST could not be formed. Aborting Christofides." << endl;
        return {{}, -1.0}; 
    }

    vector<int> odd_vertices_indices = find_odd_degree_vertices(num_cities, mst);
    vector<Edge> mwpm = find_greedy_min_weight_perfect_matching(odd_vertices_indices, problem_data.dist_matrix);
    
    vector<Edge> multigraph_all_edges = mst;
    multigraph_all_edges.insert(multigraph_all_edges.end(), mwpm.begin(), mwpm.end());

    int euler_start_node = 0;
    if (!multigraph_all_edges.empty()) {
        bool start_node_has_edges = false;
        for(const auto& edge : multigraph_all_edges) {
            if (edge.city1_idx == 0 || edge.city2_idx == 0) {
                start_node_has_edges = true;
                break;
            }
        }
        if (!start_node_has_edges && !multigraph_all_edges.empty()) { 
             euler_start_node = multigraph_all_edges[0].city1_idx;
        }
    } else if (num_cities == 1) { 
        euler_start_node = 0;
    } else if (num_cities > 1 && multigraph_all_edges.empty()){ 
         cerr << "  Error: Multigraph for Eulerian circuit is empty for N > 1. Aborting." << endl;
         return {{}, -1.0};
    }


    vector<int> eulerian_circuit_nodes = find_eulerian_circuit(num_cities, euler_start_node, multigraph_all_edges);
    if (eulerian_circuit_nodes.empty() && num_cities > 1 && !multigraph_all_edges.empty()) {
        cerr << "  Error: Eulerian circuit could not be formed. Aborting Christofides." << endl;
        return {{}, -1.0};
    }
    if (eulerian_circuit_nodes.empty() && num_cities == 1) eulerian_circuit_nodes = {0}; 


    pair<vector<int>, double> final_tour_and_cost = convert_eulerian_to_hamiltonian(
        eulerian_circuit_nodes,
        problem_data.dist_matrix
    );

    return final_tour_and_cost;
}


} // namespace ChristofidesAlgorithm