#ifndef CHRISTOFIDES_H
#define CHRISTOFIDES_H

#include "tsp_parser.h" 
#include <vector>
#include <utility> 

struct Edge {
    int city1_idx; 
    int city2_idx; 
    double weight;

    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

namespace ChristofidesAlgorithm {

    std::pair<std::vector<int>, double> solve(const TSPData& problem_data);

    std::vector<Edge> find_minimum_spanning_tree(
        int num_cities,
        const std::vector<std::vector<double>>& distance_matrix
    );

    std::vector<int> find_odd_degree_vertices(
        int num_cities,
        const std::vector<Edge>& mst_edges
    );

    std::vector<Edge> find_greedy_min_weight_perfect_matching(
        const std::vector<int>& odd_degree_city_indices,
        const std::vector<std::vector<double>>& distance_matrix
    );

    std::vector<int> find_eulerian_circuit(
        int num_cities,
        int start_city_idx,
        const std::vector<Edge>& multigraph_edges 
    );

    std::pair<std::vector<int>, double> convert_eulerian_to_hamiltonian(
        const std::vector<int>& eulerian_path_nodes,
        const std::vector<std::vector<double>>& distance_matrix
    );

} 

#endif 