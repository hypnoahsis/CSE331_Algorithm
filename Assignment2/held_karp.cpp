#include "held_karp.h"
#include <iostream>
#include <vector>
#include <algorithm>   
#include <limits>      

using namespace std; 

namespace HeldKarpAlgorithm {

    pair<vector<int>, double> solve(const TSPData& problem_data) {

        int num_cities = problem_data.dimension;
        if (num_cities == 0) {
            cerr << "Error: Number of cities is 0. Cannot solve." << endl;
            return {{}, 0.0}; 
        }

        // DP table
        vector<vector<double>> dp(1 << num_cities, vector<double>(num_cities, numeric_limits<double>::infinity()));

        // Predecessor table.  Keep track of best predecessor city along min path.
        vector<vector<int>> predecessor(1 << num_cities, vector<int>(num_cities, -1));

        // Base case: start at city 0, only visit city 0, cost is 0
        dp[1][0] = 0;
        predecessor[1][0] = 0;

        // Iterate through subsets in increasing order of size:
        int all_cities_mask = (1 << num_cities) - 1;

        for (int mask = 1; mask <= all_cities_mask; ++mask) {
            for (int end_city = 0; end_city < num_cities; ++end_city) {
                if ((mask & (1 << end_city)) != 0) {
                    for (int last_city = 0; last_city < num_cities; ++last_city) {
                        if (last_city != end_city && (mask & (1 << last_city)) != 0) {
                            int prev_mask = mask ^ (1 << end_city); 

                            if (dp[prev_mask][last_city] != numeric_limits<double>::infinity()) {
                                double current_cost = dp[prev_mask][last_city] + problem_data.dist_matrix[last_city][end_city];
                                if (current_cost < dp[mask][end_city]) {
                                    dp[mask][end_city] = current_cost;
                                    predecessor[mask][end_city] = last_city; 
                                }
                            }
                        }
                    }
                }
            }
        }

       // Find the minimum tour cost by returning from each city to the start city
        double min_tour_cost = numeric_limits<double>::infinity();
        int optimal_end_city = -1;  

        for (int end_city = 1; end_city < num_cities; ++end_city) { 
            if (dp[all_cities_mask][end_city] != numeric_limits<double>::infinity()) {
                double tour_cost = dp[all_cities_mask][end_city] + problem_data.dist_matrix[end_city][0];
                if (tour_cost < min_tour_cost) {
                    min_tour_cost = tour_cost;
                    optimal_end_city = end_city; 
                }
            }
        }
       
        vector<int> optimal_tour;
        if (optimal_end_city != -1) {
           optimal_tour = reconstruct_tour(predecessor, num_cities, all_cities_mask, optimal_end_city);
           
           if (optimal_tour.front() != 0) {
                cerr << "Warning: reconstructed tour did not start with node 0." << endl;
           }
        }
        else {
            cout << "No Optimal tour reconstruction" << endl;
             return {optimal_tour, numeric_limits<double>::infinity()}; 
        }
        return {optimal_tour, min_tour_cost};
    }

    // Reconstruct tour using the predecessor table and the calculated costs
    vector<int> reconstruct_tour(const vector<vector<int>>& predecessor,
                                    int num_cities, 
                                    int all_cities_mask, 
                                    int optimal_end_city) { 
        
        if (num_cities <= 0) return {};
        
        vector<int> tour;
        vector<int> path_stack; 
        int current_city = optimal_end_city; 
        int current_mask = all_cities_mask; 

        //Reconstruct path until at start city.
        while (current_city != 0) {

            path_stack.push_back(current_city);  

            int last_city = predecessor[current_mask][current_city];
            if (last_city == -1 || last_city == current_city) {
                cerr << "Error: Invalid tour or predecessor table. Aborting reconstruction at city "
                     << current_city << " and mask " << current_mask << endl;
                return {}; 
            }
            current_mask = current_mask ^ (1 << current_city); 
            current_city = last_city; 
        }

       //Add start and reverse the tour so it starts with 0;
        vector<int> generated_tour;
        generated_tour.push_back(0);
        while (!path_stack.empty()) {
            generated_tour.push_back(path_stack.back());
            path_stack.pop_back();
        }
        tour = generated_tour;

        return tour;
    }

} // namespace