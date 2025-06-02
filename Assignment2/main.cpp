#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <limits>
#include <chrono> 
#include <numeric> 

#include "tsp_parser.h"
#include "christofides.h"
#include "held_karp.h"
#include "novel_algorithm.h"

using namespace std;
using namespace std::chrono;


void print_algorithm_output(
    const string& algo_name,
    const string& instance_label,
    int num_cities,
    const pair<vector<int>, double>& solution, 
    double avg_runtime_ms, 
    bool print_full_path = false,
    int path_print_threshold = 50
) {
    cout << "\n--- " << algo_name << " ---" << endl;
    cout << "Instance:         " << instance_label << " (" << num_cities << " cities)" << endl;

    const vector<int>& tour = solution.first;
    double cost = solution.second;

    bool error_condition = false;
    if (cost < 0 && cost != -1.0 && cost != -2.0 && cost != -3.0 && cost != -4.0) { 
         error_condition = true;
    }
    if (cost == numeric_limits<double>::infinity()){
        error_condition = true;
    }
    if (num_cities > 1 && tour.empty() && cost >= 0 && cost != numeric_limits<double>::infinity()) {
        error_condition = true;
    }
     if (num_cities == 1 && tour.size() > 1 && cost >= 0){
        error_condition = true;
     }


    if (error_condition) {
        cout << "Status:           Failed or No Valid Tour" << endl;
        if(cost == numeric_limits<double>::infinity()) cout << "  Reason:         Cost was infinity (likely too large/timeout/error)" << endl;
        else if (cost < 0) cout << "  Reason:         Algorithm returned an error code: " << fixed << setprecision(2) << cost << endl;
        else if (tour.empty() && num_cities > 0) cout << "  Reason:         Tour path is empty" << endl;
        if (avg_runtime_ms >=0) cout << "Avg Runtime:      " << fixed << setprecision(3) << avg_runtime_ms << " ms" << endl;
        else cout << "Avg Runtime:      N/A (Not Run or Error)" << endl;
        return;
    }

     if (num_cities > 0 && tour.size() != static_cast<size_t>(num_cities) && !tour.empty()) {
        cout << "Status:           Warning - Tour size mismatch!" << endl;
        cout << "  Expected Cities: " << num_cities << ", Tour Cities: " << tour.size() << endl;
    } else if (!error_condition) {
         cout << "Status:           Success" << endl;
    }


    cout << "Tour Cost:        " << fixed << setprecision(2) << cost << endl;
    cout << "Avg Runtime:      " << fixed << setprecision(3) << avg_runtime_ms << " ms" << endl;

    if (print_full_path && num_cities <= path_print_threshold && !tour.empty()) {
        cout << "Tour Path:        ";
        for (size_t i = 0; i < tour.size(); ++i) {
            cout << tour[i] << (i == tour.size() - 1 ? "" : " -> ");
        }
        if (tour.size() > 1) {
            cout << " -> " << tour[0];
        }
        cout << endl;
    } else if (!tour.empty()) {
    } else if (num_cities == 1 && cost == 0) {
        cout << "Tour Path:        0 -> 0" << endl;
    }
}


const int NUM_RUNS_FOR_AVG = 3;

void run_all_versions_on_instance(const TSPData& data, const string& instance_file_label) {
    cout << "\n======== Processing Instance: " << instance_file_label << " (" << data.dimension << " cities) ========" << endl;

    high_resolution_clock::time_point t1, t2;
    vector<double> runtimes_ms;
    double total_duration_ms, avg_duration_ms;
    bool print_paths_for_this_instance = (data.dimension <= 30);

    pair<vector<int>, double> solution; 

    // --- Christofides ---
    runtimes_ms.clear();
    total_duration_ms = 0;
    cout << "Running Christofides " << NUM_RUNS_FOR_AVG << " times..." << endl;
    for (int i = 0; i < NUM_RUNS_FOR_AVG; ++i) {
        t1 = high_resolution_clock::now();
        pair<vector<int>, double> current_solution = ChristofidesAlgorithm::solve(data);
        t2 = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
        runtimes_ms.push_back(duration);
        total_duration_ms += duration;
        if (i == 0) solution = current_solution; 
    }
    avg_duration_ms = total_duration_ms / NUM_RUNS_FOR_AVG;
    print_algorithm_output("Christofides Algorithm", instance_file_label, data.dimension, solution, avg_duration_ms, print_paths_for_this_instance);

    // --- Held-Karp ---
    if (data.dimension > 0 && data.dimension <= 20) {
        cout << "Running Held-Karp 1 time..." << endl;
        t1 = high_resolution_clock::now();
        solution = HeldKarpAlgorithm::solve(data);
        t2 = high_resolution_clock::now();
        avg_duration_ms = duration_cast<microseconds>(t2 - t1).count() / 1000.0; 
        print_algorithm_output("Held-Karp Algorithm", instance_file_label, data.dimension, solution, avg_duration_ms, print_paths_for_this_instance);
    } else {
        cout << "\n--- Held-Karp Algorithm ---" << endl;
        cout << "Instance:         " << instance_file_label << " (" << data.dimension << " cities)" << endl;
        cout << "Status:           Skipped (N=" << data.dimension << ", N=0 or N > 20)" << endl;
        cout << "Tour Cost:        N/A" << endl;
        cout << "Avg Runtime:      N/A" << endl;
    }

    // --- Novel V0 ---
    runtimes_ms.clear();
    total_duration_ms = 0;
    cout << "\n" << endl;
    cout << "Running Novel V0 " << NUM_RUNS_FOR_AVG << " times..." << endl;
    for (int i = 0; i < NUM_RUNS_FOR_AVG; ++i) {
        t1 = high_resolution_clock::now();
        pair<vector<int>, double> current_solution = NovelAlgorithm::solve_v0_base_erosion(data);
        t2 = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
        runtimes_ms.push_back(duration);
        total_duration_ms += duration;
        if (i == 0) solution = current_solution;
    }
    avg_duration_ms = total_duration_ms / NUM_RUNS_FOR_AVG;
    print_algorithm_output("Novel V0 (Base Erosion)", instance_file_label, data.dimension, solution, avg_duration_ms, print_paths_for_this_instance);

    // --- Novel V1 ---
    runtimes_ms.clear();
    total_duration_ms = 0;
    cout << "\n" << endl;
    cout << "Running Novel V1 " << NUM_RUNS_FOR_AVG << " times..." << endl;
    for (int i = 0; i < NUM_RUNS_FOR_AVG; ++i) {
        t1 = high_resolution_clock::now();
        pair<vector<int>, double> current_solution = NovelAlgorithm::solve_v1_erosion_plus_2opt(data);
        t2 = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
        runtimes_ms.push_back(duration);
        total_duration_ms += duration;
        if (i == 0) solution = current_solution;
    }
    avg_duration_ms = total_duration_ms / NUM_RUNS_FOR_AVG;
    print_algorithm_output("Novel V1 (Erosion + 2-Opt)", instance_file_label, data.dimension, solution, avg_duration_ms, print_paths_for_this_instance);

    // --- Novel V2 ---
    runtimes_ms.clear();
    total_duration_ms = 0;
    cout << "\n" << endl;
    cout << "Running Novel V2 " << NUM_RUNS_FOR_AVG << " times..." << endl;
    for (int i = 0; i < NUM_RUNS_FOR_AVG; ++i) {
        t1 = high_resolution_clock::now();
        pair<vector<int>, double> current_solution = NovelAlgorithm::solve_v2_hull_centroid_plus_2opt(data);
        t2 = high_resolution_clock::now();
        double duration = duration_cast<microseconds>(t2 - t1).count() / 1000.0;
        runtimes_ms.push_back(duration);
        total_duration_ms += duration;
        if (i == 0) solution = current_solution;
    }
    avg_duration_ms = total_duration_ms / NUM_RUNS_FOR_AVG;
    print_algorithm_output("Novel V2 (Hull Centroid + 2-Opt)", instance_file_label, data.dimension, solution, avg_duration_ms, print_paths_for_this_instance);

    cout << "======================================================\n" << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <tsp_file_path1> [<tsp_file_path2> ...]" << endl;
        return 1;
    }

    for (int i = 1; i < argc; ++i) {
        string filename = argv[i];
        TSPData problem_instance;
        try {
            problem_instance = parse_tsp_file(filename);

            if (problem_instance.dimension == 0 && !filename.empty() && filename.find("example15") == string::npos) {
                 cerr << "Failed to parse TSP file or file is empty/invalid (dimension 0): " << filename << endl;
                 continue;
            }
             if (problem_instance.dimension > 0 && problem_instance.dist_matrix.empty() && problem_instance.coords.empty() && filename.find("example15") == string::npos){
                 cerr << "Parsed dimension for " << filename << " but no coords or matrix to proceed." << endl;
                 continue;
             }

            run_all_versions_on_instance(problem_instance, filename);

        } catch (const std::bad_alloc& ba) {
            cerr << "Memory allocation failed (std::bad_alloc) while processing " << filename 
                 << ". Algorithm likely too memory-intensive for this instance size on this system." << endl;
        }
        catch (const exception& e) {
            cerr << "Exception while processing " << filename << ": " << e.what() << endl;
        }
    }

    return 0;
}