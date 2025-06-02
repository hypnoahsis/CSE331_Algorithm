#include "tsp_parser.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath> 

using namespace std;

TSPData parse_tsp_file(const string& filename) {
    ifstream file(filename);
    TSPData data;
    string line;
    bool reading_coords = false;

    while (getline(file, line)) {
        istringstream iss(line);

        if (line.find("NODE_COORD_SECTION") != string::npos) {
            reading_coords = true;
            continue;
        }
        if (line.find("EOF") != string::npos) {
            reading_coords = false;
            break;
        }

        if (reading_coords) {
            Point p;
            if (iss >> p.id >> p.x >> p.y) {
                p.original_idx = data.coords.size(); 
                data.coords.push_back(p);
            }
        } 
    }
    build_distance_matrix(data); 

    return data;
}

void build_distance_matrix(TSPData& data) {
    data.dimension = data.coords.size(); 
    int n = 0;
    n = data.dimension;
    data.dist_matrix.resize(n, vector<double>(n, 0.0)); 

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dist = EUC_2d_distance(data.coords[i], data.coords[j]);
            data.dist_matrix[i][j] = dist;
            data.dist_matrix[j][i] = dist;
        }
    }
}

int find_point_index(const Point& p_to_find, const vector<Point>&) {
    return p_to_find.original_idx; 
}