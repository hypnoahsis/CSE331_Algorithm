#ifndef TSP_PARSER_H
#define TSP_PARSER_H

#include <vector>
#include <string>
#include <cmath> 

using namespace std;

struct Point {
    int id;
    double x, y;
    int original_idx;
};

struct TSPData {
    int dimension = 0;
    vector<Point> coords;
    vector<vector<double>> dist_matrix; 
};


inline double EUC_2d_distance(const Point& p1, const Point& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return round(sqrt(dx * dx + dy * dy));
}


TSPData parse_tsp_file(const string& filename);


void build_distance_matrix(TSPData& data);

#endif 