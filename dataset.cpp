#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <fstream>
#include <chrono>
using namespace std;

void write_to_file(const vector<int>& data, const string& filename) {
    ofstream out(filename);
    for (int x : data) out << x << " ";
    out.close();
}

vector<int> generate_sorted(int n) {
    vector<int> v(n);
    for (int i = 0; i < n; ++i) v[i] = i;
    return v;
}

vector<int> generate_reverse_sorted(int n) {
    vector<int> v(n);
    for (int i = 0; i < n; ++i) v[i] = n - i;
    return v;
}

vector<int> generate_random(int n) {
    vector<int> v(n);
    iota(v.begin(), v.end(), 0);
    shuffle(v.begin(), v.end(), default_random_engine(random_device{}()));
    return v;
}

vector<int> generate_partially_sorted(int n) {
    vector<int> v(n);
    for (int i = 0; i < n; ++i) v[i] = i;
    shuffle(v.begin() + n / 4, v.begin() + 3 * n / 4, default_random_engine(random_device{}()));
    return v;
}

int main() {
    vector<int> sizes = {1000, 10000, 100000, 500000, 1000000};
    for (int size : sizes) {
        write_to_file(generate_sorted(size), "sorted_" + to_string(size) + ".txt");
        write_to_file(generate_reverse_sorted(size), "reverse_sorted_" + to_string(size) + ".txt");
        write_to_file(generate_random(size), "random_" + to_string(size) + ".txt");
        write_to_file(generate_partially_sorted(size), "partial_" + to_string(size) + ".txt");
        cout << "Generated datasets for size: " << size << endl;
    }
    return 0;
}