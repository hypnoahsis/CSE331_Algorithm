#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <iomanip>

using namespace std;
using namespace chrono;

const int NUM_ALGOS = 12;
const char* algo_names[NUM_ALGOS] = {
    "MergeSort", "InsertionSort", "SelectionSort", "QuickSort",
    "HeapSort", "BubbleSort", "TimSort", "CocktailShakerSort",
    "CombSort", "TournamentSort", "IntroSort", "LibrarySort"
};

void MergeSort(int[], int, int);
void InsertionSort(int[], int, int);
void SelectionSort(int[], int);
void QuickSort(int[], int, int);
void HeapSort(int[], int);
void BubbleSort(int[], int);
void TimSort(int[], int);
void CocktailShakerSort(int[], int);
void CombSort(int[], int);
void TournamentSort(int[], int);
void IntroSort(int[], int);
void LibrarySort(int[], int);

bool is_sorted(int A[], int n) {
    for (int i = 1; i < n; ++i)
        if (A[i - 1] > A[i]) return false;
    return true;
}

void run_algorithm_on_file(const string& filename, const char* algoname, int trials) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cerr << "Error opening file " << filename << endl;
        return;
    }

    int* A = new int[1000005];
    int n = 0;
    while (fin >> A[n]) ++n;
    fin.close();

    cout << filename << "\n" << algoname << ": ";

    long long total_time = 0;
    bool all_passed = true;

    for (int t = 0; t < trials; ++t) {
        int* B = new int[n];
        memcpy(B, A, sizeof(int) * n);

        auto start = high_resolution_clock::now();

        if (strcmp(algoname, "MergeSort") == 0)
            MergeSort(B, 0, n - 1);
        else if (strcmp(algoname, "InsertionSort") == 0)
            InsertionSort(B, 0, n - 1);
        else if (strcmp(algoname, "SelectionSort") == 0)
            SelectionSort(B, n);
        else if (strcmp(algoname, "QuickSort") == 0)
            QuickSort(B, 0, n - 1);
        else if (strcmp(algoname, "HeapSort") == 0)
            HeapSort(B, n);
        else if (strcmp(algoname, "BubbleSort") == 0)
            BubbleSort(B, n);
        else if (strcmp(algoname, "TimSort") == 0)
            TimSort(B, n);
        else if (strcmp(algoname, "CocktailShakerSort") == 0)
            CocktailShakerSort(B, n);
        else if (strcmp(algoname, "CombSort") == 0)
            CombSort(B, n);
        else if (strcmp(algoname, "TournamentSort") == 0)
            TournamentSort(B, n);
        else if (strcmp(algoname, "IntroSort") == 0)
            IntroSort(B, n);
        else if (strcmp(algoname, "LibrarySort") == 0)
            LibrarySort(B, n);

        auto end = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(end - start).count();
        total_time += duration;

        bool correct = is_sorted(B, n);
        all_passed &= correct;

        cout << duration << (t < trials - 1 ? " " : "");

        delete[] B;
    }

    double avg_time = total_time / static_cast<double>(trials);
    cout << "\nAverage: " << fixed << setprecision(2) << avg_time << " μs "
         << (all_passed ? "✅" : "❌") << "\n\n";

    delete[] A;
}

int main() {
    srand(time(0));

    string prefixes[] = { "sorted", "reverse_sorted", "random", "partial" };
    int sizes[] = { 1000, 10000, 100000, 500000, 1000000 };

    for (const string& prefix : prefixes) {
        for (int size : sizes) {
            int trials = 1;
            if (size == 1000 || size == 10000 || size == 100000) trials = 10;
            else if (size == 500000) trials = 3;

            string filename = prefix + "_" + to_string(size) + ".txt";
            for (int i = 0; i < NUM_ALGOS; ++i)
                run_algorithm_on_file(filename, algo_names[i], trials);
        }
    }

    return 0;
}

