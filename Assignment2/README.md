# 20211093 노형진 Assignment 2: Traveling Salesman Problem

This project implements and compares several algorithms for solving the Traveling Salesman Problem (TSP). All algorithms are run sequentially on each provided dataset, and their performance (tour cost, runtime) is reported.

## Core Algorithm Implementations:

*   **`christofides.cpp` / `christofides.h` :**
    *   Source code for Christofides' heuristic algorithm.
    *   *Note:* This implementation uses a greedy approach for the Minimum Weight Perfect Matching (MWPM) step, making it an approximation of the true Christofides algorithm.
*   **`held_karp.cpp` / `held_karp.h` :**
    *   Source code for the Held-Karp algorithm, an exact dynamic programming solution for TSP.
    *   *Note:* Due to its O(N² * 2^N) complexity, this algorithm is automatically skipped in the main test harness for instances with N > 20 cities.
*   **`novel_algorithm.cpp` / `novel_algorithm.h` :**
    *   Source code for the "Centroid-Hull Erosion" novel algorithm developed for this assignment.
    *   The main test harness runs three variations:
        *   **Novel V0 (Base Erosion):** The basic erosion heuristic.
        *   **Novel V1 (Erosion + 2-Opt):** Base erosion followed by 2-Opt refinement.
        *   **Novel V2 (Hull Centroid + 2-Opt):** A variation using hull centroids with 2-Opt refinement.

## Utility and Helper Files:

*   **`tsp_parser.cpp` / `tsp_parser.h` :**
    *   Utility for parsing TSP coordinate data from `.tsp` files (TSPLIB format, specifically `NODE_COORD_SECTION`).
*   **`novel_algorithm_utils.cpp` / `novel_algorithm_utils.h` :**
    *   Helper functions specifically used by the novel algorithm (e.g., Convex Hull computation, 2-Opt refinement, geometric calculations).

## Main Program:

*   **`main.cpp` :**
    *   The main executable program (`tsp_solver`) that drives the experiments.
    *   It takes one or more TSP file paths as command-line arguments.
    *   For each TSP instance, it runs Christofides, Held-Karp (if N <= 20), and all three variations of the Novel Algorithm.
    *   Runtimes for heuristics are averaged over `NUM_RUNS_FOR_AVG` (currently 3) runs.

## Datasets:

The following datasets should be placed in a `dataset/` subdirectory relative to the executable.

*   **`dataset/a280.tsp` :** TSPLIB instance (280 cities).
*   **`dataset/xql662.tsp` :** TTD instance (662 cities).
*   **`dataset/kz9976.tsp` :** TTD instance (9976 cities, Kazakhstan).
*   **`dataset/mona-lisa100K.tsp` :** TTD instance (100,000 "cities", Mona Lisa TSP Challenge).
*   **(Optional) `dataset/a20.tsp`, `dataset/kz20.tsp`, `dataset/mona20.tsp`, `dataset/xql20.tsp` :**
    Smaller N=20 derived instances can be used for quick testing or detailed analysis of Held-Karp. (If you are not providing these, you can remove this line).

## Compilation and Execution:

1.  **Compile the project:**
    Navigate to the project's root directory in your terminal and run:
    ```bash
    make
    ```
    This will use the provided `Makefile` to compile all necessary source files and create the executable named `tsp_solver`.
    To clean compiled files, you can run `make clean`.

2.  **Run the solver:**
    Execute the compiled program, providing paths to one or more TSP dataset files as arguments. The paths should be relative to where you run the executable, or absolute paths.
    ```bash
    ./tsp_solver dataset/a280.tsp
    ```
    To run on multiple datasets:
    ```bash
    ./tsp_solver dataset/a280.tsp dataset/xql662.tsp
    ```
    For very small instances (e.g., N <= 20) to see Held-Karp run:
    ```bash
    # Assuming you have a dataset/a20.tsp file
    ./tsp_solver dataset/a20.tsp
    ```

## Expected Output:

For each algorithm and dataset, the program will output:
*   Algorithm name and instance details.
*   Status (Success, Failed, Skipped, or Warning).
*   Calculated tour cost.
*   Average runtime in milliseconds (ms).
*   The full tour path may be printed for small instances (N <= 30 by default).



