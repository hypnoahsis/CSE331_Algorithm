#ifndef CONV_SORT_H
#define CONV_SORT_H

#include <iostream>
#include <cmath>
#include <climits>

void merge(int A[], int p, int q, int r);

void MergeSort(int A[], int p, int r);

void InsertionSort(int A[], int p, int r);

void SelectionSort(int A[], int n);

int partition(int A[], int p, int r);

void QuickSort(int A[], int p, int r);

void max_heapify(int A[], int n, int i);

void build_max_heap(int A[], int n);

void HeapSort(int A[], int n);

void BubbleSort(int A[], int n);

void TimSort(int A[], int n);

void CocktailShakerSort(int A[], int n);

void CombSort(int A[], int n);

void TournamentSort(int A[], int n);

void IntroSort(int A[], int p, int r, int limit);

void IntroSort(int A[], int n);

int binary_search(int S[], int c, int key);

void rebalance(int S[], int c);

void LibrarySort(int A[], int n);

#endif