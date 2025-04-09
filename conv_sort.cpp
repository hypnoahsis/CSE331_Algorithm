#include <iostream>
#include <cmath>
#include <climits>
#include <cstring>
#include <algorithm>
#include "conv_sort.h"

using namespace std;

void merge(int A[], int p, int q, int r)
{
    int n_L = q - p + 1;
    int n_R = r - q;
    int* L = new int[n_L];
    int* R = new int[n_R];

    for (int i = 0; i < n_L; i++)
    {
        L[i] = A[p + i];
    }
    for (int j = 0; j < n_R; j++)
    {
        R[j] = A[q + j + 1];
    }

    int i = 0, j = 0, k = p;
    while (i < n_L && j < n_R)
    {
        if (L[i] <= R[j])
        {
            A[k++] = L[i++];
        }
        else A[k++] = R[j++];
    }

    while (i < n_L)
    {
        A[k++] = L[i++];
    }
    while (j < n_R)
    {
        A[k++] = R[j++];
    }

    delete[] L;
    delete[] R;
}

void MergeSort(int A[], int p, int r)
{
    if (p >= r) return;
    int q = p + (r - p) / 2;
    MergeSort(A, p, q);
    MergeSort(A, q + 1, r);
    merge(A, p, q, r);
}

void InsertionSort(int A[], int p, int r)
{
    for (int i = p + 1; i <= r; i++)
    {
        int key = A[i];
        int j = i - 1;
        while (j >= p && A[j] > key)
        {
            A[j + 1] = A[j];
            j--;
        }
        A[j + 1] = key;
    }
}

void SelectionSort(int A[], int n)
{
    for (int i = 0; i < n - 1; i++)
    {
        int min = i;
        for (int j = i + 1; j < n; j++)
        {
            if (A[j] < A[min]) min = j;
        }
        swap(A[i], A[min]);
    }
}

int partition(int A[], int p, int r)
{
    int x = A[r];
    int i = p - 1;
    for (int j = p; j < r; j++)
    {
        if (A[j] <= x)
        {
            i++;
            swap(A[i], A[j]);
        }
    }
    swap(A[i + 1], A[r]);

    return i + 1;
}

int randomized_partition(int A[], int p, int r)
{
    int i = p + rand() % (r - p + 1);
    swap(A[i], A[r]);
    return partition(A, p, r);
}

void QuickSort(int A[], int p, int r)
{
    if (p < r)
    {
        int q = randomized_partition(A, p, r);
        QuickSort(A, p, q - 1);
        QuickSort(A, q + 1, r);
    }
}

void max_heapify(int A[], int n, int i)
{
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * (i + 1);
    if (l < n && A[l] > A[i]) largest = l;
    else largest = i;
    if (r < n && A[r] > A[largest]) largest = r;
    if (largest != i)
    {
        swap(A[i], A[largest]);
        max_heapify(A, n, largest);
    }
}

void build_max_heap(int A[], int n)
{
    for (int i = n / 2 - 1; i >= 0; i--)
    {
        max_heapify(A, n, i);
    }
}

void HeapSort(int A[], int n)
{
    build_max_heap(A, n);
    for(int i = n - 1; i > 0; i--)
    {
        swap(A[0], A[i]);
        max_heapify(A, i, 0);
    }
}

void BubbleSort(int A[], int n)
{
    for (int i = 0; i < n - 1; i++)
    {
        bool swp = false;
        for (int j = 0; j < n - i - 1; j++)
        {
            if (A[j] > A[j + 1])
            {
                swap(A[j], A[j + 1]);
                swp = true;
            }
        }
        if (!swp) break;
    }
}



void TimSort(int A[], int n)
{
    int run = 32;
    for (int i = 0; i < n; i += run)
    {
        int r = min(i + run - 1, n - 1);
        InsertionSort(A, i, r);
    }

    for (int j = run; j < n; j *= 2)
    {
        for (int p = 0; p < n; p += 2 * j)
        {
            int q = min(p + j - 1, n - 1);
            int r = min(p + 2 * j - 1, n - 1);

            if (q < r) merge(A, p, q, r);

        }

    }
}

void CocktailShakerSort(int A[], int n)
{
    bool swp = true;
    int p = 0;
    int r = n - 1;
    while (swp)
    {
        swp = false;
        for (int i = p; i < r; i++)
        {
            if (A[i] > A[i + 1])
            {
                swap(A[i], A[i + 1]);
                swp = true;
            }
        }
        if (!swp) break;
        swp = false;
        r--;

        for (int i = r - 1; i >= p; i--)
        {
            if (A[i] > A[i + 1])
            {
                swap(A[i], A[i + 1]);
                swp = true;
            }
        }
        p++;
    }
}

void CombSort(int A[], int n)
{
    int gap = n;
    const double shrink = 1.3;
    bool swp = true;

    while (gap > 1 || swp)
    {
        gap = int(gap / shrink);
        if (gap < 1) gap = 1;

        swp = false;

        for (int i = 0; i + gap < n; i++)
        {
            if (A[i] > A[i + gap])
            {
                swap(A[i], A[i + gap]);
                swp = true;
            }
        }
    }
}

void TournamentSort(int A[], int n)
{
    int m = 1;
    while (m < n) m *= 2;

    int* T = new int[2 * m - 1];

    for (int i = 0; i < m; i++)
    {
        if (i < n) T[m - 1 + i] = A[i];
        else T[m - 1 + i] = INT32_MAX;
    }
    
    for (int i = m - 2; i >= 0; i--)
    {
        T[i] = min(T[2 * i + 1], T[2 * (i + 1)]);
    }

    for (int j = 0; j < n; j++)
    {
        A[j] = T[0];
        int i = m - 1;
        while (i < 2 * m - 1 && T[i] != T[0]) i++;
        T[i] = INT32_MAX;

        while (i > 0)
        {
            i = (i - 1) / 2;
            T[i] = min(T[2 * i + 1], T[2 * (i + 1)]);
        }
    }
    delete[] T;
}

void IntroSort(int A[], int p, int r, int limit)
{
    if (r - p + 1 <= 16)
    {
        InsertionSort(A, p, r);
        return;
    }

    if (limit == 0)
    {
        int n = r - p + 1;
        int* a = new int[n];
        for (int i = 0; i < n; i++)
        {
            a[i] = A[p + i];
        }
        HeapSort(a, n);
        for (int i = 0; i < n; i++)
        {
            A[p + i] = a[i];
        }
        delete[] a;
        return;
    }

    int q = partition(A, p, r);
    IntroSort(A, p, q - 1, limit - 1);
    IntroSort(A, q + 1, r, limit - 1);
}

void IntroSort(int A[], int n)
{
    int limit = 2 * ceil(log2(n));
    IntroSort(A, 0, n - 1, limit);
}

const int gap = -1;
int binary_search(int S[], int c, int key)
{
    int left = 0, right = c - 1;
    int ans = c;
    while (left <= right)
    {
        int mid = (left + right) / 2;
        int actual = 2 * mid;
        if (S[actual] == gap || S[actual] >= key)
        {
            ans = mid;
            right = mid - 1;
        }
        else
            left = mid + 1;
    }

    return 2 * ans;
}

void rebalance(int S[], int c)
{
    int r = c - 1;
    int w = 2 * (c - 1);

    while (r >= 0)
    {
        S[w + 1] = gap;
        S[w] = S[2 * r];
        r--;
        w -= 2;
    }
}

void LibrarySort(int A[], int n)
{
    int size = 2 * n * ((int)floor(log2(n)) + 2);
    int* S = new int[size];
    for (int i = 0; i < size; i++) S[i] = gap;
    S[0] = A[0];
    int c = 1;
    int phases = (int)floor(log2(n)) + 1;
    for (int i = 1; i <= phases; i++)
    {
        rebalance(S, c);
        c *= 2;
        int start = (1 << (i - 1));
        int end = min((1 << i), n);
        for (int j = start; j < end; j++)
        {
            int key = A[j];
            int idx = binary_search(S, c, key);
            int pos = idx;

            int k = pos;
            while (k < size && S[k] != gap) k += 2;
            if (k < size)
            {
                for (int m = k; m > pos; m -= 2)
                    S[m] = S[m - 2];
                S[pos] = key;
            }
            else
            {
                k = pos - 2;
                while (k >= 0 && S[k] != gap) k -= 2;
                if (k >= 0)
                {
                    for (int m = k; m < pos; m += 2)
                        S[m] = S[m + 2];
                    S[pos - 2] = key;
                }
            }
        }
    }
    int k = 0;
    for (int i = 0; i < size && k < n; i++)
        if (S[i] != gap) A[k++] = S[i];
    delete[] S;
}
