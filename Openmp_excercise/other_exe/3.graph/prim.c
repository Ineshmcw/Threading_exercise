/**
 * @file
 * @author [Timothy Maloney](https://github.com/sl1mb0)
 * @brief [Prim's algorithm](https://en.wikipedia.org/wiki/Prim%27s_algorithm)
 * implementation in C to find the MST of a weighted, connected graph.
 * @details Prim's algorithm uses a greedy approach to generate the MST of a weighted connected graph.
 * The algorithm begins at an arbitrary vertex v, and selects a next vertex u, 
 * where v and u are connected by a weighted edge whose weight is the minimum of all edges connected to v. 
 * @references Page 319 "Introduction to the Design and Analysis of Algorithms" - Anany Levitin
 *
 * To test - run './prim -test'
 * prim() will find the MST of the following adj. matrix:
 *	  
 *	  0  1  2  3
 *        1  0  4  6
 *        2  4  0  5
 *        3  6  5  0
 * 
 * The minimum spanning tree for the above weighted connected graph is given by the following adj matrix:
 *	   
 *	  0  1  2  3
 *	  1  0  0  0
 *	  2  0  0  0
 *	  3  0  0  0
 *
 *
 * The following [link](https://visualgo.net/en/mst) provides a visual representation of graphs that can be used to test/verify the algorithm for different adj
 * matrices and their weighted, connected graphs.
 */

#include <stdio.h>        /// for IO operations
#include <string.h>      /// for string comparison
#include <assert.h>     /// for assert()
#include <inttypes.h>  /// for uint16_t
#include <omp.h>        /// for OpenMP
#include <stdlib.h>     /// for malloc

#define MAX 1000
#define INF 999

/**
 * @brief Finds index of minimum element in edge list for an arbitrary vertex
 * @param arr graph row
 * @param N number of elements in arr
 * @returns index of minimum element in arr
 */
uint16_t minimum(uint16_t arr[], uint16_t N)
{
    uint16_t index = 0;
    uint16_t min = INF;

    for (uint16_t i = 0; i < N; i++)
    {
        if (arr[i] < min)
        {
            min = arr[i];
            index = i;
        }
    }
    return index;
}

/**
 * @brief Used to find MST of user-generated adj matrix G
 * @returns void
 */
void prim(uint16_t G[][MAX], uint16_t MST[][MAX], uint16_t V)
{
    uint16_t u, v;
    uint16_t E_t[MAX], path[MAX];
    uint16_t V_t[MAX], no_of_edges;

    E_t[0] = 0;  // edges for current vertex
    V_t[0] = 1;  // list of visited vertices

    for (uint16_t i = 1; i < V; i++)
    {
        E_t[i] = G[i][0];
        path[i] = 0;
        V_t[i] = 0;
    }

    no_of_edges = V - 1;

    while (no_of_edges > 0)
    {
        u = minimum(E_t, V);
        while (V_t[u] == 1)
        {
            E_t[u] = INF;
            u = minimum(E_t, V);
        }

        v = path[u];
        MST[v][u] = E_t[u];
        MST[u][v] = E_t[u];
        no_of_edges--;
        V_t[u] = 1;

        for (uint16_t i = 1; i < V; i++)
        {
            if (V_t[i] == 0 && G[u][i] < E_t[i])
            {
                E_t[i] = G[u][i];
                path[i] = v;
            }
        }
    }
}

/**
 * @brief Self-test implementations
 * @returns void
 */
static void test(uint16_t G[][MAX], uint16_t MST[][MAX], uint16_t V) {
    V = MAX;
    for (uint16_t i = 0; i < V; ++i) {
        for (uint16_t j = 0; j < V; ++j) {
            G[i][j] = (i == j) ? 0 : (rand() % 100 + 1);
        }
    }

    for (uint16_t i = 0; i < MAX; i++) {
        for (uint16_t j = 0; j < MAX; j++) {
            MST[i][j] = 0;
        }
    }

    double start_time = omp_get_wtime();
    prim(G, MST, V);
    double end_time = omp_get_wtime();

    printf("\nTime taken for Prim's algorithm: %f seconds\n", end_time - start_time);
}

int main(int argc, char const *argv[]) {   
    uint16_t G[MAX][MAX];
    uint16_t MST[MAX][MAX];
    uint16_t V = MAX;

    if (argc == 2 && strcmp(argv[1], "-test") == 0) {
        test(G, MST, V);
    } else {
        V = MAX;
        for (uint16_t i = 0; i < V; ++i) {
            for (uint16_t j = 0; j < V; ++j) {
                G[i][j] = (i == j) ? 0 : (rand() % 100 + 1);
            }
        }

        for (uint16_t i = 0; i < MAX; i++) {
            for (uint16_t j = 0; j < MAX; j++) {
                MST[i][j] = 0;
            }
        }

        double start_time = omp_get_wtime();
        prim(G, MST, V);
        double end_time = omp_get_wtime();

        printf("\nTime taken for Prim's algorithm: %f seconds\n", end_time - start_time);
    }

    return 0;
}