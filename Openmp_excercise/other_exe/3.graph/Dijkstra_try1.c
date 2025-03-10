#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>

// Structure for storing a graph
struct Graph
{
    int vertexNum;
    int **edges;
};

// Constructs a graph with V vertices and E edges
void createGraph(struct Graph *G, int V)
{
    G->vertexNum = V;
    G->edges = (int **)malloc(V * sizeof(int *));
    for (int i = 0; i < V; i++)
    {
        G->edges[i] = (int *)malloc(V * sizeof(int));
        for (int j = 0; j < V; j++)
            G->edges[i][j] = INT_MAX;
        G->edges[i][i] = 0;
    }
}

// Adds the given edge to the graph
void addEdge(struct Graph *G, int src, int dst, int weight)
{
    G->edges[src][dst] = weight;
}

// Utility function to find minimum distance vertex in mdist
int minDistance(int mdist[], int vset[], int V)
{
    int minVal = INT_MAX;
    static int minInd = -1; // remembers the previous value if not modified in the loop
    for (int i = 0; i < V; i++)
        if (vset[i] == 0 && mdist[i] < minVal)
        {
            minVal = mdist[i];
            minInd = i;
        }

    return minInd;
}

// Utility function to print distances
void print(int dist[], int V)
{
    printf("\nVertex  Distance\n");
    for (int i = 0; i < V; i++)
    {
        if (dist[i] != INT_MAX)
            printf("%d\t%d\n", i, dist[i]);
        else
            printf("%d\tINF", i);
    }
}

// The main function that finds the shortest path from given source
// to all other vertices using Dijkstra's Algorithm.It doesn't work on negative
// weights
void Dijkstra(struct Graph *graph, int src)
{
    int V = graph->vertexNum;
    int mdist[V]; // Stores updated distances to vertex
    int vset[V];  // vset[i] is true if the vertex i included
                  // in the shortest path tree

    // Initialise mdist and vset. Set distance of source as zero
    for (int i = 0; i < V; i++)
        mdist[i] = INT_MAX, vset[i] = 0;

    mdist[src] = 0;

// iterate to find shortest path
#pragma omp parallel
    {
#pragma omp single nowait
        {
            for (int count = 0; count < V - 1; count++)
            {
#pragma omp task
                {
                    int u = -1;
                    int min_dist = INT_MAX;
                    for (int i = 0; i < V; i++)
                    {
                        if (!vset[i] && mdist[i] < min_dist)
                        {
                            {
                                if (mdist[i] < min_dist)
                                {
                                    min_dist = mdist[i];
                                    u = i;
                                }
                            }
                        }
                    }

                    // if (u == -1) break;  // No reachable vertices left

                    // Relaxation step: Update distances of adjacent vertices (parallel part)
                    if (u != -1)
                    {
                        for (int v = 0; v < V; v++)
                        {
                            if (!vset[v] && graph->edges[u][v] != INT_MAX &&
                                mdist[u] + graph->edges[u][v] < mdist[v])
                            {
                                mdist[v] = mdist[u] + graph->edges[u][v];
                            }
                        }
                        vset[u] = 1;
                    }
                }
            }
        }
    }
    // print the shortest distances
    print(mdist, V);

    return;
}

// Driver Function
int main()
{
    int V = 10000;
    int E = 50000;
    int src, dst, weight;
    int gsrc = 0;

    struct Graph G;
    createGraph(&G, V);

    srand(42);

    for (int i = 0; i < E; i++)
    {
        src = rand() % V;
        dst = rand() % V;
        weight = rand() % 100 + 1;
        addEdge(&G, src, dst, weight);
    }

    double start_time = omp_get_wtime();
    Dijkstra(&G, gsrc);
    double end_time = omp_get_wtime();

    printf("\nTime taken for Dijkstra's algorithm: %f seconds\n", end_time - start_time);

    return 0;
}
