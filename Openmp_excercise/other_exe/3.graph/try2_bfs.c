#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <string.h>

#define SIZE 100000

struct node {
    int vertex;
    struct node *next;
};

struct node *createNode(int);

struct Graph {
    int numVertices;
    struct node **adjLists;
    int *visited;
};

struct Graph *createGraph(int vertices);
void addEdge(struct Graph *graph, int src, int dest);
void bfs(struct Graph *graph, int startVertex);

int main() {
    int vertices = 100000; // Number of vertices
    int edges = 500000;    // Number of edges

    // Create graph
    struct Graph *graph = createGraph(vertices);

    // Seed for reproducibility
    srand(42);

    // Generate random edges
    for (int i = 0; i < edges; i++) {
        int src = rand() % vertices;
        int dest = rand() % vertices;
        addEdge(graph, src, dest);
    }

    // Perform BFS from vertex 0
    double start_time = omp_get_wtime();
    bfs(graph, 0);
    double end_time = omp_get_wtime();

    printf("\nTime taken: %f seconds\n", end_time - start_time);

    return 0;
}

void bfs(struct Graph *graph, int startVertex) {
    int numVertices = graph->numVertices;
    int numThreads;

    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }

    // Thread-local queues
    int **localQueues = (int **)malloc(numThreads * sizeof(int *));
    int *localQueueSizes = (int *)calloc(numThreads, sizeof(int));

    for (int i = 0; i < numThreads; i++) {
        localQueues[i] = (int *)malloc(numVertices * sizeof(int));
    }

    // Initialize visited array and add start vertex
    graph->visited[startVertex] = 1;
    localQueues[0][0] = startVertex;
    localQueueSizes[0] = 1;

    // Parallel BFS loop
    int active = 1; // Flag to indicate active threads
    while (active) {
        active = 0;

        #pragma omp parallel
        {
            int threadId = omp_get_thread_num();
            int *queue = localQueues[threadId];
            int queueSize = localQueueSizes[threadId];
            int newQueueSize = 0;

            int *newQueue = (int *)malloc(numVertices * sizeof(int));

            for (int i = 0; i < queueSize; i++) {
                int currentVertex = queue[i];
                struct node *temp = graph->adjLists[currentVertex];

                while (temp) {
                    int adjVertex = temp->vertex;

                    if (__sync_bool_compare_and_swap(&graph->visited[adjVertex], 0, 1)) {
                        newQueue[newQueueSize++] = adjVertex;
                    }
                    temp = temp->next;
                }
            }

            // Update local queues
            localQueueSizes[threadId] = newQueueSize;
            memcpy(localQueues[threadId], newQueue, newQueueSize * sizeof(int));
            free(newQueue);

            #pragma omp atomic
            active += (newQueueSize > 0);
        }
    }

    // Free resources
    for (int i = 0; i < numThreads; i++) {
        free(localQueues[i]);
    }
    free(localQueues);
    free(localQueueSizes);
}

struct node *createNode(int v) {
    struct node *newNode = malloc(sizeof(struct node));
    newNode->vertex = v;
    newNode->next = NULL;
    return newNode;
}

struct Graph *createGraph(int vertices) {
    struct Graph *graph = malloc(sizeof(struct Graph));
    graph->numVertices = vertices;
    graph->adjLists = malloc(vertices * sizeof(struct node *));
    graph->visited = malloc(vertices * sizeof(int));

    for (int i = 0; i < vertices; i++) {
        graph->adjLists[i] = NULL;
        graph->visited[i] = 0;
    }

    return graph;
}

void addEdge(struct Graph *graph, int src, int dest) {
    struct node *newNode = createNode(dest);
    newNode->next = graph->adjLists[src];
    graph->adjLists[src] = newNode;

    newNode = createNode(src);
    newNode->next = graph->adjLists[dest];
    graph->adjLists[dest] = newNode;
}
