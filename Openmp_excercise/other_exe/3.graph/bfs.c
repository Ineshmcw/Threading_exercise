#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define SIZE 100000 // Increased size for large graph

struct queue
{
    int items[SIZE];
    int front;
    int rear;
};

struct queue *createQueue();
void enqueue(struct queue *q, int);
int dequeue(struct queue *q);
int isEmpty(struct queue *q);
int pollQueue(struct queue *q);

struct node
{
    int vertex;
    struct node *next;
};

struct node *createNode(int);

struct Graph
{
    int numVertices;
    struct node **adjLists;
    int *visited;
};

struct Graph *createGraph(int vertices);
void addEdge(struct Graph *graph, int src, int dest);
void bfs(struct Graph *graph, int startVertex);

int main()
{
    int vertices = 100000; // Number of vertices
    int edges = 500000;    // Number of edges

    // Create graph
    struct Graph *graph = createGraph(vertices);

    // Seed for reproducibility
    srand(42);

    // Generate random edges
    for (int i = 0; i < edges; i++)
    {
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

void bfs(struct Graph *graph, int startVertex)
{
    struct queue *q = createQueue();

    graph->visited[startVertex] = 1;
    enqueue(q, startVertex);
    // printf("Breadth First Traversal from vertex %d:\n", startVertex);

    while (!isEmpty(q))
    {
        int currentVertex = dequeue(q);
        // printf("%d ", currentVertex);

        struct node *temp = graph->adjLists[currentVertex];
        while (temp)
        {
            int adjVertex = temp->vertex;
            if (graph->visited[adjVertex] == 0)
            {
                graph->visited[adjVertex] = 1;
                enqueue(q, adjVertex);
            }
            temp = temp->next;
        }
    }
}

struct node *createNode(int v)
{
    struct node *newNode = malloc(sizeof(struct node));
    newNode->vertex = v;
    newNode->next = NULL;
    return newNode;
}

struct Graph *createGraph(int vertices)
{
    struct Graph *graph = malloc(sizeof(struct Graph));
    graph->numVertices = vertices;
    graph->adjLists = malloc(vertices * sizeof(struct node *));
    graph->visited = malloc(vertices * sizeof(int));

    for (int i = 0; i < vertices; i++)
    {
        graph->adjLists[i] = NULL;
        graph->visited[i] = 0;
    }

    return graph;
}

void addEdge(struct Graph *graph, int src, int dest)
{
    struct node *newNode = createNode(dest);
    newNode->next = graph->adjLists[src];
    graph->adjLists[src] = newNode;

    newNode = createNode(src);
    newNode->next = graph->adjLists[dest];
    graph->adjLists[dest] = newNode;
}

struct queue *createQueue()
{
    struct queue *q = malloc(sizeof(struct queue));
    q->front = -1;
    q->rear = -1;
    return q;
}

int isEmpty(struct queue *q)
{
    return q->rear == -1;
}

void enqueue(struct queue *q, int value)
{
    if (q->rear == SIZE - 1)
    {
        printf("\nQueue is Full!!\n");
        return;
    }
    if (q->front == -1)
        q->front = 0;
    q->rear++;
    q->items[q->rear] = value;
}

int dequeue(struct queue *q)
{
    if (isEmpty(q))
    {
        printf("Queue is empty\n");
        return -1;
    }
    int item = q->items[q->front];
    q->front++;
    if (q->front > q->rear)
        q->front = q->rear = -1;
    return item;
}

int pollQueue(struct queue *q)
{
    return q->items[q->front];
}
