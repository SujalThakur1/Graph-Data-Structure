Graph Algorithms in C
This repository contains the implementation of graph algorithms in C, specifically focusing on Dijkstra's algorithm using both the adjacency matrix and adjacency list representations.

Introduction
This C project aims to demonstrate the implementation of fundamental graph algorithms. The main focus is on Dijkstra's algorithm, a popular algorithm for finding the shortest paths between nodes in a graph. The project supports both the adjacency matrix and adjacency list representations.

File Structure
graphs.h: Header file containing data structures and function declarations.
graphs.c: Main source file with implementations of graph-related functions.
graph_tester.c: A tester program to demonstrate and validate the functionalities.
adjacency_matrix.txt: Sample input file for testing purposes.
Getting Started
Follow the steps below to clone and compile the project:

Open a terminal.

Clone the repository:

```bash
    git clone https://github.com/SujalThakur1/Graph-Data-Structure.git
    cd Graph-Data-Structure
```
Compile the program using the following command:

```bash
    gcc graph_tester.c graphs.c -o graph_algorithm
```
Run the compiled program:
```bash
    ./graph_algorithm
```
Graph Operations
createAdjacencyMatrix: Creates an empty adjacency matrix with a default edge value.
addEdge: Adds a new edge to the adjacency matrix.
addEdges: Adds multiple edges to the adjacency matrix.
doDepthFirstTraversal: Conducts a depth-first traversal of the graph.
loadMatrixFromFile: Loads data from a file into the adjacency matrix.
doDijkstraAlgorithm: Performs Dijkstra's algorithm on the adjacency matrix.
findShortestPathTo: Determines the shortest path between two nodes.
Dijkstra's Algorithm
Two versions of Dijkstra's algorithm are implemented:

doDijkstraAlgorithm: Uses the adjacency matrix representation.
doDijkstraAlgorithmOnAdjacencyList: Uses the adjacency list representation.
