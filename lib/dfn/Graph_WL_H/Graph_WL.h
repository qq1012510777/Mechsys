#pragma once

#include <bits/stdc++.h>
#include <cmath>
#include <ctime>
#include <iostream>
#include <queue>
#include <vector>
using namespace std;
//--------------------------------------------------------------------------
class Graph
{
public:
    //size_t NumOfFractures;
    ///< variable with the number of fractures
    size_t V;
    //std::vector<std::vector<size_t>> G;

    //std::vector<std::vector<size_t>> Adj;
    list<int> *Adj; // adjacency lists
    ///< adjacent list

    //std::vector<bool> Visited;
    ///< history: if this fracture has
    // been visited or not

public:
    Graph(const size_t NumOfFractures_1,
          const std::vector<size_t> Connections);

    void DFS(std::vector<std::vector<size_t>> & ListOfClusters);

    void DFSUtil(int s, vector<bool> &visited, vector<size_t> &onecluster);

    void addEdge(int v, int w);

    void CreateGraph_i(std::vector<std::vector<size_t>> &S);
};

inline Graph::Graph(const size_t NumOfFractures_1, const std::vector<size_t> Connections)
{
    //std::cout << "debug_0.01\n";
    Adj = new list<int>[NumOfFractures_1];
    V = NumOfFractures_1;

    //std::cout << "debug_0.03\n";

    for (size_t i = 0; i < Connections.size() / 2; ++i)
    {
        addEdge(Connections[2 * i], Connections[2 * i + 1]);
        addEdge(Connections[2 * i + 1], Connections[2 * i]);
    }
}

inline void Graph::addEdge(int v, int w)
{
    Adj[v].push_back(w); // Add w to v’s list.
}

inline void Graph::DFS(std::vector<std::vector<size_t>> & ListOfClusters)
{
    // Mark all the vertices as not visited
    vector<bool> visited(V, false);

    for (int i = 0; i < (int)V; i++)
    {
        vector<size_t> onecluster;
        if (!visited[i])
        {
            DFSUtil(i, visited, onecluster);
        }

        if (onecluster.size() > 0)
        {
            ListOfClusters.push_back(onecluster);   
        }
    }
}

inline void Graph::DFSUtil(int s, vector<bool> &visited, vector<size_t> &onecluster)
{
    // Create a stack for DFS
    stack<int> stack;

    // Push the current source node.
    stack.push(s);

    while (!stack.empty())
    {
        // Pop a vertex from stack and print it
        s = stack.top();
        stack.pop();

        // Stack may contain same vertex twice. So
        // we need to print the popped item only
        // if it is not visited.
        if (!visited[s])
        {
            //cout << s << " ";
            onecluster.push_back((size_t)s);
            visited[s] = true;
        }

        // Get all adjacent vertices of the popped vertex s
        // If a adjacent has not been visited, then push it
        // to the stack.
        for (auto i = Adj[s].begin(); i != Adj[s].end(); ++i)
        {
            if (!visited[*i])
            {
                stack.push(*i);
            };
        }
    }
}

inline void Graph::CreateGraph_i(std::vector<std::vector<size_t>> &S)
{
    DFS(S);   
}
