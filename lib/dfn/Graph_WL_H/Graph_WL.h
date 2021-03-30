#pragma once

#include <iostream>
#include <ctime>
#include <cmath>
#include <queue>
#include <vector>

class Graph
{
public:
    size_t NumOfFractures;
    ///< variable with the number of fractures

    //std::vector<std::vector<size_t>> G;

    std::vector<std::vector<size_t>> Adj;
    ///< adjacent list

    std::vector<bool> Visited;
    ///< history: if this fracture has
    // been visited or not

public:
    Graph(const size_t NumOfFractures_1,
          const std::vector<size_t> Connections);
    void DFS(size_t V, std::vector<size_t> &B,
             size_t &Tag_r);
    void CreateGraph_i(std::vector<std::vector<size_t>> &S);
};

inline Graph::Graph(const size_t NumOfFractures_1, const std::vector<size_t> Connections)
{
    //std::cout << "debug_0.01\n";
    NumOfFractures = NumOfFractures_1;

    /*G.resize(NumOfFractures);

    for (size_t i = 0; i < NumOfFractures; ++i)
        G[i].resize(NumOfFractures);

    for (size_t i = 0; i < NumOfFractures; ++i)
    {
        for (size_t j = 0; j < NumOfFractures; ++j)
            G[i][j] = 0;
    }*/
    Adj.resize(NumOfFractures_1);
    //std::cout << "debug_0.02\n";
    for (size_t i = 0; i < NumOfFractures; ++i)
    {
        Adj[i].resize(1);
        Adj[i][0] = i;
    }

    //std::cout << "debug_0.03\n";

    for (size_t i = 0; i < Connections.size() / 2; ++i)
    {
        //G[Connections[2 * i]][Connections[2 * i + 1]] = 1;
        //G[Connections[2 * i + 1]][Connections[2 * i]] = 1;
        Adj[Connections[2 * i]].push_back(Connections[2 * i + 1]);
        Adj[Connections[2 * i + 1]].push_back(Connections[2 * i]);
    }

    Visited.resize(NumOfFractures);
    for (size_t i = 0; i < NumOfFractures; ++i)
    {
        Visited[i] = false;
    }
}

inline void Graph::DFS(size_t V, std::vector<size_t> &B /*one array to record one cluster*/, size_t &Tag_r)
{
    Visited[V] = true;
    B[Tag_r] = V;
    //B.push_back(V);

    /*for (size_t i = 0; i < NumOfFractures; ++i)
    {
        if (Visited[i] == false && G[V][i] == 1)
        {
            Tag_r++;
            DFS(i, B, Tag_r);
        }
    }*/

    for (size_t i = 0; i < Adj[V].size(); i++)
    {
        if (!Visited[Adj[V][i]])
        {
            Tag_r++;
            DFS(Adj[V][i], B, Tag_r);
        }
    }
}

inline void Graph::CreateGraph_i(std::vector<std::vector<size_t>> &S)
{
    //std::cout << "graph_x\n";
    if (S.size() != 0)
    {
        std::cout << "Error! ListOfClusters should be empty array!\n";
        exit(0);
    }
    for (size_t i = 0; i < NumOfFractures; ++i)
    {
        if (!Visited[i])
        {
            size_t Tag_r = 0;
            //std::vector<size_t> A;
            std::vector<size_t> A(NumOfFractures);
            //std::cout << "debug_0.5\n";
            DFS(i, A, Tag_r); // so, at least, the first element of A will form a cluster
            //std::cout << "debug_1\n";
            //need a function to determine true size of one cluster, i.e., provided A = [0,1,2,3,0,0,0], the size of this cluster is four

            std::vector<size_t> new_A(Tag_r + 1);

            S.push_back(new_A);
            //std::cout << "debug_2\n";
            for (size_t i = 0; i < Tag_r + 1; ++i)
                S[S.size() - 1][i] = A[i];
            //std::cout << "debug_3\n";
            //S.push_back(A);
            //S.push_back(new_A); //
        }
    }
}
