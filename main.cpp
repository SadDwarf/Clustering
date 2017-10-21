#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <exception>
#include <climits>


using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::vector;


class DisjointSetUnion {
private:
    vector<int> parent_;
    vector<int> ranks_;

public:
    explicit DisjointSetUnion(size_t size)
        : parent_()
        , ranks_(size, 0)
    {
        parent_.reserve(size);
        for (size_t i = 0; i < size; ++i) {
            parent_.push_back(i);
        }
    }

    int find(int node) {
        if (parent_[node] != node) {
            parent_[node] = find(parent_[node]);
        }
        return parent_[node];
    }

    void union_sets(int first, int second) {
        int first_root = find(first);
        int second_root = find(second);
        if (first_root == second_root) {
            return;
        }

        if (ranks_[first_root] < ranks_[second_root]) {
            parent_[first_root] = second_root;
        } else if (ranks_[first_root] > ranks_[second_root]) {
            parent_[second_root] = first_root;
        } else {
            parent_[second_root] = first_root;
            ++ranks_[first_root];
        }
    }
};


struct Edge {
    size_t from;
    size_t to;
    double weight;
};


// Map arbitrary labels to 0, ..., (n-1) labels
vector<size_t> RenumerateLabels(const vector<size_t>& rawLabels) {
    vector<int> rawToNew(rawLabels.size(), -1);
    size_t indexesUsed = 0;
    vector<size_t> newLabels(rawLabels.size());
    for (size_t i = 0; i < rawLabels.size(); ++i) {
        size_t oldLabel = rawLabels[i];
        if (rawToNew[oldLabel] == -1) {
            rawToNew[oldLabel] = indexesUsed;
            ++indexesUsed;
        }
        newLabels[i] = rawToNew[oldLabel];
    }
    return newLabels;
}

bool comparator(Edge a, Edge b) {
    return a.weight < b.weight;
}

vector<size_t> ClusterGraphMST(vector<Edge> edges,  // copy for sorting
                               size_t vertexCount,
                               size_t clusterCount) {
    if (clusterCount > vertexCount)
        throw std::out_of_range("required more clusters than vertexes in MST clustering");
    vector<size_t> ans(vertexCount);
    std::sort(edges.begin(), edges.end(), comparator);
    DisjointSetUnion dsu(vertexCount);
    size_t connectedComponnents = 0;
    for (int i = 0; i < edges.size(); ++i) {
        if ((vertexCount - connectedComponnents) == clusterCount)
            break;
        if (dsu.find(edges[i].from) != dsu.find(edges[i].to)) {
            dsu.union_sets(edges[i].from, edges[i].to);
            connectedComponnents++;
        }
    }
    for (int i = 0; i < vertexCount; ++i) {
        ans[i] = dsu.find(i);
    }
    ans = RenumerateLabels(ans);
    return ans;
}



template <typename T, typename Dist>
vector<Edge> PairwiseDistances(vector<T> objects, Dist distance) {
    vector<Edge> edges;
    for (size_t i = 0; i < objects.size(); ++i) {
        for (size_t j = i + 1; j < objects.size(); ++j) {
            edges.push_back({i, j, distance(objects[i], objects[j])});
        }
    }
    return edges;
}


template <typename T, typename Dist>
vector<size_t> ClusterMST(const vector<T>& objects, Dist distance, size_t clusterCount) {
    vector<Edge> edges = PairwiseDistances(objects, distance);
    return ClusterGraphMST(edges, objects.size(), clusterCount);
}

template <typename T, typename Dist>
size_t choseFurthest(const vector<T>& objects,
                     const vector<size_t>& centers,
                     Dist distance) {
    double maxDist = -1, minDist;
    size_t localBestVert = 0;
    for (int i = 0; i < objects.size(); ++i) {
        minDist = -1;
        for (int j = 0; j < centers.size(); ++j) {
            double temp = distance(objects[i], objects[centers[j]]);
            if (minDist == -1 || temp < minDist)
                minDist = temp;
        }
        if (minDist > maxDist) {
            localBestVert = i;
            maxDist = minDist;
        }
    }
    return localBestVert;
}

template <typename T, typename Dist>
vector<size_t> findNearestCentr(const vector<T>& objects,
                     const vector<size_t>& centers,
                     Dist distance) {
    double minDist = -1;
    size_t bestVert;
    std::vector<size_t> ans(objects.size());
    for (int i = 0; i < objects.size(); ++i) {
        minDist = -1;
        for (int j = 0; j < centers.size(); ++j) {
            double temp = distance(objects[i], objects[centers[j]]);
            if (minDist == -1 || temp < minDist) {
                minDist = temp;
                bestVert = j;
            }
        }
        ans[i] = bestVert;
    }
    return ans;
}

template <typename T, typename Dist>
vector<size_t> ClusterMinDistToCenter(const vector<T>& objects,
                                      Dist distance,
                                      size_t clusterCount)
{
    if (clusterCount > objects.size())
        throw std::out_of_range("required more clusters than vertexes in MDC clustering");
    auto baseGenerator = std::default_random_engine();
    auto generateFirstVert = std::uniform_int_distribution<size_t>
                                                           (0,objects.size() - 1);
    vector<size_t> centers;
    vector<size_t> ans(objects.size());
    centers.push_back(generateFirstVert(baseGenerator));
    for (int i = 0; i < clusterCount; ++i) {
        centers.push_back(choseFurthest(objects, centers, distance));
    }
    return findNearestCentr(objects, centers, distance);
}


struct Point2D {
    double x, y;
};


double EuclidianDistance(const Point2D& first, const Point2D& second) {
    return std::sqrt((first.x - second.x) * (first.x - second.x) +
                     (first.y - second.y) * (first.y - second.y));
}


vector<Point2D> Random2DClusters(const vector<Point2D>& centers,
                                 const vector<double>& xVariances,
                                 const vector<double>& yVariances,
                                 size_t pointsCount)
{
    auto baseGenerator = std::default_random_engine();
    auto generateCluster = std::uniform_int_distribution<size_t>(0, centers.size() - 1);
    auto generateDeviation = std::normal_distribution<double>();

    vector<Point2D> results;
    for (size_t i = 0; i < pointsCount; ++i) {
        size_t c = generateCluster(baseGenerator);
        double x = centers[c].x + generateDeviation(baseGenerator) * xVariances[c];
        double y = centers[c].y + generateDeviation(baseGenerator) * yVariances[c];
        results.push_back({x, y});
    }

    return results;
}


// Generate files for plotting in gnuplot
void GNUPlotClusters2D(const vector<Point2D>& points,
                       const vector<size_t>& labels,
                       size_t clustersCount,
                       const string& outFolder)
{
    std::ofstream scriptOut(outFolder + "/script.txt");
    scriptOut << "plot ";

    for (size_t cluster = 0; cluster < clustersCount; ++cluster) {
        string filename = std::to_string(cluster) + ".dat";
        std::ofstream fileOut(outFolder + "/" + filename);
        scriptOut << "\"" << filename << "\"" << " with points, ";

        for (size_t i = 0; i < points.size(); ++i) {
            if (labels[i] == cluster) {
                fileOut << points[i].x << "\t" << points[i].y << "\n";
            }
        }
    }
}


int main() {
    auto points = Random2DClusters(
        {{1, 0}, {0, 0}, {2, 0}},
        {0.2, 0.2, 0.2},
        {0.2, 0.4, 0.4},
        1000);

    vector<size_t> labels(points.size(), 0);
    GNUPlotClusters2D(points, labels, 1, "/home/egor/CLionProjects/clustering/plot_base");

    size_t clustersCount = 3;

    labels = ClusterMST(points, EuclidianDistance, clustersCount);
    GNUPlotClusters2D(points, labels, clustersCount, "/home/egor/CLionProjects/clustering/plot_mst");

    labels = ClusterMinDistToCenter(points, EuclidianDistance, clustersCount);
    GNUPlotClusters2D(points, labels, clustersCount, "/home/egor/CLionProjects/clustering/plot_mdc");


    return 0;
}