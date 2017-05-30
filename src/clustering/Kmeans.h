#ifndef KMEANS_H
#define KMEANS_H

#include <vector>
#include <limits>
#include "DGtal/math/Statistic.h"

template<typename Container>
class KMeans {

public:
    std::vector<Container> kmeansAlgorithm(int k);

private:
    std::vector<Container> initialize2Clusters();

    void assignClasses(std::vector<Container> &clusters, const Container &centroids);

    void recomputeCentroids(Container &centroids, const std::vector<Container> &clusters);

private:
    Container myContainer;
};

template<typename Container>
std::vector<Container> KMeans<Container>::kmeansAlgorithm(int k) {
    std::vector<Container> clusters = initialize2Clusters();

    Container centroids;
    for (const Container &point : clusters) {
        centroids.push_back(point[0]);
    }

    int cpt = 0;
    bool convergence = false;
    while (!convergence && cpt < 1000) {
        Container previousCentroids = centroids;

        assignClasses(clusters, centroids);
        recomputeCentroids(centroids, clusters);
        cpt++;

        if (centroids == previousCentroids)
            convergence = true;
    }
    return clusters;
}

template<typename Container>
std::vector<Container> KMeans<Container>::initialize2Clusters() {
    typedef typename Container::value_type Scalar;
    std::vector<Container> clusters;
    Scalar minVal = *std::min_element(myContainer.begin(), myContainer.end());
    Scalar maxVal = *std::max_element(myContainer.begin(), myContainer.end());
    clusters.push_back({minVal});
    clusters.push_back({maxVal});
    return clusters;
}

template<typename Container>
void KMeans<Container>::assignClasses(std::vector<Container> &clusters, const Container &centroids) {
    typedef typename Container::value_type Scalar;

    int k = clusters.size();
    for (int i = 0; i < k; i++) {
        clusters[i].clear();
    }

    for (const Scalar &data : myContainer) {
        Scalar distance = std::numeric_limits<Scalar>::max();
        int index = 0;
        for (int i = 0; i < k; i++) {
            Scalar currentDistance = sqrt(pow((data - centroids[i]), 2));
            if (currentDistance < distance) {
                distance = currentDistance;
                index = i;
            }
        }
        clusters[index].push_back(data);
    }
}

template<typename Container>
void KMeans<Container>::recomputeCentroids(Container &centroids, const std::vector<Container> &clusters) {
    typedef typename Container::value_type Scalar;
    for (int i = 0, end = clusters.size(); i < end; i++) {
        DGtal::Statistic<Scalar> stats;
        stats.addValues(clusters[i].begin(), clusters[i].end());
        double mean = stats.mean();
        centroids[i] = mean;
    }
}

#endif
