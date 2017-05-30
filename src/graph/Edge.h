#ifndef EDGE_H
#define EDGE_H


#include <vector>
#include <algorithm>
#include "geometry/Distance.h"
#include "geometry/SetProcessor.h"
#include "DGtal/topology/DomainMetricAdjacency.h"

template<typename Container>
class Edge : public Container {
public:
    typedef typename Container::value_type Point;
public:
    Edge(const Container &aCurve) : Container(aCurve) {}

    Edge(const Edge &other) : Container(other) {}

public:
    std::vector<Edge<Container> *> neighboringEdges(const std::vector<Edge<Container> *> &edges,
                                                    const Container &branchingPoints);

};

template<typename Container>
std::vector<Edge<Container> *> Edge<Container>::neighboringEdges(const std::vector<Edge<Container> *> &edges,
                                                                 const Container &branchingPoints) {
    typedef typename Container::value_type Point;
    typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
    typedef DGtal::MetricAdjacency<Space, 3> MetricAdjacency;

    std::vector<Edge<Container> *> neighbors;
    Point branchPoint;
    for (const Point &b : branchingPoints) {
        if (std::find(this->begin(), this->end(), b) != this->end()) {
            branchPoint = b;
        }
    }

    std::vector<Point> nb;
    std::back_insert_iterator<std::vector<Point> > inserter(nb);
    MetricAdjacency::writeNeighbors(inserter, branchPoint);

    for (Edge<Container> *edge : edges) {
        Container setEdge = *edge;
        SetProcessor<Container> setProcessor(setEdge);
        if (setProcessor.sameContainer((Container) (*this))) continue;
        for (const Point &n : nb) {
            if (find(setEdge.begin(), setEdge.end(), n) != setEdge.end())
                neighbors.push_back(edge);
        }
    }

    return neighbors;
}

#endif
