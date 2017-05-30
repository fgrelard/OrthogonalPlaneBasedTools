#ifndef DIJKSTRA_ALGORITHM_H
#define DIJKSTRA_ALGORITHM_H

#include <set>
#include <vector>
#include "geometry/path/ShortestPathAlgorithm.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/graph/BreadthFirstVisitor.h"
#include "DGtal/kernel/sets/CDigitalSet.h"
#include "DGtal/topology/DigitalTopology.h"
#include "DGtal/topology/Object.h"

template<typename Point, typename Container>
class DijkstraAlgorithm : public ShortestPathAlgorithm<Point, Container> {

    BOOST_CONCEPT_ASSERT((DGtal::concepts::CDigitalSet<Container>));

public:
    typedef DGtal::PointVector<Point::dimension, double> RealVector;
public:
    DijkstraAlgorithm() {}

    DijkstraAlgorithm(const Point &source,
                      const Point &destination,
                      const Container &aSet) : ShortestPathAlgorithm<Point, Container>(source, destination, aSet) {}

    DijkstraAlgorithm(const DijkstraAlgorithm &other) : ShortestPathAlgorithm<Point, Container>(other) {}

    virtual typename LinkPointAlgorithm<Point>::Path linkPoints();

};


template<typename Point, typename Container>
typename LinkPointAlgorithm<Point>::Path
DijkstraAlgorithm<Point, Container>::linkPoints() {
    typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
    typedef DGtal::MetricAdjacency<Space, 1> Adj6;
    typedef DGtal::MetricAdjacency<Space, 3> Adj26;
    typedef DGtal::DigitalTopology<Adj26, Adj6> DT26_6;
    typedef DGtal::Object<DT26_6, Container> Graph;
    typedef DGtal::BreadthFirstVisitor<Graph, std::set<Point>> Visitor;

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    Graph graph(dt26_6, this->mySet);
    Visitor visitor(graph, this->mySource);
    std::vector<Point> points = ShortestPathAlgorithm<Point, Container>::shortestPath(visitor);
    typename LinkPointAlgorithm<Point>::Path path(points);
    return path;
}

#endif
