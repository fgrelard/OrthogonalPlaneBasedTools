#ifndef SHORTEST_PATH_ALGORITHM_H
#define SHORTEST_PATH_ALGORITHM_H

#include <set>
#include <vector>
#include <map>
#include "geometry/path/LinkPointAlgorithm.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/graph/CGraphVisitor.h"

template<typename Point, typename Container>
class ShortestPathAlgorithm : public LinkPointAlgorithm<Point> {
public:
    typedef DGtal::PointVector<Point::dimension, double> RealVector;
public:
    ShortestPathAlgorithm() {}

    ShortestPathAlgorithm(const Point &source,
                          const Point &destination,
                          const Container &aSet) : LinkPointAlgorithm<Point>(source, destination),
                                                   mySet(aSet) {}

    ShortestPathAlgorithm(const ShortestPathAlgorithm &other) : LinkPointAlgorithm<Point>(other),
                                                                mySet(other.mySet) {}

    virtual typename LinkPointAlgorithm<Point>::Path linkPoints() = 0;

protected:
    template<typename Visitor>
    std::vector<Point> shortestPath(Visitor &visitor);

    std::vector<Point> reconstructPath(const std::map<Point, Point> &aMapPrevious);


protected:
    Container mySet;
};

template<typename Point, typename Container>
template<typename Visitor>
std::vector<Point>
ShortestPathAlgorithm<Point, Container>::shortestPath(Visitor &visitor) {

    BOOST_CONCEPT_ASSERT((DGtal::concepts::CGraphVisitor<Visitor>));

    typedef typename Visitor::Node MyNode;

    MyNode node;
    std::vector<Point> thePath;
    std::map<Point, Point> aMapPrevious;
    while (!visitor.finished()) {
        node = visitor.current();
        if (node.first == this->myDestination) {
            return reconstructPath(aMapPrevious);
        }
        std::vector<Point> neighbors;
        std::back_insert_iterator<std::vector<Point>> iter(neighbors);
        visitor.graph().writeNeighbors(iter, node.first);
        for (auto it = neighbors.begin(), ite = neighbors.end(); it != ite; ++it) {
            auto itMapExisting = aMapPrevious.find(*it);
            if (itMapExisting == aMapPrevious.end()) {
                aMapPrevious[*it] = node.first;
            }
        }
        visitor.expand();
    }
    return thePath; //path not found = empty vector
}

template<typename Point, typename Container>
std::vector<Point>
ShortestPathAlgorithm<Point, Container>::reconstructPath(const std::map<Point, Point> &aMapPrevious) {
    std::vector<Point> path;
    Point aPoint = this->myDestination;

    while (aPoint != this->mySource) {
        path.push_back(aPoint);
        aPoint = aMapPrevious.at(aPoint);
    }


    path.push_back(this->mySource);
    return path;
}

#endif
