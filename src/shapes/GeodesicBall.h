#ifndef GEODESIC_BALL_H
#define GEODESIC_BALL_H

#include <vector>
#include <set>
#include "DGtal/base/Common.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "geometry/DistanceToPointFunctor.h"

template<typename Distance, typename Point>
class GeodesicBall {
public:
    GeodesicBall() : myCenter(), myRadius(0) {}

    GeodesicBall(const Distance &distance,
                 const Point &center,
                 double radius) : myDistance(distance),
                                  myCenter(center),
                                  myRadius(radius) {}

    template<typename Graph>
    std::vector<Point> surfaceIntersection(const Graph &setSurface);

private:
    Distance myDistance;
    Point myCenter;
    double myRadius;
};

template<typename Distance, typename Point>
template<typename Graph>
std::vector<Point> GeodesicBall<Distance, Point>::surfaceIntersection(const Graph &setSurface) {
    typedef DistanceToPointFunctor<Distance> DistanceFunctor;
    typedef DGtal::DistanceBreadthFirstVisitor<Graph, DistanceFunctor, std::set<Point>> Visitor;
    typedef typename Visitor::Node MyNode;

    DistanceFunctor distance(myDistance, myCenter);
    Visitor visitor(setSurface, distance, myCenter);
    MyNode node;
    std::vector<Point> intersection;
    while (!visitor.finished()) {
        node = visitor.current();
        if (node.second >= myRadius) break;
        intersection.push_back(node.first);
        visitor.expand();
    }
    return intersection;
}

#endif
