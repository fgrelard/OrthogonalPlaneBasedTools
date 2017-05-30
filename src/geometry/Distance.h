#ifndef DISTANCE_H
#define DISTANCE_H

#include <math.h>
#include "DGtal/base/Common.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/topology/DigitalTopology.h"
#include "DGtal/topology/Object.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"

namespace Distance {
    inline double euclideanDistance(float x1, float y1, float x2, float y2) {
        return sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
    }

    template<typename Point>
    inline double euclideanDistance(Point p, Point other) {
        int dim = Point::dimension;
        double sum = 0;
        for (int i = 0; i < dim; i++) {
            sum += pow((p[i] - other[i]), 2);
        }
        return sqrt(sum);
    }

    template<typename Container>
    inline double hausdorffDistance(const Container &first, const Container &second) {
        typedef typename Container::value_type Point;
        typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
        typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;

        L2Metric l2Metric;
        double distanceMax = 0;
        for (auto it = first.begin(), ite = first.end(); it != ite; ++it) {
            Point firstP = *it;
            Point closestPointInTheoretical = *std::min_element(second.begin(), second.end(),
                                                           [&](const Point &one, const Point &two) {
                                                               return l2Metric(one, firstP) < l2Metric(two, firstP);
                                                           });
            double distance = l2Metric(closestPointInTheoretical, firstP);
            if (distance > distanceMax)
                distanceMax = distance;
        }
        return distanceMax;
    }

    template<typename Container>
    inline double distanceSet(const Container &first, const Container &second) {
        typedef typename Container::value_type Point;
        typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
        typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;

        L2Metric l2Metric;
        double distanceMax = std::numeric_limits<double>::max();
        if (first.size() == 0 || second.size() == 0) return distanceMax;
        for (auto it = first.begin(), ite = first.end(); it != ite; ++it) {
            Point firstP = *it;
            Point closestPointInTheoretical = *std::min_element(second.begin(), second.end(),
                                                           [&](const Point &one, const Point &two) {
                                                               return l2Metric(one, firstP) < l2Metric(two, firstP);
                                                           });
            double distance = l2Metric(closestPointInTheoretical, firstP);
            if (distance < distanceMax)
                distanceMax = distance;
        }
        return distanceMax;
    }


    template<typename Point, typename Container>
    inline double geodesicDistance(const Point &first, const Point &second, const Container &object) {
        typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
        typedef DGtal::MetricAdjacency<Space, 1> Adj6;
        typedef DGtal::MetricAdjacency<Space, 3> Adj26;
        typedef DGtal::DigitalTopology<Adj26, Adj6> DT26_6;
        typedef DGtal::Object<DT26_6, Container> ObjectType;
        typedef DGtal::BreadthFirstVisitor<ObjectType, std::set<Point> > Visitor;
        typedef typename Visitor::Node MyNode;

        Adj26 adj26;
        Adj6 adj6;
        DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
        ObjectType graph(dt26_6, object);
        Visitor visitor(graph, first);
        MyNode node;
        double distance = std::numeric_limits<double>::max();
        while (!visitor.finished()) {
            node = visitor.current();
            if (node.first == second) {
                distance = node.second;
                break;
            }
            visitor.expand();
        }
        return distance;
    }
};

#endif
