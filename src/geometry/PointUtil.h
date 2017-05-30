#ifndef POINT_UTIL_H
#define POINT_UTIL_H


#include <vector>
#include <set>
#include <limits>
#include "DGtal/base/Common.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "geometry/path/BresenhamAlgorithm.h"

namespace PointUtil {


    template<typename Domain, typename Container>
    Domain computeBoundingBox(const Container &points);

    template<typename Domain>
    typename Domain::Point box(const typename Domain::Point &p,
                               const Domain &domain);

    template<typename Container, typename Point, typename RealVector>
    Point trackPoint(const Container &container, const Point &start, const RealVector &vector);

    template<typename Container, typename Point, typename RealVector>
    Container traversedLineTracking(const Container &container, const Point &start, const RealVector &vector);

    template<typename Container, typename Point, typename RealVector>
    std::pair<Point, Point> twoClosestPointsTrackingWithNormal(const Container &container,
                                                               const Point &ref, const RealVector &vRef,
                                                               const Point &other, const RealVector &vOther);
}


template<typename Domain, typename Container>
Domain PointUtil::computeBoundingBox(const Container &points) {
    typedef typename Container::value_type Point;
    int min = std::numeric_limits<int>::max();
    int max = std::numeric_limits<int>::min();

    Point low = Point::diagonal(min);
    Point up = Point::diagonal(max);
    for (const auto &point : points) {
        for (int i = 0; i < Point::dimension; i++) {
            low[i] = (point[i] < low[i]) ? point[i] : low[i];
            up[i] = (point[i] > up[i]) ? point[i] : up[i];
        }
    }
    Domain domain(low, up);
    return domain;
}

template<typename Domain>
typename Domain::Point
PointUtil::
box(const typename Domain::Point &p,
    const Domain &domain) {
    typedef typename Domain::Point Point;
    Point lower = domain.lowerBound();
    Point upper = domain.upperBound();
    return p.sup(lower).inf(upper);
}

template<typename Container, typename Point, typename RealVector>
Point
PointUtil::
trackPoint(const Container &container, const Point &start, const RealVector &vector) {
    Point point(start);
    Point previous(start);
    int scalar = 1;
    while (container.find(point) != container.end()) {
        previous = point;
        point = start + vector * scalar;
        scalar++;
    }
    return previous;
}

template<typename Container, typename Point, typename RealVector>
Container
PointUtil::
traversedLineTracking(const Container &volume, const Point &start, const RealVector &dir) {
    Container container(volume.domain());
    RealVector vector = dir.getNormalized();
    Point tracked = trackPoint(volume, start, vector);
    BresenhamAlgorithm<Point> bresenham(start, tracked);
    std::vector<Point> path = bresenham.linkPoints();
    container.insert(path.begin(), path.end());
    return container;
}


template<typename Container, typename Point, typename RealVector>
std::pair<Point, Point>
PointUtil::
twoClosestPointsTrackingWithNormal(const Container &container, const Point &reference, const RealVector &dirRef,
                                   const Point &other, const RealVector &dirOther) {

    if (dirRef == RealVector::zero || dirOther == RealVector::zero) return std::make_pair(other, reference);
    RealVector normalRef = dirRef.getNormalized();
    RealVector normalOther = dirOther.getNormalized();
    RealVector dirVectorReference = (other - reference).getNormalized();
    RealVector dirVectorOther = (reference - other).getNormalized();
    if (normalRef.dot(dirVectorReference) < 0)
        normalRef = -normalRef;
    if (normalOther.dot(dirVectorOther) < 0)
        normalOther = -normalOther;


    Container traversedCurrent = PointUtil::traversedLineTracking(container, reference, normalRef);
    Container traversedReference = PointUtil::traversedLineTracking(container, other, normalOther);
    if (traversedCurrent.size() == 0 || traversedReference.size() == 0)
        return std::make_pair(other, reference);

    double distanceCR = std::numeric_limits<double>::max();
    Point closest1, closest2;
    DGtal::ExactPredicateLpSeparableMetric<typename Container::Space, 2> l2Metric;

    for (auto it = traversedCurrent.begin(), ite = traversedCurrent.end(); it != ite; ++it) {
        Point current = *it;
        Point nearest = *std::min_element(traversedReference.begin(), traversedReference.end(),
                                          [&](const Point &one, const Point &two) {
                                              return (l2Metric(one, current) < l2Metric(two, current)
                                                      && l2Metric(one, current) > sqrt(3));
                                          });
        double currentDistance = l2Metric(nearest, current);
        if (currentDistance < distanceCR) {
            distanceCR = currentDistance;
            closest1 = *it;
            closest2 = nearest;
        }
    }

    return std::make_pair(closest1, closest2);
}


#endif
