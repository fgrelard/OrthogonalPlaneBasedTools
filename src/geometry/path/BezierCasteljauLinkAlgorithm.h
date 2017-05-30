#ifndef BEZIER_CASTELJAU_LINK_ALGORITHM_H
#define BEZIER_CASTELJAU_LINK_ALGORITHM_H

#include <set>
#include <vector>
#include <map>
#include "geometry/path/LinkPointAlgorithm.h"
#include "geometry/path/BresenhamAlgorithm.h"
#include "geometry/path/BezierLinkAlgorithm.h"
#include "shapes/Polygon.h"
#include "DGtal/kernel/PointVector.h"

template<typename Point>
class BezierCasteljauLinkAlgorithm : public BezierLinkAlgorithm<Point> {
public:
    typedef DGtal::PointVector<Point::dimension, double> RealVector;
public:
    BezierCasteljauLinkAlgorithm() {};

    BezierCasteljauLinkAlgorithm(const Point &source,
                                 const Point &destination,
                                 const Point &controlPoint1,
                                 const Point &controlPoint2) : BezierLinkAlgorithm<Point>(source, destination,
                                                                                          controlPoint1,
                                                                                          controlPoint2) {}

    BezierCasteljauLinkAlgorithm(const BezierCasteljauLinkAlgorithm &other) : BezierCasteljauLinkAlgorithm<Point>(
            other) {}

    virtual typename LinkPointAlgorithm<Point>::Path linkPoints();

private:
    std::vector<Polygon<Point> > createPolygonSubdivision(const Polygon<Point> &polygon);

    void recursivePolygonSubdivisionBezier(const Polygon<Point> &polygon, std::map<Point, Point> &points);

    bool stopRecursivePolygonDecomposition(const Polygon<Point> &polygon);

    int computeDiameter(const std::set<int> &equation);

    std::vector<std::pair<Point, Point> > orderPoints(const std::map<Point, Point> &points);


};

template<typename Point>
std::vector<Polygon<Point> >
BezierCasteljauLinkAlgorithm<Point>::createPolygonSubdivision(const Polygon<Point> &polygon) {
    typedef std::pair<Point, Point> Segment;

    Point source = polygon.getPolygon()[0];
    Point controlPoint1 = polygon.getPolygon()[1];
    Point controlPoint2 = polygon.getPolygon()[2];
    Point destination = polygon.getPolygon()[3];

    std::vector<Polygon<Point> > polygons;
    Segment s1(source, controlPoint1);
    Segment s2(controlPoint1, controlPoint2);
    Segment s3(controlPoint2, destination);

    std::vector<Segment> segments = {s1, s2, s3};
    std::vector<Point> newPoints;
    while (segments.size() > 0) {
        std::vector<Point> currentPoints;
        for (const Segment &s : segments) {
            Point first = s.first;
            Point second = s.second;
            Point middle = (first + second) * 1 / 2.;
            newPoints.push_back(middle);
            currentPoints.push_back(middle);
        }

        segments.clear();
        for (int i = 0; i < currentPoints.size() - 1; i++) {
            Segment s(currentPoints[i], currentPoints[i + 1]);
            segments.push_back(s);
        }
    }
    Polygon<Point> polygon1{source, newPoints[0], newPoints[3], newPoints[5]};
    Polygon<Point> polygon2{newPoints[5], newPoints[4], newPoints[2], destination};
    polygons.push_back(polygon1);
    polygons.push_back(polygon2);
    return polygons;
}


template<typename Point>
int BezierCasteljauLinkAlgorithm<Point>::computeDiameter(const std::set<int> &equation) {
    return *equation.rbegin() - *equation.begin();
}

template<typename Point>
bool BezierCasteljauLinkAlgorithm<Point>::stopRecursivePolygonDecomposition(const Polygon<Point> &polygon) {
    std::vector<Point> points = polygon.getPolygon();
    Point dirVector = points[3] - points[0];
    std::set<int> firstEquation, secondEquation;
    for (const Point &point : points) {
        int value = -dirVector[2] * point[0] + dirVector[0] * point[2];
        int secondValue = -dirVector[2] * point[1] + dirVector[1] * point[2];
        firstEquation.insert(value);
        secondEquation.insert(secondValue);
    }
    double diameter1 = computeDiameter(firstEquation);
    double diameter2 = computeDiameter(secondEquation);

    double valueToCheck = std::max(std::abs(dirVector[0]), std::max(std::abs(dirVector[1]), std::abs(dirVector[2])));
    valueToCheck *= 4 / 3.;
    return (diameter1 <= valueToCheck && diameter2 <= valueToCheck);
}

template<typename Point>
void BezierCasteljauLinkAlgorithm<Point>::recursivePolygonSubdivisionBezier(const Polygon<Point> &polygon,
                                                                            std::map<Point, Point> &points) {
    if (!stopRecursivePolygonDecomposition(polygon)) {
        std::vector<Polygon<Point> > polygons = createPolygonSubdivision(polygon);
        this->recursivePolygonSubdivisionBezier(polygons[0], points);
        this->recursivePolygonSubdivisionBezier(polygons[1], points);
    } else {
        points[polygon.getPolygon()[0]] = polygon.getPolygon()[3];
    }
}

template<typename Point>
std::vector<std::pair<Point, Point> >
BezierCasteljauLinkAlgorithm<Point>::orderPoints(const std::map<Point, Point> &points) {
    Point next = this->mySource;
    std::vector<std::pair<Point, Point>> orderedPoints;
    while (next != this->myDestination) {
        Point current = next;
        auto iterator = points.find(next);
        if (iterator == points.end()) break;
        next = iterator->second;
        orderedPoints.push_back(make_pair(current, next));
    }
    return orderedPoints;
}


template<typename Point>
typename LinkPointAlgorithm<Point>::Path
BezierCasteljauLinkAlgorithm<Point>::linkPoints() {
    Polygon<Point> polygon{this->mySource, this->myControlPoint1, this->myControlPoint2, this->myDestination};
    std::map<Point, Point> points;
    recursivePolygonSubdivisionBezier(polygon, points);
    std::vector<std::pair<Point, Point> > orderedPoints = this->orderPoints(points);
    std::vector<Point> linkPoints;
    for (const auto &pair : orderedPoints) {
        BresenhamAlgorithm<Point> bresenhamAlgo(pair.first, pair.second);
        typename LinkPointAlgorithm<Point>::Path curve = bresenhamAlgo.linkPoints();
        if (pair.first != this->mySource && pair.first != this->myDestination)
            linkPoints.push_back(pair.first);
        linkPoints.insert(linkPoints.end(), curve.begin(), curve.end());
        if (pair.second != this->myDestination && pair.second != this->mySource)
            linkPoints.push_back(pair.second);
    }
    typename LinkPointAlgorithm<Point>::Path path(linkPoints);
    return path;
}

#endif
