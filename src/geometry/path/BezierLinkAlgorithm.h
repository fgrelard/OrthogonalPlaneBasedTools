#ifndef BEZIER_LINK_ALGORITHM_H
#define BEZIER_LINK_ALGORITHM_H

#include <set>
#include <vector>
#include "geometry/path/LinkPointAlgorithm.h"
#include "DGtal/kernel/PointVector.h"

template<typename Point>
class BezierLinkAlgorithm : public LinkPointAlgorithm<Point> {
public:
    typedef DGtal::PointVector<Point::dimension, double> RealVector;
public:
    BezierLinkAlgorithm() {}

    BezierLinkAlgorithm(const Point &source,
                        const Point &destination,
                        const Point &controlPoint1,
                        const Point &controlPoint2) : LinkPointAlgorithm<Point>(source, destination),
                                                      myControlPoint1(controlPoint1),
                                                      myControlPoint2(controlPoint2) {}

    BezierLinkAlgorithm(const BezierLinkAlgorithm &other) : LinkPointAlgorithm<Point>(other),
                                                            myControlPoint1(other.myControlPoint1),
                                                            myControlPoint2(other.myControlPoint2) {}

    virtual typename LinkPointAlgorithm<Point>::Path linkPoints();

protected:
    Point myControlPoint1;
    Point myControlPoint2;
};

template<typename Point>
typename LinkPointAlgorithm<Point>::Path
BezierLinkAlgorithm<Point>::linkPoints() {
    std::vector<Point> bezier;
    for (double t = 0; t < 1; t += 0.01) {
        Point current;
        for (int k = 0; k < Point::dimension; k++) {
            double coordinateAtK = pow((1 - t), 3) * this->mySource[k] +
                                   3 * pow((1 - t), 2) * t * this->myControlPoint1[k] +
                                   3 * (1 - t) * pow(t, 2) * this->myControlPoint2[k] +
                                   pow(t, 3) * this->myDestination[k];
            current[k] = std::round(coordinateAtK);
        }
        bezier.push_back(current);
    }
    typename LinkPointAlgorithm<Point>::Path path(bezier);
    return path;
}

#endif
