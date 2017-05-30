#ifndef BRESENHAM_ALGORITHM_H
#define BRESENHAM_ALGORITHM_H

#include <set>
#include <vector>
#include "geometry/path/LinkPointAlgorithm.h"
#include "DGtal/kernel/PointVector.h"

template<typename Point>
class BresenhamAlgorithm : public LinkPointAlgorithm<Point> {
public:
    typedef DGtal::PointVector<Point::dimension, double> RealVector;
public:
    BresenhamAlgorithm() {}

    BresenhamAlgorithm(const Point &source,
                       const Point &destination) : LinkPointAlgorithm<Point>(source, destination) {}

    BresenhamAlgorithm(const BresenhamAlgorithm &other) : LinkPointAlgorithm<Point>(other) {}

    virtual typename LinkPointAlgorithm<Point>::Path linkPoints();

};

template<typename Point>
typename LinkPointAlgorithm<Point>::Path
BresenhamAlgorithm<Point>::linkPoints() {
    Point source = this->mySource, destination = this->myDestination;

    std::vector<Point> pointsBetween;
    RealVector slope = this->myDestination - this->mySource;
    std::set<int> coordinates = {0, 1, 2};
    double maxValue = 0;
    int indexForMaxValue = 0, indexForSecond = 0, indexForThird = 0;
    for (auto it = coordinates.begin(), ite = coordinates.end(); it != ite; ++it) {
        double value = std::abs(slope[*it]);
        if (value > maxValue) {
            maxValue = value;
            indexForMaxValue = *it;
        }
    }

    maxValue = -1;
    for (auto it = coordinates.begin(), ite = coordinates.end(); it != ite; ++it) {
        if (*it == indexForMaxValue) continue;
        double value = std::abs(slope[*it]);
        if (value > maxValue) {
            maxValue = value;
            indexForSecond = *it;
        }
    }

    for (auto it = coordinates.begin(), ite = coordinates.end(); it != ite; ++it) {
        if (*it == indexForMaxValue || *it == indexForSecond) continue;
        indexForThird = *it;
    }


    bool toReverse = false;
    if (this->mySource[indexForMaxValue] > this->myDestination[indexForMaxValue]) {
        source = this->myDestination;
        destination = this->mySource;
        slope = -slope;
        toReverse = true;
    }

    //y-init
    double deltaX = slope[indexForMaxValue];
    int signX = (slope[indexForMaxValue] > 0) - (slope[indexForMaxValue] < 0);

    double deltaY = std::abs(slope[indexForSecond]);
    int signY = (slope[indexForSecond] > 0) - (slope[indexForSecond] < 0);
    double ey = 2 * deltaY - deltaX;
    double yinc1 = 2 * deltaY;
    double yinc2 = 2 * (deltaY - deltaX);

    //z-init
    double deltaZ = std::abs(slope[indexForThird]);
    int signZ = (slope[indexForThird] > 0) - (slope[indexForThird] < 0);
    double ez = 2 * deltaZ - deltaX;
    double zinc1 = 2 * deltaZ;
    double zinc2 = 2 * (deltaZ - deltaX);


    Point current = source, next = source;

    Point offsetX, offsetY, offsetZ;
    offsetX[indexForMaxValue] = signX;
    offsetY[indexForSecond] = signY;
    offsetZ[indexForThird] = signZ;

    int i = current[indexForMaxValue];
    int imax = destination[indexForMaxValue];
    while (i < imax) {
        i++;
        if (current != source && current != destination)
            pointsBetween.push_back(current);
        current += offsetX;
        if (ey < 0) {
            ey += yinc1;

        } else {
            ey += yinc2;
            current += offsetY;
        }

        if (ez < 0) {
            ez += zinc1;
        } else {
            ez += zinc2;
            current += offsetZ;
        }
    }
    if (toReverse)
        std::reverse(pointsBetween.begin(), pointsBetween.end());
    typename LinkPointAlgorithm<Point>::Path path(pointsBetween);
    return path;
}

#endif
