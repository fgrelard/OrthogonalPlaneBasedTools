#ifndef ELLIPSE_H
#define ELLIPSE_H

#include <vector>
#include <set>
#include <limits>
#include <algorithm>
#include "geometry/Distance.h"
#include "shapes/Ball.h"


template<typename Point>
class Ellipse {
public:
    Ellipse(const std::vector<Point> &points, const Point &center);

    bool computeFocalPoints(Point &, Point &);

    bool distanceConsistentEllipse(const std::vector<Point> &points);

public:
    double myMajorAxis = 0.0;
    double myMinorAxis = 0.0;
    Point myCenter;
    std::pair<Point, Point> myPointsForMajor;
    std::pair<Point, Point> myPointsForMinor;
private:
    void computeParameters(const std::vector<Point> &points);
};


template<typename Point>
Ellipse<Point>::Ellipse(const std::vector<Point> &points, const Point &center) {
    myCenter = center;
    computeParameters(points);

}

template<typename Point>
void Ellipse<Point>::computeParameters(const std::vector<Point> &points) {
    //need 5 points for ellipse computation
    if (points.size() < 6) return;
    double maximum = 0.0;
    double minimum = std::numeric_limits<double>::max();

    Point maxi = *(std::max_element(points.begin(), points.end(), [&](const Point &one, const Point &two) {
        if (euclideanDistance(two, myCenter) > maximum) {
            maximum = euclideanDistance(two, myCenter);
            return true;
        } else return false;
    }));
    Point mini = *(std::max_element(points.begin(), points.end(), [&](const Point &one, const Point &two) {
        if (euclideanDistance(two, myCenter) < minimum) {
            minimum = euclideanDistance(two, myCenter);
            return true;
        } else return false;
    }));

    Point maxiToCenter = myCenter - maxi;
    Point miniToCenter = myCenter - mini;

    Point symmetricMaxi = myCenter + maxiToCenter;
    Point symmetricMini = myCenter + miniToCenter;

    myPointsForMajor = std::make_pair(maxi, symmetricMaxi);
    myPointsForMinor = std::make_pair(mini, symmetricMini);

    myMinorAxis = minimum * 2;
    myMajorAxis = maximum * 2;
}


template<typename Point>
bool Ellipse<Point>::computeFocalPoints(Point &focal, Point &otherFocal) {
    Point bigA = myPointsForMajor.first;
    Point bigAprime = myPointsForMajor.second;
    Point directionVector = bigAprime - bigA;

    std::set<Point> lineMajorAxis;
    for (float i = -myMajorAxis; i <= myMajorAxis; i += 0.1) {
        int x = bigA[0] + directionVector[0] * i;
        int y = bigA[1] + directionVector[1] * i;
        int z = bigA[2] + directionVector[2] * i;
        lineMajorAxis.insert({x, y, z});
    }

    Ball<Point> ball(myPointsForMinor.first, myMajorAxis / 2);
    std::set<Point> surfaceBall = ball.pointsSurfaceBall();
    std::set<Point> intersection;
    set_intersection(lineMajorAxis.begin(), lineMajorAxis.end(), surfaceBall.begin(), surfaceBall.end(),
                     inserter(intersection, intersection.begin()));
    if (intersection.size() == 2) {
        focal = *(intersection.begin());
        otherFocal = *(++(intersection.begin()));
        return true;
    }
    return false;
}

template<typename Point>
bool Ellipse<Point>::distanceConsistentEllipse(const std::vector<Point> &points) {
    Point focalOne;
    Point focalTwo;

    bool foundFocalPoints = computeFocalPoints(focalOne, focalTwo);
    if (focalOne == focalTwo) return false;
    if (foundFocalPoints) {
        for (auto it = points.begin(), itE = points.end(); it != itE; ++it) {
            int sumCeil = ceil(euclideanDistance(*it, focalOne) + euclideanDistance(*it, focalTwo));
            int sumFloor = floor(euclideanDistance(*it, focalOne) + euclideanDistance(*it, focalTwo));
            if (sumCeil != (int) myMajorAxis && sumFloor != (int) myMajorAxis) {
                return false;
            }
        }
        return true;
    }
    return false;
}

#endif
