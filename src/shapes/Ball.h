#ifndef BALL_H
#define BALL_H


#include <vector>
#include <set>
#include <iostream>
#include "DGtal/base/Common.h"
#include "geometry/Distance.h"
#include "geometry/PointUtil.h"
#include "shapes/DigitalPlane.h"

template<typename Point>
class Ball {
public:
    typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
    typedef DGtal::HyperRectDomain<Space> Domain;
    typedef typename DGtal::DigitalSetSelector<Domain, DGtal::BIG_DS + DGtal::HIGH_BEL_DS>::Type DigitalSet;
    typedef typename Space::RealVector RealVector;
    typedef typename Space::Dimension Dimension;
public:
    Ball() : myRadius(0.0), myCenter({0, 0, 0}) {}

    Ball(const Point &center, double radius) : myCenter(center), myRadius(radius) {}

    Ball(const Ball &other) : myCenter(other.myCenter), myRadius(other.myRadius) {}

public:
    inline bool contains(const Point &point) const { return Distance::euclideanDistance(point, myCenter) <= myRadius; }


    template<typename Image>
    Image intersection(const Image &image);

    template<typename Image>
    Image surfaceIntersection(const Image &image);

    DigitalSet intersection(const DigitalSet &setPoint);

    DigitalSet surfaceIntersection(const DigitalSet &setSurface);

    DigitalSet pointSet() const;

    DigitalSet pointsInHalfBall(const RealVector &normal = RealVector(0, 1, 0)) const;

    DigitalSet pointsSurfaceBall() const;

public:
    Point getCenter() const { return myCenter; }

    double getRadius() const { return myRadius; }

public:
    bool operator!=(const Ball &other) const { return (myCenter != other.myCenter || myRadius != other.myRadius); }

private:
    Point myCenter;
    double myRadius;
};

template<typename Point>
template<typename Image>
Image
Ball<Point>::
intersection(const Image &image) {
    DigitalSet points = pointSet();
    Domain domainImage = image.domain();
    Domain domainPointSet = points.domain();
    Point lowerS = domainPointSet.lowerBound();
    Point upperS = domainPointSet.upperBound();
    Domain domain(PointUtil::box(lowerS, domainImage),
                  PointUtil::box(upperS, domainImage));
    Image other(domain);
    for (auto it = domain.begin(), ite = domain.end();
         it != ite; ++it) {
        Point p = *it;
        if (contains(p))
            other.setValue(p, image(p));
        else
            other.setValue(p, 0);
    }
    return other;
}

template<typename Point>
template<typename Image>
Image
Ball<Point>::
surfaceIntersection(const Image &image) {
    DigitalSet points = pointSet();
    Domain domainImage = image.domain();
    Domain domainPointSet = points.domain();
    Point lowerS = domainPointSet.lowerBound();
    Point upperS = domainPointSet.upperBound();
    Domain domain(PointUtil::box(lowerS, domainImage),
                  PointUtil::box(upperS, domainImage));
    Image other(domain);
    DigitalSet surface = pointsSurfaceBall();
    for (auto it = domain.begin(), ite = domain.end();
         it != ite; ++it) {
        Point p = *it;
        if (surface.find(p) != surface.end())
            other.setValue(p, image(p));
        else
            other.setValue(p, 0);
    }
    return other;
}

template<typename Point>
typename Ball<Point>::DigitalSet Ball<Point>::intersection(const DigitalSet &setPoint) {
    DigitalSet intersection(setPoint.domain());
    for (const Point &p : setPoint) {
        double distance = Distance::euclideanDistance(p, myCenter);
        if (distance <= myRadius) {
            intersection.insert(p);
        }
    }
    return intersection;
}

template<typename Point>
typename Ball<Point>::DigitalSet Ball<Point>::surfaceIntersection(const DigitalSet &setSurface) {
    DigitalSet intersection(setSurface.domain());
    for (auto it = setSurface.begin(), ite = setSurface.end(); it != ite; ++it) {
        double distance = Distance::euclideanDistance(*it, myCenter);
        if (distance >= myRadius - 1 && distance <= myRadius) {
            intersection.insert(*it);
        }
    }
    return intersection;
}

template<typename Point>
typename Ball<Point>::DigitalSet Ball<Point>::pointSet() const {
    Point lower = myCenter + Point::diagonal(-myRadius);
    Point upper = myCenter + Point::diagonal(myRadius + 1);
    Domain domain(lower, upper);
    DigitalSet points(domain);
    for (const Point &p : domain) {
        if (contains(p))
            points.insert(p);
    }

    return points;
}


template<typename Point>
typename Ball<Point>::DigitalSet Ball<Point>::pointsInHalfBall(const RealVector &normal) const {
    DigitalSet points = pointSet();
    DigitalSet pointsInHalf(points.domain());
    DigitalPlane<Space> plane(myCenter, normal);
    for (const Point &p : points) {
        if (plane.isPointAbove(p))
            pointsInHalf.insert(p);
    }
    return pointsInHalf;
}

template<typename Point>
typename Ball<Point>::DigitalSet Ball<Point>::pointsSurfaceBall() const {
    Point lower = myCenter + Point::diagonal(-myRadius);
    Point upper = myCenter + Point::diagonal(myRadius + 1);
    Domain domain(lower, upper);
    DigitalSet points(domain);
    for (const Point &p : domain) {
        double distance = Distance::euclideanDistance(p, myCenter);
        if (distance >= myRadius - 1 && distance <= myRadius)
            points.insert(p);
    }


    return points;
}

#endif
