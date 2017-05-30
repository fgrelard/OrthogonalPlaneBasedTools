#ifndef CONE_H
#define CONE_H

#include <vector>
#include <set>
#include <cmath>
#include "shapes/DigitalPlane.h"
#include "shapes/Ball.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/base/Common.h"
#include "geometry/PointUtil.h"

template<typename Point, typename Vector>
class Cone {
public:
    typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
    typedef DGtal::HyperRectDomain<Space> Domain;
    typedef typename DGtal::DigitalSetSelector<Domain, DGtal::BIG_DS + DGtal::HIGH_BEL_DS>::Type DigitalSet;

public:
    Cone(const Point &aVertex, const Vector &anOrientation,
         double aRadius, double aHeight) : myVertex(aVertex), myOrientation(anOrientation),
                                           myRadius(aRadius), myHeight(aHeight) {}

    bool contains(const Point &point);

    DigitalSet pointSet();

private:
    Point myVertex;
    Vector myOrientation;
    double myRadius;
    double myHeight;
};


template<typename Point, typename Vector>
bool Cone<Point, Vector>::contains(const Point &point) {
    if (point == myVertex) return true;
    double halfAperture = atan(myRadius / myHeight);

    Vector vertexToPoint = (point - myVertex);
    Vector vertexToCenter = myOrientation * myHeight;
    double x = vertexToPoint.getNormalized().dot(myOrientation);
    if (x < -1.0) x = -1.0;
    else if (x > 1.0) x = 1.0;
    bool inCone = (acos(x) <= halfAperture);
    bool isProjection = (vertexToPoint.dot(vertexToCenter) / vertexToCenter.norm()) < vertexToCenter.norm();

    return (inCone &&
            isProjection);
}

template<typename Point, typename Vector>
typename Cone<Point, Vector>::DigitalSet Cone<Point, Vector>::pointSet() {
    Point l1 = myVertex - Point::diagonal(myRadius);
    Point l2 = myVertex + Point::diagonal(myRadius);

    Point center = myVertex + myOrientation * myHeight;
    Point l3 = center - Point::diagonal(myRadius);
    Point l4 = center + Point::diagonal(myRadius);

    Point lower = std::min(l1, std::min(l2, std::min(l3, l4)));
    Point upper = std::max(l1, std::max(l2, std::max(l3, l4)));
    std::vector<Point> points;
    Domain domain(lower, upper);
    for (const auto &p : domain) {
        if (contains(p))
            points.push_back(p);
    }
    DigitalSet aSet(domain);
    aSet.insert(points.begin(), points.end());
    return aSet;
}

#endif
