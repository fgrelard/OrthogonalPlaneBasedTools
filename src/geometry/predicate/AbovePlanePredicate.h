#ifndef ABOVE_PLANE_PREDICATE_H
#define ABOVE_PLANE_PREDICATE_H

#include "shapes/DigitalPlane.h"

template<typename Space>
class AbovePlanePredicate {
public:
    typedef DigitalPlane<Space> Plane;
    typedef typename Space::Point Point;
public:
    AbovePlanePredicate() {}

    AbovePlanePredicate(const std::vector<Plane> &planes) : myPlanes(planes) {}

public:
    bool operator()(const Point &reference, const Point &pointToTest) {
        if (myPlanes.size() == 0) return true;
        bool add = false;
        auto iterator = std::find_if(myPlanes.begin(), myPlanes.end(),
                                     [&](const Plane &plane) {
                                         return (plane.getCenter() == reference);
                                     });
        if (iterator == myPlanes.end()) return true;

        Plane plane = *iterator;
        return plane.isPointAbove(pointToTest);
    }

private:
    std::vector<Plane> myPlanes;
};

#endif
