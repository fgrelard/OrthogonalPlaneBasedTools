#ifndef TRUE_PREDICATE_H
#define TRUE_PREDICATE_H

#include <vector>
#include "shapes/DigitalPlane.h"

template<typename Space>
class TruePredicate {
public:
    typedef DigitalPlane<Space> Plane;
    typedef typename Space::Point Point;
public:
    TruePredicate() {}

    TruePredicate(const std::vector<Plane> &planes) {}

public:
    bool operator()(const Point &one, const Point &second) {
        return true;
    }
};

#endif
