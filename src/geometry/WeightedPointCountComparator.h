#ifndef WEIGHTEDPOINTCOUNTCOMPARATOR_H
#define WEIGHTEDPOINTCOUNTCOMPARATOR_H

#include "WeightedPointCount.h"

template<typename Point>
struct WeightedPointCountComparator {
    bool operator()(const Point *lhs, const Point *rhs) const {
        return lhs->myWeight >= rhs->myWeight;
    }
};

#endif
