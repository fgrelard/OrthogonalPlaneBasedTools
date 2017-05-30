#ifndef WEIGHTED_POINT_COUNT_H
#define WEIGHTED_POINT_COUNT_H

#include "geometry/WeightedPoint.h"

template<typename Point>
class WeightedPointCount : public WeightedPoint<Point> {
    typedef WeightedPoint<Point> Base;
public:
    WeightedPointCount() {}

    WeightedPointCount(const Point &aPoint, double aWeight, int aCount) : Base(aPoint, aWeight), myCount(aCount) {}

    WeightedPointCount(const WeightedPointCount &other) : Base(other), myCount(other.myCount) {}

    using Base::Base;

    double getWeightedValue() const { return Base::myWeight / myCount; }

    friend bool operator<(const WeightedPointCount &it, const WeightedPointCount &other) {
        return (it.myWeight >= other.myWeight);
    }

public:
    int myCount = 0;
};

#endif
