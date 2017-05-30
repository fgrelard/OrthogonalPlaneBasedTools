#ifndef DISTANCE_TO_POINT_FUNCTOR_H
#define DISTANCE_TO_POINT_FUNCTOR_H

#include "DGtal/base/Clone.h"

template<typename Distance>
struct DistanceToPointFunctor {
    typedef typename Distance::Space Space;
    typedef typename Distance::Value Value;
    typedef typename Space::Point Point;


    DistanceToPointFunctor(DGtal::Clone<Distance> distance,
                           const Point &aP)
            : myDistance(distance), p(aP) {}

    Value operator()(const Point &q) const {
        return myDistance(p, q);
    }

    Distance myDistance;
    Point p;
};

#endif
