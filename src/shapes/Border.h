#ifndef BORDER_H
#define BORDER_H

#include "DGtal/kernel/sets/CDigitalSet.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"

template<typename Container>
class Border {
    BOOST_CONCEPT_ASSERT((DGtal::concepts::CDigitalSet<Container>));

protected:
    typedef typename Container::Point Point;
    typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
    typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;

public:
    Border() = delete;

    Border(const Container &aSet) : myBorder(convertToBorder(aSet)) {}

    Border(const Border &other) : myBorder(other.myBorder) {}


    Container pointSet() const { return myBorder; }

private:
    Container convertToBorder(const Container &aSet) {
        typedef DGtal::DistanceTransformation<Space, Container, L2Metric> DTL2;

        L2Metric l2Metric;
        Container surfacePoints(aSet.domain());
        DTL2 dt(&aSet.domain(), &aSet, &l2Metric);

        for (const Point &p : aSet) {
            if (dt(p) > 0 && dt(p) <= 1)
                surfacePoints.insert(p);
        }
        return surfacePoints;
    }

private:
    Container myBorder;
};

#endif
