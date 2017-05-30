#ifndef MEDIAL_AXIS_H
#define MEDIAL_AXIS_H

#include <vector>
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/kernel/sets/CDigitalSet.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/topology/Object.h"
#include "DGtal/topology/DigitalTopology.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "geometry/DistanceToPointFunctor.h"


template<typename Container>
class MedialAxis {

    BOOST_CONCEPT_ASSERT((DGtal::concepts::CDigitalSet<Container>));

protected:
    typedef typename Container::value_type Point;
    typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
    typedef DGtal::MetricAdjacency<typename Container::Space, 1> Adj6;
    typedef DGtal::MetricAdjacency<typename Container::Space, 3> Adj26;
    typedef DGtal::DigitalTopology<Adj26, Adj6> DT26_6;
    typedef DGtal::Object<DT26_6, Container> ObjectType;
    typedef typename Container::Space::RealVector RealVector;

public:
    MedialAxis() = delete;

    MedialAxis(const Container &aSet) : mySet(compute(aSet)) {}

    MedialAxis(const MedialAxis &other) : mySet(other.mySet) {}

private:
    Container compute(const Container &aSet);

public:
    Container pointSet() const { return mySet; }

private:
    Container mySet;
};

template<typename Container>
Container
MedialAxis<Container>::compute(const Container &aSet) {
    typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;
    typedef DGtal::DistanceTransformation<Space, Container, L2Metric> DTL2;

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    ObjectType obj(dt26_6, aSet);

    L2Metric l2Metric;
    DTL2 imageFct(&aSet.domain(), &aSet, &l2Metric);
    Container medialAxis(aSet.domain());

    for (const Point &p : aSet) {
        bool add = true;
        double d = imageFct(p);
        std::vector<Point> neighbors;
        std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
        obj.writeNeighbors(inserter, p);
        for (const Point &n : neighbors) {
            float distanceToBoundary = imageFct(n);
            float minEnclosingRadius = sqrt(1 + pow(d, 2));
            if (d <= 1 || distanceToBoundary > d) {
                add = false;
            }
        }
        if (add)
            medialAxis.insert(p);
    }
    return medialAxis;
}

#endif
