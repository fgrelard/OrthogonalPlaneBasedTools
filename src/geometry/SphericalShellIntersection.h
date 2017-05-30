#ifndef SPHERICAL_SHELL_INTERSECTION_H
#define SPHERICAL_SHELL_INTERSECTION_H

#include "DGtal/kernel/sets/CDigitalSet.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/topology/DigitalTopology.h"
#include "DGtal/topology/Object.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"

#include "geometry/DistanceToPointFunctor.h"
#include "shapes/DigitalPlane.h"

template<typename Container>
class SphericalShellIntersection {
    BOOST_CONCEPT_ASSERT((DGtal::concepts::CDigitalSet<Container>));

public:
    typedef typename Container::Point Point;
    typedef typename Container::Space Space;
    typedef typename Space::RealVector RealVector;
    typedef typename DGtal::HyperRectDomain<Space> Domain;
    typedef typename DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;
    typedef DistanceToPointFunctor<L2Metric> DistanceFunctor;
    typedef DGtal::MetricAdjacency<Space, 1> Adj6;
    typedef DGtal::MetricAdjacency<Space, 3> Adj26;
    typedef DGtal::DigitalTopology<Adj26, Adj6> DT26_6;
    typedef DGtal::Object<DT26_6, Container> ObjectType;

public:
    SphericalShellIntersection() = delete;

    SphericalShellIntersection(const Container &setVolume, const Point &center, double radiusInnerBall,
                               const RealVector &dirVector = RealVector::zero);

    SphericalShellIntersection(const Container &setVolume, const Point &center, double radiusInnerBall,
                               double radiusOuterBall, const RealVector &dirVector = RealVector::zero);

    SphericalShellIntersection(const SphericalShellIntersection &other);

    ~SphericalShellIntersection();

public:
    Container shell();

    unsigned int degree(const Container &shell, double minCCSize = 0);

    unsigned int degree(double minCCSize = 0);

private:
    Container *myContainer;
    Point myCenter;
    double myRadiusInner;
    double myRadiusOuter;
    RealVector myDirectionVector;
};

template<typename Container>
SphericalShellIntersection<Container>::
SphericalShellIntersection(const Container &container,
                           const Point &center,
                           double radiusInner,
                           const RealVector &dirVector) {
    myContainer = new Container(container);
    myCenter = center;
    myRadiusInner = radiusInner;
    myRadiusOuter = myRadiusInner + 1;
    myDirectionVector = dirVector;
}

template<typename Container>
SphericalShellIntersection<Container>::
SphericalShellIntersection(const Container &container,
                           const Point &center,
                           double radiusInner,
                           double radiusOuter,
                           const RealVector &dirVector) {
    myContainer = new Container(container);
    myCenter = center;
    myRadiusInner = radiusInner;
    myRadiusOuter = radiusOuter;
    myDirectionVector = dirVector;
}

template<typename Container>
SphericalShellIntersection<Container>::
SphericalShellIntersection(const SphericalShellIntersection &other) {
    myContainer = new Container(*other.myContainer);
    myCenter = other.myCenter;
    myRadiusInner = other.myRadiusInner;
    myRadiusOuter = other.myRadiusOuter;
    myDirectionVector = other.myDirectionVector;

}

template<typename Container>
SphericalShellIntersection<Container>::
~SphericalShellIntersection() {
    if (myContainer != 0) {
        delete myContainer;
        myContainer = 0;
    }

}

template<typename Container>
Container
SphericalShellIntersection<Container>::
shell() {
    typedef DGtal::BreadthFirstVisitor<ObjectType, std::set<Point> > Visitor;
    typedef typename Visitor::Node Node;

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    ObjectType obj(dt26_6, *myContainer);
    L2Metric l2Metric;
    DistanceFunctor functor(l2Metric, myCenter);
    Visitor visitor(obj, myCenter);
    Container shell(myContainer->domain());
    int radius = (int) myRadiusInner;
    while (!visitor.finished()) {
        Node node = visitor.current();
        if (node.second > (int) myRadiusOuter) break;
        bool add = (node.second >= myRadiusInner);
        if (myDirectionVector != RealVector::zero) {
            DigitalPlane<Space> digPlane(myCenter, myDirectionVector);
            add &= digPlane.isPointAbove(node.first);
        }
        if (add)
            shell.insert(node.first);
        visitor.expand();

    }
    return shell;
}

template<typename Container>
unsigned int
SphericalShellIntersection<Container>::
degree(const Container &shell, double minCCSize) {
    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    ObjectType objectImage(dt26_6, shell);
    std::vector<ObjectType> objects;
    std::back_insert_iterator<std::vector<ObjectType> > inserter(objects);
    unsigned int nbConnectedComponents = objectImage.writeComponents(inserter);
    unsigned int cpt = 0;

    for (const auto &obj : objects) {
        if (obj.size() > minCCSize) //noise
            cpt++;
    }
    return cpt;
}

template<typename Container>
unsigned int
SphericalShellIntersection<Container>::
degree(double minCCSize) {
    Container s = shell();
    return degree(s, minCCSize);
}

#endif
