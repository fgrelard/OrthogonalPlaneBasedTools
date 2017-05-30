#ifndef CUTTING_PLANE_ESTIMATOR_H
#define CUTTING_PLANE_ESTIMATOR_H

#include <vector>
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "geometry/CurveProcessor.h"

template<typename Plane>
class CuttingPlaneEstimator {
public:
    typedef typename Plane::DigitalSet Container;
    typedef typename Plane::Point Point;
    typedef typename Plane::Vector RealVector;
    typedef typename Plane::Domain Domain;
    typedef typename Plane::Space Space;
    typedef typename DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;
    typedef DGtal::DistanceTransformation<Space, Container, L2Metric> DTL2;

public:
    CuttingPlaneEstimator();

    CuttingPlaneEstimator(const std::vector<Container> &branches,
                          const Container &setVolume);

    CuttingPlaneEstimator(const CuttingPlaneEstimator &other);

    ~CuttingPlaneEstimator();

public:
    template<typename PlaneEstimator>
    std::vector<Plane> cuttingPlanes(PlaneEstimator &planeEstimator,
                                     const DTL2 &dt);

    std::vector<Plane> cuttingPlanes(const std::vector<Plane> &planes);

public:

    template<typename PlaneEstimator>
    Plane referencePlane(const Point &referencePoint,
                         PlaneEstimator &planeEstimator,
                         const DTL2 &dt);

    Plane referencePlane(const Point &referencePoint,
                         const std::vector<Plane> &planes);

    template<typename PlaneEstimator>
    std::vector<Plane> computePlanes(const std::vector<Point> &orientedEdge,
                                     PlaneEstimator &planeEstimator,
                                     const DTL2 &dt);

    std::vector<Plane> extractPlanes(const std::vector<Point> &orientedEdge,
                                     const std::vector<Plane> &planes);

    Plane cuttingPlane(const std::vector<Plane> &planes);


    std::vector<Plane> filteredPlanes(const std::vector<Plane> &cuttingPlanes,
                                      const std::vector<Plane> &endPlanes,
                                      const Plane &referencePlane);

    Point extractBranchingPoint();

private:
    std::vector<Container> *myBranches;
    Container *myVolume;

};

template<typename Plane>
CuttingPlaneEstimator<Plane>::
CuttingPlaneEstimator() {
    myBranches = 0;
    myVolume = 0;
}

template<typename Plane>
CuttingPlaneEstimator<Plane>::
CuttingPlaneEstimator(const std::vector<Container> &branches,
                      const Container &setVolume) {
    myBranches = new std::vector<Container>(branches);
    myVolume = new Container(setVolume);
}

template<typename Plane>
CuttingPlaneEstimator<Plane>::
CuttingPlaneEstimator(const CuttingPlaneEstimator &other) {
    myBranches = new std::vector<Container>(*other.myBranches);
    myVolume = new Container(*other.myContainer);
}

template<typename Plane>
CuttingPlaneEstimator<Plane>::
~CuttingPlaneEstimator() {
    if (myBranches) {
        delete myBranches;
        myBranches = 0;
    }
    if (myVolume) {
        delete myVolume;
        myVolume = 0;
    }

}

template<typename Plane>
template<typename PlaneEstimator>
std::vector<Plane>
CuttingPlaneEstimator<Plane>::
cuttingPlanes(PlaneEstimator &planeEstimator, const DTL2 &dt) {
    std::vector<Plane> cuttingPlanes;
    std::vector<Plane> endPlanes;
    Point b = extractBranchingPoint();

    for (const Container &branch : *myBranches) {
        std::vector<Point> orderedBranch = CurveProcessor<Container>(branch).convertToOrderedCurve(b);
        std::vector<Plane> planes = computePlanes(orderedBranch, planeEstimator, dt);
        Plane cutting = cuttingPlane(planes);
        Plane end = *(planes.rbegin());

        cuttingPlanes.push_back(cutting);
        endPlanes.push_back(end);
    }
    Plane planeReference = referencePlane(b, planeEstimator, dt);
    cuttingPlanes = filteredPlanes(cuttingPlanes, endPlanes, planeReference);
    return cuttingPlanes;

}

template<typename Plane>
std::vector<Plane>
CuttingPlaneEstimator<Plane>::
cuttingPlanes(const std::vector<Plane> &branchesToPlanes) {
    std::vector<Plane> cuttingPlanes;
    std::vector<Plane> endPlanes;
    Point b = extractBranchingPoint();

    for (const Container &branch : *myBranches) {
        std::vector<Point> orderedBranch = CurveProcessor<Container>(branch).convertToOrderedCurve(b);
        std::vector<Plane> planes = extractPlanes(orderedBranch, branchesToPlanes);
        Plane cutting = cuttingPlane(planes);
        Plane end = *(planes.rbegin());

        cuttingPlanes.push_back(cutting);
        endPlanes.push_back(end);
    }
    Plane planeReference = referencePlane(b, branchesToPlanes);
    cuttingPlanes = filteredPlanes(cuttingPlanes, endPlanes, planeReference);
    return cuttingPlanes;

}

template<typename Plane>
template<typename PlaneEstimator>
Plane
CuttingPlaneEstimator<Plane>::
referencePlane(const Point &referencePoint,
               PlaneEstimator &planeEstimator,
               const DTL2 &dt) {
    double radius = dt(referencePoint) + 2;
    planeEstimator.setRadius(radius);
    return planeEstimator.convergentPlaneAt(referencePoint, *myVolume, radius * 5);
}

template<typename Plane>
Plane
CuttingPlaneEstimator<Plane>::
referencePlane(const Point &referencePoint,
               const std::vector<Plane> &branchesToPlanes) {
    auto iterator = std::find_if(branchesToPlanes.begin(),
                                 branchesToPlanes.end(),
                                 [&](const Plane &plane) {
                                     return (plane.getCenter() == referencePoint);
                                 });
    if (iterator != branchesToPlanes.end()) {
        return *iterator;
    }
    return Plane();
}


template<typename Plane>
template<typename PlaneEstimator>
std::vector<Plane>
CuttingPlaneEstimator<Plane>::
computePlanes(const std::vector<Point> &orientedEdge,
              PlaneEstimator &planeEstimator,
              const DTL2 &dt) {
    std::vector<Plane> planes;
    for (const Point &p : orientedEdge) {
        double radius = dt(p) + 2;
        planeEstimator.setRadius(radius);
        Plane plane = planeEstimator.convergentPlaneAt(p, *myVolume, radius * 5);
        planes.push_back(plane);

    }
    return planes;
}


template<typename Plane>
std::vector<Plane>
CuttingPlaneEstimator<Plane>::
extractPlanes(const std::vector<Point> &orientedEdge,
              const std::vector<Plane> &branchesToPlane) {
    std::vector<Plane> planes;
    for (const Point &p : orientedEdge) {
        Plane plane = referencePlane(p, branchesToPlane);
        //if (plane == Plane())
        planes.push_back(plane);
    }
    return planes;
}

template<typename Plane>
Plane
CuttingPlaneEstimator<Plane>::
cuttingPlane(const std::vector<Plane> &planes) {
    typedef DGtal::MetricAdjacency<Space, 3> MAdj;

    Plane candidate;
    double previousFactor = 0;
    for (const Plane &plane : planes) {
        Point p = plane.getCenter();
        Container current = plane.intersectionWithSetOneCC(*myVolume);
        double currentValue = current.size();
        std::vector<Point> neighbors;
        std::back_insert_iterator<std::vector<Point>> inserter(neighbors);
        MAdj::writeNeighbors(inserter, p);
        for (const Point &n : neighbors) {
            auto nInMap = find_if(planes.begin(), planes.end(), [&](const Plane &planeN) {
                return planeN.getCenter() == n;
            });
            if (nInMap != planes.end()) {
                double valueNeighbor = nInMap->intersectionWithSetOneCC(*myVolume).size();
                double factor = valueNeighbor / currentValue;
                if (factor > previousFactor) {
                    previousFactor = factor;
                    candidate = plane;
                }
            }
        }
    }
    return candidate;

}


template<typename Plane>
std::vector<Plane>
CuttingPlaneEstimator<Plane>::
filteredPlanes(const std::vector<Plane> &cuttingPlanes, const std::vector<Plane> &endPlanes,
               const Plane &referencePlane) {

    std::vector<Plane> first, second;

    for (int i = 0, e = endPlanes.size(); i < e; i++) {
        Plane endPlane = endPlanes[i];
        Plane initPlane = cuttingPlanes[i];
        if (referencePlane.contains(endPlane.getCenter())) continue;
        if (referencePlane.isPointAbove(endPlane.getCenter())) {
            first.push_back(initPlane);
        } else {
            second.push_back(initPlane);
        }
    }
    return (second.size() > first.size()) ? second : first;
}

template<typename Plane>
typename CuttingPlaneEstimator<Plane>::Point
CuttingPlaneEstimator<Plane>::
extractBranchingPoint() {
    for (int i = 0; i < myBranches->size(); i++) {
        Container setCurrent = (*myBranches)[i];
        for (int j = i + 1; j < myBranches->size(); j++) {
            Container setOther = (*myBranches)[j];
            Container intersection = SetProcessor<Container>(setCurrent).intersection(setOther);
            if (intersection.size() == 1)
                return *intersection.begin();
        }
    }
    return Point::zero;
}


#endif
