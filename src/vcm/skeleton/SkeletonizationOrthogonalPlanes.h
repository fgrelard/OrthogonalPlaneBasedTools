#ifndef SKELETONIZATION_ORTHOGONAL_PLANES_H
#define SKELETONIZATION_ORTHOGONAL_PLANES_H

/**
 * @file SkeletonizationOrthogonalPlanes.h
 * @author Florent Grelard (florent.grelard@labri.fr)
 * LaBRI, Bordeaux University
 *
 * @date 2017/01/20
 *
 */

#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "geometry/DistanceToPointFunctor.h"
#include "geometry/MedialAxis.h"
#include "shapes/DigitalPlaneSet.h"
#include "geometry/Distance.h"
#include "geometry/SetProcessor.h"
#include "geometry/CurveProcessor.h"
#include "geometry/SphericalShellIntersection.h"
#include "vcm/OrthogonalPlaneEstimator.h"
#include "geometry/junction/SSIJunctionDetection.h"
#include "vcm/skeleton/post/JunctionProcessingSkeleton.h"
#include "geometry/predicate/AbovePlanePredicate.h"
#include "ShapeDescriptor.h"

/**
   * Description of template class 'SkeletonizationOrthogonalPlanes' <p>
   * \brief Aim: This class aims at computing the skeleton of a
   * given tubular volume through Orthogonal plane estimation.
   * The skeleton points are given as the centers of mass of the
   * intersection between the orthogonal planes and the volume.
   * A tracking procedure ensures reasonable computation times.
   *
   * You may obtain the skeleton by \ref skeletonize.
   *
   *
   * @tparam Container type of Digital Set (model of CDigitalSet).
   *
   * @tparam JunctionDetection type allowing to detect junctions at
   * tracked points. For instance, SSIJunctionDetection and
   * NoJunctionDetection are models of this type.
   *
   * @tparam PostProcessing type allowing to post process the skeleton
   * Typically, necessary for junctions. For instance,
   * JunctionProcessingSkeleton or NoPostProcessingSkeleton are
   * models of this type
   *
   */

template<typename Container,
        typename JunctionDetection = SSIJunctionDetection<Container>,
        typename PostProcessing = JunctionProcessingSkeleton<Container, AbovePlanePredicate<typename Container::Space> > >
class SkeletonizationOrthogonalPlanes {
public:
    typedef typename Container::Space Space;
    typedef typename Space::Point Point;
    typedef typename Space::RealPoint RealPoint;
    typedef typename Space::RealVector RealVector;
    typedef typename DGtal::functors::BallConstantPointFunction<Point, double> KernelFunction;
    typedef OrthogonalPlaneEstimator<Container, KernelFunction> PlaneEstimator;
    typedef DGtal::MetricAdjacency<Space, 1> Adj6;
    typedef DGtal::MetricAdjacency<Space, 3> Adj26;
    typedef DGtal::DigitalTopology<Adj26, Adj6> DT26_6;
    typedef DGtal::Object<DT26_6, Container> ObjectType;
    typedef DigitalPlane<Space> Plane;
    typedef DigitalPlaneSet<Space> PlaneSet;
    typedef typename PlaneEstimator::L2Metric L2Metric;
    typedef DGtal::DistanceTransformation<Space, Container, L2Metric> DTL2;

public:
    SkeletonizationOrthogonalPlanes() = delete;

    /**
     * Constructor.
     *
     * @param setVolume the input volume on which to compute the
     * skeleton
     *
     * @param junctionDetection the junction detector
     *
     * @param R the offset radius for the set of points. Voronoi cells
     * are intersected with this offset. The unit corresponds to a step in the digital space.
     */
    SkeletonizationOrthogonalPlanes(const Container &setVolume,
                                    const JunctionDetection &junctionDetection,
                                    double R = 10);

    SkeletonizationOrthogonalPlanes(const SkeletonizationOrthogonalPlanes &other);

    ~SkeletonizationOrthogonalPlanes();

    // ----------------------- Interface  --------------------------------------
public:
    /// @return the skeleton
    Container skeletonize();

    // ----------------------- Internal methods --------------------------------------
public:
    Point trackNextPoint(const PlaneSet &plane);

    void markPoints(const Point &point);

    void markPoints(const Container &points);

    void markPointsBetweenPlanes(const PlaneSet &currentPlane,
                                 const PlaneSet &previousPlane,
                                 double distanceMax = std::numeric_limits<double>::max());

    Container filterIsolatedPoints(int minSize = 1);

    std::vector<Plane> orientEndPoints();

private:
    Container restrictPlaneSet(const PlaneSet &planeSet, double radius);

    Plane orientPlane(const Plane &undirectedPlane,
                      const Container &endPoints);

    bool isInJunction(const PlaneSet &p, double radius);

// ------------------------- Private Datas --------------------------------

private:
    Container *myVolume;
    PlaneEstimator *myPlaneEstimator;
    JunctionDetection *myJunctionDetection;
    double myBigR;


// ------------------------- Internals --------------------------------

private:
    Container *mySkeleton;
    Container *myMarkedVertices;
    DTL2 *myDT;
    Container *myMedialAxis;
    std::vector<Plane> myPlanes;
};


template<typename Container, typename JunctionDetection, typename PostProcessing>
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
SkeletonizationOrthogonalPlanes(const Container &setVolume,
                                const JunctionDetection &junctionDetection,
                                double R) {
    L2Metric l2Metric;
    myVolume = new Container(setVolume);
    myJunctionDetection = new JunctionDetection(junctionDetection);
    double r = 5;
    myBigR = R;
    KernelFunction chi(1.0, r);
    myPlaneEstimator = new PlaneEstimator(*myVolume, chi, myBigR, r);

    mySkeleton = new Container(myVolume->domain());
    myMarkedVertices = new Container(setVolume.domain());
    myDT = new DTL2(myVolume->domain(), *myVolume, l2Metric);
    MedialAxis<Container> ma(*myVolume);
    myMedialAxis = new Container(ma.pointSet());
}

template<typename Container, typename JunctionDetection, typename PostProcessing>
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
SkeletonizationOrthogonalPlanes(const SkeletonizationOrthogonalPlanes &other) {
    myVolume = new Container(*other.myVolume);
    myPlaneEstimator = new PlaneEstimator(*other.myPlaneEstimator);
    myBigR = other.myBigR;
    myJunctionDetection = new JunctionDetection(*other.myJunctionDetection);

    mySkeleton = new Container(*other.mySkeleton);
    myMarkedVertices = new Container(*other.myMarkedVertices);
    myDT = new DTL2(*other.myDT);
    myMedialAxis = new Container(*other.myMedialAxis);
    myPlanes = other.myPlanes;
}

template<typename Container, typename JunctionDetection, typename PostProcessing>
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
~SkeletonizationOrthogonalPlanes() {
    if (myVolume) {
        delete myVolume;
        myVolume = 0;
    }
    if (myPlaneEstimator) {
        delete myPlaneEstimator;
        myPlaneEstimator = 0;
    }
    if (myJunctionDetection) {
        delete myJunctionDetection;
        myJunctionDetection = 0;
    }
    if (mySkeleton) {
        delete mySkeleton;
        mySkeleton = 0;
    }
    if (myMarkedVertices) {
        delete myMarkedVertices;
        myMarkedVertices = 0;
    }
    if (myDT) {
        delete myDT;
        myDT = 0;
    }
    if (myMedialAxis) {
        delete myMedialAxis;
        myMedialAxis = 0;
    }
}

template<typename Container, typename JunctionDetection, typename PostProcessing>
Container
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
skeletonize() {
    PlaneSet previous;
    L2Metric l2Metric;
    size_t nbPoints = myVolume->size();
    SetProcessor<Container> procMA(*myMedialAxis);
    Container a3ShellPoints(myVolume->domain());
    Point p = *(max_element(myVolume->begin(), myVolume->end(), [&](const Point &one,
                                                                    const Point &two) {
        return ((*myDT)(one) < (*myDT)(two));
    }));
    double distanceMax = (*myDT)(p) + 2.0;

    DGtal::trace.beginBlock("Computing skeleton");
    while (myMarkedVertices->size() < nbPoints) {
        DGtal::trace.progressBar(myMarkedVertices->size(),
                                 nbPoints);
        Point closestPoint = procMA.closestPointAt(p);
        double radius = (*myDT)(closestPoint) + 2.0;
        myPlaneEstimator->setRadius(radius);
        Plane plane = myPlaneEstimator->convergentPlaneAt(p, *myVolume, distanceMax);
        radius = myPlaneEstimator->getRadius();
        Container planePoints = plane.intersectionWithSetOneCC(*myVolume);
        PlaneSet planeSet(plane, planePoints);
        PlaneSet planeSetG = previous;
        markPoints(p);
        markPoints(restrictPlaneSet(planeSet, radius));

        RealPoint centerOfMass = ShapeDescriptor<Container>(planePoints).extractCenterOfMass();

        if (centerOfMass != RealPoint::zero && Distance::euclideanDistance(centerOfMass, RealPoint(p)) <= sqrt(3)) {
            Point g = SetProcessor<Container>(planePoints).closestPointAt(centerOfMass);
            RealVector normal = plane.getPlaneEquation().normal();
            int connexity = plane.getConnexity();
            Plane planeG(g, normal, connexity);
            planeSetG = PlaneSet(planeG, planePoints);

            if (previous.isDefined() &&
                l2Metric(plane.getCenter(), previous.digitalPlane().getCenter()) <= 2 * sqrt(3)) {
                markPointsBetweenPlanes(planeSetG, previous, radius);
            }
            if (isInJunction(planeSetG, radius))
                a3ShellPoints.insert(g);

            else if (CurveProcessor<Container>(*mySkeleton).isPointThin(g)) {
                mySkeleton->insert(g);
                myPlanes.push_back(planeG);
                previous = planeSetG;
            }

        }
        p = trackNextPoint(planeSetG);
    }
    DGtal::trace.endBlock();

    DGtal::trace.beginBlock("PostProcessing");
    Container fillHoles = CurveProcessor<Container>(*mySkeleton).fillHoles(sqrt(3), 2 * sqrt(3));
    delete mySkeleton;
    mySkeleton = new Container(fillHoles);

    Container filteredSkeleton = filterIsolatedPoints();
    filteredSkeleton = CurveProcessor<Container>(filteredSkeleton).fillHolesNotInSet(a3ShellPoints, *myVolume);
    delete mySkeleton;
    mySkeleton = new Container(filteredSkeleton);

    myPlanes = orientEndPoints();

    PostProcessing algo(*mySkeleton, a3ShellPoints, *myVolume, myPlanes);
    Container postProcessedSkeleton = algo.postProcess();
    postProcessedSkeleton = CurveProcessor<Container>(postProcessedSkeleton).fillHoles(*myVolume);

    delete mySkeleton;
    mySkeleton = new Container(postProcessedSkeleton);
    DGtal::trace.endBlock();
    return *mySkeleton;
}

template<typename Container, typename JunctionDetection, typename PostProcessing>
typename SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::Point
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
trackNextPoint(const PlaneSet &planeSet) {
    const Point center = planeSet.digitalPlane().getCenter();
    const RealVector direction = planeSet.digitalPlane().getPlaneEquation().normal();
    const Container set = planeSet.pointSet();

    Point current = center;
    double scalar = 1.0;
    while (current == center ||
           set.find(current) != set.end()) {
        current = Point(RealPoint(center) + direction * scalar);
        scalar += 0.5;
    }
    if (myVolume->find(current) == myVolume->end() ||
        myMarkedVertices->find(current) != myMarkedVertices->end()) {
        scalar = 1.0;
        current = center;
        while (current == center ||
               set.find(current) != set.end()) {
            current = Point(RealPoint(center - direction * scalar));
            scalar += 0.5;
        }
        if (myVolume->find(current) == myVolume->end() ||
            myMarkedVertices->find(current) != myMarkedVertices->end()) {
            double distance = 0;
            for (const Point &v : *myVolume) {
                if (myMarkedVertices->find(v) == myMarkedVertices->end()) {
                    double currentD = (*myDT)(v);
                    if (currentD > distance) {
                        distance = currentD;
                        current = v;
                    }
                }
            }
        }
    }
    return current;
}

template<typename Container, typename JunctionDetection, typename PostProcessing>
void
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
markPoints(const Point &point) {
    myMarkedVertices->insert(point);
}

template<typename Container, typename JunctionDetection, typename PostProcessing>
void
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
markPoints(const Container &points) {
    myMarkedVertices->insert(points.begin(), points.end());
}

template<typename Container, typename JunctionDetection, typename PostProcessing>
void
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
markPointsBetweenPlanes(const PlaneSet &current,
                        const PlaneSet &previous,
                        double distanceMax) {

    L2Metric l2Metric;
    Container difference(myVolume->domain());
    Plane currentPlane = current.digitalPlane();
    Plane previousPlane = previous.digitalPlane();
    RealVector dirCurrent = currentPlane.getPlaneEquation().normal();
    RealVector dirPrevious = previousPlane.getPlaneEquation().normal();
    Point pCurrent = currentPlane.getCenter();
    Point pPrevious = previousPlane.getCenter();

    if (dirPrevious.dot(dirCurrent) > 0) {
        dirCurrent = -dirCurrent;
    }

    currentPlane = Plane(pCurrent, dirCurrent, currentPlane.getConnexity());

    for (const Point &p : *myVolume) {
        if (l2Metric(p, pCurrent) > distanceMax) continue;
        if (currentPlane.isPointAbove(p) &&
            previousPlane.isPointAbove(p))
            difference.insert(p);
    }
    myMarkedVertices->insert(difference.begin(), difference.end());
}

template<typename Container, typename JunctionDetection, typename PostProcessing>
Container
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
filterIsolatedPoints(int minSize) {
    Container filteredSkeleton(mySkeleton->domain());
    SetProcessor<Container> setProc(*mySkeleton);
    std::vector<Container> components = setProc.toConnectedComponents();
    for (const Container &cc : components) {
        if (cc.size() > minSize) {
            filteredSkeleton.insert(cc.begin(), cc.end());
        }
    }
    return filteredSkeleton;
}

template<typename Container, typename JunctionDetection, typename PostProcessing>
std::vector<typename SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::Plane>
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
orientEndPoints() {
    std::vector<Plane> planeEndPoints;
    SetProcessor<Container> setProc(*mySkeleton);
    std::vector<Container> components = setProc.toConnectedComponents();
    L2Metric l2Metric;
    for (const Container &cc : components) {
        if (cc.size() < 2) continue;
        CurveProcessor<Container> curveProc(cc);
        Container ccConnectivity = curveProc.ensureConnectivity();
        CurveProcessor<Container> connProc(ccConnectivity);
        Container branching = connProc.branchingPoints();
        Container localEnd = connProc.endPoints();
        SetProcessor<Container> setProc(cc);
        for (const Point &p : localEnd) {
            if (setProc.intersectionNeighborhoodAt(p, branching).size() > 0) continue;
            Plane plane = *(std::min_element(
                    myPlanes.begin(),
                    myPlanes.end(), [&](const Plane &plane, const Plane &otherPlane) {
                        return l2Metric(plane.getCenter(), p) < l2Metric(otherPlane.getCenter(), p);
                    }));
            RealVector normal = plane.getPlaneEquation().normal();
            plane = Plane(p, normal, plane.getConnexity());
            plane = orientPlane(plane, localEnd);
            planeEndPoints.push_back(plane);
        }
    }
    return planeEndPoints;
}

template<typename Container, typename JunctionDetection, typename PostProcessing>
Container
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
restrictPlaneSet(const PlaneSet &planeSet, double radius) {
    L2Metric l2Metric;
    Container set = planeSet.pointSet();
    Point center = planeSet.digitalPlane().getCenter();
    Container toMark(set.domain());
    for (const Point &p : set) {
        if (l2Metric(p, center) <= radius)
            toMark.insert(p);
    }
    return toMark;
}

template<typename Container, typename JunctionDetection, typename PostProcessing>
typename SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::Plane
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
orientPlane(const Plane &undirectedPlane, const Container &endPoints) {

    L2Metric l2Metric;
    Point center = undirectedPlane.getCenter();
    RealVector normal = undirectedPlane.getPlaneEquation().normal();
    Point candidate = *(std::max_element(endPoints.begin(), endPoints.end(), [&](const Point &one, const Point &two) {
        return (l2Metric(one, center) < l2Metric(two, center));
    }));
    RealVector dir = (center - candidate).getNormalized();
    normal = (normal.dot(dir) < 0) ? -normal : normal;

    for (typename RealVector::Dimension i = 0;
         i < RealVector::dimension; ++i) {
        normal[i] = (std::fabs(normal[i] - 0.0) < std::numeric_limits<double>::epsilon()) ? 0.0 : normal[i];
    }
    Plane plane(center, normal, undirectedPlane.getConnexity());
    return plane;
}


template<typename Container, typename JunctionDetection, typename PostProcessing>
bool
SkeletonizationOrthogonalPlanes<Container, JunctionDetection, PostProcessing>::
isInJunction(const PlaneSet &planeSet, double radius) {

    Point center = planeSet.digitalPlane().getCenter();
    RealVector normal = planeSet.digitalPlane().getPlaneEquation().normal();
    int connexity = planeSet.digitalPlane().getConnexity();
    Container set = planeSet.pointSet();
    Container minusSet = Plane(center, -normal, connexity).intersectionWithSetOneCC(*myVolume);
    double radiusCurrent = ShapeDescriptor<Container>(set).lengthMajorAxis() + 2.0;
    double radiusCurrentMinus = ShapeDescriptor<Container>(minusSet).lengthMajorAxis() + 2.0;
    double radiusShell = std::max(4.0, std::max(radiusCurrent, radiusCurrentMinus));
    radiusShell *= 1.2;
    double noise = radiusShell / 2.0;
    if (myJunctionDetection->isInJunction(center, radiusShell, noise)) {
        return true;
    }
    return false;
}

#endif
