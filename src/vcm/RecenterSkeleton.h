#ifndef RECENTER_SKELETON_H
#define RECENTER_SKELETON_H


/**
 * @file RecenterSkeleton.h
 * @author Florent Grelard (florent.grelard@labri.fr)
 * LaBRI, Bordeaux University
 *
 * @date 2017/01/20
 *
 */

#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/kernel/Point2ScalarFunctors.h"
#include "vcm/OrthogonalPlaneEstimator.h"
#include "vcm/CuttingPlaneEstimator.h"
#include "geometry/CurveDecomposition.h"
#include "ShapeDescriptor.h"

/**
   * Description of template class 'RecenterSkeleton' <p>
   * \brief Aim: This class aims at recentering points
   * in an existing skeleton. It is based on cutting plane estimation
   * around branching points of the skeleton. Cutting planes are then
   * realigned to delineate subvolumes.
   *
   * You may obtain the recentered skeleton by using the
   * method \ref recenter.
   *
   *
   * @tparam Container type of Digital Set (model of CDigitalSet).
   *
   */
template<typename Container>
class RecenterSkeleton {
public:
    typedef typename Container::Point Point;
    typedef typename Container::Space Space;
    typedef typename Space::RealVector RealVector;
    typedef typename DGtal::functors::BallConstantPointFunction<Point, double> KernelFunction;
    typedef typename DGtal::HyperRectDomain<Space> Domain;
    typedef OrthogonalPlaneEstimator<Container, KernelFunction> PlaneEstimator;
    typedef typename PlaneEstimator::Plane Plane;
    typedef typename PlaneEstimator::L2Metric L2Metric;
    typedef DGtal::DistanceTransformation<Space, Container, L2Metric> DTL2;
    typedef Edge<Container> GraphEdge;
    typedef DGtal::MetricAdjacency<Space, 1> Adj6;
    typedef DGtal::MetricAdjacency<Space, 3> Adj26;
    typedef DGtal::DigitalTopology<Adj26, Adj6> DT26_6;
    typedef DGtal::Object<DT26_6, Container> ObjectType;
    typedef CuttingPlaneEstimator<Plane> CuttingPlane;


public:
    RecenterSkeleton() = delete;

    /**
     * Constructor.
     *
     * @param skeleton the input non centered skeleton
     *
     * @param volume the corresponding volume
     *
     */
    RecenterSkeleton(const Container &skeleton,
                     const Container &volume);

    ~RecenterSkeleton();

    RecenterSkeleton(const RecenterSkeleton &other);


// ----------------------- Interface  --------------------------------------
public:
    /// @return the recentered skeleton
    Container recenter();


    // ----------------------- Internal methods --------------------------------------
public:

    std::vector<Plane> computePlanes();

    std::vector<Container> adjacentBranchesToPoint(const Point &p,
                                                   const std::vector<GraphEdge> &graph);

    std::vector<Container> restrictBranches(const std::vector<Container> &branches, const Point &p);

    Container referenceBranch(const std::vector<Container> &branches, const std::vector<Container> &subSetBranches);

    Container planeToBranch(const Plane &plane, const std::vector<Container> &branches);

    std::vector<Container> planesToBranches(const std::vector<Plane> &planes, const std::vector<Container> &branches);

    std::vector<Plane> alignPlanes(const Plane &currentPlane,
                                   const std::vector<Plane> &planes);

    std::vector<Plane> orientNormalPlanes(const std::vector<Plane> &planes,
                                          const Point &p);


    std::vector<Point> orderedBranchToRecenter(const Container &currentBranch,
                                               const Container &referenceBranch,
                                               const Container &restrictedVolume,
                                               const Point &branchingPoint);

    Container subVolume(const Container &volume, const Plane &plane);

    Container subVolume(const Container &volume, const std::vector<Plane> &planes);

    Container relevantCC(const Container &subVolume, const Container &constraint);

    template<typename OtherContainer>
    Container recenterSkeletonPoints(const Container &subVolume,
                                     const OtherContainer &existingSkeleton,
                                     const Container &computedSkeleton);

    Container postProcess(const Container &existingSkeleton,
                          const Container &processedEdges);

    // ------------------------- Private Datas --------------------------------
private:
    Container *mySkeleton;
    Container *myVolume;


    // ----------------------- Internals  --------------------------------------
private:
    PlaneEstimator *myPlaneEstimator;
    DTL2 *myDT;

};

template<typename Container>
RecenterSkeleton<Container>::
RecenterSkeleton(const Container &skeleton,
                 const Container &volume) {
    mySkeleton = new Container(skeleton);
    myVolume = new Container(volume);
    L2Metric l2Metric;
    myDT = new DTL2(myVolume->domain(), *myVolume, l2Metric);
    double r = 10;
    double R = 10;
    KernelFunction chi(1.0, r);
    int connexity = 6;
    myPlaneEstimator = new PlaneEstimator(volume, chi, R, r, connexity);
}

template<typename Container>
RecenterSkeleton<Container>::
~RecenterSkeleton() {
    if (mySkeleton) {
        delete mySkeleton;
        mySkeleton = 0;
    }
    if (myVolume) {
        delete myVolume;
        myVolume = 0;
    }
    if (myDT) {
        delete myDT;
        myDT = 0;
    }
    if (myPlaneEstimator) {
        delete myPlaneEstimator;
        myPlaneEstimator = 0;
    }
}

template<typename Container>
RecenterSkeleton<Container>::
RecenterSkeleton(const RecenterSkeleton &other) {
    mySkeleton = new Container(*other.mySkeleton);
    myVolume = new Container(*other.myVolume);
    myPlaneEstimator = new PlaneEstimator(*other.myPlaneEstimator);
    myDT = new DTL2(*other.myDT);
}

template<typename Container>
Container
RecenterSkeleton<Container>::
recenter() {

    Container branching = CurveProcessor<Container>(*mySkeleton).branchingPoints();
    std::vector<GraphEdge> graph = CurveDecomposition<Container>(*mySkeleton, branching).branchDecomposition();

    DGtal::trace.beginBlock("Recentering");
    std::vector<Plane> planes = computePlanes();
    Container skeletonPoints(mySkeleton->domain());
    Container processedEdges(mySkeleton->domain());
//#pragma omp parallel for schedule(dynamic) shared(skeletonPoints, processedEdges)
    for (size_t i = 0; i < branching.size(); i++) {
        auto begin = branching.begin();
        std::advance(begin, i);
        Point b = *begin;
        std::vector<Container> adjacentEdges = adjacentBranchesToPoint(b, graph);
        double radius = max_element(adjacentEdges.begin(), adjacentEdges.end(), [&](const Container &e1,
                                                                                    const Container &e2) {
            return e1.size() < e2.size();
        })->size() * 0.5;

        std::vector<Container> restricted = restrictBranches(adjacentEdges, b);
        CuttingPlane cuttingPE(restricted, *myVolume);
        std::vector<Plane> cuttingPlanes = cuttingPE.cuttingPlanes(planes);

        std::vector<Plane> orientedPlanes = orientNormalPlanes(cuttingPlanes, b);
        if (cuttingPlanes.size() < 2) continue;
        std::vector<Container> branchesRecentering = planesToBranches(cuttingPlanes, restricted);
        Container refBranch = referenceBranch(restricted, branchesRecentering);
        refBranch = SetProcessor<Container>(refBranch).subSet(b, refBranch.size() * 0.6);
        Container restrictedVolume = SetProcessor<Container>(*myVolume).subSet(b, radius * 1.5);
        size_t j = 0;
        for (const Plane &plane : orientedPlanes) {
            std::vector<Plane> alignedPlanes = alignPlanes(plane, orientedPlanes);
            alignedPlanes = orientNormalPlanes(alignedPlanes, b);
            Plane currentPlane(plane.getCenter(),
                               -plane.getPlaneEquation().normal(),
                               6);
            Container subVolume1 = subVolume(restrictedVolume, currentPlane);
            Container subVolume2 = subVolume(restrictedVolume, alignedPlanes);

            subVolume1.insert(subVolume2.begin(), subVolume2.end());
            subVolume1 = relevantCC(subVolume1, refBranch);

            Container correspondingBranch = branchesRecentering[j];
            correspondingBranch = SetProcessor<Container>(correspondingBranch).subSet(b,
                                                                                      correspondingBranch.size() * 0.6);

            if (refBranch.size() == 0 || correspondingBranch.size() == 0) break;
            std::vector<Point> pointsToRecenter = orderedBranchToRecenter(correspondingBranch, refBranch,
                                                                          restrictedVolume, b);
            Container recenteredPoints = recenterSkeletonPoints(subVolume1, pointsToRecenter, skeletonPoints);

            skeletonPoints.insert(recenteredPoints.begin(), recenteredPoints.end());
            processedEdges.insert(correspondingBranch.begin(), correspondingBranch.end());
            processedEdges.insert(refBranch.begin(), refBranch.end());
            j++;

        }
    }
    Container post = postProcess(*mySkeleton, processedEdges);
    skeletonPoints.insert(post.begin(), post.end());
    skeletonPoints = CurveProcessor<Container>(skeletonPoints).fillHoles(sqrt(3), 2 * sqrt(3));
    skeletonPoints = CurveProcessor<Container>(skeletonPoints).fillHoles(*myVolume);
    DGtal::trace.endBlock();
    return skeletonPoints;
}

template<typename Container>
std::vector<typename RecenterSkeleton<Container>::Plane>
RecenterSkeleton<Container>::
computePlanes() {
    std::vector<Plane> planes;
    for (const Point &p : *mySkeleton) {
        double radius = (*myDT)(p) + 2;
        myPlaneEstimator->setRadius(radius);
        Plane plane = myPlaneEstimator->convergentPlaneAt(p, *myVolume, radius * 5);
        planes.push_back(plane);

    }
    return planes;
}


template<typename Container>
std::vector<Container>
RecenterSkeleton<Container>::
adjacentBranchesToPoint(const Point &p,
                        const std::vector<GraphEdge> &graph) {
    std::vector<Container> adjacent;
    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    ObjectType objSkel(dt26_6, *mySkeleton);
    std::vector<Point> neigh;
    std::back_insert_iterator<std::vector<Point>> inserter(neigh);
    objSkel.writeNeighbors(inserter, p);

    for (const Container &edge : graph) {
        for (const Point &n : neigh) {
            if (edge.find(n) != edge.end()) {
                adjacent.push_back(edge);
                break;
            }
        }
    }
    return adjacent;
}


template<typename Container>
Container
RecenterSkeleton<Container>::
referenceBranch(const std::vector<Container> &branches, const std::vector<Container> &subSetBranches) {
    Container restrictEdgeB(myVolume->domain());
    for (const Container &branch : branches) {
        unsigned int cpt = 0;
        for (const Container &branchSub : subSetBranches) {
            if (!SetProcessor<Container>(branch).sameContainer(branchSub)) {
                cpt++;
            }
        }
        if (cpt == subSetBranches.size())
            restrictEdgeB = branch;
    }
    return restrictEdgeB;
}

template<typename Container>
Container
RecenterSkeleton<Container>::
planeToBranch(const Plane &plane,
              const std::vector<Container> &branches) {
    auto iterator = std::find_if(branches.begin(), branches.end(), [&](const Container &c) {
        return (c.find(plane.getCenter()) != c.end());
    });
    if (iterator != branches.end())
        return *iterator;
    return Container(Domain(Point::zero, Point::zero));
}

template<typename Container>
std::vector<Container>
RecenterSkeleton<Container>::
planesToBranches(const std::vector<Plane> &planes,
                 const std::vector<Container> &branches) {
    std::vector<Container> subBranches;
    for (const Plane &plane : planes) {
        Container branch = planeToBranch(plane, branches);
        subBranches.push_back(branch);
    }
    return subBranches;
}

template<typename Container>
std::vector<Container>
RecenterSkeleton<Container>::
restrictBranches(const std::vector<Container> &branches,
                 const Point &p) {
    std::vector<Container> restrictedBranches;
    for (const Container &branch : branches) {
        Container restricted = SetProcessor<Container>(branch).subSet(p, branch.size());
        if (restricted.size() > 0)
            restrictedBranches.push_back(restricted);
    }
    return restrictedBranches;
}

template<typename Container>
std::vector<typename RecenterSkeleton<Container>::Plane>
RecenterSkeleton<Container>::
alignPlanes(const Plane &currentPlane,
            const std::vector<Plane> &planes) {

    Point currentPoint = currentPlane.getCenter();
    Container currentSet = currentPlane.intersectionWithSetOneCC(*myVolume);

    std::vector<Container> otherEdges;
    std::vector<Plane> planesRotated;
    for (size_t j = 0, jend = planes.size(); j < jend; j++) {
        if (currentPlane == planes[j]) continue;
        Plane otherPlane = planes[j];
        Point otherPoint = otherPlane.getCenter();
        if (otherPoint == Point::zero || currentPoint == Point::zero) continue;
        Container otherSet = otherPlane.intersectionWithSetOneCC(*myVolume);

        std::pair<Point, Point> closestPointsInter = SetProcessor<Container>(currentSet).twoClosestPoints(otherSet);
        Point current = closestPointsInter.first;
        Point current2 = closestPointsInter.second;

        RealVector normal1 = currentPlane.getPlaneEquation().normal();
        RealVector normal2 = otherPlane.getPlaneEquation().normal();

        RealVector normalRot = normal1.crossProduct(normal2);
        RealVector normalPlane2 = normalRot.crossProduct(normal1);
        Plane newPlane2 = Plane(current2, normalPlane2, 6);
        planesRotated.push_back(newPlane2);
    }
    return planesRotated;
}


template<typename Container>
std::vector<typename RecenterSkeleton<Container>::Plane>
RecenterSkeleton<Container>::
orientNormalPlanes(const std::vector<Plane> &planes,
                   const Point &p) {
    std::vector<Plane> oriented;
    for (const Plane &plane : planes) {
        Point center = plane.getCenter();
        RealVector n = plane.getPlaneEquation().normal();
        RealVector dir = (p - center).getNormalized();
        n = (dir.dot(n) < 0) ? -n : n;
        Plane newp(center, n, 6);
        oriented.push_back(newp);
    }
    return oriented;
}

template<typename Container>
std::vector<typename RecenterSkeleton<Container>::Point>
RecenterSkeleton<Container>::orderedBranchToRecenter(const Container &currentBranch,
                                                     const Container &referenceBranch,
                                                     const Container &restrictedVolume,
                                                     const Point &branchingPoint) {

    Container eEdge = CurveProcessor<Container>(currentBranch).endPoints();
    Point candEdge = branchingPoint;
    for (const Point &p : eEdge) {
        if (p != branchingPoint)
            candEdge = p;
    }

    Container currentBranchInVol(myVolume->domain()), referenceBranchInVol(myVolume->domain());
    for (const Point &pResVol : restrictedVolume) {
        if (currentBranch.find(pResVol) != currentBranch.end())
            currentBranchInVol.insert(pResVol);
        if (referenceBranch.find(pResVol) != referenceBranch.end())
            referenceBranchInVol.insert(pResVol);
    }

    std::vector<Point> currentBranchOriented = CurveProcessor<Container>(currentBranchInVol).convertToOrderedCurve(
            candEdge);
    std::vector<Point> referenceBranchOriented = CurveProcessor<Container>(referenceBranchInVol).convertToOrderedCurve(
            branchingPoint);

    std::vector<Point> edgesVolume;
    edgesVolume.insert(edgesVolume.end(), currentBranchOriented.begin(), currentBranchOriented.end());
    edgesVolume.insert(edgesVolume.end(), referenceBranchOriented.begin(), referenceBranchOriented.end());
    return edgesVolume;
}

template<typename Container>
Container
RecenterSkeleton<Container>::
subVolume(const Container &restrictedVolume,
          const Plane &plane) {
    Container subVolume(restrictedVolume.domain());
    for (const Point &p : restrictedVolume) {
        if (plane.isPointAbove(p))
            subVolume.insert(p);
    }
    return subVolume;
}

template<typename Container>
Container
RecenterSkeleton<Container>::
subVolume(const Container &restrictedVolume,
          const std::vector<Plane> &cuttingPlanes) {
    Container subVolume(restrictedVolume.domain());
    for (const Point &p : restrictedVolume) {
        bool add = true;
        for (const Plane &digitalPlane : cuttingPlanes) {
            add &= digitalPlane.isPointAbove(p);
        }
        if (add && cuttingPlanes.size() > 0)
            subVolume.insert(p);
    }
    return subVolume;
}

template<typename Container>
Container
RecenterSkeleton<Container>::
relevantCC(const Container &subVolume,
           const Container &constraint) {
    Container vol(subVolume.domain());
    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    ObjectType obj(dt26_6, subVolume);
    std::vector<ObjectType> objects;
    std::back_insert_iterator<std::vector<ObjectType> > inserter(objects);
    unsigned int nbCC = obj.writeComponents(inserter);
    if (nbCC > 1) {
        Container candidate = subVolume;
        double number = 0;
        for (const ObjectType &o : objects) {
            Container set = o.pointSet();
            Container intersect = SetProcessor<Container>(set).intersection(constraint);
            if (intersect.size() > number) {
                number = intersect.size();
                candidate = set;
            }
        }
        return candidate;
    }
    return subVolume;
}

template<typename Container>
template<typename OtherContainer>
Container
RecenterSkeleton<Container>::
recenterSkeletonPoints(const Container &subVolume,
                       const OtherContainer &existingSkeleton,
                       const Container &computedSkeleton) {

    typedef typename Space::RealPoint RealPoint;
    Container setVCM(subVolume.domain());
    setVCM.insert(existingSkeleton.begin(), existingSkeleton.end());
    Container smoothSkeleton(subVolume.domain());
    if (existingSkeleton.size() <= 2) {
        smoothSkeleton.insert(existingSkeleton.begin(), existingSkeleton.end());
        return smoothSkeleton;
    }

    int i = 0;
    bool add = false;

    double radius = existingSkeleton.size() * 0.4;
    radius = (radius < 2) ? 2 : radius;
    KernelFunction chi(1.0, radius);
    PlaneEstimator planeEstimator(setVCM, chi, 20, radius, 6);
    L2Metric l2Metric;

    for (const Point &cp : existingSkeleton) {

        Point currentPoint = SetProcessor<Container>(subVolume).closestPointAt(cp);
        Plane plane = planeEstimator.planeAt(currentPoint);
        Container planeSet = plane.intersectionWithSetOneCC(subVolume);
        RealPoint realCenter = ShapeDescriptor<Container>(planeSet).extractCenterOfMass();
        Point centerOfMass = SetProcessor<Container>(planeSet).closestPointAt(realCenter);
        bool stop = false;
        for (const Point &p : computedSkeleton)
            if (l2Metric(centerOfMass, p) <= sqrt(3))
                stop = true;
        if (stop && smoothSkeleton.size() > 0) break;
        if (realCenter != RealPoint::zero && !stop) {
            smoothSkeleton.insert(centerOfMass);
        }
        i++;
    }
    return smoothSkeleton;
}

template<typename Container>
Container
RecenterSkeleton<Container>::
postProcess(const Container &existingSkeleton,
            const Container &processedEdges) {

    Container skeletonPoints(existingSkeleton.domain());
    Container notProcessed(existingSkeleton.domain());
    for (const Point &s : existingSkeleton) {
        if (processedEdges.find(s) == processedEdges.end())
            notProcessed.insert(s);
    }

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    ObjectType objNotProcessed(dt26_6, notProcessed);
    std::vector<ObjectType> objects;
    std::back_insert_iterator<std::vector<ObjectType> > inserter(objects);
    objNotProcessed.writeComponents(inserter);
    for (const ObjectType &o : objects) {
        Container currentSet = o.pointSet();
        Container smoothedSkeleton = recenterSkeletonPoints(*myVolume, currentSet, Container(myVolume->domain()));
        skeletonPoints.insert(smoothedSkeleton.begin(), smoothedSkeleton.end());
    }
    return skeletonPoints;
}


#endif
