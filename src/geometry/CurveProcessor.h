#ifndef CURVE_PROCESSOR_H
#define CURVE_PROCESSOR_H

//Forward declaration due to dependent inclusion
template<typename Space>
class ConnectedComponentMerger;

#include <vector>
#include <queue>

#include "DGtal/base/Common.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/topology/DigitalTopology.h"
#include "DGtal/topology/Object.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "shapes/Ball.h"
#include "geometry/path/AStarAlgorithm.h"
#include "ConnectedComponentMerger.h"


template<typename Container>
class CurveProcessor {

    BOOST_CONCEPT_ASSERT((DGtal::concepts::CDigitalSet<Container>));


public:
    typedef typename Container::value_type Point;
    typedef typename Container::Space Space;
    typedef DGtal::MetricAdjacency<Space, 1> Adj6;
    typedef DGtal::MetricAdjacency<Space, 3> Adj26;
    typedef DGtal::DigitalTopology<Adj26, Adj6> DT26_6;
    typedef DGtal::Object<DT26_6, Container> ObjectType;
    typedef typename Container::Space::RealVector RealVector;
    typedef typename ConnectedComponentMerger<Space>::Domain Domain;

public:
    CurveProcessor(const Container &container) : myCurve(container) {}

public:

    Container ensureConnectivity();

    Container endPoints();

    Container branchingPoints();

    template<typename DTL2>
    Container subCurve(const DTL2 &dt,
                       const Container &constraintInSet);


    Container fillHoles(double lowerBound = std::numeric_limits<double>::min(),
                        double upperBound = std::numeric_limits<double>::max());

    Container fillHoles(const Container &container);

    Container fillHolesNotInSet(const Container &set, const Container &setVolume);

    bool isPointThin(const Point &point);

    std::vector<Point> convertToOrderedCurve();

    std::vector<Point> convertToOrderedCurve(const Point &startingPoint);

    Container intersectionNeighborhood(const Container &otherCurve);

private:
    Container myCurve;
};


template<typename Container>
Container CurveProcessor<Container>::ensureConnectivity() {

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    Container cleanSet(myCurve.domain());
    ObjectType obj(dt26_6, myCurve);
    Container &S = obj.pointSet();

    Container endSet = endPoints();
    cleanSet = S;
    for (auto it = S.begin(), ite = S.end(); it != ite; ++it) {
        ObjectType obj(dt26_6, cleanSet);
        if (obj.isSimple(*it) &&
            endSet.find(*it) == endSet.end()) {
            cleanSet.erase(*it);
        }
    }

    return cleanSet;
}

template<typename Container>
Container CurveProcessor<Container>::endPoints() {

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);

    Container set = myCurve;
    ObjectType objectSet(dt26_6, set);
    Container endPoints(set.domain());
    for (auto it = set.begin(), ite = set.end(); it != ite; ++it) {
        Point p = *it;
        std::vector<Point> neighbors;
        std::back_insert_iterator<std::vector<Point>> inserter(neighbors);
        objectSet.writeNeighbors(inserter, p);
        if (neighbors.size() <= 1)
            endPoints.insert(p);

        //Is it in same quadrant: case connectivity != 26
        else {
            RealVector previous;
            bool isEndPoint = true;
            std::vector<RealVector> vectors;
            for (const Point &n : neighbors) {
                RealVector dir = (n - p).getNormalized();
                vectors.push_back(dir);
            }
            //Min angle (resp max dot product) determined by two points with one varying coordinate
            for (int i = 0; i < vectors.size(); i++) {
                for (int j = i + 1; j < vectors.size(); j++) {
                    if (vectors[i].dot(vectors[j]) <= (1 / (1 + sqrt(2))))
                        isEndPoint = false;
                }
            }
            if (isEndPoint)
                endPoints.insert(p);
        }
    }
    return endPoints;
}


template<typename Container>
Container CurveProcessor<Container>::branchingPoints() {

    typedef DGtal::DepthFirstVisitor<ObjectType, std::set<Point> > Visitor;
    typedef typename Visitor::Node MyNode;
    typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    L2Metric l2Metric;

    ObjectType graph(dt26_6, myCurve);
    Container branchingPoints(myCurve.domain());
    Container end = endPoints();
    if (end.size() == 0) return branchingPoints;
    Visitor visitor(graph, *end.begin());
    MyNode node;

    std::vector<Point> existingSkeletonOrdered;
    std::pair<Point, double> previous;
    while (!visitor.finished()) {
        node = visitor.current();
        if (node.second != 0 && ((int) node.second - previous.second) <= 0) {
            std::vector<Point> neighbors;
            std::back_insert_iterator<std::vector<Point>> inserter(neighbors);
            graph.writeNeighbors(inserter, node.first);
            double minDistance = std::numeric_limits<double>::max();
            Point cand;
            for (const Point &n : neighbors) {
                if (find(existingSkeletonOrdered.begin(), existingSkeletonOrdered.end(), n) !=
                    existingSkeletonOrdered.end()) {
                    double currentDistance = l2Metric(n, node.first);
                    if (currentDistance < minDistance) {
                        minDistance = currentDistance;
                        cand = n;
                    }
                }
            }
            branchingPoints.insert(cand);
            existingSkeletonOrdered.push_back(cand);

        }
        previous = node;
        existingSkeletonOrdered.push_back(node.first);
        visitor.expand();
    }
    return branchingPoints;
}


template<typename Container>
template<typename DTL2>
Container
CurveProcessor<Container>::
subCurve(const DTL2 &dt, const Container &constraintInSet) {
    Container ep = endPoints();
    Point e;
    for (const Point &p : ep) {
        if (constraintInSet.find(p) != constraintInSet.end())
            e = p;
    }
    double ballRadius = (dt.domain().isInside(e)) ? dt(e) : 2;
    ballRadius = (ballRadius < 2) ? 2 : ballRadius;
    Ball<Point> ball(e, ballRadius);
    Container restrictedEdge(myCurve.domain());
    for (const Point &e : myCurve) {
        if (!ball.contains(e))
            restrictedEdge.insert(e);
    }
    if (restrictedEdge.size() < 2) restrictedEdge = myCurve;
    return restrictedEdge;
}


template<typename Container>
Container
CurveProcessor<Container>::
fillHoles(double lowerBound, double upperBound) {
    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);

    ObjectType graph(dt26_6, Container(Domain(Point::zero, Point::zero)));

    ObjectType objectImage(dt26_6, myCurve);

    std::vector<ObjectType> skeletonCC;
    std::back_insert_iterator<std::vector<ObjectType> > inserterCC(skeletonCC);
    objectImage.writeComponents(inserterCC);

    Container myOneCCCurve = myCurve;

    bool shouldStop = false;
    while (!shouldStop) {
        int indexRef = 0;
        ConnectedComponentMerger<Space> vPair(skeletonCC, graph, lowerBound, upperBound);
        if (!vPair.isUndefined()) {
            Point first = vPair.first();
            Point second = vPair.second();
            LinkPointAlgorithm<Point> *linkAlgo;
            linkAlgo = new BresenhamAlgorithm<Point>(first, second);
            std::vector<Point> link = linkAlgo->linkPoints();
            delete linkAlgo;
            vPair.mergeObjects(skeletonCC, link);
            myOneCCCurve.insert(link.begin(), link.end());
        } else
            shouldStop = true;
    }
    return myOneCCCurve;

}


template<typename Container>
Container
CurveProcessor<Container>::
fillHoles(const Container &setVolume) {
    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);

    ObjectType graph(dt26_6, setVolume);

    ObjectType objectImage(dt26_6, myCurve);

    std::vector<ObjectType> skeletonCC;
    std::back_insert_iterator<std::vector<ObjectType> > inserterCC(skeletonCC);
    objectImage.writeComponents(inserterCC);
    sort(skeletonCC.begin(), skeletonCC.end(), [&](const ObjectType &one,
                                                   const ObjectType &two) {
             return one.size() > two.size();
         });
    Container myOneCCCurve = myCurve;

    bool shouldStop = false;
    while (!shouldStop) {
        ConnectedComponentMerger<Space> vPair(skeletonCC, graph, std::numeric_limits<double>::min(),
                                              std::numeric_limits<double>::max());
        if (!vPair.isUndefined()) {
            Point first = vPair.first();
            Point second = vPair.second();
            LinkPointAlgorithm<Point> *linkAlgo;
            linkAlgo = new AStarAlgorithm<Point, Container>(first, second, setVolume);
            std::vector<Point> link = linkAlgo->linkPoints();
            delete linkAlgo;
            vPair.mergeObjects(skeletonCC, link);
            myOneCCCurve.insert(link.begin(), link.end());
        } else
            shouldStop = true;
    }
    return myOneCCCurve;

}

template<typename Container>
Container
CurveProcessor<Container>::
fillHolesNotInSet(const Container &set, const Container &setVolume) {
    typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    L2Metric l2Metric;
    ObjectType graph(dt26_6, setVolume);

    ObjectType objectImage(dt26_6, myCurve);

    std::vector<ObjectType> skeletonCC;
    std::back_insert_iterator<std::vector<ObjectType> > inserterCC(skeletonCC);
    objectImage.writeComponents(inserterCC);

    Container myOneCCCurve = myCurve;

    bool shouldStop = false;
    while (!shouldStop) {
        ConnectedComponentMerger<Space> vPair(skeletonCC);
        if (!vPair.isUndefined()) {
            Point first = vPair.first();
            Point second = vPair.second();

            AStarAlgorithm<Point, Container> linkAstar(first, second, setVolume);
            std::vector<Point> linkA = linkAstar.linkPoints();

            BresenhamAlgorithm<Point> linkBresenham(first, second);
            std::vector<Point> linkB = linkBresenham.linkPoints();
            if (linkB.size() * 5 < linkA.size())
                vPair = ConnectedComponentMerger<Space>(skeletonCC, graph);

            first = vPair.first();
            second = vPair.second();

            RealVector dirFirst = (second - first).getNormalized();
            RealVector dirSecond = (first - second).getNormalized();
            double radius = l2Metric(first, second) + 1;
            Ball<Point> ballFirst(first, radius);
            Ball<Point> ballSecond(second, radius);

            Container traversedFirst = ballFirst.pointsInHalfBall(dirFirst);
            Container traversedSecond = ballSecond.pointsInHalfBall(dirSecond);
            bool add = true;
            for (const Point &p : set) {
                if (traversedFirst.find(p) != traversedFirst.end() ||
                    traversedSecond.find(p) != traversedSecond.end())
                    add = false;
            }
            std::vector<Point> link;
            if (add) {
                AStarAlgorithm<Point, Container> linkAlgo(first, second, setVolume);
                link = linkAlgo.linkPoints();
            }
            vPair.mergeObjects(skeletonCC, link);
            myOneCCCurve.insert(link.begin(), link.end());
        } else
            shouldStop = true;
    }
    return myOneCCCurve;

}


template<typename Container>
bool
CurveProcessor<Container>::
isPointThin(const Point &point) {
    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    ObjectType obj(dt26_6, myCurve);
    std::vector<Point> neighbors;
    std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
    obj.writeNeighbors(inserter, point);
    if (neighbors.size() < 2) {
        return true;
    }
    return false;
}


template<typename Container>
std::vector<typename CurveProcessor<Container>::Point> CurveProcessor<Container>::convertToOrderedCurve() {
    std::vector<Point> orientedEdge;

    Container set = myCurve;
    if (set.size() == 0) return orientedEdge;

    Container e = endPoints();
    Point start = *(e.begin());
    return convertToOrderedCurve(start);

}


template<typename Container>
std::vector<typename CurveProcessor<Container>::Point>
CurveProcessor<Container>::convertToOrderedCurve(const typename CurveProcessor<Container>::Point &startingPoint) {
    std::vector<Point> orientedEdge;
    Container edge = myCurve;
    if (edge.size() == 0) return orientedEdge;

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);

    ObjectType objEdge(dt26_6, edge);
    Point start = startingPoint;
    orientedEdge.push_back(start);
    bool toAdd = true;
    while (toAdd) {
        std::vector<Point> neighbors;
        std::back_insert_iterator<std::vector<Point>> inserter(neighbors);
        objEdge.writeNeighbors(inserter, start);
        unsigned int cpt = 0;
        for (const Point &n : neighbors) {
            if (std::find(orientedEdge.begin(), orientedEdge.end(), n) == orientedEdge.end()) {
                orientedEdge.push_back(n);
                start = n;
            } else
                cpt++;
        }
        if (cpt == neighbors.size())
            toAdd = false;
    }
    return orientedEdge;
}

template<typename Container>
Container
CurveProcessor<Container>::
intersectionNeighborhood(const Container &otherCurve) {
    Container container(myCurve.domain());
    for (const Point &p : myCurve) {
        std::vector<Point> neighbors;
        std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
        Adj26::writeNeighbors(inserter, p);
        for (const Point &o : otherCurve) {
            auto iterator = find(neighbors.begin(), neighbors.end(), o);
            if (iterator != neighbors.end())
                container.insert(p);
        }
    }
    return container;
}

#endif
