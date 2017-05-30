#ifndef CONNECTED_COMPONENT_MERGER_H
#define CONNECTED_COMPONENT_MERGER_H

template<typename Container>
class CurveProcessor;

template<typename Container>
class SetProcessor;


#include <vector>
#include "DGtal/kernel/sets/CDigitalSet.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/topology/DigitalTopology.h"
#include "DGtal/topology/Object.h"
#include "geometry/CurveProcessor.h"
#include "geometry/SetProcessor.h"

template<typename Space>
class ConnectedComponentMerger {
    BOOST_CONCEPT_ASSERT((DGtal::concepts::CSpace<Space>));

public:
    typedef typename Space::Point Scalar;
    typedef DGtal::MetricAdjacency<Space, 1> Adj6;
    typedef DGtal::MetricAdjacency<Space, 3> Adj26;
    typedef DGtal::DigitalTopology<Adj26, Adj6> DT26_6;
    typedef typename DGtal::HyperRectDomain<Space> Domain;
    typedef typename DGtal::DigitalSetSelector<Domain, DGtal::BIG_DS + DGtal::HIGH_BEL_DS>::Type Container;
    typedef DGtal::Object<DT26_6, Container> ObjectType;

public:
    ConnectedComponentMerger(const std::vector<ObjectType> &cc,
                             const ObjectType &surroundingSetGraph,
                             double lowerBound = std::numeric_limits<double>::min(),
                             double upperBound = std::numeric_limits<double>::max(),
                             int indexRef = -1) {
        if (indexRef == -1) {
            if (surroundingSetGraph.size() > 0)
                findObjectToMerge(cc, surroundingSetGraph, lowerBound, upperBound);
            else
                findObjectToMerge(cc, lowerBound, upperBound);
        } else
            findObjectToMerge(cc, surroundingSetGraph, lowerBound, upperBound, indexRef);

    }

    ConnectedComponentMerger(const std::vector<ObjectType> &cc,
                             double lowerBound = std::numeric_limits<double>::min(),
                             double upperBound = std::numeric_limits<double>::max()) {
        findObjectToMerge(cc, lowerBound, upperBound);
    }

public:
    Scalar first() const { return myFirstPoint; }

    Scalar second() const { return mySecondPoint; }

    int firstPosition() const { return myFirstIndex; }

    int secondPosition() const { return mySecondIndex; }

    bool isUndefined() const { return (myFirstIndex == -1 || mySecondIndex == -1); }

public:

    void findObjectToMerge(const std::vector<ObjectType> &cc,
                           const ObjectType &surroundingSetGraph,
                           double lowerBound,
                           double upperBound,
                           int index);

    void findObjectToMerge(const std::vector<ObjectType> &cc,
                           double lowerBound,
                           double upperBound);

    void findObjectToMerge(const std::vector<ObjectType> &cc,
                           const ObjectType &surroundingSetGraph,
                           double lowerBound,
                           double upperBound);

    void mergeObjects(std::vector<ObjectType> &cc,
                      const std::vector<Scalar> &linkBetweenObjects);

private:
    Scalar myFirstPoint;
    Scalar mySecondPoint;
    int myFirstIndex = -1;
    int mySecondIndex = -1;
};

template<typename Space>
void ConnectedComponentMerger<Space>::findObjectToMerge(const std::vector<ObjectType> &objects,
                                                        const ObjectType &graph,
                                                        double lowerBound,
                                                        double upperBound,
                                                        int index) {
    typedef typename DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;

    typedef DistanceToPointFunctor<L2Metric> DistanceFunctor;

    typedef DGtal::DistanceBreadthFirstVisitor<ObjectType, DistanceFunctor, std::set<Scalar> > Visitor;
    typedef typename Visitor::Node MyNode;

    myFirstIndex = index;

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    double distance = std::numeric_limits<double>::max();

    Container currentCC = objects[myFirstIndex].pointSet();
    CurveProcessor<Container> curveProc(currentCC);
    Container endPointCurrent = curveProc.endPoints();
    Container merged(graph.domain());
    for (int i = 0; i < objects.size(); i++) {
        if (i != index) {
            Container set = objects[i].pointSet();
            merged.insert(set.begin(), set.end());
        }
    }

    for (const Scalar &e : endPointCurrent) {
        L2Metric l2Metric;
        DistanceFunctor functor(l2Metric, e);
        Visitor visitor(graph, functor, e);
        MyNode node;
        double currentDistance = std::numeric_limits<double>::max();
        Scalar candidate;
        while (!visitor.finished()) {
            node = visitor.current();
            if (node.second > upperBound || node.second > distance) break;
            if (node.second > lowerBound) {
                auto iteratorObj = merged.find(node.first);
                if (iteratorObj != merged.end()) {
                    currentDistance = node.second;
                    candidate = node.first;
                    break;
                }
            }
            visitor.expand();
        }
        if (currentDistance < distance) {
            myFirstPoint = e;
            mySecondPoint = candidate;
            distance = currentDistance;
        }
    }

    for (auto itCC = objects.begin(), itCCe = objects.end(); itCC != itCCe; ++itCC) {
        Container oSet = itCC->pointSet();
        auto iterator = oSet.find(mySecondPoint);
        if (iterator != oSet.end()) {
            mySecondIndex = itCC - objects.begin();
            break;
        }
    }
}


template<typename Space>
void
ConnectedComponentMerger<Space>::findObjectToMerge(const std::vector<ObjectType> &objects,
                                                   double lowerBound,
                                                   double upperBound) {
    typedef typename DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;
    L2Metric l2Metric;
    double distance = std::numeric_limits<double>::max();
    for (size_t i = 0; i < objects.size(); i++) {
        Container currentCC = objects[i].pointSet();
        CurveProcessor<Container> curveProc(currentCC);
        Container endPointCurrent = curveProc.endPoints();
        for (size_t j = i + 1; j < objects.size(); j++) {
            Container otherCC = objects[j].pointSet();
            for (const Scalar &e : endPointCurrent) {
                Scalar closest = SetProcessor<Container>(otherCC).closestPointAt(e);
                double currentDistance = l2Metric(closest, e);
                if (currentDistance > lowerBound &&
                    currentDistance <= upperBound &&
                    currentDistance < distance) {
                    myFirstPoint = e;
                    mySecondPoint = closest;
                    distance = currentDistance;
                    myFirstIndex = i;
                    mySecondIndex = j;
                }
            }

        }
    }
}

template<typename Space>
void
ConnectedComponentMerger<Space>::findObjectToMerge(const std::vector<ObjectType> &objects,
                                                   const ObjectType &graph, double lowerBound, double upperBound) {
    typedef typename DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;

    typedef DistanceToPointFunctor<L2Metric> DistanceFunctor;

    typedef DGtal::DistanceBreadthFirstVisitor<ObjectType, DistanceFunctor, std::set<Scalar> > Visitor;
    typedef typename Visitor::Node MyNode;

    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);
    double distance = std::numeric_limits<double>::max();


    for (size_t i = 0; i < objects.size(); i++) {
        Container currentCC = objects[i].pointSet();
        CurveProcessor<Container> curveProc(currentCC);
        Container endPointCurrent = curveProc.endPoints();
        Container merged(graph.domain());
        for (int index = 0; index < objects.size(); index++) {
            if (index != i) {
                Container set = objects[index].pointSet();
                merged.insert(set.begin(), set.end());
            }
        }
        for (const Scalar &e : endPointCurrent) {
            L2Metric l2Metric;
            DistanceFunctor functor(l2Metric, e);
            Visitor visitor(graph, functor, e);
            MyNode node;
            double currentDistance = std::numeric_limits<double>::max();
            Scalar candidate;
            while (!visitor.finished()) {
                node = visitor.current();
                if (node.second > upperBound || node.second > distance) break;
                if (node.second > lowerBound) {
                    auto iteratorObj = merged.find(node.first);
                    if (iteratorObj != merged.end()) {
                        currentDistance = node.second;
                        candidate = node.first;
                        break;
                    }
                }
                visitor.expand();
            }
            if (currentDistance < distance) {
                myFirstPoint = e;
                mySecondPoint = candidate;
                distance = currentDistance;
                myFirstIndex = i;
            }
        }
    }
    for (auto itCC = objects.begin(), itCCe = objects.end(); itCC != itCCe; ++itCC) {
        Container oSet = itCC->pointSet();
        auto iterator = oSet.find(mySecondPoint);
        if (iterator != oSet.end()) {
            mySecondIndex = itCC - objects.begin();
            break;
        }
    }
}

template<typename Space>
void
ConnectedComponentMerger<Space>::
mergeObjects(std::vector<ObjectType> &objects, const std::vector<Scalar> &link) {


    Adj26 adj26;
    Adj6 adj6;
    DT26_6 dt26_6(adj26, adj6, DGtal::JORDAN_DT);

    Container toKeep = objects[myFirstIndex].pointSet();
    Container toDelete = objects[mySecondIndex].pointSet();
    toKeep.insert(link.begin(), link.end());
    toKeep.insert(toDelete.begin(),
                  toDelete.end());
    ObjectType toAdd(dt26_6, toKeep);
    if (myFirstIndex == objects.size() - 1) {
        objects.pop_back();
        objects[mySecondIndex] = objects[myFirstIndex - 1];
        objects.pop_back();
    } else if (myFirstIndex == objects.size() - 2) {
        objects[mySecondIndex] = objects[myFirstIndex + 1];
        objects.pop_back();
        objects.pop_back();
    } else if (mySecondIndex == objects.size() - 1) {
        objects.pop_back();
        objects[myFirstIndex] = objects[mySecondIndex - 1];
        objects.pop_back();
    } else if (mySecondIndex == objects.size() - 2) {
        objects[myFirstIndex] = objects[mySecondIndex + 1];
        objects.pop_back();
        objects.pop_back();
    } else {
        objects[myFirstIndex] = objects[objects.size() - 1];
        objects[mySecondIndex] = objects[objects.size() - 2];
        objects.pop_back();
        objects.pop_back();
    }
    objects.insert(objects.begin(), toAdd);
}


#endif
