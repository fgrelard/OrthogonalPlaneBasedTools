#ifndef MULTI_PATH_THINNER_H
#define MULTI_PATH_THINNER_H

#include "DGtal/kernel/sets/CDigitalSet.h"
#include "geometry/path/LinkPointAlgorithm.h"
#include "geometry/PointUtil.h"
#include "geometry/path/BezierLinkAlgorithm.h"
#include "geometry/CurveProcessor.h"

template<typename Container>
class MultiPathThinner {
    BOOST_CONCEPT_ASSERT((DGtal::concepts::CDigitalSet<Container>));
public:
    typedef typename Container::Space Space;
    typedef typename Space::Point Point;
    typedef typename Space::RealVector RealVector;
    typedef std::pair<Point, RealVector> PointToTangent;
    typedef typename LinkPointAlgorithm<Point>::Path Path;

public:
    MultiPathThinner() = delete;

    MultiPathThinner(const Container &container,
                     const std::vector<PointToTangent> &pointsToLink,
                     const PointToTangent &referencePoint) : myReference(referencePoint),
                                                             myPoints(pointsToLink) {
        myContainer = new Container(container);
    }

    MultiPathThinner(const MultiPathThinner &other) : myReference(other.myReference),
                                                      myPoints(other.myPoints) {
        myContainer = new Container(*other.myContainer);
    }

    ~MultiPathThinner() {
        if (myContainer) {
            delete myContainer;
            myContainer = 0;
        }
    }

public:
    Container linkPoints();

    Container linkPointsThin();

private:
    std::vector<Path> multiPaths();

    Container ensureThinness(const std::vector<Path> &points);

    Point firstIntersectingPoint(const Path &first, const Path &second);

    void orderPaths(std::vector<Path> &otherPaths);

    bool isInVolume(const Path &path);

private:
    Container *myContainer;
    PointToTangent myReference;
    std::vector<PointToTangent> myPoints;
};

template<typename Container>
Container
MultiPathThinner<Container>::linkPoints() {
    Container container(myContainer->domain());
    std::vector<Path> paths = multiPaths();
    for (const Path &p : paths) {
        container.insert(p.begin(), p.end());
    }
    return container;
}


template<typename Container>
Container
MultiPathThinner<Container>::linkPointsThin() {
    std::vector<Path> paths = multiPaths();
    Container container = ensureThinness(paths);
    return container;
}


template<typename Container>
std::vector<typename MultiPathThinner<Container>::Path>
MultiPathThinner<Container>::multiPaths() {
    std::vector<Path> paths;
    for (const PointToTangent &pToT : myPoints) {
        Point p = pToT.first;
        RealVector n = pToT.second;
        std::pair<Point, Point> controlPoints = PointUtil::twoClosestPointsTrackingWithNormal(*myContainer, p, n,
                                                                                              myReference.first,
                                                                                              myReference.second);
        Path path = BezierLinkAlgorithm<Point>(p, myReference.first, controlPoints.first,
                                               controlPoints.second).linkPoints();
        if (isInVolume(path))
            paths.push_back(path);
    }
    return paths;
}

template<typename Container>
typename MultiPathThinner<Container>::Point
MultiPathThinner<Container>::
firstIntersectingPoint(const Path &first, const Path &second) {
    Container setFirst(myContainer->domain());
    Container setSecond(myContainer->domain());
    setFirst.insert(first.begin(), first.end());
    setSecond.insert(second.begin(), second.end());

    CurveProcessor<Container> curveProc(setFirst);
    Container intersection = curveProc.intersectionNeighborhood(setSecond);

    if (intersection.size() > 0) {
        Point candidate = (*intersection.begin());
        for (const Point &p : first) {
            if (intersection.find(p) != intersection.end())
                return p;
        }
    }
    return *first.rbegin();
}

template<typename Container>
void
MultiPathThinner<Container>::
orderPaths(std::vector<Path> &paths) {
    typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;
    L2Metric l2;
    std::vector<int> indices(paths.size(), -1);
    for (int i = 0; i < paths.size(); i++) {
        Path first = paths[i];
        double distance = std::numeric_limits<double>::max();
        for (int j = 0; j < paths.size(); j++) {
            if (i == j) continue;
            Path second = paths[j];
            Point candidate = firstIntersectingPoint(first, second);
            double currentDistance = l2(candidate, *first.begin());
            if (currentDistance < distance) {
                distance = currentDistance;
                indices[i] = j;
            }
        }
    }
    std::vector<Path> tmpPath = paths;
    for (int i = 0; i < paths.size(); i++) {
        if (i == indices[indices[i]]) {
            int index1 = i;
            int index2 = indices[i];
            tmpPath[0] = paths[index1];
            tmpPath[1] = paths[index2];
            int indexJ = 2;
            for (int j = 0; j < paths.size(); j++) {
                if (j == index1 || j == index2) continue;
                tmpPath[indexJ] = paths[j];
                indexJ++;
            }
        }
    }
    paths = tmpPath;

}

template<typename Container>
Container
MultiPathThinner<Container>::
ensureThinness(const std::vector<Path> &paths) {
    typedef DGtal::MetricAdjacency<Space, 3> Adj26;

    std::vector<Path> linkThin = paths;
    Container resultingThin(myContainer->domain());

    while (linkThin.size() >= 2) {
        orderPaths(linkThin);
        Path first = linkThin[0];
        Path second = linkThin[1];
        Container setFirst(myContainer->domain());
        setFirst.insert(first.begin(), first.end());
        Container setSecond(myContainer->domain());
        setSecond.insert(second.begin(), second.end());

        CurveProcessor<Container> curveProc(setFirst);
        Container intersection = curveProc.intersectionNeighborhood(setSecond);

        if (intersection.size() > 0) {
            Point candidate = firstIntersectingPoint(first, second);
            auto iteratorFirst = find(first.begin(), first.end(), candidate);
            resultingThin.insert(first.begin(), iteratorFirst);

            std::vector<Point> neighbors;
            std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
            Adj26::writeNeighbors(inserter, candidate);
            auto iteratorSecond = find_if(second.begin(), second.end(), [&](const Point &p) {
                return (find(neighbors.begin(), neighbors.end(), p) != neighbors.end());
            });

            if (iteratorSecond != second.end()) {
                ++iteratorSecond;
            }
            resultingThin.insert(second.begin(), iteratorSecond);
            RealVector dir = -myReference.second;
            Point newPoint = *iteratorFirst;
            std::pair<Point, Point> controlPoints = PointUtil::twoClosestPointsTrackingWithNormal(*myContainer,
                                                                                                  newPoint, dir,
                                                                                                  myReference.first,
                                                                                                  myReference.second);
            BezierLinkAlgorithm<Point> bezierAlgo(newPoint, myReference.first, controlPoints.first,
                                                  controlPoints.second);
            Path path = bezierAlgo.linkPoints();
            if (isInVolume(path))
                linkThin.push_back(path);
        } else {
            resultingThin.insert(first.begin(), first.end());
            resultingThin.insert(second.begin(), second.end());
        }
        std::vector<Path>(linkThin.begin() + 2, linkThin.end()).swap(linkThin);
    }

    if (linkThin.size() == 1) {
        resultingThin.insert(linkThin[0].begin(), linkThin[0].end());
    }
    return resultingThin;
}

template<typename Container>
bool
MultiPathThinner<Container>::
isInVolume(const Path &path) {
    bool add = true;
    for (const Point &p : path) {
        add &= (myContainer->find(p) != myContainer->end());
    }
    return add;
}

#endif
