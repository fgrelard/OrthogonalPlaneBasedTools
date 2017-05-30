#ifndef WATERSHED_H
#define WATERSHED_H

#include "geometry/WeightedPointCount.h"
#include <vector>
#include <queue>
#include <map>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/MetricAdjacency.h"
#include <boost/bimap.hpp>

#define MASK -3
#define WATERSHED -2
#define UNLABELLED -1
#define FIRST_LABEL 0


template<typename Point>
class Watershed {


    class WatershedPoint : public WeightedPointCount<Point> {
        typedef WeightedPointCount<Point> Base;
    public:
        using Base::Base;
    public:
        WatershedPoint() {}

        WatershedPoint(const Point &aPoint, double aWeight, int aCount, double aDistance) : Base(aPoint, aWeight,
                                                                                                 aCount),
                                                                                            myDistance(aDistance) {}

        WatershedPoint(const WatershedPoint &other) : Base(other), myDistance(other.myDistance) {}

        friend bool operator<(const WatershedPoint &it, const WatershedPoint &other) {
            return (it.myWeight < other.myWeight);
        }

    public:
        double myDistance;
    };


    class WatershedInformation {
    public:
        WatershedInformation() : myLabel(UNLABELLED), myDistance(0.0), myValue(0.0) {}

        WatershedInformation(double aValue, int aIdentifier) : myLabel(UNLABELLED), myDistance(0.0), myValue(aValue),
                                                               myIdentifier(aIdentifier) {}

        void setLabel(int aLabel) { myLabel = aLabel; }

        int getLabel() const { return myLabel; }

        void setValue(double aValue) { myValue = aValue; }

        double getValue() const { return myValue; }

        void setDistance(double aDistance) { myDistance = aDistance; }

        double getDistance() const { return myDistance; }

        int getIdentifier() { return myIdentifier; }

        friend bool operator<(WatershedInformation it, WatershedInformation other) {
            return (it.getValue() < other.getValue() ||
                    (it.getValue() == other.getValue() && it.getIdentifier() < other.getIdentifier()));
        }

    private:
        int myLabel;
        double myDistance;
        double myValue;
        int myIdentifier;
    };


public:
    typedef WatershedInformation WatershedInformation;

public:
    template<typename Container>
    Watershed(const Container &container, double aEpsilon);

    void compute();

    void pixelsAtSameAltitude(std::queue<Point> &fifo, Point &current, int &currentI, double currentAltitude);

    void extendBasins(std::queue<Point> &fifo);

    void detectMinima(std::queue<Point> &fifo, int &label, double altitude);

    std::map<Point, WatershedInformation> getWatershed();

    int getBins();

    void closeWatershed();

private:
    std::map<Point, WatershedInformation> myImageWatershed;
    std::map<WatershedInformation, Point> sortedMap;
    double myEpsilon;
    Point myFictitious;
};


template<typename Point>
template<typename Container>
Watershed<Point>::Watershed(const Container &container, double aEpsilon) {
    myEpsilon = aEpsilon;
    double uniqueLabel = 0;
    for (const auto &p : container) {
        WatershedInformation wi(p.second, uniqueLabel++);
        myImageWatershed[p.first] = wi;
        sortedMap[wi] = p.first;
        if (p.first < myFictitious)
            myFictitious = p.first - Point::diagonal(100);
    }
}


template<typename Point>
void Watershed<Point>::compute() {
    int label = UNLABELLED;
    std::queue<Point> fifo;
    Point current = sortedMap.begin()->second;
    int currentI = 0;
    DGtal::trace.beginBlock("Watershed");
    double altitude = 0;
    while (currentI < myImageWatershed.size()) {
        DGtal::trace.progressBar(currentI, myImageWatershed.size());
        altitude = myImageWatershed[current].getValue();
        pixelsAtSameAltitude(fifo, current, currentI, altitude);
        extendBasins(fifo);
        detectMinima(fifo, label, altitude);
    }
    closeWatershed();
    DGtal::trace.endBlock();
}


template<typename Point>
void
Watershed<Point>::pixelsAtSameAltitude(std::queue<Point> &fifo, Point &current, int &currentI, double currentAltitude) {
    DGtal::MetricAdjacency<DGtal::Z3i::Space, 3> adj;
    for (auto it = std::next(sortedMap.begin(), currentI), ite = sortedMap.end(); it != ite; ++it) {
        //Once we found the right position in vector, i is not incremented anymore
        std::pair<WatershedInformation, Point> pairWatershed = *it;
        Point p = pairWatershed.second;
        WatershedInformation wi = pairWatershed.first;
        current = p;
        if (wi.getValue() > currentAltitude + myEpsilon) break;

        wi.setLabel(MASK);
        myImageWatershed[p] = wi;
        std::vector<Point> neighbors;
        std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
        adj.writeNeighbors(inserter, p);
        for (const Point &n : neighbors) {
            if (myImageWatershed.find(n) != myImageWatershed.end()) {
                WatershedInformation wn = myImageWatershed.at(n);
                if (wn.getLabel() > 0 || wn.getLabel() == WATERSHED) {
                    wi.setDistance(1);
                    myImageWatershed[p] = wi;
                    fifo.push(p);
                    break;
                }
            }
        }
        currentI++;

    }

    //To break from main loop: all the points have been processed
    // if (i == sortedMap.size())
    //         currentI = sortedMap.size();
}

template<typename Point>
void Watershed<Point>::extendBasins(std::queue<Point> &fifo) {
    // The first pixels treated are the closest (d=1), then d=2...
    int d_cur = 1;
    fifo.push(myFictitious);
    DGtal::MetricAdjacency<DGtal::Z3i::Space, 3> adj;

    while (true) {
        Point p = fifo.front();
        fifo.pop();
        if (p == myFictitious) {
            if (fifo.empty()) // this altitude is processed
                break;
            else {
                fifo.push(myFictitious);
                d_cur++;
                p = fifo.front();
                fifo.pop();
            }
        }
        // Labelling p by inspecting neighbours
        std::vector<Point> neighbors;
        std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
        adj.writeNeighbors(inserter, p);

        WatershedInformation wp = myImageWatershed.at(p);
        for (const Point &n : neighbors) {
            if (myImageWatershed.find(n) != myImageWatershed.end()) {
                WatershedInformation wn = myImageWatershed.at(n);
                if (wn.getDistance() <= d_cur &&
                    (wn.getLabel() == WATERSHED || wn.getLabel() > 0)) {
                    if (wn.getLabel() > 0) {
                        if (wp.getLabel() == MASK) {
                            wp.setLabel(wn.getLabel());
                            myImageWatershed[p] = wp;
                        } else if (wp.getLabel() != wn.getLabel()) {
                            wp.setLabel(WATERSHED);
                            myImageWatershed[p] = wp;
                        }
                    } else if (wp.getLabel() == MASK) {
                        wp.setLabel(WATERSHED);
                        myImageWatershed[p] = wp;
                    }
                } else if (wn.getLabel() == MASK && wn.getDistance() == 0) {
                    wn.setDistance(d_cur + 1);
                    myImageWatershed[n] = wn;
                    fifo.push(n);
                }
            }
        }

    } // End : while ( true )
}


template<typename Point>
void Watershed<Point>::detectMinima(std::queue<Point> &fifo, int &label, double altitude) {
    DGtal::MetricAdjacency<DGtal::Z3i::Space, 3> adj;
    for (const std::pair<WatershedInformation, Point> &pairWatershed : sortedMap) {
        Point p = pairWatershed.second;
        WatershedInformation wp = myImageWatershed.at(p);
        if (wp.getValue() >= altitude && wp.getValue() <= altitude + myEpsilon) {
            wp.setDistance(0);
            myImageWatershed[p] = wp;
            if (wp.getLabel() == MASK) {
                label++;
                fifo.push(p);
                wp.setLabel(label);
                myImageWatershed[p] = wp;
                while (!fifo.empty()) {
                    Point q = fifo.front();
                    fifo.pop();
                    // Labelling p by inspecting neighbours
                    std::vector<Point> neighbors;
                    std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
                    adj.writeNeighbors(inserter, q);
                    for (const Point &n : neighbors) {
                        if (myImageWatershed.find(n) != myImageWatershed.end()) {
                            WatershedInformation wn = myImageWatershed.at(n);
                            if (wn.getLabel() == MASK) {
                                fifo.push(n);
                                wn.setLabel(label);
                                myImageWatershed[n] = wn;
                            }
                        }
                    }

                }
            }
        }
    }
}

template<typename Point>
void Watershed<Point>::closeWatershed() {
    DGtal::MetricAdjacency<DGtal::Z3i::Space, 1> adj;
    for (const std::pair<Point, WatershedInformation> &pairWatershed : myImageWatershed) {
        Point p = pairWatershed.first;
        WatershedInformation wp = pairWatershed.second;
        std::vector<Point> neighbors;
        std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
        adj.writeNeighbors(inserter, p);
        for (const Point &n : neighbors) {
            if (myImageWatershed.find(n) != myImageWatershed.end()) {
                WatershedInformation wn = myImageWatershed.at(n);
                if (wp.getLabel() != wn.getLabel() &&
                    wn.getLabel() != WATERSHED &&
                    wp.getLabel() != WATERSHED) {
                    wp.setLabel(WATERSHED);
                    myImageWatershed[p] = wp;
                }
            }
        }
    }
}

template<typename Point>
std::map<Point, typename Watershed<Point>::WatershedInformation> Watershed<Point>::getWatershed() {
    return myImageWatershed;
}

template<typename Point>
int Watershed<Point>::getBins() {
    return (*max_element(myImageWatershed.begin(), myImageWatershed.end(),
                         [](const std::pair<Point, WatershedInformation> &one,
                            const std::pair<Point, WatershedInformation> &two) {
                             return one.second.getLabel() < two.second.getLabel();
                         })).second.getLabel();
}


#endif
