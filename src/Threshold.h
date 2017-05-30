#ifndef THRESHOLD_H
#define THRESHOLD_H

#include <vector>
#include <algorithm>

#include "DGtal/base/Common.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include <Eigen/Dense>
#include "DGtal/math/Histogram.h"
#include "DGtal/math/Statistic.h"
#include "geometry/Distance.h"
#include "DGtal/images/CImage.h"

template<typename Container>
class Threshold {
    BOOST_CONCEPT_ASSERT((DGtal::concepts::CCommutativeRing<typename Container::value_type>));
public:
    Threshold() = delete;

    Threshold(const Container &aData) : myData(aData) {}

    double otsuThreshold();

    double unimodalThreshold();

private:
    Container myData;
};


template<typename Container>
double Threshold<Container>::otsuThreshold() {
    using namespace DGtal;
    double proba = 0;                // first order cumulative
    double mu = 0;                // second order cumulative
    double mean = 0;               // total mean level
    double threshold = 0;        // optimal threshold value
    double max = 0.0;

    Statistic<double> stats;
    stats.addValues(myData.begin(), myData.end());
    stats.terminate(); // stats are computed.

    Histogram<double> *hist = new Histogram<double>();
    hist->init(Histogram<double>::SquareRoot, stats);
    hist->addValues(myData.begin(), myData.end());
    hist->terminate();
    double myWidth = (stats.max() - stats.min()) / hist->size();
    double myBin = stats.min();
    for (int i = 0; i < hist->size(); i++) {
        myBin += myWidth;
//		std::cout << myBin << " " << hist->pdf(i) << endl;
        mean += ((double) i / hist->size()) * hist->pdf(i);
    }
    for (int i = 0; i < hist->size(); i++) {
        proba += hist->pdf(i);
        mu += ((double) i / hist->size()) * hist->pdf(i);
        double currentValue = pow((mean * proba - mu), 2) * proba * (1 - proba);
        if (currentValue > max) {
            max = currentValue;
            threshold = ((double) i / hist->size());
        }

    }

    return threshold;
}

template<typename Container>
double Threshold<Container>::unimodalThreshold() {

    typedef DGtal::SpaceND<2, DGtal::int32_t> Space;
    typedef typename Space::RealPoint RealPoint;
    typedef typename Space::RealVector RealVector;

    DGtal::Statistic<double> stats;
    stats.addValues(myData.begin(), myData.end());
    stats.terminate(); // stats are computed.

    DGtal::Histogram<double> *hist = new DGtal::Histogram<double>();
    hist->init(DGtal::Histogram<double>::SquareRoot, stats);
    hist->addValues(myData.begin(), myData.end());
    hist->terminate();
    double myWidth = (stats.max() - stats.min()) / hist->size();
    RealPoint maxPeak(0, 0);
    for (int i = 1; i < hist->size(); i++) {
//		cout << i*myWidth+stats.min() << " " << hist->pdf(i) << endl;
        if (hist->pdf(i) > maxPeak[1])
            maxPeak = RealPoint(i * myWidth, hist->pdf(i));
    }
    RealPoint tail(stats.max(), hist->pdf(hist->size() - 1));
    RealVector directionLine = (tail - maxPeak).getNormalized();
    double maxDistanceOrthogonal = 0.0;
    double threshold = 0.0;

    //Start from maxPeak (origin)
    int begin = maxPeak[0] / myWidth;
    for (int i = begin + 1; i < hist->size(); i++) {
        RealPoint currentPoint(i * myWidth, hist->pdf(i));
        RealVector v = currentPoint - maxPeak;
        RealPoint orthogonalProjection = ((v.dot(directionLine)) / (directionLine.dot(directionLine))) * directionLine;

        //Need to change basis (go back to true origin)
        orthogonalProjection += maxPeak;
        double currentOrthogonalDistance = Distance::euclideanDistance(orthogonalProjection, currentPoint);
        if (currentOrthogonalDistance > maxDistanceOrthogonal) {
            maxDistanceOrthogonal = currentOrthogonalDistance;
            threshold = currentPoint[0];
        }
    }
    threshold = threshold + (threshold - maxPeak[0]);
    return threshold;
}


#endif
