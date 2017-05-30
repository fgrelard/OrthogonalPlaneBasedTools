#ifndef CONCATENATION_H
#define CONCATENATION_H

#include <vector>
#include <functional>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

class Concatenation {
public:
    Concatenation(const std::vector<DGtal::Z3i::DigitalSet> &edges, int level) : myEdges(edges), myLevel(level) {}

    double computeSumFunction(const std::function<double(const DGtal::Z3i::DigitalSet &aSet)> &func,
                              const std::function<bool(const DGtal::Z3i::DigitalSet &aSet)> &pred = {}) const {
        double sumValue = 0;
        int cpt = 0;
        for (const DGtal::Z3i::DigitalSet &edge : myEdges) {
            bool checkPred = true;
            if (pred) checkPred = pred(edge);
            if (checkPred) {
                sumValue += func(edge);
                cpt++;
            }
        }
        if (cpt == 0) return 0;
        return sumValue;
    }

    double computeAverageFunction(const std::function<double(const DGtal::Z3i::DigitalSet &aSet)> &func,
                                  const std::function<bool(const DGtal::Z3i::DigitalSet &aSet)> &pred = {}) const {
        double sumValue = 0;
        int cpt = 0;
        for (const DGtal::Z3i::DigitalSet &edge : myEdges) {
            bool checkPred = true;
            if (pred) checkPred = pred(edge);
            if (checkPred) {
                sumValue += func(edge);
                cpt++;
            }
        }
        if (cpt == 0) return 0;
        return sumValue / cpt;
    }


public:
    std::vector<DGtal::Z3i::DigitalSet> myEdges;
    int myLevel;
};

#endif
