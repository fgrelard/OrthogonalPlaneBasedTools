#ifndef LEVEL_CONCATENATION_H
#define LEVEL_CONCATENATION_H

#include <vector>
#include <functional>
#include "Concatenation.h"
#include "DGtal/base/Common.h"

class LevelConcatenation {
public:
    LevelConcatenation(const std::vector<Concatenation> &concats, int level) : myConcatenations(concats),
                                                                               myLevel(level) {}

    double computeAverageLevelFunction(const std::function<double(const DGtal::Z3i::DigitalSet &aSet)> &func,
                                       const std::function<bool(const DGtal::Z3i::DigitalSet &aSet)> &pred = {}) const {
        double sumValue = 0;
        for (const Concatenation &concat : myConcatenations) {
            sumValue += concat.computeSumFunction(func, pred);
        }
        return sumValue / myConcatenations.size();
    }

public:
    std::vector<Concatenation> myConcatenations;
    int myLevel;
};

#endif
