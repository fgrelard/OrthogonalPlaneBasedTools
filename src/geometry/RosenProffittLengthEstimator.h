#ifndef ROSEN_PROFFITT_LE
#define ROSEN_PROFFITT_LE

#include <math.h>
#include <iostream>

template<typename SurfelContainer>
class RosenProffittLengthEstimator {
public:
    typedef typename SurfelContainer::iterator Iterator;

public:
    RosenProffittLengthEstimator() : myFactorDirect{M_PI * ((sqrt(2) + 1) / 8.)},
                                     myFactorDiagonal{M_PI * ((sqrt(2) + 2) / 16.)} {}

    RosenProffittLengthEstimator(double factorDirect, double factorDiagonal) : myFactorDirect{factorDirect},
                                                                               myFactorDiagonal{factorDiagonal} {}

    double eval(Iterator itb, Iterator ite);


private:
    double computeLength(int, int);

private:
    double myFactorDirect;
    double myFactorDiagonal;
};

template<typename SurfelContainer>
double RosenProffittLengthEstimator<SurfelContainer>::
eval(Iterator itb, Iterator ite) {
    int numberDiagonal = 0;
    int numberDirect = 0;

    Iterator nextIt = itb;
    if (ite != itb)
        nextIt++;
    else
        return 0;
    while (nextIt != ite) {
        typename SurfelContainer::value_type::Point difference = nextIt->myCoordinates - itb->myCoordinates;
        bool nextDirect = abs(difference[0]) > 1 || abs(difference[1]) > 1 || abs(difference[2]) > 1;
        if (nextDirect)
            numberDirect++;
        else
            numberDiagonal++;
        ++itb;
        ++nextIt;
    }
    std::cout << numberDirect << " " << numberDiagonal << std::endl;
    return computeLength(numberDiagonal, numberDirect);
}

template<typename SurfelContainer>
double RosenProffittLengthEstimator<SurfelContainer>::
computeLength(int numberDiagonal, int numberDirect) {
    return numberDiagonal * myFactorDiagonal + numberDirect * myFactorDirect;
}

#endif
