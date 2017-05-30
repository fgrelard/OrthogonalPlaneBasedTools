#ifndef SSI_JUNCTION_DETECTION_H
#define SSI_JUNCTION_DETECTION_H

#include "geometry/SphericalShellIntersection.h"

template<typename Container>
class SSIJunctionDetection {
    BOOST_CONCEPT_ASSERT((DGtal::concepts::CDigitalSet<Container>));
public:
    typedef SphericalShellIntersection<Container> SSI;
    typedef typename Container::Point Point;
public:
    SSIJunctionDetection() = delete;

    SSIJunctionDetection(const Container &container) {
        myContainer = new Container(container);
    }

    SSIJunctionDetection(const SSIJunctionDetection &other) {
        myContainer = new Container(*other.myContainer);
    }

    ~SSIJunctionDetection() {
        if (myContainer) {
            delete myContainer;
            myContainer = 0;
        }

    }

public:
    bool isInJunction(const Point &p, double radius, double noise = 0);

private:
    Container *myContainer;
};

template<typename Container>
bool
SSIJunctionDetection<Container>::
isInJunction(const Point &p, double radius, double noise) {
    SSI ssi(*myContainer, p, radius);
    unsigned int degree = ssi.degree(noise);
    return (degree >= 3);
}


#endif
