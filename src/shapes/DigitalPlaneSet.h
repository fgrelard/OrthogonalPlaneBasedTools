#ifndef DIGITAL_PLANE_SET_H
#define DIGITAL_PLANE_SET_H

#include "DigitalPlane.h"

template<typename TSpace>
class DigitalPlaneSet {
public:
    typedef DigitalPlane<TSpace> Plane;
    typedef typename Plane::DigitalSet DigitalSet;
    typedef typename Plane::Domain Domain;
    typedef typename Plane::Point Point;
public:
    DigitalPlaneSet() : myDigitalSet(Domain(Point::zero, Point::zero)) {}

    DigitalPlaneSet(const Plane &digitalPlane,
                    const DigitalSet &aDigitalSet) : myDigitalPlane(digitalPlane),
                                                     myDigitalSet(aDigitalSet) {}

    DigitalPlaneSet(const DigitalPlaneSet &other) : myDigitalPlane(other.myDigitalPlane),
                                                    myDigitalSet(other.myDigitalSet) {}

    DigitalSet pointSet() const { return myDigitalSet; }

    Plane digitalPlane() const { return myDigitalPlane; }

    bool isDefined() { return myDigitalSet.size() != 0; }

private:
    Plane myDigitalPlane;
    DigitalSet myDigitalSet;
};

#endif
