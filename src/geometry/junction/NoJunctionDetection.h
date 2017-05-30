#ifndef NO_JUNCTION_DETECTION_H
#define NO_JUNCTION_DETECTION_H

class NoJunctionDetection {

public:
    NoJunctionDetection() {}

public:
    template<typename Point>
    bool isInJunction(const Point &p, double radius, double noise = 0);

};

template<typename Point>
bool
NoJunctionDetection::
isInJunction(const Point &p, double radius, double noise) {
    return false;
}


#endif
