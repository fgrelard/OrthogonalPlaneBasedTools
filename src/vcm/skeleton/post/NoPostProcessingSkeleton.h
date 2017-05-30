#ifndef NO_POST_PROCESSING_SKELETON_H
#define NO_POST_PROCESSING_SKELETON_H

#include <vector>
#include "shapes/DigitalPlane.h"

template<typename Container>
class NoPostProcessingSkeleton {
public:
    typedef typename Container::Space Space;
    typedef DigitalPlane<Space> Plane;
public:
    NoPostProcessingSkeleton() = delete;

    NoPostProcessingSkeleton(const Container &skeletonPoints,
                             const Container &a3ShellPoints,
                             const Container &setVolume,
                             const std::vector<Plane> &planesEndPoints) {
        mySkeleton = new Container(skeletonPoints);
    }

    NoPostProcessingSkeleton(const NoPostProcessingSkeleton &other) {
        mySkeleton = new Container(*other.mySkeleton);
    }

    ~NoPostProcessingSkeleton() {
        if (mySkeleton != 0) {
            delete mySkeleton;
            mySkeleton = 0;
        }
    }

public:
    Container postProcess() { return *mySkeleton; }

private:
    Container *mySkeleton;
};

#endif
