#ifndef LINK_POINT_ALGORITHM_H
#define LINK_POINT_ALGORITHM_H

#include "DGtal/base/Common.h"
#include "DGtal/base/CBidirectionalRange.h"

template<typename Point>
class LinkPointAlgorithm {
    BOOST_CONCEPT_ASSERT((DGtal::concepts::CBidirectionalRange<Point>));

public:
    typedef std::vector<Point> Path;
public:
    LinkPointAlgorithm() {}

    LinkPointAlgorithm(const Point &aSource, const Point &aDestination) : mySource(aSource),

                                                                          myDestination(aDestination) {}

    LinkPointAlgorithm(const LinkPointAlgorithm &other) : mySource(other.mySource),
                                                          myDestination(other.myDestination) {}

    virtual ~LinkPointAlgorithm() = default;

    virtual Path linkPoints() = 0;


protected:
    Point mySource;
    Point myDestination;
};

#endif
