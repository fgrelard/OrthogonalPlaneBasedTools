#ifndef __MORPHOMATHS__
#define __MORPHOMATHS__

#include <limits>
#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/kernel/domains/DomainPredicate.h"
#include "DGtal/images/CImage.h"

template<typename Image>
class Morphomaths {

    BOOST_CONCEPT_ASSERT((DGtal::concepts::CImage<Image>));

public:
    typedef typename Image::Point Point;
    typedef typename Image::Domain Domain;
    typedef typename Image::Value Value;

public:
    Morphomaths() = delete;

    Morphomaths(const Image &image, int aSize = 1) : myImage(image), mySize(aSize) {}

    Morphomaths(const Morphomaths &other) : myImage(other.myImage), mySize(other.mySize) {}

public:
    bool process(Point p, Value& valueMin, Value& valueMax);

    Image constructOnePxBorderImage();

    Image erosion();

    Image dilation();

    Image open();

    Image close();

protected:
    Image myImage;
    int mySize;

};


template<typename Image>
bool Morphomaths<Image>::process(Point p, Value& valueMin, Value& valueMax) {
    Domain rect(p - Point::diagonal(mySize), p + Point::diagonal(mySize));
    for (const Point& r : rect) {
        if (myImage.domain().isInside(r)) {
            Value value = myImage(r);
            if (r != p) {
                if (value < valueMin) valueMin = value;
                if (value > valueMax) valueMax = value;
            }
        }
    }
    return (valueMin != valueMax);
}

template<typename Image>
Image Morphomaths<Image>::constructOnePxBorderImage() {
    typedef typename Image::Domain Domain;
    typedef typename Image::Point Point;

    Domain domain = myImage.domain();

    Image toReturn(Domain(domain.lowerBound() - Point::diagonal(), domain.upperBound() + Point::diagonal()));
    for (auto it = toReturn.domain().begin(), ite = toReturn.domain().end(); it != ite; ++it) {
        Point p = *it;
        if (domain.isInside(p))
            toReturn.setValue(p, myImage(p));
        else
            toReturn.setValue(p, 0);
    }

    return toReturn;
}

template<typename Image>
Image Morphomaths<Image>::erosion() {
    typedef typename Image::Domain Domain;
    typedef typename Image::Point Point;

    Domain domain = myImage.domain();
    Image toReturn = constructOnePxBorderImage();
    myImage = toReturn;

    for (const Point& p : domain) {
        Value valueMin = std::numeric_limits<int>::max();
        Value valueMax = 0;
        bool shouldBeEroded = process(p, valueMin, valueMax);
        if (shouldBeEroded) toReturn.setValue(p, valueMin);
    }
    Image out(domain);
    for (auto it = domain.begin(), ite = domain.end(); it != ite; ++it) {
        out.setValue(*it, toReturn(*it));
    }
    return out;
}

template<typename Image>
Image Morphomaths<Image>::dilation() {
    typedef typename Image::Domain Domain;
    typedef typename Image::Point Point;
    typedef DGtal::functors::NotPointPredicate<Image> BackgroundPredicate;

    Domain domain = myImage.domain();
    Image toReturn = constructOnePxBorderImage();
    myImage = toReturn;
    BackgroundPredicate backgroundPredicate(myImage);


   for (const Point& p : domain) {
        Value valueMin = std::numeric_limits<int>::max();
        Value valueMax = 0;
        bool shouldBeEroded = process(p, valueMin, valueMax);
        if (shouldBeEroded) toReturn.setValue(p, valueMax);
    }
    Image out(domain);
    for (auto it = domain.begin(), ite = domain.end(); it != ite; ++it) {
        out.setValue(*it, toReturn(*it));
    }
    return out;
}

template<typename Image>
Image Morphomaths<Image>::open() {
    myImage = erosion();
    Image dilat = dilation();
    return dilat;
}

template<typename Image>
Image Morphomaths<Image>::close() {
    myImage = dilation();
    Image eros = erosion();
    return eros;
}

#endif
