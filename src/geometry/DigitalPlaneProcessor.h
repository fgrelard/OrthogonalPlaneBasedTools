#ifndef DIGITAL_PLANE_PROCESSOR_H
#define DIGITAL_PLANE_PROCESSOR_H

#include "shapes/DigitalPlane.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/images/ConstImageAdapter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "geometry/SetProcessor.h"

template<typename TSpace>
class DigitalPlaneProcessor {
    typedef DigitalPlane<TSpace> DigPlane;
    typedef DGtal::ExactPredicateLpSeparableMetric<TSpace, 2> L2Metric;
    typedef typename TSpace::RealVector RealVector;
    typedef typename TSpace::Point Point;
    typedef DGtal::HyperRectDomain<TSpace> Domain;
    typedef typename DGtal::DigitalSetSelector<Domain, DGtal::BIG_DS + DGtal::HIGH_BEL_DS>::Type DigitalSet;
    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image;
    typedef typename Image::Value Value;

    typedef DGtal::SpaceND<Point::dimension - 1, DGtal::int32_t> SubSpace;
    typedef DGtal::HyperRectDomain<SubSpace> SubDomain;
    typedef DGtal::ImageContainerBySTLVector<SubDomain, unsigned char> SubImage;
    typedef typename SubSpace::Point SubPoint;

    typedef DGtal::ConstImageAdapter<Image, SubDomain, DGtal::functors::Point2DEmbedderIn3D<Domain>, Value, DGtal::functors::Identity> ImageAdapter;

public:
    DigitalPlaneProcessor() : myDigitalPlane(DigPlane()) {}

    DigitalPlaneProcessor(const DigPlane &aDigitalPlane) : myDigitalPlane(aDigitalPlane) {}

    DigitalPlaneProcessor(const DigitalPlaneProcessor &other) : myDigitalPlane(other.myDigitalPlane) {}

public:
    std::vector<RealVector> planeToQuadrangle();

    SubImage sliceFromPlane(const Image &image, int patch_width);

private:
    DigPlane myDigitalPlane;
};

template<typename TSpace>
std::vector<typename DigitalPlaneProcessor<TSpace>::RealVector>
DigitalPlaneProcessor<TSpace>::planeToQuadrangle() {
    Point origin = myDigitalPlane.getCenter();
    RealVector normal = myDigitalPlane.getPlaneEquation().normal();

    double a = normal[0];
    double b = normal[1];
    double c = normal[2];
    double d = -a * origin[0] - b * origin[1] - c * origin[2];

    RealVector pRefOrigin;
    if (a != 0) {
        pRefOrigin[0] = -d / a;
        pRefOrigin[1] = 0.0;
        pRefOrigin[2] = 0.0;
        if (pRefOrigin == origin) {
            pRefOrigin[1] = -1.0;
        }
    } else if (b != 0) {
        pRefOrigin[0] = 0.0;
        pRefOrigin[1] = -d / b;
        pRefOrigin[2] = 0.0;
        if (pRefOrigin == origin) {
            pRefOrigin[0] = -1.0;
        }
    } else if (c != 0) {
        pRefOrigin[0] = 0.0;
        pRefOrigin[1] = 0.0;
        pRefOrigin[2] = -d / c;
        if (pRefOrigin == origin) {
            pRefOrigin[0] = -1.0;
        }
    }
    RealVector uDir1;
    uDir1 = (pRefOrigin - origin) / ((pRefOrigin - origin).norm());
    RealVector uDir2;
    uDir2[0] = uDir1[1] * c - uDir1[2] * b;
    uDir2[1] = uDir1[2] * a - uDir1[0] * c;
    uDir2[2] = uDir1[0] * b - uDir1[1] * a;

    uDir2 /= uDir2.norm();


    RealVector myFirstAxisEmbeddedDirection = -uDir1;
    RealVector mySecondAxisEmbeddedDirection = -uDir2;

    std::vector<RealVector> fourPointsForPlane;
    RealVector p1, p2, p3, p4;
    p1 = {(RealVector) -myFirstAxisEmbeddedDirection - mySecondAxisEmbeddedDirection};
    p2 = {(RealVector) -myFirstAxisEmbeddedDirection + mySecondAxisEmbeddedDirection};
    p3 = {(RealVector) myFirstAxisEmbeddedDirection + mySecondAxisEmbeddedDirection};
    p4 = {(RealVector) myFirstAxisEmbeddedDirection - mySecondAxisEmbeddedDirection};

    fourPointsForPlane = {p1, p2, p3, p4};
    return fourPointsForPlane;
}

template<typename TSpace>
typename DigitalPlaneProcessor<TSpace>::SubImage
DigitalPlaneProcessor<TSpace>::sliceFromPlane(const typename DigitalPlaneProcessor<TSpace>::Image &image,
                                              int patch_width) {

    Point origin = myDigitalPlane.getCenter();
    RealVector normal = myDigitalPlane.getPlaneEquation().normal();
    Domain domain3Dyup(image.domain().lowerBound() + Point(-patch_width, -patch_width, -patch_width),
                       image.domain().upperBound() + Point(patch_width, patch_width, patch_width));
    SubDomain domainImage2D(SubPoint(0, 0),
                            SubPoint(patch_width, patch_width));
    DGtal::functors::Identity idV;
    DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain> embedder(domain3Dyup, origin, normal, patch_width,
                                                                      domain3Dyup.lowerBound());

    ImageAdapter extractedImage(image, domainImage2D, embedder, idV);
    SubImage outImage(extractedImage.domain());

    for (const typename ImageAdapter::Point &p : extractedImage.domain())
        outImage.setValue(p, extractedImage(p));
    return outImage;
}


#endif
