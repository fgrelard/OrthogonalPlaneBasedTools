#ifndef SHAPE_DESCRIPTOR_H
#define SHAPE_DESCRIPTOR_H

#include <vector>
#include <algorithm>

#include "DGtal/base/Common.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include <Eigen/Dense>
#include "DGtal/math/Histogram.h"
#include "DGtal/math/Statistic.h"
#include "geometry/Distance.h"
#include "DGtal/images/CImage.h"


template<typename Container>
class ShapeDescriptor {
public:
    typedef typename Container::Space Space;
    typedef typename Space::Point Point;
    typedef typename Space::RealVector RealVector;
    typedef typename Space::RealPoint RealPoint;
    typedef typename DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;
    typedef Eigen::MatrixXd Matrix;
public:
    ShapeDescriptor() = delete;

    ShapeDescriptor(const Container &aData) : myData(aData) {}

public:

    std::pair<Point, Point> majorAxis();

    double lengthMajorAxis();

    template<typename Image>
    RealPoint centerOfMass(const Image &image);

    RealPoint extractCenterOfMass();

    RealVector computeNormalFromLinearRegression();

    RealVector computeNormalFromCovarianceMatrix();

    RealVector extractEigenVector(const Matrix &m, int colNumber);

    double extractEigenValue(const Matrix &m, int colNumber);

    Matrix computeCovarianceMatrix();

    template<typename Image2D>
    Matrix computeCovarianceMatrixImage(const Image2D &image);

private:
    Container myData;
};


template<typename Container>
std::pair<typename ShapeDescriptor<Container>::Point, typename ShapeDescriptor<Container>::Point>
ShapeDescriptor<Container>::majorAxis() {
    Point p1, p2;
    L2Metric l2Metric;
    double distanceFarthestPoint = 0;
    for (const Point &p : myData) {
        for (const Point &o : myData) {
            double currentDistance = l2Metric(p, o);
            if (l2Metric(p, o) > distanceFarthestPoint) {
                distanceFarthestPoint = currentDistance;
                p1 = p;
                p2 = o;
            }
        }
    }
    std::pair<Point, Point> pair = std::make_pair(p1, p2);
    return pair;
}


template<typename Container>
double
ShapeDescriptor<Container>::lengthMajorAxis() {
    L2Metric l2Metric;
    std::pair<Point, Point> farthest = majorAxis();
    double distanceFarthestPoint = l2Metric(farthest.first, farthest.second);
    double radius = distanceFarthestPoint / 2.0;
    return radius;
}


template<typename Container>
template<typename Image>
typename Container::Space::RealPoint ShapeDescriptor<Container>::centerOfMass(const Image &image) {

    BOOST_CONCEPT_ASSERT((DGtal::concepts::CImage<Image>));

    double m000 = 0.0;
    std::vector<double> masses(Container::Space::dimension, 0.0);

    for (auto it = image.domain().begin(), ite = image.domain().end(); it != ite; ++it) {
        typename Container::Space::Point current = *it;
        m000 += image(current);
        for (typename Container::Space::Dimension i = 0; i < Container::Space::dimension; i++) {
            masses[i] += current[i] * image(current);
        }
    }

    if (m000 != 0.0) {
        typename Container::Space::RealPoint v;
        for (typename Container::Space::Dimension i = 0; i < Container::Space::dimension; i++)
            v[i] = masses[i] * 1.0 / m000;
        return v;
    }
    return Container::Space::RealPoint::zero;
}

template<typename Container>
typename Container::Space::RealPoint ShapeDescriptor<Container>::
extractCenterOfMass() {
    BOOST_CONCEPT_ASSERT((DGtal::concepts::CDigitalSet<Container>));

    if (myData.size() != 0) {
        typedef typename DGtal::ImageSelector<typename Container::Domain, unsigned char>::Type Image;
        Image image = DGtal::ImageFromSet<Image>::create(myData, 150);
        typename Container::Space::RealPoint g = centerOfMass(image);
        return g;
    }
    return Container::Space::RealPoint::zero;
}


/**
 * Computes the normal of a plane from a set of points
 * Method : linear regression
 */
template<typename Container>
typename Container::Space::RealVector ShapeDescriptor<Container>::computeNormalFromLinearRegression() {
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXi;
    unsigned int size = myData.size();
    Matrix A(size, 3);
    VectorXi b = VectorXi::Zero(size, 1);

    for (int i = 0; i < size; i++) {
        A(i, 0) = (double) myData[i][0] * 1.0;
        A(i, 1) = (double) myData[i][1] * 1.0;
        A(i, 2) = 1.0;
        b(i, 0) = (double) myData[i][2] * 1.0;
    }
    Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
    typename Container::Space::RealVector normal;
    normal[0] = x(0, 0);
    normal[1] = x(1, 0);
    normal[2] = -1;
    return normal.getNormalized();
}

/**
 * Computes the normal of a plane from a set of points
 * Method : covariance matrix
 */
template<typename Container>
typename Container::Space::RealVector ShapeDescriptor<Container>::computeNormalFromCovarianceMatrix() {

    unsigned int size = myData.size();
    if (size < 2) return (Container::Space::RealVector::zero);

    Matrix A(size, 3);
    for (int i = 0; i < size; i++) {
        A(i, 0) = (double) myData[i][0] * 1.0;
        A(i, 1) = (double) myData[i][1] * 1.0;
        A(i, 2) = (double) myData[i][2] * 1.0;
    }
    Matrix centered = A.rowwise() - A.colwise().mean();
    Matrix cov = (centered.adjoint() * centered) / double(A.rows() - 1);
    Eigen::SelfAdjointEigenSolver<Matrix> eig(cov);
    typename Container::Space::RealVector normal;
    auto veigen = eig.eigenvectors().col(0);
    normal[0] = veigen[0];
    normal[1] = veigen[1];
    normal[2] = veigen[2];
    return normal;
}

template<typename Container>
typename ShapeDescriptor<Container>::Matrix
ShapeDescriptor<Container>::computeCovarianceMatrix() {


    typedef typename Container::ConstIterator ConstIterator;
    typedef typename Container::Domain Domain;
    typedef typename Domain::Point Point;

    int dimens = Point::dimension;
    int size = myData.size();
    Matrix A(size, dimens);
    if (size < dimens) return Matrix(0, 0);

    int i = 0;
    for (ConstIterator it = myData.begin(), ite = myData.end();
         it != ite; ++it) {
        Point point = *it;
        for (int j = 0; j < dimens; j++)
            A(i, j) = (double) point[j] * 1.0;
        i++;
    }
    Matrix centered = A.rowwise() - A.colwise().mean();
    Matrix cov = (centered.adjoint() * centered) / double(A.rows() - 1);
    return cov;
}

template<typename Container>
template<typename Image>
typename ShapeDescriptor<Container>::Matrix
ShapeDescriptor<Container>::computeCovarianceMatrixImage(const Image &image) {
    BOOST_CONCEPT_ASSERT((DGtal::concepts::CImage<Image>));

    typedef typename Image::Domain Domain;
    typedef typename Domain::Point Point;
    int size = 0;
    Container aSet(image.domain());
    for (typename Domain::ConstIterator it = image.domain().begin(), ite = image.domain().end();
         it != ite; ++it) {
        Point point = *it;
        if (image(*it) > 0) {
            size++;
            aSet.insert(*it);
        }
    }
    myData = aSet;
    return computeCovarianceMatrix();
}

template<typename Container>
typename Container::Space::RealVector ShapeDescriptor<Container>::extractEigenVector(const Matrix &m, int colNumber) {


    Eigen::SelfAdjointEigenSolver<Matrix> eig(m);
    typename Container::Space::RealVector vector;
    auto veigen = eig.eigenvectors().col(colNumber);
    for (typename Container::Space::RealVector::Dimension i = 0; i < Container::Space::RealVector::dimension; i++) {
        vector[i] = veigen[i];
    }
    return vector;
}

template<typename Container>
double ShapeDescriptor<Container>::extractEigenValue(const Matrix &m, int columnNumber) {


    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(m);
    typename Container::Space::RealVector vector;
    auto veigen = eig.eigenvalues().col(0);
    if (columnNumber < Space::dimension)
        return veigen[columnNumber];
    return 0.0;
}


#endif
