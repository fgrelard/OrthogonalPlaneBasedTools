#ifndef VORONOI_COVARIANCE_MEASURE_ADJUSTABLE_RADIUS_H
#define VORONOI_COVARIANCE_MEASURE_ADJUSTABLE_RADIUS_H

// Inclusions
#include <cmath>
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/geometry/volumes/estimation/VoronoiCovarianceMeasure.h"
#include "DGtal/kernel/sets/DigitalSetSelector.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "shapes/Ball.h"

template<typename TSpace, typename TSeparableMetric>
class VCMAdjustableRadius : public DGtal::VoronoiCovarianceMeasure<TSpace, TSeparableMetric> {
public:
    typedef DGtal::VoronoiCovarianceMeasure<TSpace, TSeparableMetric> Base;
    typedef typename Base::Space Space;
    typedef typename Space::RealVector Vector;
    typedef typename Base::MatrixNN MatrixNN;
    typedef typename Base::Point Point;
    typedef typename Base::Metric Metric;
    typedef typename Base::Scalar Scalar;
    typedef typename Base::Domain Domain;
    typedef typename DGtal::DigitalSetSelector<Domain, DGtal::BIG_DS + DGtal::HIGH_BEL_DS>::Type DigitalSet;
    // ----------------------- Standard services ------------------------------
public:

    /**
     * Constructor.
     *
     * @param _R the offset radius for the set of points. Voronoi cells
     * are intersected with this offset. The unit corresponds to a step in the digital space.
     *
     * @param _r (an upper bound of) the radius of the support of
     * forthcoming kernel functions (\f$ \chi_r \f$). The unit
     * corresponds to a step in the digital space. This parameter is
     * used for preparing the data structure that answers to proximity
     * queries.
     *
     * @param aMetric an instance of the metric.
     * @param verbose if 'true' displays information on ongoing computation.
     */
    VCMAdjustableRadius(double _R, double _r, Metric aMetric = Metric(), bool verbose = false) : Base(_R, _r, aMetric,
                                                                                                      verbose),
                                                                                                 myContainer(
                                                                                                         Domain(Point::zero,
                                                                                                                Point::zero)),
                                                                                                 mySmallRAdjustable(
                                                                                                         _r) {};

    /**
     * Destructor.
     */
    ~VCMAdjustableRadius();

    template<typename PointInputIterator>
    void init(PointInputIterator itb, PointInputIterator ite);

    void setMySmallR(double r) { this->mySmallRAdjustable = r; }

    Scalar r() const { return mySmallRAdjustable; }

    DigitalSet pointSet() const { return myContainer; }

    template<typename Point2ScalarFunction>
    MatrixNN measure(const Point2ScalarFunction &chi_r, const Point &p) const;

    template<typename Point2ScalarFunction>
    MatrixNN measure(const DigitalSet &neighbors, Point2ScalarFunction chi_r, Point p) const;

    template<typename Point2ScalarFunction>
    MatrixNN measureJunction(const Vector &dirVector, Point2ScalarFunction chi_r, Point p) const;

    // ----------------------- Interface --------------------------------------

    // ------------------------- Protected Datas ------------------------------
    // ------------------------- Private Datas --------------------------------
protected:
    DigitalSet myContainer;
    double mySmallRAdjustable;

    // ------------------------- Hidden services ------------------------------
protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    VCMAdjustableRadius() : myContainer(Domain(Point::zero, Point::zero)) {}

public:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    VCMAdjustableRadius(const VCMAdjustableRadius &other) : Base(other), myContainer(other.myContainer) {};

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    VCMAdjustableRadius &operator=(const VCMAdjustableRadius &other);

    // ------------------------- Internals ------------------------------------

};

template<typename TSpace, typename TSeparableMetric>
inline
VCMAdjustableRadius<TSpace, TSeparableMetric>::
~VCMAdjustableRadius() {
    this->clean();
}


template<typename TSpace, typename TSeparableMetric>
inline
VCMAdjustableRadius<TSpace, TSeparableMetric> &
VCMAdjustableRadius<TSpace, TSeparableMetric>::
operator=(const VCMAdjustableRadius &other) {
    Base::operator=(other);
    myContainer = other.myContainer;
    return *this;
}

template<typename TSpace, typename TSeparableMetric>
template<typename PointInputIterator>
void
VCMAdjustableRadius<TSpace, TSeparableMetric>::init(PointInputIterator itb,
                                                    PointInputIterator ite) {
    Base::init(itb, ite);
    this->myContainer.insert(itb, ite);
}


template<typename TSpace, typename TSeparableMetric>
template<typename Point2ScalarFunction>
typename VCMAdjustableRadius<TSpace, TSeparableMetric>::MatrixNN
VCMAdjustableRadius<TSpace, TSeparableMetric>::measure(const Point2ScalarFunction &chi_r,
                                                       const typename VCMAdjustableRadius<TSpace, TSeparableMetric>::Point &p) const {

    Ball<Point> ball(p, r());
    DigitalSet neighbors = ball.intersection(this->myContainer);
    MatrixNN vcm;
    // std::cout << *it << " has " << neighbors.size() << " neighbors." << std::endl;
    for (auto it_neighbors = neighbors.begin(),
                 it_neighbors_end = neighbors.end(); it_neighbors != it_neighbors_end; ++it_neighbors) {
        Point q = *it_neighbors;
        Scalar coef = chi_r(q - p);
        if (coef > 0.0) {
            MatrixNN vcm_q = this->vcmMap().at(q);
            vcm_q *= coef;
            vcm += vcm_q;
        }
    }
    return vcm;
}

//-----------------------------------------------------------------------------
template<typename TSpace, typename TSeparableMetric>
template<typename Point2ScalarFunction>
typename VCMAdjustableRadius<TSpace, TSeparableMetric>::MatrixNN
VCMAdjustableRadius<TSpace, TSeparableMetric>::
measure(const DigitalSet &neighbors,
        Point2ScalarFunction chi_r,
        typename VCMAdjustableRadius<TSpace, TSeparableMetric>::Point p) const {

    MatrixNN vcm;
    // std::cout << *it << " has " << neighbors.size() << " neighbors." << std::endl;
    for (auto it_neighbors = neighbors.begin(),
                 it_neighbors_end = neighbors.end(); it_neighbors != it_neighbors_end; ++it_neighbors) {
        Point q = *it_neighbors;
        Scalar coef = chi_r(q - p);
        // Scalar coef = 1.0;
        if (coef > 0.0) {
            typename std::map<Point, MatrixNN>::const_iterator it = this->vcmMap().find(q);
            if (it != this->vcmMap().end()) {
                MatrixNN vcm_q = it->second;
                vcm_q *= coef;
                vcm += vcm_q;
            }
        }
    }
    return vcm;
}


template<typename TSpace, typename TSeparableMetric>
template<typename Point2ScalarFunction>
typename VCMAdjustableRadius<TSpace, TSeparableMetric>::MatrixNN
VCMAdjustableRadius<TSpace, TSeparableMetric>::
measureJunction(const Vector &dirVector, Point2ScalarFunction chi_r, Point p) const {

    typedef typename TSpace::RealVector Vector;

    Ball<Point> ball(p, r());
    DigitalSet neighbors = ball.pointsInHalfBall(dirVector);

    MatrixNN vcm, evec, null;

    Vector eval;
    std::map<Point, Vector> mapPoint2Normal;
    for (auto it_neighbors = neighbors.begin(),
                 it_neighbors_end = neighbors.end(); it_neighbors != it_neighbors_end; ++it_neighbors) {
        Point q = *it_neighbors;
        Scalar coef = chi_r(q - p);
        if (coef > 0.0) {
            typename std::map<Point, MatrixNN>::const_iterator it = this->vcmMap().find(q);
            if (it != this->vcmMap().end()) {
                MatrixNN vcm_q = it->second;
                vcm_q *= coef;
                vcm += vcm_q;
            }
        }
    }
    return vcm;
}


#endif
