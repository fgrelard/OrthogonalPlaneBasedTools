#ifndef ORTHOGONAL_PLANE_ESTIMATOR_H
#define ORTHOGONAL_PLANE_ESTIMATOR_H

/**
 * @file OrthogonalPlaneEstimator.h
 * @author Florent Grelard (florent.grelard@labri.fr)
 * LaBRI, Bordeaux University
 *
 * @date 2017/01/20
 *
 */

#include "vcm/VCMAdjustableRadius.h"
#include "shapes/DigitalPlane.h"

/**
   * Description of template class 'OrthogonalPlaneEstimator' <p>
   * \brief Aim: This class aims at estimating orthogonal planes
   * using VoronoiCovarianceMeasure. Estimation can be performed
   * directly on a set of voxels.
   *
   * You may obtain the planes at a given point by using the
   * method \ref planeAt or \ref convergencePlaneAt.
   *
   *
   * @tparam Container type of Digital Set (model of CDigitalSet).
   *
   * @tparam KernelFunction  the type of a functor
   * Point->Scalar. For instance functors::HatPointFunction and
   * functors::BallConstantPointFunction are models of this type.
   *
   */

template<typename TContainer, typename KernelFunction>
class OrthogonalPlaneEstimator {
private:
    static TContainer emptyContainer;
public:
    typedef TContainer Container;
    typedef typename Container::value_type Point;
    typedef DGtal::SpaceND<Point::dimension, DGtal::int32_t> Space;
    typedef typename Space::RealVector RealVector;
    typedef typename Space::RealPoint RealPoint;
    typedef DGtal::HyperRectDomain<Space> Domain;
    typedef DGtal::ExactPredicateLpSeparableMetric<Space, 2> L2Metric;
    typedef VCMAdjustableRadius<Space, L2Metric> VCM;
    typedef DigitalPlane<Space> Plane;
    typedef DGtal::EigenDecomposition<Point::dimension, double> LinearAlgebraTool;

public:
    OrthogonalPlaneEstimator() = delete;

    /**
    * Constructor.
    *
    * @param container the input volume on which to estimate
    * orthogonal planes
    *
    * @param chi the functor Point->Scalar
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
    * @param aConnexity Digital plane connexity (6, 18, 26)
    *
    * @param aMetric an instance of the metric.
    * @param verbose if 'true' displays information on ongoing computation.
    *
    */
    OrthogonalPlaneEstimator(const Container &container, const KernelFunction &chi, double R, double r,
                             int aConnexity = 26, const L2Metric &l2Metric = L2Metric(), bool verbose = false);

    OrthogonalPlaneEstimator(const OrthogonalPlaneEstimator &other);

    ~OrthogonalPlaneEstimator();


// ----------------------- Interface  --------------------------------------
public:
    /** @return Orthogonal plane at point
     * @param point the point where to estimate the plane
     *
     * @param dirVector the direction of a half ball for point
     * integration (no vector by default = full ball)
     *
     * @param points points for integration (no points by default)
     */
    Plane planeAt(const Point &point,
                  const RealVector &dirVector = RealVector::zero,
                  const Container &points = emptyContainer) const;

    /** @return Orthogonal plane at point
     * @param point the point where to estimate the plane
     *
     * @param volume corresponding volume where plane intersect
     *
     * @param maxRadius maximum radius for convergence estimation
     *
     * @param dirVector the direction of a half ball for point
     * integration (no vector by default = full ball)
     *
     * @param points points for integration (no points by default)
     */
    Plane convergentPlaneAt(const Point &point,
                            const Container &volume,
                            double maxRadius,
                            const RealVector &dirVector = RealVector::zero,
                            const Container &points = emptyContainer);

    void setRadius(double radius);

    double getRadius() const { return myVCM->r(); }

public:
    OrthogonalPlaneEstimator &operator=(const OrthogonalPlaneEstimator &other);

    // ----------------------- Internal methods --------------------------------------
private:
    bool isConvergenceReached(const Container &volume,
                              const Plane &digPlane,
                              double currentRadius,
                              const RealVector &dirVector = RealVector::zero);

    // ------------------------- Private Datas --------------------------------

protected:
    VCM *myVCM;
    KernelFunction myChi;
    int myConnexity;
};

template<typename Container, typename KernelFunction>
Container OrthogonalPlaneEstimator<Container, KernelFunction>::emptyContainer = Container(
        Domain(Point(0, 0, 0), Point(0, 0, 0)));


template<typename Container, typename KernelFunction>
OrthogonalPlaneEstimator<Container, KernelFunction>::
OrthogonalPlaneEstimator(const Container &container, const KernelFunction &chi, double R, double r, int aConnexity,
                         const L2Metric &l2Metric, bool verbose) : myChi(chi), myConnexity(aConnexity) {
    myVCM = new VCM(R, r, l2Metric, verbose);
    myVCM->init(container.begin(), container.end());
}

template<typename Container, typename KernelFunction>
OrthogonalPlaneEstimator<Container, KernelFunction>::
OrthogonalPlaneEstimator(const OrthogonalPlaneEstimator &other) : myChi(other.myChi) {
    L2Metric l2Metric;
    myVCM = new VCM(other.myVCM->R(), other.myVCM->r(), l2Metric, false);
    Container set = other.myVCM->pointSet();
    myVCM->init(set.begin(), set.end());
    myConnexity = other.myConnexity;
}

template<typename Container, typename KernelFunction>
OrthogonalPlaneEstimator<Container, KernelFunction>::
~OrthogonalPlaneEstimator() {
    if (myVCM != 0) {
        delete myVCM;
        myVCM = 0;
    }
}

template<typename Container, typename KernelFunction>
OrthogonalPlaneEstimator<Container, KernelFunction> &
OrthogonalPlaneEstimator<Container, KernelFunction>::
operator=(const OrthogonalPlaneEstimator &other) {
    L2Metric l2Metric;
    myVCM = new VCM(other.myVCM->R(), other.myVCM->r(), l2Metric, false);
    Container set = myVCM->pointSet();
    myChi = other.myChi;
    myConnexity = other.myConnexity;
    return *this;
}

template<typename Container, typename KernelFunction>
typename OrthogonalPlaneEstimator<Container, KernelFunction>::Plane
OrthogonalPlaneEstimator<Container, KernelFunction>::
planeAt(const Point &point, const RealVector &dirVector, const Container &points) const {

    typename LinearAlgebraTool::Matrix vcm_r, evec;
    RealVector eval;
    // Compute VCM and diagonalize it.
    if (dirVector != RealVector::zero)
        vcm_r = myVCM->measureJunction(dirVector, myChi, point);
    else if (!points.empty())
        vcm_r = myVCM->measure(points, myChi, point);
    else
        vcm_r = myVCM->measure(myChi, point);
    LinearAlgebraTool::getEigenDecomposition(vcm_r, evec, eval);
    // Display normal
    RealVector normal = evec.column(0);
    Plane plane(point, normal, myConnexity);
    return plane;
}

template<typename Container, typename KernelFunction>
bool
OrthogonalPlaneEstimator<Container, KernelFunction>::
isConvergenceReached(const Container &volume,
                     const Plane &plane,
                     double currentRadius,
                     const RealVector &dirVector) {

    L2Metric l2Metric;
    bool alright = true;
    typename Plane::DigitalSet intersection = plane.intersectionWithSetOneCC(volume);
    if (dirVector == RealVector::zero) {
        for (const Point &p : intersection) {
            if (l2Metric(p, plane.getCenter()) >= currentRadius) {
                alright = false;
            }
        }
    } else {
        Point current = plane.getCenter();
        int scalar = 1;
        auto itSetVolume = std::find(volume.begin(), volume.end(), current);
        while (itSetVolume != volume.end()) {
            current = Point(RealPoint(plane.getCenter()) + dirVector * scalar);
            scalar++;
            itSetVolume = std::find(volume.begin(), volume.end(), current);
        }
        double distance = l2Metric(current, plane.getCenter());
        alright = (distance < currentRadius);
    }
    return alright;
}

template<typename Container, typename KernelFunction>
typename OrthogonalPlaneEstimator<Container, KernelFunction>::Plane
OrthogonalPlaneEstimator<Container, KernelFunction>::
convergentPlaneAt(const Point &point,
                  const Container &volume,
                  double maxRadius,
                  const RealVector &dirVector,
                  const Container &points) {

    bool isConvergent = false;
    double currentRadius = myVCM->r();
    Plane convergentPlane;
    do {
        currentRadius++;
        setRadius(currentRadius);
        convergentPlane = planeAt(point, dirVector, points);
        isConvergent = isConvergenceReached(volume,
                                            convergentPlane,
                                            currentRadius,
                                            dirVector);

    } while (!isConvergent && currentRadius < maxRadius);

    return convergentPlane;
}

template<typename Container, typename KernelFunction>
void
OrthogonalPlaneEstimator<Container, KernelFunction>::
setRadius(double radius) {
    myVCM->setMySmallR(radius);
    myChi = KernelFunction(1.0, radius);
}

#endif
