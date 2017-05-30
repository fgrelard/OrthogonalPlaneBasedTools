#ifndef VCM_ON_DIGITAL_SURFACE_ADJUSTABLE_RADIUS_H
#define VCM_ON_DIGITAL_SURFACE_ADJUSTABLE_RADIUS_H

// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/base/CountedConstPtrOrConstPtr.h"
#include "DGtal/kernel/Point2ScalarFunctors.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/topology/CDigitalSurfaceContainer.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/geometry/volumes/distance/CSeparableMetric.h"
#include "vcm/VCMAdjustableRadius.h"
#include "geometry/MedialAxis.h"
#include "DGtal/math/ScalarFunctors.h"
#include "DGtal/geometry/surfaces/estimation/LocalEstimatorFromSurfelFunctorAdapter.h"
#include "DGtal/geometry/surfaces/estimation/estimationFunctors/ElementaryConvolutionNormalVectorEstimator.h"

/// Possible embeddings for surfel as point(s)
enum Surfel2PointEmbedding {
    Pointels = 0, InnerSpel = 1, OuterSpel = 2
};

/////////////////////////////////////////////////////////////////////////////
// template class VCMOnDigitalSurfaceAdjustableRadius
/**
 * Description of template class
 * 'VCMOnDigitalSurfaceAdjustableRadius' <p> \brief Aim: This
 * class specializes the Voronoi covariance measure for digital
 * surfaces. It adds notably the embedding of surface elements, the
 * diagonalisation of the VCM, and the orientation of the first VCM
 * eigenvector toward the interior of the surface.
 *
 * @note Documentation in \ref moduleVCM_sec3_1.
 *
 * @see VoronoiCovarianceMeasure
 *
 * @tparam TDigitalSurfaceContainer the type of digital surface
 * container (model of CDigitalSurfaceContainer).
 *
 * @tparam TSeparableMetric a model of CSeparableMetric used for
 * computing the Voronoi map (e.g. Euclidean metric is
 * DGtal::ExactPredicateLpSeparableMetric<TSpace, 2> )
 *
 * @tparam TKernelFunction the type of the kernel function chi_r used
 * for integrating the VCM, a map: Point -> Scalar.
 */
template<typename TDigitalSurfaceContainer, typename TSeparableMetric,
        typename TKernelFunction>
class VCMOnDigitalSurfaceAdjustableRadius {
    BOOST_CONCEPT_ASSERT((DGtal::concepts::CDigitalSurfaceContainer<TDigitalSurfaceContainer>));
    BOOST_CONCEPT_ASSERT((DGtal::concepts::CSeparableMetric<TSeparableMetric>));
    // ----------------------- public types ------------------------------
public:
    typedef TDigitalSurfaceContainer DigitalSurfaceContainer; ///< the chosen container
    typedef TSeparableMetric Metric;  ///< the chosen metric
    typedef TKernelFunction KernelFunction;  ///< the kernel function
    typedef DGtal::DigitalSurface<DigitalSurfaceContainer> Surface;  ///< the chosen digital surface
    typedef typename DigitalSurfaceContainer::KSpace KSpace;  ///< the cellular space
    typedef typename DigitalSurfaceContainer::Surfel Surfel;  ///< the n-1 cells
    typedef typename KSpace::SCell SCell;  ///< the signed cells
    typedef typename KSpace::Space Space;  ///< the digital space
    typedef typename KSpace::Point Point;  ///< the digital points
    typedef VCMAdjustableRadius<Space, Metric> VCM;  ///< the Voronoi Covariance Measure
    typedef typename VCM::Scalar Scalar;  ///< the "real number" type
    typedef typename Surface::ConstIterator ConstIterator;  ///< the iterator for traversing the surface
    typedef DGtal::EigenDecomposition<KSpace::dimension, Scalar> LinearAlgebraTool;  ///< diagonalizer (nD).
    typedef typename VCM::VectorN VectorN;  ///< n-dimensional R-vector
    typedef typename VCM::MatrixNN MatrixNN;  ///< nxn R-matrix
    typedef DGtal::HyperRectDomain<Space> Domain;
    typedef typename DGtal::DigitalSetSelector<Domain, DGtal::BIG_DS + DGtal::HIGH_BEL_DS>::Type DigitalSet;


    BOOST_CONCEPT_ASSERT((DGtal::concepts::CUnaryFunctor<KernelFunction, Point, Scalar>));

    /// Structure to hold a diagonalized matrix.
    struct EigenStructure {
        VectorN values;   ///< eigenvalues from the smallest to the biggest
        MatrixNN vectors; ///< corresponding eigenvectors
    };
    /// Structure to hold the normals for each surfel (the VCM one and the trivial one).
    struct Normals {
        VectorN vcmNormal;
        VectorN trivialNormal;
    };
    typedef std::map<Point, EigenStructure> Point2EigenStructure;  ///< the map Point -> EigenStructure
    typedef std::map<Surfel, Normals> Surfel2Normals;    ///< the map Surfel -> Normals
    typedef std::map<Point, double> Point2Radius;
    // ----------------------- Standard services ------------------------------
public:

    /**
     * Destructor.
     */
    ~VCMOnDigitalSurfaceAdjustableRadius();

    /**
     * Constructor. Computes the VCM of the given \a surface.
     *
     * @param _surface the digital surface that is aliased in this. The
     * user can \b secure the aliasing by passing a
     * CountedConstPtrOrConstPtr.
     *
     * @param _surfelEmbedding the chosen embedding for surfels.
     *
     * @param _R the offset radius for the set of points. Voronoi cells
     * are intersected with this offset. The unit corresponds to a step in the digital space.
     *
     * @param _r (an upper bound of) the radius of the support of the
     * kernel function \a chi_r (note \f$\chi_r\f$ in the VCM
     * paper). The unit corresponds to a step in the digital
     * space. This parameter is used for preparing the data structure
     * that answers to proximity queries.
     *
     * @param chi_r the kernel function whose support has radius less
     * or equal to \a r.
     *
     * @param t the radius for the trivial normal estimator, which is
     * used for finding the correct orientation inside/outside for the
     * VCM.
     *
     * @param aMetric an instance of the metric.
     *
     * @param verbose if 'true' displays information on ongoing computation.
     */
    VCMOnDigitalSurfaceAdjustableRadius(DGtal::ConstAlias<Surface> _surface,
                                        Surfel2PointEmbedding _surfelEmbedding,
                                        Scalar _R, Scalar _r,
                                        KernelFunction chi_r,
                                        const MedialAxis<DigitalSet> &aMedialAxis,
                                        Scalar t = 2.5, Metric aMetric = Metric(),
                                        bool verbose = false);

    /// the const-aliased digital surface.
    DGtal::CountedConstPtrOrConstPtr<Surface> surface() const;

    /// the chosen embedding Surfel -> Point(s)
    Surfel2PointEmbedding surfelEmbedding() const;

    /// @return the parameter R in the VCM, i.e. the offset radius for
    /// the compact set K.
    Scalar R() const;

    /// @return the parameter r in VCM(chi_r), i.e. an upper bound for
    /// the diameter of the support of kernel functions.
    Scalar r() const;

    /// @return the radius for the trivial normal estimator, which is
    /// used for finding the correct orientation inside/outside for
    /// the VCM.
    Scalar radiusTrivial() const;

    VCM &vcm();

    /**
       @param[in] outIt an output iterator on Point to write the point(s) associated to surfel \a s.
       @param[in] s the surfel that is embedded in the digital space according to mySurfelEmbedding.
       @return the (modified) output iterator after the write operations.
    */
    template<typename PointOutputIterator>
    PointOutputIterator getPoints(PointOutputIterator outIt, Surfel s) const;

    /// @return a const-reference to the map Surfel -> Normals (vcm and trivial normal).
    const Surfel2Normals &mapSurfel2Normals() const;

    /// @return a const-reference to the map Point ->
    /// EigenStructure of the chi_r VCM (eigenvalues and
    /// eigenvectors).
    const Point2EigenStructure &mapPoint2ChiVCM() const;

    /**
       Gets the eigenvalues of the chi_r VCM at surfel \a s sorted from lowest to highest.
       @param[out] values the eigenvalues of the chi_r VCM at \a s.
       @param[in] s the surfel
       @return 'true' is the surfel \a s was valid.
    */
    bool getChiVCMEigenvalues(VectorN &values, Surfel s) const;

    /**
       Gets the eigen decomposition of the chi_r VCM at surfel \a s.
       @param[out] values the eigenvalues of the chi_r VCM at \a s sorted from lowest to highest..
       @param[out] vectors the eigenvectors of the chi_r VCM at \a s associated to \a values.
       @param[in] s the surfel
       @return 'true' is the surfel \a s was valid.
    */
    bool getChiVCMEigenStructure(VectorN &values, MatrixNN &vectors, Surfel s) const;

    // ----------------------- Interface --------------------------------------
public:


    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay(std::ostream &out) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
protected:
    /// (possibly secure) alias of the digital surface
    DGtal::CountedConstPtrOrConstPtr<Surface> mySurface;
    /// The chosen embedding for the surfels.
    Surfel2PointEmbedding mySurfelEmbedding;
    /// The kernel function chi_r
    KernelFunction myChi;
    /// Stores the voronoi covariance measure of the point embedding of the surface.
    VCM myVCM;
    /// Stores the radius for the trivial normal estimator, which is
    /// used for finding the correct orientation inside/outside for
    /// the VCM.
    Scalar myRadiusTrivial;
    /// Stores for each point p its convolved VCM, i.e. VCM( chi_r( p ) )
    Point2EigenStructure myPt2EigenStructure;
    /// Stores for each surfel its vcm normal and its trivial normal.
    Surfel2Normals mySurfel2Normals;
    /// Medial axis
    MedialAxis<DigitalSet> myMedialAxis;


    // ------------------------- Private Datas --------------------------------
private:

    // ------------------------- Hidden services ------------------------------
protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    VCMOnDigitalSurfaceAdjustableRadius();

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    VCMOnDigitalSurfaceAdjustableRadius(const VCMOnDigitalSurfaceAdjustableRadius &other);

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    VCMOnDigitalSurfaceAdjustableRadius &operator=(const VCMOnDigitalSurfaceAdjustableRadius &other);

    // ------------------------- Internals ------------------------------------
private:

}; // end of class VCMOnDigitalSurfaceAdjustableRadius


/**
 * Overloads 'operator<<' for displaying objects of class 'VCMOnDigitalSurfaceAdjustableRadius'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'VCMOnDigitalSurfaceAdjustableRadius' to write.
 * @return the output stream after the writing.
 */
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
std::ostream &
operator<<(std::ostream &out,
           const VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction> &object);

template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
~VCMOnDigitalSurfaceAdjustableRadius() {
}


template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
VCMOnDigitalSurfaceAdjustableRadius(const VCMOnDigitalSurfaceAdjustableRadius &other) : myMedialAxis(
        other.myMedialAxis) {
}


template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction> &
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
operator=(const VCMOnDigitalSurfaceAdjustableRadius &other) {
    myMedialAxis = other.myMedialAxis;
    return *this;
}

//-----------------------------------------------------------------------------
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
VCMOnDigitalSurfaceAdjustableRadius(DGtal::ConstAlias<Surface> _surface,
                                    Surfel2PointEmbedding _surfelEmbedding,
                                    Scalar _R, Scalar _r,
                                    KernelFunction chi_r,
                                    const MedialAxis<DigitalSet> &aMedialAxis,
                                    Scalar delta, Metric aMetric, bool verbose)
        : mySurface(_surface), mySurfelEmbedding(_surfelEmbedding), myChi(chi_r),
          myVCM(_R, _r, aMetric, verbose), myRadiusTrivial(delta), myMedialAxis(aMedialAxis) {
    const KSpace &ks = this->mySurface->container().space();
    std::vector<Point> vectPoints;

// Get points.
    if (verbose) DGtal::trace.beginBlock("Getting points.");
    std::set<Point> pointSet;
    for (ConstIterator it = this->mySurface->begin(), itE = this->mySurface->end(); it != itE; ++it)
        getPoints(std::inserter(pointSet, pointSet.begin()), *it);
    vectPoints.resize(pointSet.size());
    std::copy(pointSet.begin(), pointSet.end(), vectPoints.begin());
    if (verbose) DGtal::trace.endBlock();

// Compute Voronoi Covariance Matrix for all points.
    myVCM.init(vectPoints.begin(), vectPoints.end());

// Compute VCM( chi_r ) for each point.
    if (verbose) DGtal::trace.beginBlock("Integrating VCM( chi_r(p) ) for each point.");
    int i = 0;
    DigitalSet medialAxisSet = myMedialAxis.pointSet();
    for (auto it = vectPoints.begin(), itE = vectPoints.end();
         it != itE; ++it) {
        Point p = *it;
        if (verbose) DGtal::trace.progressBar(++i, vectPoints.size());
        Point closestPointToCurrent;
        double distance = std::numeric_limits<double>::max();
        for (const Point &m : medialAxisSet) {
            if (aMetric(m, closestPointToCurrent) < distance) {
                distance = aMetric(m, closestPointToCurrent);
                closestPointToCurrent = m;
            }
        }
        std::vector<DGtal::Z3i::RealPoint> normalsToPoint{DGtal::Z3i::RealPoint(0, 0, 1),
                                                          DGtal::Z3i::RealPoint(0, 1, 0),
                                                          DGtal::Z3i::RealPoint(1, 0, 0)};
        double radius = aMetric(closestPointToCurrent, p) + 2.0;
        this->myChi = KernelFunction(1.0, radius);
        myVCM.setMySmallR(radius);


        MatrixNN measure = myVCM.measure(this->myChi, p);
// On diagonalise le résultat.
        EigenStructure &evcm = this->myPt2EigenStructure[p];
        LinearAlgebraTool::getEigenDecomposition(measure, evcm.vectors, evcm.values);
    }
    myVCM.clean(); // free some memory.
    if (verbose) DGtal::trace.endBlock();

    if (verbose) DGtal::trace.beginBlock("Computing average orientation for each surfel.");
    typedef DGtal::functors::HatFunction<Scalar> Functor;
    Functor fct(1.0, this->myRadiusTrivial);
    typedef DGtal::functors::ElementaryConvolutionNormalVectorEstimator<Surfel, DGtal::CanonicSCellEmbedder<KSpace> >
            SurfelFunctor;
    typedef DGtal::LocalEstimatorFromSurfelFunctorAdapter<DigitalSurfaceContainer, Metric, SurfelFunctor, Functor>
            NormalEstimator;

    DGtal::CanonicSCellEmbedder<KSpace> canonic_embedder(ks);
    SurfelFunctor surfelFct(canonic_embedder, 1.0);
    NormalEstimator estimator;
    estimator.attach(*(this->mySurface));
    estimator.setParams(aMetric, surfelFct, fct, this->myRadiusTrivial);
    estimator.init(1.0, this->mySurface->begin(), this->mySurface->end());
    i = 0;
    std::vector<Point> pts;
    int surf_size = this->mySurface->size();
    for (ConstIterator it = this->mySurface->begin(), itE = this->mySurface->end(); it != itE; ++it) {
        if (verbose) DGtal::trace.progressBar(++i, surf_size);
        Surfel s = *it;
        Normals &normals = this->mySurfel2Normals[s];
// get rough estimation of normal
        normals.trivialNormal = estimator.eval(it);
// get points associated with surfel s
        getPoints(std::back_inserter(pts), s);
        for (typename std::vector<Point>::const_iterator itPts = pts.begin(), itPtsE = pts.end();
             itPts != itPtsE; ++itPts) {
            Point p = *itPts;
            const EigenStructure &evcm = this->myPt2EigenStructure[p];
            VectorN n = evcm.vectors.column(2);
            if (n.dot(normals.trivialNormal) < 0) normals.vcmNormal -= n;
            else normals.vcmNormal += n;
        }
        if (pts.size() > 1) normals.vcmNormal /= pts.size();
        pts.clear();
    }
    if (verbose) DGtal::trace.endBlock();


}

//-----------------------------------------------------------------------------
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
DGtal::CountedConstPtrOrConstPtr<typename VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::Surface>
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
surface() const {
    return mySurface;
}

//-----------------------------------------------------------------------------
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
Surfel2PointEmbedding
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
surfelEmbedding() const {
    return mySurfelEmbedding;
}

//-----------------------------------------------------------------------------
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
typename VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::Scalar
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
R() const {
    return myVCM.R();
}

//-----------------------------------------------------------------------------
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
typename VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::Scalar
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
r() const {
    return myVCM.r();
}

template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
typename VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::VCM &
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
vcm() {
    return myVCM;
}

//-----------------------------------------------------------------------------
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
typename VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::Scalar
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
radiusTrivial() const {
    return myRadiusTrivial;
}

//-----------------------------------------------------------------------------
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
const typename VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::Surfel2Normals &
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
mapSurfel2Normals() const {
    return mySurfel2Normals;
}

//-----------------------------------------------------------------------------
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
const typename VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::Point2EigenStructure &
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
mapPoint2ChiVCM() const {
    return myPt2EigenStructure;
}


//-----------------------------------------------------------------------------
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
bool
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
getChiVCMEigenvalues(VectorN &values, Surfel s) const {
    std::vector<Point> pts;
    getPoints(std::back_inserter(pts), s);
    bool ok = true;
    int i = 0;
    values = VectorN(); // Setting values to 0 before averaging.
    for (typename std::vector<Point>::const_iterator itPts = pts.begin(), itPtsE = pts.end();
         itPts != itPtsE; ++itPts, ++i) {
        Point p = *itPts;
        typename Point2EigenStructure::const_iterator itEigen = myPt2EigenStructure.find(p);
        if (itEigen == myPt2EigenStructure.end()) {
            ok = false;
            break;
        }
        const EigenStructure &evcm = itEigen->second;
        values += evcm.values;
    }
    if (i > 1) values /= i;
    return ok;
}

//-----------------------------------------------------------------------------
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
bool
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
getChiVCMEigenStructure(VectorN &values, MatrixNN &vectors, Surfel s) const {
    std::vector<Point> pts;
    getPoints(std::back_inserter(pts), s);
    bool ok = true;
    int i = 0;
    values = VectorN();   // Setting values to 0 before averaging.
    vectors = MatrixNN(); // Setting values to 0 before averaging.
    for (typename std::vector<Point>::const_iterator itPts = pts.begin(), itPtsE = pts.end();
         itPts != itPtsE; ++itPts, ++i) {
        Point p = *itPts;
        typename Point2EigenStructure::const_iterator itEigen = myPt2EigenStructure.find(p);
        if (itEigen == myPt2EigenStructure.end()) {
            ok = false;
            break;
        }
        const EigenStructure &evcm = itEigen->second;
        values += evcm.values;
        vectors += evcm.vectors;
    }
    if (i > 1) {
        values /= i;
        vectors /= i;
    }
    return ok;
}

//-----------------------------------------------------------------------------
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
template<typename PointOutputIterator>
inline
PointOutputIterator
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
getPoints(PointOutputIterator outIt, Surfel s) const {
    BOOST_CONCEPT_ASSERT((boost::OutputIterator<PointOutputIterator, Point>));
    const KSpace &ks = mySurface->container().space();
    DGtal::Dimension k = ks.sOrthDir(s);
    switch (mySurfelEmbedding) {
        case Pointels: {
            typename KSpace::Cells faces = ks.uFaces(ks.unsigns(s));
            for (typename KSpace::Cells::const_iterator it = faces.begin(), itE = faces.end();
                 it != itE; ++it) {
                if (ks.uDim(*it) == 0)  // get only pointels (cell of dim 0)
                    *outIt++ = ks.uCoords(*it);
            }
            // Dimension i = (k+1)%3;
            // Dimension j = (i+1)%3;
            // SCell l1 = ks.sIncident( s, i, true );
            // SCell l2 = ks.sIncident( s, i, false );
            // *outIt++ = ks.sCoords( ks.sIncident( l1, j, true ) );
            // *outIt++ = ks.sCoords( ks.sIncident( l1, j, false ) );
            // *outIt++ = ks.sCoords( ks.sIncident( l2, j, true ) );
            // *outIt++ = ks.sCoords( ks.sIncident( l2, j, false ) );
        }
            break;
        case InnerSpel:
            *outIt++ = ks.sCoords(ks.sDirectIncident(s, k));
            break;
        case OuterSpel:
            *outIt++ = ks.sCoords(ks.sIndirectIncident(s, k));
            break;
    }
    return outIt;
}

///////////////////////////////////////////////////////////////////////////////
// Interface - public :

/**
 * Writes/Displays the object on an output stream.
 * @param out the output stream where the object is written.
 */
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
void
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
selfDisplay(std::ostream &out) const {
    out << "[VCMOnDigitalSurfaceAdjustableRadius"
        << " #pts=" << myPt2EigenStructure.size()
        << " #surf=" << mySurfel2Normals.size()
        << "]";
}

/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
bool
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
isValid() const {
    return true;
}



///////////////////////////////////////////////////////////////////////////////
// Implementation of inline itktools                                        //

template<typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
inline
std::ostream &
operator<<(std::ostream &out,
           const VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction> &object) {
    object.selfDisplay(out);
    return out;
}


#endif
