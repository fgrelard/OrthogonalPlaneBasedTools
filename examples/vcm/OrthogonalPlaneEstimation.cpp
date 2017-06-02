#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "geometry/MedialAxis.h"
#include "vcm/OrthogonalPlaneEstimator.h"
#include "shapes/DigitalPlane.h"
#include "geometry/DigitalPlaneProcessor.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


Z3i::DigitalSet
markPointsBetweenPlanes(const Z3i::DigitalSet& setVolume,
                        const DigitalPlane<Z3i::Space> &currentPlane,
                        const DigitalPlane<Z3i::Space> &previousPlane,
                        double distanceMax) {

    typedef DigitalPlane<Z3i::Space> Plane;
    Z3i::L2Metric l2Metric;
    Z3i::DigitalSet difference(setVolume.domain());
    Z3i::RealVector dirCurrent = currentPlane.getPlaneEquation().normal();
    Z3i::RealVector dirPrevious = previousPlane.getPlaneEquation().normal();
    Z3i::Point pCurrent = currentPlane.getCenter();
    Z3i::Point pPrevious = previousPlane.getCenter();

    if (dirPrevious.dot(dirCurrent) >= 0) {
        dirCurrent = -dirCurrent;
    }

    Plane currentPlane2(pCurrent, dirCurrent, currentPlane.getConnexity());

    for (const Z3i::Point &p : setVolume) {
        if (l2Metric(p, pCurrent) > distanceMax) continue;
        if (currentPlane2.isPointAbove(p) &&
            previousPlane.isPointAbove(p))
            difference.insert(p);
    }
    return difference;
}


int main( int  argc, char**  argv )
{
    typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image;
    typedef functors::BallConstantPointFunction<Z3i::Point,double> KernelFunction;
    typedef OrthogonalPlaneEstimator<Z3i::DigitalSet, KernelFunction> OrthoPlaneEstimator;

    typedef DGtal::ExactPredicateLpSeparableMetric<Z3i::Space, 2> L2Metric;
    typedef DGtal::DistanceTransformation<Z3i::Space, Z3i::DigitalSet, L2Metric> DTL2;
    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
        ("thresholdMin,m", po::value<int>()->default_value(1), "minimum threshold for binarization")
        ("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
        ("radiusInside,R", po::value<double>()->default_value(7), "radius of the ball inside voronoi cell")
        ("radiusNeighbour,r", po::value<double>()->default_value(15), "radius of the ball for the neighbourhood")
        ;

    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    } catch(const std::exception& ex){
        parseOK=false;
        trace.info()<< "Error checking program options: "<< ex.what()<< endl;
    }
    po::notify(vm);
    if( !parseOK || vm.count("help")||argc<=1)
    {
        std::cout << "Usage: " << argv[0] << " [input]\n"
                  << "Display volume file as a voxel set by using QGLviewer"<< endl
                  << general_opt << "\n";
        return 0;
    }
    if(!vm.count("input"))
    {
        trace.error() << " The file name was not defined" << endl;
        return 0;
    }

    string inputFilename = vm["input"].as<std::string>();
    int thresholdMin = vm["thresholdMin"].as<int>();
    int thresholdMax = vm["thresholdMax"].as<int>();
    double R = vm["radiusInside"].as<double>();
    double r = vm["radiusNeighbour"].as<double>();

    Image volume = VolReader<Image>::importVol(inputFilename);
    Z3i::Domain domainVolume = volume.domain();
    Z3i::DigitalSet setVolume(domainVolume);
    SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume,
                                                  thresholdMin-1, thresholdMax);
    QApplication application(argc,argv);
    Viewer3D<> viewer;
    viewer.show();

    double radiusVCM = r;
    KernelFunction chi(1.0, radiusVCM);
    OrthoPlaneEstimator orthogonalPlaneEstimator(setVolume, chi, R, r);
    orthogonalPlaneEstimator.setRadius(radiusVCM);
    int sliceNumber = 0;
    int moduloFactor = setVolume.size() / 1000;
    L2Metric l2Metric;
    double distanceB = 8;
    Z3i::DigitalSet setB(setVolume.domain());
    DTL2 dt(&setVolume.domain(), &setVolume, &l2Metric);
    for (auto it = setVolume.begin(), ite = setVolume.end();
         it != ite; ++it) {
        if (setB.find(*it) != setB.end()) continue;
        sliceNumber++;
        // Compute VCM and diagonalize it.
        viewer.setFillColor(Color::Gray);
        viewer.setFillTransparency(255);

        Z3i::Point current= *it; //it->getPoint();
        double radius = dt(current) + 2;
        orthogonalPlaneEstimator.setRadius(radius);
        DigitalPlane<Z3i::Space> plane = orthogonalPlaneEstimator.convergentPlaneAt(current, setVolume, 100);
        Z3i::RealVector normal = plane.getPlaneEquation().normal();
        DigitalPlane<Z3i::Space> planeBefore(current - distanceB * normal, normal);
        DigitalPlane<Z3i::Space> planeAfter(current + distanceB * normal, -normal);
        Z3i::DigitalSet diff = markPointsBetweenPlanes(setVolume, planeBefore, planeAfter, 100);
        setB.insert(diff.begin(), diff.end());
        DigitalPlaneProcessor<Z3i::Space> planeProc(plane);
        std::vector<Z3i::RealVector> points = planeProc.planeToQuadrangle();
        double f = 20.0;

        viewer.setLineColor(Color::Blue);
        viewer.setFillColor(Color::Blue);
        viewer.setFillTransparency(150);
        viewer.addQuad(current + points[0] * f,
                       current + points[1] * f,
                       current + points[2] * f,
                       current + points[3] * f);


    }

    for (auto it = setVolume.begin(), ite = setVolume.end(); it != ite; ++it) {
        if (volume(*it) >= thresholdMin)
            viewer << CustomColors3D(Color(0,0,255,20), Color(0,0,255,20))<<*it;
    }

    viewer << CustomColors3D(Color(220,220,220,20), Color(220,220,220,20)) << setVolume;
    viewer << Viewer3D<>::updateDisplay;
    application.exec();
    return 0;


}
