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
#include "vcm/skeleton/SkeletonizationOrthogonalPlanes.h"
#include "vcm/skeleton/post/NoPostProcessingSkeleton.h"
#include "geometry/MedialAxis.h"
#include "vcm/OrthogonalPlaneEstimator.h"
#include "geometry/predicate/AbovePlanePredicate.h"
#include "vcm/skeleton/post/JunctionProcessingSkeleton.h"


using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


int main( int  argc, char**  argv )
{
        typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image;
        typedef functors::BallConstantPointFunction<Z3i::Point,double> KernelFunction;
        typedef OrthogonalPlaneEstimator<Z3i::DigitalSet, KernelFunction> OrthoPlaneEstimator;
        typedef AbovePlanePredicate<Z3i::Space> Predicate;
        typedef JunctionProcessingSkeleton<Z3i::DigitalSet, Predicate> PostProcessing;
        typedef SSIJunctionDetection<Z3i::DigitalSet> JunctionDetection;
//typedef NoPostProcessingSkeleton<Z3i::DigitalSet> PostProcessing;

        po::options_description general_opt("Allowed options are: ");
        general_opt.add_options()
                ("help,h", "display this message")
                ("input,i", po::value<std::string>(), "vol file (corresponding volume)")
                ("output,o", po::value<std::string>(), "output skeleton filename")
                 ("radiusInside,R", po::value<double>()->default_value(7), "Big R (radius to trim Voronoi cells)")
                ("thresholdMin,m", po::value<int>()->default_value(1), "minimum threshold for binarization")
                ("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
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
        string outFilename  = vm["output"].as<std::string>();
        int thresholdMin = vm["thresholdMin"].as<int>();
        int thresholdMax = vm["thresholdMax"].as<int>();
        double R = vm["radiusInside"].as<double>();

        Image volume = VolReader<Image>::importVol(inputFilename);
        Z3i::Domain domainVolume = volume.domain();
        Z3i::DigitalSet setVolume(domainVolume);
        SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume,
                                                      thresholdMin-1, thresholdMax);

        QApplication application(argc,argv);
        Viewer3D<> viewer;
        viewer.show();

        JunctionDetection jd(setVolume);
        SkeletonizationOrthogonalPlanes<Z3i::DigitalSet> skeletonization(setVolume, jd, R);
        Z3i::DigitalSet skeleton = skeletonization.skeletonize();

        viewer << CustomColors3D(Color::Red, Color::Red) << skeleton;
        viewer << CustomColors3D(Color(210,210,210,20), Color(210,210,210,20)) << setVolume;
        for (auto it = setVolume.begin(), ite = setVolume.end(); it != ite; ++it) {
                if (volume(*it) >= thresholdMin)
                        viewer << CustomColors3D(Color(0,0,255,20), Color(0,0,255,20))<<*it;
        }
        viewer << CustomColors3D(Color(220,220,220,20), Color(220,220,220,20)) << setVolume;

        Image outImage(volume.domain());

        DGtal::imageFromRangeAndValue(skeleton.begin(), skeleton.end(), outImage, 10);
        VolWriter<Image>::exportVol(outFilename, outImage);

        viewer << Viewer3D<>::updateDisplay;
        application.exec();
        return 0;


}
