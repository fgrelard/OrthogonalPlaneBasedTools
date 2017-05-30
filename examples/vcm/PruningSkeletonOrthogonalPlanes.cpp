#include <iostream>
#include <QtGui/qapplication.h>
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
#include "vcm/PruningOrthogonalPlanes.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


int main( int  argc, char**  argv )
{
        typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image;
        typedef functors::BallConstantPointFunction<Z3i::Point,double> KernelFunction;
        typedef OrthogonalPlaneEstimator<Z3i::DigitalSet, KernelFunction> OrthoPlaneEstimator;
        po::options_description general_opt("Allowed options are: ");
        general_opt.add_options()
                ("help,h", "display this message")
                ("curve,c", po::value<std::string>(), "vol file (curve)")
                ("input,i", po::value<std::string>(), "input vol file (corresponding volume)")
                ("output,o", po::value<std::string>(), "output vol file (corresponding volume")
                ("thresholdPruning,t", po::value<double>()->default_value(25), "threshold for pruning (angle in degrees)")
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

        string curveFilename = vm["curve"].as<std::string>();
        string inputFilename = vm["input"].as<std::string>();
        string outFilename = vm["output"].as<std::string>();
        int thresholdMin = vm["thresholdMin"].as<int>();
        int thresholdMax = vm["thresholdMax"].as<int>();
        double threshold = vm["thresholdPruning"].as<double>();


        Image volume = VolReader<Image>::importVol(inputFilename);
        Z3i::Domain domainVolume = volume.domain();
        Z3i::DigitalSet setVolume(domainVolume);
        SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume,
                                                      thresholdMin-1, thresholdMax);

        Image curve = VolReader<Image>::importVol(curveFilename);
        Z3i::Domain domainCurve = curve.domain();
        Z3i::DigitalSet setCurve(domainCurve);
        SetFromImage<Z3i::DigitalSet>::append<Image> (setCurve, curve,
                                                      thresholdMin-1, thresholdMax);


        QApplication application(argc,argv);
        Viewer3D<> viewer;
        viewer.show();

        PruningOrthogonalPlanes<Z3i::DigitalSet> pruning(setCurve, setVolume, threshold);
        Z3i::DigitalSet skeletonPruned = pruning.prune();

        viewer << CustomColors3D(Color::Red, Color::Red) << skeletonPruned;
        viewer << CustomColors3D(Color::Yellow, Color::Yellow) << setCurve;
        for (auto it = setVolume.begin(), ite = setVolume.end(); it != ite; ++it) {
                if (volume(*it) >= thresholdMin)
                        viewer << CustomColors3D(Color(0,0,255,20), Color(0,0,255,20))<<*it;
        }

        viewer << CustomColors3D(Color(220,220,220,20), Color(220,220,220,20)) << setVolume;

        Image outImage(volume.domain());

        DGtal::imageFromRangeAndValue(skeletonPruned.begin(), skeletonPruned.end(), outImage, 10);
        VolWriter<Image>::exportVol(outFilename, outImage);

        viewer << Viewer3D<>::updateDisplay;
        application.exec();
        return 0;


}
