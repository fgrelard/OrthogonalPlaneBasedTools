#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/readers/ITKReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/writers/ITKWriter.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

int main(int argc, char** argv) {

	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
		("output,o",  po::value<std::string>(), "output itk file" ) ; 

	bool parseOK=true;
	po::variables_map vm;
	try{
		po::store(po::parse_command_line(argc, argv, general_opt), vm);  
	}catch(const std::exception& ex){
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
	string outputFilename = vm["output"].as<std::string>();
	typedef Z3i::Space Space;
	typedef Z3i::KSpace KSpace;
	typedef HyperRectDomain<Space> Domain;
	typedef ImageSelector<Domain, unsigned char>::Type Image;
	Z3i::Point translationVector(0, 0, 0);
	Image image = ITKReader<Image>::importITK(inputFilename);
	Domain aDomain(image.domain().lowerBound() + translationVector, image.domain().upperBound() + translationVector);
	Image out(aDomain);
	
	for (auto it = aDomain.begin(), ite = aDomain.end(); it != ite; ++it) {
		out.setValue(*it, image(*it));
	}
	
	VolWriter<Image>::exportVol(outputFilename, out);
//	ITKWriter<Image>::exportITK(outputFilename, image);


// Z3i::Point upperBound = image.domain().upperBound();
	// Domain domain(Z3i::Point(0,0,0), Z3i::Point(512,512,(upperBound[2]-1)*2));
	// Image output(domain);
	// for (auto it = domain.begin(), ite = domain.end(); it != ite; ++it) {
	// 	Z3i::Point current = *it;
	// 	Z3i::Point inImage(current[0], current[1], current[2] / 2);
	// 	output.setValue(current, image(inImage));
	// }
	return 0;
}
