#include <CompactNSearch>

#include <iostream>
#include <fstream>
#include "json.hpp"
#include "helpers.hpp"
#include "sphsolver.hpp"
#include "particleAttrib.hpp"
#include "filewriter.hpp"
#include "poissonsample.hpp"

using json = nlohmann::json;
using namespace CompactNSearch;

int main ( void ){

	std::cout << ">>> Starting Program." << std::endl;
	std::map<std::string,ParticleAttributes*> pData;

	json simDataInput, fluidDataInput, boundaryDataInput;
	readJSON(simDataInput,      "../inputs/constants.json");
	readJSON(fluidDataInput,    "../inputs/fluid.json");
	readJSON(boundaryDataInput, "../inputs/boundary.json");

	if( simDataInput["operation"] == "SPHSimulation" ){
		
		pData["fluid"] = new ParticleAttributes(simDataInput, fluidDataInput);
		pData["boundary"] = new ParticleAttributes(simDataInput, boundaryDataInput);


		FileWriter writer(simDataInput["outputFile"].dump());

		SPHSolver solver(simDataInput, pData);
			for(int t=0;t<(Uint)simDataInput["steps"];t++){
				solver.neighborSearch();
				solver.marchTime(t);

				if (t % (Uint)simDataInput["outper"] == 0)
					writer.write(pData);
				
				if (t % (Uint)simDataInput["newfileper"] == 0)
					writer.openNextFile();

				std::cout << "End of timestep : " << t << std::endl;
			}

		writer.close();

	} else if ( simDataInput["operation"] == "GeometryGeneration" ){

		Real r = simDataInput["dropletRadius"];
		Real3 center = Real3{(Real)simDataInput["dropletCenter"][0],(Real)simDataInput["dropletCenter"][1],(Real)simDataInput["dropletCenter"][2]};
		Real3 x_min = sub(center,mult(r,Real3{1.0,1.0,1.0}));
		Real3 x_max = add(center,mult(r,Real3{1.0,1.0,1.0}));
		std::cout << "... Generating Initial Configuration..." << std::endl;

		std::vector<Real3> samples = thinks::poissonDiskSampling(0.8 * (float)simDataInput["dx"], x_min, x_max);
		std::vector<Real> x,y,z;

		int totalParticles = 0;
		for (int i = 0; i < samples.size(); i++ ){
			if(length(sub(samples[i],center)) > r) continue;
			x.push_back( samples[i][0] );
			y.push_back( samples[i][1] );
			z.push_back( samples[i][2] );
			totalParticles ++;
		}
		std::string fname = simDataInput["geomOutputFile"].dump();
		fname = fname.substr(1, fname.size()); fname = fname.substr(0, fname.size()-1);
		H5PartFile* fileWriter = H5PartOpenFile(fname.c_str(),H5PART_WRITE);

		H5PartSetStep(fileWriter,0);
		H5PartSetNumParticles(fileWriter,totalParticles);

		H5PartWriteDataFloat64(fileWriter,"x",&x[0]);
		H5PartWriteDataFloat64(fileWriter,"y",&y[0]);
		H5PartWriteDataFloat64(fileWriter,"z",&z[0]);
		H5PartCloseFile(fileWriter);
		std::cout << ">>> Geometry Generation Complete. Rejection Sampled to " << x.size() << " samples." << std::endl;
	}

	std::cout << ">>> Exiting Program." << std::endl;
	return 0;

}


