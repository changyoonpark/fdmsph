#include "particleAttrib.hpp"
#include <cstdlib>
#include <H5Part.h>
#include <assert.h>


using std::ofstream;


ParticleAttributes::ParticleAttributes(const json& _simDataIn, const json& _parDataIn)
 : simDataIn(_simDataIn), parDataIn(_parDataIn)
{
	initParticles();
}



void ParticleAttributes::readInitialPlacement(std::string fileName){

	std::cout << "... Reading Particle Placements From " << fileName << std::endl;
	H5PartFile *fileReader = H5PartOpenFile(fileName.c_str(),H5PART_READ);
	H5PartSetStep(fileReader,0);
	Uint nparticles = H5PartGetNumParticles(fileReader);
	Uint nsteps     = H5PartGetNumSteps(fileReader);
	std::cout << "... input file (.h5part) has " << nsteps << " timesteps with " << nparticles << " particles." << std::endl;

	std::vector<double> x,y,z;
	x.resize(nparticles);
	y.resize(nparticles);
	z.resize(nparticles);

	H5PartReadDataFloat64(fileReader,"x",&x[0]);
	H5PartReadDataFloat64(fileReader,"y",&y[0]);
	H5PartReadDataFloat64(fileReader,"z",&z[0]);

	int n = 0;
  	smoothingLength = simDataIn["smoothingLength"];
	numParticles = 0;
	for(int i=0;i<x.size();i++){

		Real3 posToAdd{x[i],y[i],z[i]};
		addDefaultFluidParticleAtPosition(posToAdd,getT0());
		n+=1;
	}

	H5PartCloseFile(fileReader);
	std::cout << "--- Particle Data Loaded : " << n << " Particles from input file." << std::endl;
}

void ParticleAttributes::boundaryInit(){

	numParticles = 0;
	smoothingLength = simDataIn["smoothingLength"];
	dx              = simDataIn["dx"];


std::cout << "... Computing Volume that should be used (2D) " << std::endl;

	Real kernelSum = 0.0, dx = (Real)simDataIn["dx"];		
	int neighs = 0;
	for(int i = -20; i <= 20; i++){
		for(int j = -20; j <= 20; j++){
			Real xTemp = ((Real) i) * dx;
			Real yTemp = ((Real) j) * dx;
			Real3 relPos{xTemp, yTemp, 0};
			Real dist = length(relPos);
			// std::cout << dist << std::endl;
			if ( dist < ((Real) simDataIn["smoothingLength"]) ){
				kernelSum += W_Wendland_2D_2h(dist, (Real) simDataIn["smoothingLength"]);
				neighs ++;
			}
		}
	}	

	std::cout << "... Number of neighbors for the center particle : " << neighs << std::endl;
	defaultVolume = 1.0 / kernelSum;	
	defaultMass = defaultVolume * (Real) parDataIn["rho0"];

	std::cout << "... The default volume & mass used for the boundary is v_o = " << defaultVolume << ", m_o = " << defaultMass << std::endl;

	for(auto boundaryEntry : parDataIn["objects"]){
		int n = 0;


		// Initialize Frame Type Solid Boundary.
		if (boundaryEntry["type"] == "frame"){

			int offset = 0;

			if (boundaryEntry["includeBoundary"] == "True"){
				offset = 0;
			} else{
				offset = 1;
			}

			Real x0 = boundaryEntry["begin"][0];
			Real y0 = boundaryEntry["begin"][1];
			Real z0 = boundaryEntry["begin"][2];
			Real x1 = boundaryEntry["end"][0];
			Real y1 = boundaryEntry["end"][1];
			Real z1 = boundaryEntry["end"][2];

			Real wallt = boundaryEntry["thickness"];
			Real rho0 = 0, T0 = 0;
			if (parDataIn.find("rho0") != parDataIn.end() ){
				rho0 = parDataIn["rho0"];
			}
			if (boundaryEntry.find("T0") != boundaryEntry.end() ){
				T0 = boundaryEntry["T0"];
			}

			if (simDataIn["dimensions"] == 3){

				for (int i = offset; i <= (int)((x1 - x0)/dx) - offset; i ++ ){
				for (int j = offset; j <= (int)((y1 - y0)/dx) - offset; j ++ ){
				for (int k = offset; k <= (int)((z1 - z0)/dx) - offset; k ++ ){

					Real3 posToAdd{x0 + dx * i, y0 + dx * j, z0 + dx * k};
					Real3 zeroVector{0,0,0};
					Real3x3 zeromat{zeroVector,zeroVector,zeroVector};


					if ((posToAdd[0] > x0 + wallt + offset * dx || posToAdd[0] < x1 - wallt - offset * dx) &&
					    (posToAdd[1] > y0 + wallt + offset * dx || posToAdd[1] < y1 - wallt - offset * dx) &&
					    (posToAdd[2] > z0 + wallt + offset * dx || posToAdd[2] < z1 - wallt - offset * dx) ) continue;

					addDefaultFluidParticleAtPosition(posToAdd,getT0());
					n+=1;

				}}}

			}
			else if (simDataIn["dimensions"] == 2){

				for (int i = offset; i <= (int)((x1 - x0)/dx) - offset; i ++ ){
				for (int k = offset; k <= (int)((z1 - z0)/dx) - offset; k ++ ){

					Real3 posToAdd{x0 + dx * i, 0, z0 + dx * k};

					if ( (posToAdd[0] >= x0 + wallt - 0.99 * offset * dx && posToAdd[0] <= x1 - wallt + 0.99 * offset *  dx ) &&
					     (posToAdd[2] >= z0 + wallt - 0.99 * offset * dx && posToAdd[2] <= z1 - wallt + 0.99 * offset *  dx ) ) continue;

					addDefaultFluidParticleAtPosition(posToAdd,getT0());
					n+=1;

				}}

			} else{
				std::cout << "Dimensions Not Specified. Check Input File." << std::endl;
				assert(0);
			}
		}

		// Initialize Block Type Solid Boundary.
		else if(boundaryEntry["type"] == "block"){

			int offset = 0, offset2;
			if (boundaryEntry["includeBoundary"] == "True"){
				offset = 0;
			} else{
				offset = 1;
			}


			Real x0 = (Real)boundaryEntry["begin"][0];
			Real y0 = (Real)boundaryEntry["begin"][1];
			Real z0 = (Real)boundaryEntry["begin"][2];
			Real x1 = (Real)boundaryEntry["end"][0];
			Real y1 = (Real)boundaryEntry["end"][1];
			Real z1 = (Real)boundaryEntry["end"][2];
			Real rho0 = 0, T0 = 0;

			if (parDataIn.find("rho0") != parDataIn.end() ){
				rho0 = parDataIn["rho0"];
			}
			if (boundaryEntry.find("T0") != boundaryEntry.end() ){
				T0 = boundaryEntry["T0"];
			}

      std::mt19937 rng;
      rng.seed(std::random_device()());
      std::uniform_int_distribution<std::mt19937::result_type> dist6(0,100000);

			if (simDataIn["dimensions"] == 3){
				std::cout << "... Generating 3D Block Boundary" << std::endl;
				for (int i = offset; i <= (int)((x1 - x0)/dx) - offset; i ++ ){
				for (int j = offset; j <= (int)((y1 - y0)/dx) - offset; j ++ ){
				for (int k = offset; k <= (int)((z1 - z0)/dx) - offset; k ++ ){
					Real3 posToAdd{x0 + dx * (Real)i, y0 + dx * (Real)j, z0 + dx * (Real)k};

          Real3 _perturb{((Real)dist6(rng)/100000.0-0.5) / 10.0 * dx,
                        ((Real)dist6(rng)/100000.0-0.5) / 10.0 * dx,
                        ((Real)dist6(rng)/100000.0-0.5) / 10.0 * dx};

          			posToAdd = add(posToAdd, _perturb);

					addDefaultFluidParticleAtPosition(posToAdd,getT0());
					n+=1;
				}}}

			} else if (simDataIn["dimensions"] == 2){
				std::cout << "... Generating 2D Block Boundary" << std::endl;
				for (int i = offset; i <= (int)((x1 - x0)/dx) - offset; i ++ ){
				for (int k = offset; k <= (int)((y1 - y0)/dx) - offset; k ++ ){

					Real3 posToAdd{x0 + dx * i, y0 + dx * k, 0};
					addDefaultFluidParticleAtPosition(posToAdd,getT0());
					n+=1;

				}}

			} else{
				std::cout << "Dimensions Not Specified. Check Input File." << std::endl;
				assert(0);

			}

		}

		std::cout << "Initialized Boundary Type : " << boundaryEntry["type"] << " with " << n << " Particles. " << std::endl;

	}
}

void ParticleAttributes::fluidInit(){

	std::cout << "--- Initializing Fluid Inlet " << std::endl;

	inletWidth    = (Real)parDataIn["inlet"]["width"];


	if (simDataIn["dimensions"] == 2){

		std::cout << "... Assigning the inlet properties " << std::endl;

		inletCenter   = Real3 { (Real) parDataIn["inlet"]["center"][0],
							    (Real) parDataIn["inlet"]["center"][1], 0 };

		inletTangent  = Real3 {((Real) parDataIn["inlet"]["tangent"][0]),
							    (Real) parDataIn["inlet"]["tangent"][1], 0 };				

		inletVelocity = Real3 { (Real) parDataIn["inlet"]["vel"] * (Real) parDataIn["inlet"]["normal"][0],
		 						(Real) parDataIn["inlet"]["vel"] * (Real) parDataIn["inlet"]["normal"][1], 
								0};

		inletSpeed = length(inletVelocity);

		std::cout << "... inlet speed is " << inletSpeed << " [m/s]" << std::endl; 

		std::cout << "... Computing Volume that should be used (2D) " << std::endl;

		Real kernelSum = 0.0, dx = (Real)simDataIn["dx"];		
		int neighs = 0;
		for(int i = -20; i <= 20; i++){
			for(int j = -20; j <= 20; j++){
				Real xTemp = ((Real) i) * dx;
				Real yTemp = ((Real) j) * dx;
				Real3 relPos{xTemp, yTemp, 0};
				Real dist = length(relPos);
				// std::cout << dist << std::endl;
				if ( dist < ((Real) simDataIn["smoothingLength"]) ){
					kernelSum += W_Wendland_2D_2h(dist, (Real) simDataIn["smoothingLength"]);
					neighs ++;
				}
			}
		}	

		std::cout << "... Number of neighbors for the center particle : " << neighs << std::endl;
		defaultVolume = 1.0 / kernelSum;	
		defaultMass = defaultVolume * (Real) parDataIn["rho0"];

		std::cout << "... The default volume & mass is v_o = " << defaultVolume << ", m_o = " << defaultMass << std::endl;
	} 

	std::cout << "--- Initializing Fluid Particles" << std::endl;

	if(parDataIn["geometry"]["type"] == "block"){

		std::cout << "... Defining Block of Fluid" << std::endl;


		Real x0 = (Real)parDataIn["geometry"]["begin"][0];
		Real y0 = (Real)parDataIn["geometry"]["begin"][1];
		Real z0 = (Real)parDataIn["geometry"]["begin"][2];

		Real x1 = (Real)parDataIn["geometry"]["end"][0];
		Real y1 = (Real)parDataIn["geometry"]["end"][1];
		Real z1 = (Real)parDataIn["geometry"]["end"][2];

		std::cout << "... Reading Smoothing Length, Stencil" << std::endl;
		smoothingLength = simDataIn["smoothingLength"];
		dx              = simDataIn["dx"];

		int offset = 0, jstart, jend, kstart, kend;
		if (parDataIn["geometry"]["includeBoundary"] == "True"){
			offset = 0;
		} else{
			offset = 1;
		}

		jstart = offset;
		jend = (int)((y1 - y0)/dx) - offset;

		kstart = offset;
		kend = (int)((z1 - z0)/dx) - offset;
		if (simDataIn["dimensions"] == 2){
			kstart = 0; kend = 0;
		} else if (simDataIn["dimensions"] == 1){
			jstart = 0; jend = 0;
			kstart = 0; kend = 0;
		}

    	std::mt19937 rng;
    	rng.seed(std::random_device()());
    	std::uniform_int_distribution<std::mt19937::result_type> dist6(0,1000);

		int n = 0;
		numParticles = 0;
		std::cout << "... Generating Particles" << std::endl;
		for (int i = offset; i <= (int)((x1 - x0)/dx) - offset; i ++ ){
		for (int j = jstart; j <= jend; j ++ ){
		for (int k = kstart; k <= kend; k ++ ){

			Real3 zeroVector{0,0,0};
			Real3x3 zeromat{zeroVector,zeroVector,zeroVector};
			Real3x3x3 zeroijk{zeromat,zeromat,zeromat};
			Real3 _perturb{0,0,0};

			Real3 posToAdd{x0 + dx * (Real)i,
                 	       y0 + dx * (Real)j,
                     	   z0 + dx * (Real)k};

		    posToAdd = add(posToAdd,_perturb);
			addDefaultFluidParticleAtPosition(posToAdd,getT0());
			n += 1;

		}}}
		// pos[0] = Real3{0,0,0}; temp[0] = 2000.0;
		// pos[pos.size()-1] = Real3{50.0E-6,0,0}; temp[temp.size()-1] = 0;


		std::cout << "Initialized Block of " << n << " Particles from input file.\n\n" << std::endl;
	} else if (parDataIn["geometry"]["type"] == "pointCloud"){
		std::string fileName = parDataIn["geometry"]["file"];
		readInitialPlacement(fileName);
	} else if (parDataIn["geometry"]["type"] == "JSONpointCloud"){
		std::string fileName = parDataIn["geometry"]["file"];
		json pointCloud;
		readJSON(pointCloud, fileName);
		Uint n = 0;
		numParticles = 0;		
		for(int i = 0; i < pointCloud.size(); i++){
			std::cout << pointCloud[i] << std::endl;
			Real3 posToAdd{pointCloud[i]["x"],pointCloud[i]["y"],pointCloud[i]["z"]};
			addDefaultFluidParticleAtPosition(posToAdd,getT0());			
			numParticles+=1;
			n+=1;			
		}
		std::cout << "Initialized " << n << " Particles from JSON input file.\n\n" << std::endl;		
	}

}

void ParticleAttributes::addParticlesToFluid(){	

	if ( simDataIn["dimensions"] == 2 ){		

		addDefaultFluidParticleAtPosition(inletCenter);

		for ( int i = 1; i <= (int)( inletWidth / (2.0 * dx) ) + 1; i++ ){

			Real3 posToAdd = add( inletCenter, mult(  ((Real) i) * dx, inletTangent ));

			addDefaultFluidParticleAtPosition( posToAdd );

			posToAdd = sub( inletCenter, mult(  ((Real) i) * dx, inletTangent ));
			addDefaultFluidParticleAtPosition( posToAdd );

		}
		
	}
}

void ParticleAttributes::addDefaultFluidParticleAtPosition(Real3& posToAdd){

	mass.push_back(defaultMass);
	vol.push_back(defaultVolume);
	particleDensity.push_back(0);
	particleDensityGrad.push_back(zeroVector);
	perturb.push_back(zeroVector);
	shift.push_back(zeroVector);
	pos.push_back(posToAdd);
	vel.push_back(inletVelocity);
	acc.push_back(zeroVector);
	dens.push_back(getRho0());
	densdot.push_back(0);
	densGrad.push_back(zeroVector);
	tempGrad.push_back(zeroVector);
	normalVec.push_back(zeroVector);
	temp.push_back(getT0());
	isFS.push_back(false);
	conditionNumber.push_back(0);
	curvature.push_back(0);
	L.push_back(zeromat);
	L2.push_back(zeromat);	
	velGrad.push_back(zeromat);
	tau.push_back(zeromat);
	tauDot.push_back(zeromat);
	tauGrad.push_back(zeroijk);

	enthalpy.push_back(0);
	enthalpydot.push_back(0);
	type.push_back("Fluid");
	force.push_back(zeroVector);
	isSensor.push_back(false);
	heatSensed.push_back(0);
	forceSensed.push_back(zeroVector);		
	
	numParticles += 1;

}

void ParticleAttributes::addDefaultFluidParticleAtPosition(Real3& posToAdd, Real temperature){

	mass.push_back(defaultMass);
	vol.push_back(defaultVolume);
	particleDensity.push_back(0);
	particleDensityGrad.push_back(zeroVector);
	perturb.push_back(zeroVector);
	shift.push_back(zeroVector);
	pos.push_back(posToAdd);
	vel.push_back(inletVelocity);
	acc.push_back(zeroVector);
	dens.push_back(getRho0());
	densdot.push_back(0);
	densGrad.push_back(zeroVector);
	tempGrad.push_back(zeroVector);
	normalVec.push_back(zeroVector);
	temp.push_back(temperature);
	isFS.push_back(false);
	conditionNumber.push_back(0);
	curvature.push_back(0);
	L.push_back(zeromat);
	L2.push_back(zeromat);	
	velGrad.push_back(zeromat);
	tau.push_back(zeromat);
	tauDot.push_back(zeromat);
	tauGrad.push_back(zeroijk);

	enthalpy.push_back(0);
	enthalpydot.push_back(0);
	type.push_back("Fluid");
	force.push_back(zeroVector);
	isSensor.push_back(false);
	heatSensed.push_back(0);
	forceSensed.push_back(zeroVector);		
	
	numParticles += 1;

}

void ParticleAttributes::initParticles(){

	if      ( parDataIn["type"] == "Fluid"    ) fluidInit();
	else if ( parDataIn["type"] == "Boundary" ) boundaryInit();
	else{
			std::cout << "Input File Error. Please Check input." << std::endl;
			assert(0);
	}

}
