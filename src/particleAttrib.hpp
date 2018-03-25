#ifndef __PARTICLEATTRIB__
#define __PARTICLEATTRIB__

#include <iostream>
#include <fstream>
#include <random>
#include <omp.h>
#include "helpers.hpp"
#include "csv.hpp"

using namespace CompactNSearch;
using namespace RealOps;

using Uint = unsigned int;

class ParticleAttributes{

public:
	unsigned int numParticles;
	unsigned int setID;
	Real dx, smoothingLength;
	Real oneOverLambda, eta;

	Real defaultVolume, defaultMass;
	Real inletWidth;
	Real inletSpeed;
	Real3 inletVelocity;
	Real3 inletTangent;
	Real3 inletCenter;

	std::vector<Real>  vol;
	std::vector<Real>  particleDensity;
	std::vector<Real3> particleDensityGrad;
	std::vector<Real> mass;
	std::vector<Real3> pos;
	std::vector<Real3> perturb;
	std::vector<Real3> shift;
	std::vector<Real3> vel;
	std::vector<Real3> acc;
	std::vector<Real3> force;
	std::vector<Real>  dens;
	std::vector<Real3> densGrad;
	std::vector<Real3> normalVec;
	std::vector<Real>  curvature;
	std::vector<Real3> tempGrad;
	std::vector<Real3x3> L;
	std::vector<Real> conditionNumber;
	std::vector<Real3x3> L2;

	std::vector<Real3x3> tau;
	std::vector<Real3x3> tauDot;
	std::vector<Real3x3x3> tauGrad;
	std::vector<Real3x3> velGrad;


	std::vector<Real>  temp;
	std::vector<Real>  enthalpy;
	std::vector<Real>  enthalpydot;
	std::vector<Real>  densdot;
	std::vector<bool>  isSensor;
	std::vector<bool>  isFS;
	std::vector<Real>  heatSensed;
	std::vector<Real3> forceSensed;
	std::vector<std::string> type;

	ParticleAttributes(const json& inputData, const json& particleData);

	void addParticlesToFluid();

//  This should only work for fluid particles
	inline Real getT0(){ 		  		 return parDataIn["T0"];}
	inline Real3 getv0(){                return Real3{parDataIn["v0"][0],parDataIn["v0"][1],parDataIn["v0"][2]};}
	inline Real getSpecificHeat(){ 		 return parDataIn["specificHeat"];}
	inline Real getRho0(){        		 return parDataIn["rho0"];}
	inline Real getMaterialDensity(){    return parDataIn["rho_material"];}
	inline Real getMu(){				 return parDataIn["mu"];}
	inline Real getSurfaceTensionCoeff(){ return parDataIn["surfaceTensionCoeff"];}
	inline Real getSoundSpeed() {		 return parDataIn["soundSpeed"];}
	inline Real getOneOverLambda(){      return 1./(Real)parDataIn["lambda"];}
	inline Real getEta(){                return (Real)parDataIn["eta"];}
	inline Real getThermalExpansion() {  return parDataIn["thermalExpansion"];}
	inline Real getViscosity(){          return parDataIn["viscosity"];}
	inline std::string getType(){        return parDataIn["type"];}
private:

	std::vector<Real> x,y,z,vx,vy,vz;
	const json& simDataIn, parDataIn;

	const Real3     zeroVector{0,0,0};
	const Real3x3   zeromat{Real3{0,0,0},Real3{0,0,0},Real3{0,0,0}};
	const Real3x3x3 zeroijk{ Real3x3{Real3{0,0,0},Real3{0,0,0},Real3{0,0,0}},
							 Real3x3{Real3{0,0,0},Real3{0,0,0},Real3{0,0,0}},
							 Real3x3{Real3{0,0,0},Real3{0,0,0},Real3{0,0,0}}} ;

	void initParticles();
	void fluidInit();
	void boundaryInit();
	void addDefaultFluidParticleAtPosition(Real3& posToAdd);
	void addDefaultFluidParticleAtPosition(Real3& posToAdd, Real temp);
	void readInitialPlacement(std::string fileName);
};

#define __PARTICLEATTRIBDEF_

#endif
