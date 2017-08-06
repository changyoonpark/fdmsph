#include "sphsolver.hpp"

#define ZERO6X6 {{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}}

SPHSolver::SPHSolver(const json& _simData, std::map<std::string,ParticleAttributes*>& _pData)
 : simData(_simData), pData(_pData)
{

// Generate Particle Atrribution Arrays for the fluid
	nsearch = new NeighborhoodSearch(simData["smoothingLength"],true);

// load the particlesets onto the nsearcher
	std::cout << " ***** Initializing SPH Solver *****" << std::endl;

	for (const auto& pDatEntry : pData){
		ids[pDatEntry.first] = nsearch->add_point_set( pDatEntry.second->pos.front().data(), pDatEntry.second->pos.size(), true, true);
		std::cout << "... Particle set \"" << pDatEntry.first << "\" with " << pDatEntry.second->numParticles << " Particles Set Loaded onto CompactNSearch." << std::endl;
		if (pDatEntry.second->pos.size() != pDatEntry.second->numParticles)	assert(!(pDatEntry.second->pos.size() == pDatEntry.second->numParticles));
		setNames.push_back(pDatEntry.first);
	}

	totParticles = 0;
	for (const auto& pSet : nsearch->point_sets()){
		totParticles += pSet.n_points();
	}


	// Set the SPH model.
	std::cout << "... defining models for SPH" << std::endl;
	setEOS();
	setDiffusiveTerm();
	setKernels();
	setBodyForce();
	setPressureGradientFormulation();
	setViscosityFormulation();
	setViscosityConstantFormulation();
	setEOS();
	setThermalConductivityModel();
	setHeatEquationDiscretization();
	setTemperatureEnthalpyRelation();

	//Set the sensor particles.
	setSensorParticles();

	// Initialize the mass / volume
	neighborSearch();
	initializeMass();

	std::cout << "--- total number of particle sets : " << nsearch->n_point_sets() << std::endl;
	std::cout << "--- total number of particles     : " << totParticles << std::endl;

}

void SPHSolver::initializeMass(){

	std::cout << "--- Initializing particle volume / mass." << std::endl;
	const Real smoothingLength = (Real) simData["smoothingLength"];

	for (const auto& setName_i : setNames){

		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);

		// Real totmass = 0;
		for (int i = 0; i < ps_i.n_points(); ++i){
			Real kernelSum = 0.;

			for (const auto& setName_j : setNames){

				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);

				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){
					Uint const j = ps_i.neighbor(setID_j, i, _j);

					Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
					Real     dist = length(relpos);

					if ( dist > smoothingLength ) continue;

					Real3  relvel = sub(pData[setName_i]->vel[i],pData[setName_j]->vel[j]);
					Real3  reldir = div(relpos,dist);

					Real      Wij =  W_ij(dist, smoothingLength);
					kernelSum += Wij;

				}

			}
			kernelSum += W_ij(0,smoothingLength);

			pData[setName_i]->vol[i]  = (1./kernelSum);
			pData[setName_i]->mass[i] =  pData[setName_i]->dens[i] * pData[setName_i]->vol[i];

		}


	}

}

void SPHSolver::setSensorParticles(){
	for(auto& sensor : simData["sensors"]){
		if ( sensor["type"] == "sensorZPlane" ){
			Real sensorLocation = (Real)sensor["coord"];
			std::cout << "--- Creating Sensor Plane (Z = " << sensor["coord"] << ")" << std::endl;
			for (const auto& setName : setNames){
				const int setID = ids[setName];
				const auto& ps = nsearch->point_set(setID);
				if (setName == "boundary"){
					#pragma omp parallel for num_threads(NUMTHREADS)
					for (int i = 0; i < ps.n_points(); ++i){
						if ( abs(pData[setName]->pos[i][2] - sensorLocation) < simData["smoothingLength"] ){
							pData[setName]->isSensor[i] = true;
						}
					}
				}
			}

	}}
}

void SPHSolver::setThermalConductivityModel(){
	if( simData["thermalConductivity"] == "Constant"){
		std::cout << "--- Thermal Conductivity Model : Constant conductivity" << std::endl;
		thermalConductivity = constantConductivity;
	} else{
		assert(1 &&& "!!! Unimplemented Thermal Conductivity Model.");
	}
}

void SPHSolver::setHeatEquationDiscretization(){

	if (simData["heatEqDiscretization"] == "Cleary"){
		std::cout << "--- Heat Equation Discretization Model : Monaghan-Cleary Formulation" << std::endl;
		// heatTransfer = cleary;
		assert(0 && "!!! Unimplemented Heat Equation Discretization.");
	} else if (simData["heatEqDiscretization"] == "Consistent"){
		std::cout << "--- Heat Equation Discretization Model : Consistent Renormalized Laplacian Formulation" << std::endl;
		heatTransfer = consistentHeatTransfer;
	} else{
		assert(1 && "!!! Unimplemented Heat Equation Discretization.");
	}
}


void SPHSolver::setTemperatureEnthalpyRelation(){
	if (simData["temperatureEnthalpyRelation"] == "T=H"){
		std::cout << "--- Temperature-Enthalpy Relation : T=H (Debugging ONLY)" << std::endl;
		TvsH = equal;
	} else{
		assert(1 && "!!! Unimplemented Temperature / Enthalpy Relation.");
	}
}

void SPHSolver::setEOS(){
	if (simData["EOS"] == "Linear"){
		std::cout << "--- EOS : Linear (P-V-T Variant)" << std::endl;
		EOS = linearEOS;
	} else if (simData["EOS"] == "Tait"){
		std::cout << "--- EOS : Tait (P-V Variant)" << std::endl;
		EOS = taitEOS;
	} else{
		assert(1 && "!!! Unimplemented Equation of state.");
	}
}


void SPHSolver::setDiffusiveTerm(){
	if (simData["diffusionModel"] == "DeltaSPH"){
		std::cout << "--- Diffusion Model : DeltaSPH" << std::endl;
		diffusiveTerm_ij = delta_SPH;
	} else if (simData["diffusionModel"] == "None"){
		std::cout << "--- Diffusion Model : None" << std::endl;
		diffusiveTerm_ij = [](Real  rho_i,      Real rho_j,
 						  	  Real  vol_j,        Real delta,
 						  	  Real  soundSpeed, Real smoothingLength,
 						  	  Real  dist,
 						  	  Real3 relpos,    	Real3 gWij,
 						  	  Real3 densGrad_i, Real3 densGrad_j){return 0;};
	} else{
		assert(0 && "!!! Unimplemented Diffusion Model.");
	}
}

void SPHSolver::setKernels(){
	if (simData["kernelType"] == "Wendland"){
		if (simData["dimensions"] == 2){
			std::cout << "--- Kernel Type : Wendland (2D)" << std::endl;
			W_ij  =  W_Wendland_2D;
			gW_ij = gW_Wendland_2D;
		} else{
			std::cout << "--- Kernel Type : Wendland (3D)" << std::endl;
			W_ij  =  W_Wendland_3D_2h;
			gW_ij = gW_Wendland_3D_2h;
		}
	} else if (simData["kernelType"] == "Wendland_2h"){
		if (simData["dimensions"] == 2){
			std::cout << "--- Kernel Type : Wendland (2D, 2h version)" << std::endl;
			W_ij  =  W_Wendland_2D_2h;
			gW_ij = gW_Wendland_2D_2h;
		} else{
			std::cout << "--- Kernel Type : Wendland (3D, 2h version)" << std::endl;
			W_ij  =  W_Wendland_3D_2h;
			gW_ij = gW_Wendland_3D_2h;
		}
	} else{
		assert(0 && "!!! Unimplemented Kernel Type.");
	}
}

void SPHSolver::setBodyForce(){
	if (simData["bodyForce"] == "Gravity"){
		std::cout << "--- Body Force Type : Gravity" << std::endl;
		bodyForceAcc_i = gravity_acc;
	} else if (simData["bodyForce"] == "None"){
		std::cout << "--- Body Force Type : None" << std::endl;
		bodyForceAcc_i = [](){return Real3{0,0,0};};
	} else{
		assert(1 && "!!! Unimplemented BodyForce Type.");
	}
}

void SPHSolver::setViscosityConstantFormulation(){
	if (simData["viscosityConstantFormulation"] == "temperatureDependent"){
	} else if (simData["viscosityConstantFormulation"] == "Fixed"){
		std::cout << "--- Viscosity Constant Formulation : Fixed Viscosity" << std::endl;
		viscosityConstant = fixedViscosity;
	}
}

void SPHSolver::setViscosityFormulation(){
	if (simData["viscosityFormulation"] == "Shao"){
		std::cout << "--- Viscosity Formulation : Shao" << std::endl;
		viscosityAcc_ij = viscosity_acc_ij_Shao;
	} else{
		assert(1 && "!!! Unimplemented viscosity formualtion.");
	}
}

void SPHSolver::setPressureGradientFormulation(){

	if ((int)simData["pressureFormulation"] == 1){
		std::cout << "--- Pressure Gradient Formulation : Type 1" << std::endl;
		pressureAcc_ij = pressure_acc_ij_1;
	} else if ((int)simData["pressureFormulation"] == 2){
		std::cout << "--- Pressure Gradient Formulation : Type 2" << std::endl;
		pressureAcc_ij = pressure_acc_ij_2;
	} else if ((int)simData["pressureFormulation"] == 3){
		std::cout << "--- Pressure Gradient Formulation : Type 3" << std::endl;
		pressureAcc_ij = pressure_acc_ij_3;
	} else{
		assert(1 && "!!! Unimplemented pressure gradient formualtion.");
	}

}

void SPHSolver::neighborSearch(){

	std::cout << "... Searching Neighbors" << std::endl;
	auto t0 = std::chrono::high_resolution_clock::now();
	nsearch->find_neighbors();
	std::cout << "--- Neighborhood search took " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count() << "ms" << std::endl;


	// I have no idea what this does. Documentation is unclear.
	// nsearch->z_sort();
	// for (const auto& pDatEntry : pData){
	// 	auto& ps = nsearch->point_set(ids[pDatEntry.first]);
	// 	ps.sort_field(pDatEntry.second-> pos.data());
	// 	ps.sort_field(pDatEntry.second-> vel.data());
	// 	ps.sort_field(pDatEntry.second-> acc.data());
	// 	ps.sort_field(pDatEntry.second->dens.data());
	// 	ps.sort_field(pDatEntry.second->temp.data());
	// 	ps.sort_field(pDatEntry.second->type.data());
	// }

}



void SPHSolver::computeInteractions(){
	using namespace RealOps;
	// The fluid.
	const Real dx = (Real) simData["dx"];
	const Uint dims = simData["dimensions"];
	const Real smoothingLength = (Real) simData["smoothingLength"];
	VectorXd neg_delta_mn(6);
	neg_delta_mn(0) = -1.0;	neg_delta_mn(1) = -1.0;	neg_delta_mn(2) = -1.0;
	neg_delta_mn(3) = 0.0; neg_delta_mn(4) = 0.0; neg_delta_mn(5) = 0.0;

	const Real3 zerovec{0,0,0};
	const Real3x3 zeromat{zerovec,zerovec,zerovec};

	std::cout << "		|------ ... Performing Fluid - (Fluid / Boundary) Interactions"<< std::endl;
	for (const auto& setName_i : setNames){

		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);

		#pragma omp parallel for num_threads(NUMTHREADS)
		for (int i = 0; i < ps_i.n_points(); ++i){

			Real    m_i  = pData[setName_i]->mass[i];
			Real  rho_i  = pData[setName_i]->dens[i];

			// Renormalization Tensor (First Derivative)
			Real3x3& L_i = pData[setName_i]->L[i];
			L_i = zeromat;
			// Renormalization Tensor (Second Derivative)
			SymTensor3 A_i_kmn;
			MatrixXd G(6,6); G = MatrixXd::Zero(6,6);
			Real3x3& L2_i = pData[setName_i]->L2[i];
			L2_i = zeromat;


			// Normalized Temperature Gradient
			pData[setName_i]->tempGrad[i] = zerovec;
			// Renormalized Density Gradient.
			pData[setName_i]->densGrad[i] = zerovec;
      // Renormalized Color Gradient.
      pData[setName_i]->normalVec[i] = zerovec;

			// Clear the time derivatives.ff
			pData[setName_i]->acc[i] = Real3{0,0,0};
			pData[setName_i]->densdot[i] = 0;
			pData[setName_i]->enthalpydot[i] = 0;

			// Initialize the temperature wrt to the enthalpy.
			// pData[setName_i]->temp[i] = TvsH(pData[setName_i]->enthalpy[i]);


			// FIRST PHASE.

			// Compute the Renormalization Matrix / density gradient using the renormalization matrix.
			// L_i is the renormalization matrix for the first derivative.
			// A_i_knm is the third order tensor used for calculating the laplacian normalizer.
			for (const auto& setName_j : setNames){

				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);

				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){

					// Index of neighbor particle in setID_j
					Uint const j = ps_i.neighbor(setID_j, i, _j);

					// Define all the ingedients for the particle interaction
					Real    rho_j = pData[setName_j]->dens[j];
					Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
					Real     dist = length(relpos);

					if ( dist > smoothingLength ) continue;

					Real3  reldir = div(relpos,dist);
					Real    vol_j = pData[setName_j]->vol[j];
          Real      Wij = W_ij(dist, smoothingLength);
					Real3    gWij = gW_ij(dist, reldir, smoothingLength);

					// Compute First derivative Renormalization Tensor
					L_i = add(L_i,mult(vol_j,tensorProduct(gWij,relpos)));
					// Compute Renormalized Density Gradient
					pData[setName_i]->densGrad[i]  = add(pData[setName_i]->densGrad[i], mult((rho_i - rho_j) * vol_j,gWij));
					pData[setName_i]->tempGrad[i]  = add(pData[setName_i]->tempGrad[i], mult((pData[setName_i]->temp[i] - pData[setName_j]->temp[j]) * vol_j,gWij));
          pData[setName_i]->normalVec[i] = add(pData[setName_i]->normalVec[i],mult(vol_j,gWij));
				}
			}



			// Compute the inverse of the renormalization tensor.
			if (dims == 2){
				L_i[1][0] = 0; L_i[1][1] = 1.0; L_i[1][2] = 0; L_i[2][1] = 0;
			}
      checkSingularity(L_i);
			L_i = inv(L_i);
      MatrixXd _L_i(3,3); assign(_L_i,L_i);
      VectorXcd eigs_im = _L_i.eigenvalues(); VectorXd eigs = - eigs_im.real();
      Real lambda = maxEig(eigs);


			// This is <\nabla rho_i>, used for delta-SPH
			pData[setName_i]->densGrad[i]  = mult(L_i,pData[setName_i]->densGrad[i]);
			// This is <\nabla T_i>, used for normalized heat transfer
			pData[setName_i]->tempGrad[i]  = mult(L_i,pData[setName_i]->tempGrad[i]);
      // This is the normal vector.

      if (lambda > 1.5){
        pData[setName_i]->isFS[i] = 1.0;
        pData[setName_i]->normalVec[i] = mult(L_i,pData[setName_i]->normalVec[i]);
        pData[setName_i]->normalVec[i] = dir(pData[setName_i]->normalVec[i]);
      } else{
        pData[setName_i]->normalVec[i] = zerovec;
      }

			// SECOND PHASE.
			// Compute the Renormalization Tensor for the Laplacian.
			// Referenced from "Error Estimation in Smoothed Particle Hydrodynamics and a New Scheme for Second Derivatives", Fatehi et al.
			for (const auto& setName_j : setNames){

				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);

				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){

					// Index of neighbor particle in setID_j
					Uint const j = ps_i.neighbor(setID_j, i, _j);

					// Define all the ingedients for the particle interaction
					Real    rho_j = pData[setName_j]->dens[j];
					Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
					Real     dist = length(relpos);
					Real3x3&  L_i = pData[setName_i]->L[i];
					if ( dist > smoothingLength ) continue;

					Real3  reldir = div(relpos,dist);
					Real    vol_j = pData[setName_j]->vol[j];
					Real3    gWij = gW_ij(dist, reldir, smoothingLength);

					// Compute Second derivative Renormalization Tensor
					for(int k = 0; k < 3; k ++) for(int m = 0; m < 3; m ++) for(int n = 0; n < 3; n ++){
						for(int q = 0; q < 3; q++){
							A_i_kmn.add(k, m, n, L_i[k][q] * relpos[m] * relpos[n] * gWij[q] * vol_j);
						}
					}

				}
			}


			// THIRD PHASE.
			// Compute the matrix for inversion for the second derivative renormalization tensor.
			// Referenced from "Error Estimation in Smoothed Particle Hydrodynamics and a New Scheme for Second Derivatives", Fatehi et al.
			for (const auto& setName_j : setNames){

				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);

				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){

					// Index of neighbor particle in setID_j
					Uint const j = ps_i.neighbor(setID_j, i, _j);

					// Define all the ingedients for the particle interaction
					Real    rho_j = pData[setName_j]->dens[j];
					Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
					Real     dist = length(relpos);
					Real3x3&  L_i = pData[setName_i]->L[i];
					if ( dist > smoothingLength ) continue;

					Real3  reldir = div(relpos,dist);
					Real    vol_j = pData[setName_j]->vol[j];
					Real3    gWij = gW_ij(dist, reldir, smoothingLength);

					// Compute Second derivative Renormalization Tensor
					for(int I = 0; I < 6; I ++){
						Uint m = IDXPAIR[I][0], n = IDXPAIR[I][1];
						for(int J = 0; J < 6; J ++){
							Uint o = IDXPAIR[J][0], p = IDXPAIR[J][1];
							Real A_kI_ek = 0;
							for(int k = 0 ;k < 3; k++) A_kI_ek += A_i_kmn(k,m,n) * reldir[k];
							if (o == p){
								G(I,J) +=  (A_kI_ek + relpos[m] * reldir[n]) * (reldir[o] * gWij[p]) * vol_j;
							} else{
								G(I,J) +=  (A_kI_ek + relpos[m] * reldir[n]) * (reldir[o] * gWij[p] + reldir[p] * gWij[o]) * vol_j;
							}

						}
					}

				}
			}
			// std::cout << G << std::endl;
			// std::cout << "---------\n\n" << std::endl;
			L2_i = toReal3x3(G.colPivHouseholderQr().solve(neg_delta_mn));

		}
	}



	for (const auto& setName_i : setNames){

		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);

		#pragma omp parallel for num_threads(NUMTHREADS)
		for (int i = 0; i < ps_i.n_points(); ++i){
			const Real&   m_i = pData[setName_i]->mass[i];
			const Real& rho_i = pData[setName_i]->dens[i];
			const Real&   T_i = pData[setName_i]->temp[i];
      const Real& vol_i = pData[setName_i]->vol[i];
			const Real&   k_i = thermalConductivity(pData[setName_i]->getType(),
				  							 T_i);

			const Real&   P_i = EOS(T_i,   pData[setName_i]->getT0(),
							 rho_i, pData[setName_i]->getRho0(),
							 	    pData[setName_i]->getSoundSpeed(),
							 	    pData[setName_i]->getThermalExpansion());
			const Real3& tempGrad_i = pData[setName_i]->tempGrad[i];
			const Real3& densGrad_i = pData[setName_i]->densGrad[i];
      const Real3& normal_i =  pData[setName_i]->normalVec[i];
			const Real3x3& L2_i = pData[setName_i]->L2[i];


      Real& curvature = pData[setName_i]->curvature[i];
      curvature = 0;
      Real curvatureCorrection = (vol_i) * W_ij(0, smoothingLength);


			bool isSensor_i   = (pData[setName_i]->isSensor[i]);
      bool isFS_i = (pData[setName_i]->isFS[i]);
			bool isBoundary_i = (setName_i == "boundary") ? true : false;

			for (const auto& setName_j : setNames){
				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);
				const bool isBoundary_j = (setName_j == "boundary") ? true : false;

				// Loop over the neighbors of particl i in the indexed particle set.
				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){

					// Index of neighbor particle in setID_j
					Uint const j = ps_i.neighbor(setID_j, i, _j);
					// Define all the ingedients for the particle interaction
					Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
					Real3 relposj = mult(-1.0,relpos);

					Real     dist = length(relpos);
					if ( dist > smoothingLength ) continue;

					const Real&      m_j = pData[setName_j]->mass[j];
					const Real&    vol_j = pData[setName_j]->vol[j];
          const Real3& normal_j =  pData[setName_j]->normalVec[j];
          Real   rho_i_temp, rho_j_temp;

					Real3  relvel = sub(pData[setName_i]->vel[i],pData[setName_j]->vel[j]);
					Real3  reldir = div(relpos,dist);

					Real      Wij =  W_ij(dist,         smoothingLength);
					Real3    gWij = gW_ij(dist, reldir, smoothingLength);
					Real    rho_j = pData[setName_j]->dens[j];
					Real3 densGrad_j = pData[setName_j]->densGrad[j];
					Real   	  T_j = pData[setName_j]->temp[j];
					Real     mu_j = viscosityConstant(pData[setName_j]->getMu(),T_j);

					Real   	  P_j = EOS(   T_j, pData[setName_j]->getT0(),
										 rho_j, pData[setName_j]->getRho0(),
									 	    	pData[setName_j]->getSoundSpeed(),
									 	    	pData[setName_j]->getThermalExpansion());

					Real   	  k_j = thermalConductivity(pData[setName_j]->getType(),
														T_j);

          bool isFS_j = (pData[setName_j]->isFS[j]);

					// Continuity : update the density for the fluid particles.
					// Continuity : update the density for the boundary particles.

					// Delta SPH Diffusion
					pData[setName_i]->densdot[i] += dot(mult(rho_i, relvel), gWij) * vol_j;
					pData[setName_i]->densdot[i] += diffusiveTerm_ij(rho_i,rho_j,
																		 vol_j, (Real) simData["delta"],
																		 pData[setName_j]->getSoundSpeed(), smoothingLength,
																		 dist,
																		 relpos,gWij,
																		 densGrad_i,densGrad_j);

					// Heat Transfer between particles. Note that the boundary densities must be set to the
					// Actual density of the boundary material, to account for the correct thermal diffusivity.
					rho_i_temp = isBoundary_i ? pData[setName_i]->getMaterialDensity() : rho_i ;
					rho_j_temp = isBoundary_j ? pData[setName_j]->getMaterialDensity() : rho_j ;

					Real heat = heatTransfer(L2_i,tempGrad_i,
											 T_i, T_j, m_j, rho_i_temp, rho_j_temp, k_i, k_j,
											 relpos, reldir,
											 dist, gWij, vol_j);


					pData[setName_i]->enthalpydot[i] += heat;


					// If the particle is not a boundary particle, accumulate the momentum contribution.
					if (!isBoundary_i){
					// Acceleration Due to pressure gradient.
						Real3 a = add(pressureAcc_ij(m_j,gWij,P_i,P_j,rho_i,rho_j,vol_j), viscosityAcc_ij(m_j, relpos, gWij, rho_i, rho_j,dist, relvel, mu_j));
						pData[setName_i]->acc[i]      = add(pData[setName_i]->acc[i], a);
            // If particle i and particle j are both free-surface particles, add the surface tension contribution.
            if(isFS_i && isFS_j){
              curvature += (vol_j) * dot(sub(normal_j, normal_i), gWij);
              curvatureCorrection += (vol_j) * Wij;
            }

					}

					// // If the boundary particle is a sensor, accumulate the data.
					if(isSensor_i && !isBoundary_j){
						Real3 f = add(pressureAcc_ij(m_j,gWij,P_i,P_j,rho_i,rho_j,vol_j), viscosityAcc_ij(m_j, relpos, gWij, rho_i, rho_j,dist, relvel, mu_j));
						f = mult(m_i,f);
						pData[setName_i]->heatSensed[i] += heat;
						pData[setName_i]->forceSensed[i] = add(pData[setName_i]->forceSensed[i], f);
					}


				}
			}

      curvature = curvature / curvatureCorrection;
      // Acceleration Due to Surface Tension. Apply only to Free-surface particles.
      if (isFS_i && isBoundary_i){
        pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i],
          mult(- (pData[setName_i]->getSurfaceTensionCoeff() / m_i) * curvature,normal_i)
        );
      }
			// Acceleration Due to Body Force. Apply only to non-boundary particles.
			if (!isBoundary_i){
				pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i],bodyForceAcc_i());
			} else{
				pData[setName_i]->acc[i] = zerovec;
				pData[setName_i]->vel[i] = zerovec;
			}

		}


	}


}



void SPHSolver::marchTime(Uint t){
	using namespace RealOps;

	std::cout << "-------------------- Timestep #" << t << std::endl;
	const Real dt = (Real) simData["dt"];
	const Real dx = (Real) simData["dx"];
	for (const auto& setName : setNames){

		const int setID = ids[setName];
		const auto& ps = nsearch->point_set(setID);

		if( setName == "fluid"){
		// Fluid Particles    : march the position, velocity, density.
			#pragma omp parallel for num_threads(NUMTHREADS)
			for (int i = 0; i < ps.n_points(); ++i){
				pData[setName]->vel[i]  = add(pData[setName]->vel[i], mult(dt * 0.5,pData[setName]->acc[i]));
				pData[setName]->pos[i]  = add(pData[setName]->pos[i], mult(dt,pData[setName]->vel[i]));

				pData[setName]->dens[i] = pData[setName]->dens[i] + dt * 0.5 * pData[setName]->densdot[i];
				pData[setName]->temp[i] = pData[setName]->temp[i] + dt * 0.5 * pData[setName]->enthalpydot[i] / pData[setName]->getSpecificHeat();
			}
		} else if (setName == "boundary"){
		// Boundary Particles : march the density only.
			#pragma omp parallel for num_threads(NUMTHREADS)
			for (int i = 0; i < ps.n_points(); ++i){
				pData[setName]->dens[i] = pData[setName]->dens[i]  + dt * 0.5 * pData[setName]->densdot[i];
				pData[setName]->temp[i] = pData[setName]->temp[i] + dt * 0.5 * pData[setName]->enthalpydot[i] / pData[setName]->getSpecificHeat();
			}
		}

	}

//  Forward Euler. (For Debugging)
	// for (const auto& setName : setNames){
	// 	const int setID = ids[setName];
	// 	const auto& ps = nsearch->point_set(setID);
	// 	if( setName == "fluid"){
	// 		#pragma omp parallel for num_threads(NUMTHREADS)
	// 		for (int i = 0; i < ps.n_points(); ++i){
	// 			pData[setName]->temp[i] = pData[setName]->temp[i] + dt * pData[setName]->enthalpydot[i];
	// 		}
	// 	}
	// }

	std::cout << "	|--- Calculating Interatctions" << std::endl;
	computeInteractions();
	std::cout << "	|--- Updating Variables" << std::endl;


	for (const auto& setName : setNames){

		const int setID = ids[setName];
		const auto& ps = nsearch->point_set(setID);

		if( setName == "fluid"){
		// Fluid Particles    : march the position, velocity, density.
			#pragma omp parallel for num_threads(NUMTHREADS)
			for (int i = 0; i < ps.n_points(); ++i){
				pData[setName]->vel[i]  = add(pData[setName]->vel[i], mult(0.5 * dt,pData[setName]->acc[i]));

				pData[setName]->dens[i] = pData[setName]->dens[i] + dt * 0.5 * pData[setName]->densdot[i];
				pData[setName]->temp[i] = pData[setName]->temp[i] + dt * 0.5 * pData[setName]->enthalpydot[i] / pData[setName]->getSpecificHeat();
			}
		} else if (setName == "boundary"){
		// Boundary Particles : march the density only.
			#pragma omp parallel for num_threads(NUMTHREADS)
			for (int i = 0; i < ps.n_points(); ++i){
				pData[setName]->dens[i] = pData[setName]->dens[i] + dt * 0.5 * pData[setName]->densdot[i];
				pData[setName]->temp[i] = pData[setName]->temp[i] + dt * 0.5 * pData[setName]->enthalpydot[i] / pData[setName]->getSpecificHeat();

			}
		}

	}





}
