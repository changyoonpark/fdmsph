#include "sphsolver.hpp"


SPHSolver::SPHSolver(const json& _simData, std::map<std::string,ParticleAttributes*>& _pData)
 : simData(_simData), pData(_pData)
{

// Generate Particle Atrribution Arrays for the fluid
	nsearch = new NeighborhoodSearch((Real)simData["smoothingLength"] + EPSL_SMALL,true);
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
					Real3  reldir = divide(relpos,dist);

					Real      Wij =  W_ij(dist, smoothingLength);
					kernelSum += Wij;
				}

			}
			kernelSum += W_ij(0.0,smoothingLength);
			pData[setName_i]->normalVec[i][0] = (Real)ps_i.n_neighbors(ids["fluid"],i);			
			pData[setName_i]->vol[i]  = (1.0/kernelSum);
			// pData[setName_i]->vol[i]  =  9.391435E-17;
			pData[setName_i]->mass[i] =  pData[setName_i]->dens[i] * pData[setName_i]->vol[i];


      // For Debugging.
      // pData[setName_i]->temp[i] = pData[setName_i]->pos[i][2] * pData[setName_i]->pos[i][2] * 0.5;

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
		} else if (simData["dimensions"] == 3){
			std::cout << "--- Kernel Type : Wendland (3D)" << std::endl;
			W_ij  =  W_Wendland_3D_2h;
			gW_ij = gW_Wendland_3D_2h;
		}
	} else if (simData["kernelType"] == "Wendland_2h"){
		if (simData["dimensions"] == 2){
			std::cout << "--- Kernel Type : Wendland (2D, 2h version)" << std::endl;
			W_ij  =  W_Wendland_2D_2h;
			gW_ij = gW_Wendland_2D_2h;
			gW_ij_nd = gW_Wendland_2D_2h_Nondim; 			
		} else if (simData["dimensions"] == 3){
			std::cout << "--- Kernel Type : Wendland (3D, 2h version)" << std::endl;
			W_ij  =  W_Wendland_3D_2h;
			gW_ij = gW_Wendland_3D_2h;
			gW_ij_nd = gW_Wendland_3D_2h_Nondim; 			
		} else{
			std::cout << "--- Kernel Type : Wendland (1D, 2h version)" << std::endl;
			W_ij  =  W_Wendland_1D_2h;
			gW_ij = gW_Wendland_1D_2h;
			gW_ij_nd = gW_Wendland_1D_2h_Nondim; 						
		}
	} else if (simData["kernelType"] == "Quintic_Spline"){
		if (simData["dimensions"] == 3){
			std:: cout << "--- Kernel Type : Quintic Spline (3D, h version)" << std::endl;
			W_ij  = W_QuinticSpline_3D;
			gW_ij = gW_QuinticSpline_3D;
		}
	} 
	else{
		assert(0 && "!!! Unimplemented Kernel Type.");
	}
}

void SPHSolver::setBodyForce(){
	if (simData["bodyForce"] == "Gravity"){
		std::cout << "--- Body Force Type : Gravity" << std::endl;
		bodyForceAcc_i = gravity_acc;
	} else if (simData["bodyForce"] == "None"){
		std::cout << "--- Body Force Type : None" << std::endl;
		bodyForceAcc_i = [](){return Real3{0.0,0.0,0.0};};
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
	neg_delta_mn(3) =  0.0; neg_delta_mn(4) =  0.0; neg_delta_mn(5) =  0.0;

	Real3 neg_delta_mn_2D = Real3{-1.0,  // (0,0)
		-1.0,  // (2,2)
		   0}; // (0,2)


	const Real3 zerovec{0.0,0.0,0.0};
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
			MatrixXd _L_i(3,3); _L_i = MatrixXd::Zero(3,3);
			Real3x3& L_i = pData[setName_i]->L[i];
			L_i = zeromat;
			
			// Renormalization Tensor (Second Derivative)
			MatrixXd _Q_i(6,6); _Q_i = MatrixXd::Zero(6,6);
			Real3x3& L2_i = pData[setName_i]->L2[i];
			L2_i = zeromat;

			// Renormalization Tensor (Second Derivative)
			SymTensor3 A_i_kmn;
			MatrixXd G(6,6); G = MatrixXd::Zero(6,6);
			// Real3x3 G = zeromat;

			// Normalized Temperature Gradient
			pData[setName_i]->tempGrad[i] = zerovec;
			// Renormalized Density Gradient.
			pData[setName_i]->densGrad[i] = zerovec;

			// Clear the time derivatives.ff
			pData[setName_i]->acc[i] = Real3{0.0,0.0,0.0};
			pData[setName_i]->densdot[i] = 0.0;
			pData[setName_i]->enthalpydot[i] = 0.0;

			// FIRST PHASE.

			// Compute the Renormalization Matrix / density gradient using the renormalization matrix.
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
					Real3  reldir = divide(relpos,dist);
					Real    vol_j = pData[setName_j]->vol[j];
          			Real      Wij = W_ij(dist, smoothingLength);
					Real3    gWij = gW_ij(dist, reldir, smoothingLength);
					

          			Real3x3 v_j_gWij_rij = mult(- vol_j,tensorProduct(gWij,relpos));
					// Compute First derivative Renormalization Tensor
					L_i = add(L_i, v_j_gWij_rij);
					// Compute Renormalized Density Gradient
					pData[setName_i]->densGrad[i]  = add(pData[setName_i]->densGrad[i], mult((rho_j - rho_i) * vol_j,gWij));
					pData[setName_i]->tempGrad[i]  = add(pData[setName_i]->tempGrad[i], mult( (pData[setName_j]->temp[j] - pData[setName_i]->temp[i]) * vol_j, gWij));

				}
			}


			// Compute the inverse of the renormalization tensor.
			setDims(L_i,dims);

			toMatrix3d(L_i,_L_i);
			JacobiSVD<MatrixXd> svd_L_i(_L_i);
			Real cond_first_i = svd_L_i.singularValues()(0) / svd_L_i.singularValues()(svd_L_i.singularValues().size()-1);		
			pData[setName_i]->isFS[i] = cond_first_i;
			_L_i = _L_i.inverse();
			toReal3x3(_L_i,L_i);
			
			// if(cond_first_i < 100.0){
			pData[setName_i]->densGrad[i]  = mult(L_i,pData[setName_i]->densGrad[i]);
			pData[setName_i]->tempGrad[i]  = mult(L_i,pData[setName_i]->tempGrad[i]);
			// }




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

					Real3  reldir = divide(relpos,dist);
					Real    vol_j = pData[setName_j]->vol[j];
					Real3    gWij = gW_ij(dist, reldir, smoothingLength);

					if ( dist < EPSL_SMALL2){
						reldir = Real3{0,0,0};
						gWij = Real3{0,0,0};
					}
					// Compute Second derivative Renormalization Tensor
					for(int k = 0; k < 3; k ++) for(int m = 0; m < 3; m ++) for(int n = 0; n < 3; n ++){
						for(int q = 0; q < 3; q++){
							A_i_kmn.add(k, m, n, L_i[k][q] * relpos[m] * relpos[n] * gWij[q] * vol_j);
						}
					}

				}

			}































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

					Real3  reldir = divide(relpos,dist);
					Real    vol_j = pData[setName_j]->vol[j];
					Real3    gWij = gW_ij(dist, reldir, smoothingLength);

					for ( Uint I = 0; I < 6; I++){
						
						Uint m = IDXPAIR[I][0], n = IDXPAIR[I][1];
						Real Nij_mn = (A_i_kmn(0,m,n) * reldir[0] + A_i_kmn(1,m,n) * reldir[1] + A_i_kmn(2,m,n) * reldir[2]) + relpos[m] * reldir[n];
						// truncate1(Nij_mn);
						for ( Uint J = 0; J < 6; J++){
							
							Uint o = IDXPAIR[J][0], p = IDXPAIR[J][1];
							Real Mij_op = reldir[o] * gWij[p];
							Real Mij_po = reldir[p] * gWij[o];

							if (o == p){
								G(I,J) = G(I,J) + Nij_mn * Mij_op * vol_j;
							} else{
								G(I,J) = G(I,J) + Nij_mn * (Mij_op + Mij_po) * vol_j;
							}	

						}
					}
					
				}
			}

			// const Uint IDXPAIR[6][2] = {{0,0},{1,1},{2,2},{0,1},{1,2},{0,2}};
			if(dims == 2){
				G(2,2) = 1.0; G(4,4) = 1.0; G(5,5) = 1.0;
			} else if(dims == 1){
				G(1,1) = 1.0; G(2,2) = 1.0;
				G(3,3) = 1.0; G(4,4) = 1.0;
				G(5,5) = 1.0;
			}

			VectorXd L2_i_vec = G.fullPivLu().solve(neg_delta_mn);
			L2_i = toReal3x3From6(L2_i_vec,dims);
			setDims2(L2_i,dims);
											





















































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

			const Real&   P_i = EOS(T_i, pData[setName_i]->getT0(),
									     rho_i, pData[setName_i]->getRho0(),
									     pData[setName_i]->getSoundSpeed(),
										 pData[setName_i]->getThermalExpansion());
											
			const Real3& tempGrad_i = pData[setName_i]->tempGrad[i];
			const Real3& densGrad_i = pData[setName_i]->densGrad[i];
	        const Real3& normal_i =  pData[setName_i]->normalVec[i];
			const Real3x3& L2_i = pData[setName_i]->L2[i];
			const Real& cond_i = pData[setName_i]->isFS[i];

			bool isSensor_i   = (pData[setName_i]->isSensor[i]);
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
					if ( dist > smoothingLength) continue;
					
					const Real&      m_j = pData[setName_j]->mass[j];
					const Real&    vol_j = pData[setName_j]->vol[j];
					const Real3& normal_j =  pData[setName_j]->normalVec[j];
					const Real3& tempGrad_j = pData[setName_j]->tempGrad[j];

					Real   rho_i_temp, rho_j_temp;

					Real3  relvel = sub(pData[setName_i]->vel[i],pData[setName_j]->vel[j]);
					Real3  reldir = divide(relpos,dist);

					Real      Wij =  W_ij(dist,         smoothingLength);
					Real3    gWij = gW_ij(dist, reldir, smoothingLength);
					Real3    gWij_nd = gW_ij_nd(dist, reldir, smoothingLength);
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
					

					// If the laplacian corrector cannot be defined, use the conventional operator.
					Real heat;
					// if (cond_i < 100.0){
						heat = consistentHeatTransfer(L2_i,tempGrad_i,
													  T_i, T_j, m_j, rho_i_temp, rho_j_temp, k_i, k_j,
													  relpos, reldir,
													  dist, gWij, vol_j);
					// } else{
						// heat = inconsistentHeatTransfer(tempGrad_i, tempGrad_j,
						// 				    			T_i, T_j, m_j, rho_i_temp, rho_j_temp, k_i, k_j,
						// 				    			relpos, reldir,
						// 				 				dist, gWij, vol_j);		
					// }
					pData[setName_i]->enthalpydot[i] += heat;


					// If the particle is not a boundary particle, accumulate the momentum contribution.
					if (!isBoundary_i){
					  // Acceleration Due to pressure gradient.
						Real3 a = add(pressureAcc_ij(m_j,gWij,P_i,P_j,rho_i,rho_j,vol_j), viscosityAcc_ij(m_j, relpos, gWij, rho_i, rho_j,dist, relvel, mu_j));
            			// Interparticle forces, (IFF surface tension model).
            			a = add(a, iif_acc(simData["IIFCoeff"], reldir, dist, smoothingLength, m_i, m_j));

						pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i], a);

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


			// if( dims != 1 ){
				// pData[setName_i]->enthalpydot[i] = pData[setName_i]->enthalpydot[i] * ((7.0/(4.0*M_PI))/(smoothingLength*smoothingLength*smoothingLength));
			// } else {
				// pData[setName_i]->enthalpydot[i] = pData[setName_i]->enthalpydot[i] * (-0.625 / (smoothingLength*smoothingLength)); 
			// }

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
				// pData[setName]->vel[i]  = add(pData[setName]->vel[i], mult(dt * 0.5,pData[setName]->acc[i]));
				// pData[setName]->pos[i]  = add(pData[setName]->pos[i], mult(dt,pData[setName]->vel[i]));

				// pData[setName]->dens[i] = pData[setName]->dens[i] + dt * 0.5 * pData[setName]->densdot[i];

				if(pData[setName]->pos[i][0] > 1.0E-8 && pData[setName]->pos[i][0] < 50.0E-6 - 1.0E-8){
				// if(pData[setName]->pos[i][0] > 1.0E-8){
					pData[setName]->temp[i] = pData[setName]->temp[i] + dt * 0.5 * pData[setName]->enthalpydot[i] * 15.0 / (7900.0 * 450.0);					
					// pData[setName]->temp[i] = pData[setName]->temp[i] + dt * 0.5 * pData[setName]->enthalpydot[i];					
				} 

			}
		} else if (setName == "boundary"){
		// Boundary Particles : march the density only.
			#pragma omp parallel for num_threads(NUMTHREADS)
			for (int i = 0; i < ps.n_points(); ++i){
				pData[setName]->dens[i] = pData[setName]->dens[i]  + dt * 0.5 * pData[setName]->densdot[i];
				pData[setName]->temp[i] = pData[setName]->temp[i] + dt * 0.5 * pData[setName]->enthalpydot[i] / (pData[setName]->getSpecificHeat());
			}
		}

	}

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
				// pData[setName]->vel[i]  = add(pData[setName]->vel[i], mult(0.5 * dt,pData[setName]->acc[i]));

				// pData[setName]->dens[i] = pData[setName]->dens[i] + dt * 0.5 * pData[setName]->densdot[i];
				// pData[setName]->temp[i] = pData[setName]->temp[i] + dt * 0.5 * pData[setName]->enthalpydot[i] / pData[setName]->getSpecificHeat();

				if(pData[setName]->pos[i][0] > 1.0E-8 && pData[setName]->pos[i][0] < 50.0E-6 - 1.0E-8){
				// if(pData[setName]->pos[i][0] > 1.0E-8){					
					pData[setName]->temp[i] = pData[setName]->temp[i] + dt * 0.5 * pData[setName]->enthalpydot[i] * 15.0 / (7900.0 * 450.0);
					// pData[setName]->temp[i] = pData[setName]->temp[i] + dt * 0.5 * pData[setName]->enthalpydot[i];					
				} 

				// if(pData[setName]->temp[i] > pData[setName]->getT0()) pData[setName]->temp[i] = pData[setName]->getT0();
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




























































		// 	Unused code, formulation from Fatehi et al. Very sensitive with numerical precision.
		// 	Tensor4 Q_i, R_i;
		// 	Tensor3 S_i, P_i;
		// 	for (const auto& setName_j : setNames){

		// 		const int setID_j = ids[setName_j];
		// 		const auto& ps_j = nsearch->point_set(setID_j);

		// 		for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){

		// 			// Index of neighbor particle in setID_j
		// 			Uint const j = ps_i.neighbor(setID_j, i, _j);

		// 			// Define all the ingedients for the particle interaction
		// 			Real    rho_j = pData[setName_j]->dens[j];
		// 			Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
		// 			Real     dist = length(relpos);

		// 			if ( dist > smoothingLength) continue;
		// 			Real3  reldir = divide(relpos,dist);
		// 			Real    vol_j = pData[setName_j]->vol[j];
        //   			Real      Wij = W_ij(dist, smoothingLength);
		// 			Real3    gWij = gW_ij(dist, reldir, smoothingLength);

		// 			for (int m=0;m<3;m++) for(int n=0;n<3;n++) for(int o=0;o<3;o++) for(int p=0;p<3;p++)
		// 				R_i(m,n,o,p) = R_i(m,n,o,p) + vol_j * relpos[m] * reldir[n] * reldir[o] * gWij[p];

		// 			for (int m=0;m<3;m++) for(int n=0;n<3;n++) for(int k=0;k<3;k++)
		// 				S_i(m,n,k) = S_i(m,n,k) + vol_j * reldir[m] * reldir[n] * gWij[k];
					
		// 			for (int l=0;l<3;l++) for(int o=0;o<3;o++) for(int p=0;p<3;p++)
		// 				P_i(l,o,p) = P_i(l,o,p) + vol_j * relpos[l] * relpos[o] * gWij[p];
		// 		}
		// 	}


		// 	contract_323_to_4(S_i,L_i,P_i,Q_i);


		// 	Q_i.add(R_i);


			

		// 	for(int I=0;I<6;I++){
		// 		Uint m = IDXPAIR[I][0]; Uint n = IDXPAIR[I][1];
		// 		for(int J=0;J<6;J++){
		// 			Uint o = IDXPAIR[J][0]; Uint p = IDXPAIR[J][1];
		// 			_Q_i(I,J) = Q_i(m,n,o,p);
		// 		}
		// 	}
		// 	if(i == 1300){
		// 		std::cout << "before" << std::endl;
		// 		std::cout << _Q_i << std::endl;
		// 	}

		// 	if(dims == 2){
		// 		_Q_i(2,2) = 1.0; _Q_i(4,4) = 1.0; _Q_i(5,5) = 1.0;
		// 	} else if (dims == 1){
		// 		_Q_i(1,1) = 1.0; _Q_i(2,2) = 1.0; _Q_i(3,3) = 1.0; _Q_i(4,4) = 1.0; _Q_i(5,5) = 1.0;		
		// 	}

		// 	JacobiSVD<MatrixXd> svd_L2_i(_Q_i);
		// 	Real cond_second_i = svd_L2_i.singularValues()(0) / svd_L2_i.singularValues()(svd_L2_i.singularValues().size()-1);		
		// 	pData[setName_i]->isFS[i] = cond_second_i;			

		// 	L2_i = toReal3x3(_Q_i.fullPivLu().solve(neg_delta_mn));
			