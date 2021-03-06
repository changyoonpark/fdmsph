#include "sphsolver.hpp"
#define CALC_HEAT 1

SPHSolver::SPHSolver(json& _simData, std::map<std::string,ParticleAttributes*>& _pData)
// SPHSolver::SPHSolver(const json& _simData, std::map<std::string,ParticleAttributes*>& _pData)
: simData(_simData), pData(_pData)
{

// Generate Particle Atrribution Arrays for the fluid
	nsearch = new NeighborhoodSearch((double)simData["smoothingLength"]*1.1,true);
// load the particlesets onto the nsearcher
	std::cout << " ***** Initializing SPH Solver *****" << std::endl;
	currentTime = 0.0;

	for (const auto& pDatEntry : pData){
		ids[pDatEntry.first] = nsearch->add_point_set( pDatEntry.second->pos.front().data(), pDatEntry.second->pos.size(), true, true);
		std::cout << "... Particle set \"" << pDatEntry.first << "\" with " << pDatEntry.second->numParticles << " Particles Set Loaded onto CompactNSearch." << std::endl;
		if (pDatEntry.second->pos.size() != pDatEntry.second->numParticles)	assert(!(pDatEntry.second->pos.size() == pDatEntry.second->numParticles));
		setNames.push_back(pDatEntry.first);
	}

	// Uint min = 99999;
	// for (int i=0;i<pData["fluid"]->pos.size();i++){
	// 	Uint neighs = 0;
	// 	for (int j=0;j<pData["fluid"]->pos.size();j++){
	// 		if (length(sub(pData["fluid"]->pos[i],pData["fluid"]->pos[j])) <= (Real)simData["smoothingLength"] && i != j){
	// 			neighs += 1;
	// 		}			
	// 	}
	// 	if (neighs < min){
	// 		std::cout << "new min neighs : " << min << std::endl;
	// 		min = neighs;
	// 	}
	// }

	totParticles = 0;
	Uint numSets = 0;

	for (const auto& pSet : nsearch->point_sets()){
		totParticles += pSet.n_points();
		numSets += 1;
	}

	std::cout << "... Loaded " << totParticles << " particles, with " << numSets << " point sets." << std::endl;

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
	setInitialConfigurationNeighbors();

	if( simData["useVariableVolume"] == "Yes" ) initializeMass();


	std::cout << "--- total number of particle sets : " << nsearch->n_point_sets() << std::endl;
	std::cout << "--- total number of particles     : " << totParticles << std::endl;

	// setInitialDeformation();
}

void SPHSolver::addFluidInletParticles(int t){

	if ( t % (int)( ((Real) (pData["fluid"]->dx)) / ((Real) (pData["fluid"]->inletSpeed) * (Real) simData["dt"] ) ) != 0 ) return;

	std::cout << "... Adding inlet particles" << std::endl;
	pData["fluid"]->addParticlesToFluid();

	std::cout << "... Freeing nsearch instance" << std::endl;
	delete nsearch;
	nsearch = new NeighborhoodSearch((double)simData["smoothingLength"],true);

	for (const auto& pDatEntry : pData){
		ids[pDatEntry.first] = nsearch->add_point_set( pDatEntry.second->pos.front().data(), pDatEntry.second->pos.size(), true, true);
		std::cout << "... Particle set \"" << pDatEntry.first << "\" with " << pDatEntry.second->numParticles << " Particles Set Loaded onto CompactNSearch." << std::endl;
		if (pDatEntry.second->pos.size() != pDatEntry.second->numParticles)	assert(!(pDatEntry.second->pos.size() == pDatEntry.second->numParticles));
	}

	// There is a bug in the nsearch code that freezes the code at this part.
	// nsearch->resize_point_set(ids["fluid"], 
	// 							  pData["fluid"]->pos.front().data(),
	// 							  pData["fluid"]->pos.size());
								
	// if (pData["fluid"]->pos.size() != pData["fluid"]->numParticles)	assert(!(pData["fluid"]->pos.size() == pData["fluid"]->numParticles));

	std::cout << "... Resized searching algorithm for the inlet particles" << std::endl;

}


void SPHSolver::initializeMass(){

	std::cout << "--- Initializing particle volume / mass." << std::endl;
	const Real smoothingLength = (Real) simData["smoothingLength"];
	Real totVol = 0.0;

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

			pData[setName_i]->vol[i]  = (1.0/kernelSum);
			totVol += pData[setName_i]->vol[i];
			pData[setName_i]->mass[i] =  pData[setName_i]->dens[i] * pData[setName_i]->vol[i];

		}
  }
  simData["totalVolume"] = (Real)totVol;  

}

void SPHSolver::setInitialConfigurationNeighbors(){
	std::cout << "--- Initializing Initial Configuration Neighbor Map" << std::endl;
	const Real smoothingLength = (Real) simData["smoothingLength"];
	Real totVol = 0.0;
	for (const auto& setName_i : setNames){
		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);
		for (int i = 0; i < ps_i.n_points(); ++i){
			for (const auto& setName_j : setNames){

				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);

				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){
					Uint const j = ps_i.neighbor(setID_j, i, _j);					
					pData[setName_i]->nMap[i].push_back(std::tuple<std::string,Uint>{setName_j,j});
				}
			}			
		}
	}
}

void SPHSolver::setInitialDeformation(){
	for (const auto& setName_i : setNames){

		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);

		for (int i = 0; i < ps_i.n_points(); ++i){
				if(setName_i == "fluid" && pData[setName_i]->isSolid[i]){
					pData[setName_i]->pos[i] = add(pData[setName_i]-> originPos[i] , Real3{((0) * pData[setName_i]->originPos[i][0]  + 0.0  * pData[setName_i]->originPos[i][1]),
																							(0.01 * pData[setName_i]->originPos[i][0]   + 0.0  * pData[setName_i]->originPos[i][1]),0} );				
				}
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
 						  	  Real  vol_j,      Real delta,
 						  	  Real  soundSpeed, Real smoothingLength,
 						  	  Real  mass_j,     Real  dist,
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
	} 
	else if (simData["kernelType"] == "Wendland6"){
		if (simData["dimensions"] == 2){
			std::cout << "--- Kernel Type : Wendland 6 (2D)" << std::endl;
			W_ij  =  W_Wendland6_2D_h;
			gW_ij = gW_Wendland6_2D_h;
		} else{
			assert(0 && "!!! Unimplemented Kernel Type.");
		}
	}
	else if (simData["kernelType"] == "Wendland_2h"){
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

	// for (const auto& setName_i : setNames){
		
	// 			const int setID_i = ids[setName_i];
	// 			const auto& ps_i = nsearch->point_set(setID_i);
		
	// 			for (int i = 0; i < ps_i.n_points(); ++i){		
	// 				Real3 pos_i  = pData[setName_i]->pos[i];
	// 				std::cout << pos_i[0] << ", " << pos_i[1] << ", " << pos_i[2] << std::endl;
	// 			}
	// }

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



void SPHSolver::polarDecompose(	MatrixXd& F, MatrixXd& R, MatrixXd& U){
	MatrixXd W = MatrixXd::Zero(3,3), SIGM = MatrixXd::Zero(3,3), V = MatrixXd::Zero(3,3);
	JacobiSVD<MatrixXd> svd_F(F, ComputeFullU | ComputeFullV);
	SIGM = svd_F.singularValues().asDiagonal();
	W    = svd_F.matrixU();
	V    = svd_F.matrixV();
	R    = W * V.transpose();
	U	 = V * SIGM * V.transpose();
}

void SPHSolver::getStretch(MatrixXd& U, Real3& lambdas, Real3x3& dirs, Real3x3& dirs_before){

	EigenSolver<MatrixXd> es(U);

	Real3x3 basis{zeromat};
	Real3   _lambdas{zerovec};

	toReal3(es.eigenvalues().real(), _lambdas);
	toReal3(es.eigenvectors().col(0).real(), basis[0]);
	toReal3(es.eigenvectors().col(1).real(), basis[1]);
	toReal3(es.eigenvectors().col(2).real(), basis[2]);

	//ugly stupid routine to sort the eigvecs
	Real a = std::abs(dot(basis[0], dirs_before[0])), b = std::abs(dot(basis[0], dirs_before[1])), c = std::abs(dot(basis[0], dirs_before[2]));
	if( a > b && a > c){
		// now[0] = before[0]
		dirs[0] = basis[0];
		lambdas[0] = _lambdas[0];
		Real d = std::abs(dot(basis[1], dirs_before[1])), e = std::abs(dot(basis[1], dirs_before[2]));
		if ( d > e ){
			// now[1] = before[1], now[2] = before[2]
			dirs[1] = basis[1]; dirs[2] = basis[2];
			lambdas[1] = _lambdas[1]; lambdas[2] = _lambdas[2];
		} else{
			// now[1] = before[2], now[2] = before[1]
			dirs[2] = basis[1]; dirs[1] = basis[2];
			lambdas[2] = _lambdas[1]; lambdas[1] = _lambdas[2];
		}
	} else if ( b > a && b > c){
		// now[0] = before[1]
		dirs[1] = basis[0];
		lambdas[1] = _lambdas[0];
		Real d = std::abs(dot(basis[1], dirs_before[0])), e = std::abs(dot(basis[1], dirs_before[2]));
		if ( d > e ){
			// now[1] = before[0], now[2] = before[2]
			dirs[0] = basis[1]; dirs[2] = basis[2];
			lambdas[0] = _lambdas[1]; lambdas[2] = _lambdas[2];
		} else{
			// now[1] = before[2], now[2] = before[0]
			dirs[2] = basis[1]; dirs[0] = basis[2];
			lambdas[2] = _lambdas[1]; lambdas[0] = _lambdas[2];
		}		
	} else{
		// now[0] = before[2]
		dirs[2] = basis[0];
		lambdas[2] = _lambdas[0];
		Real d = std::abs(dot(basis[1], dirs_before[0])), e = std::abs(dot(basis[1], dirs_before[1]));
		if ( d > e ){
			// now[1] = before[0], now[2] = before[1]
			dirs[0] = basis[1]; dirs[1] = basis[2];
			lambdas[0] = _lambdas[1]; lambdas[1] = _lambdas[2];
		} else{
			// now[1] = before[1], now[2] = before[0]
			dirs[1] = basis[1]; dirs[0] = basis[2];
			lambdas[1] = _lambdas[1]; lambdas[0] = _lambdas[2];
		}		
	}
}



void SPHSolver::computeDeformationGradient(Uint t){
	using namespace RealOps;
	const Real dx = (Real) simData["dx"];
	const Real dt = (Real) simData["dt"];
	const Uint dims = simData["dimensions"];
	const Real smoothingLength = (Real) simData["smoothingLength"];
	
	for (const auto& setName_i : setNames){


		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);
		
		// #pragma omp parallel for num_threads(NUMTHREADS)
		for (int i = 0; i < ps_i.n_points(); ++i){


			// This should only be computed for a solidified fluid particle.
			if ( !pData[setName_i]->isSolid[i] ) continue;		


			Real3& pos_i = pData[setName_i]->pos[i];
			// Renormalization Tensor (First Derivative)
			// This can claimed to be "inconsistent", but it is something that we accept
			MatrixXd _L_i(3,3); _L_i = MatrixXd::Zero(3,3);
			Real3x3& L_i = pData[setName_i]->L[i];
			L_i = zeromat;

			// Initialize Deformation Gradient and the cauchy stress tensor
			if( t <= 1000){
				pData[setName_i]->temp[i] = 10.0 * pData[setName_i]->originPos[i][1] * pData[setName_i]->originPos[i][1] * (Real)t * 0.2;
			} else{
				pData[setName_i]->temp[i] = std::max(
											(Real)(10.0 * pData[setName_i]->originPos[i][1] * pData[setName_i]->originPos[i][1] * 1000.0 * 0.2 - 
											10.0 * pData[setName_i]->originPos[i][1] * pData[setName_i]->originPos[i][1] * ((Real)t-1000.0) * 0.2),
											(Real)0.0);
			}

			Real thermalDeform = 1.0 + 0.0005 * pData[setName_i]->temp[i];
			pData[setName_i]->defoGrad[i] = zeromat;
			if ( pData[setName_i]->isSolid[i] && setName_i == "fluid" ){
				pData[setName_i]->defoGrad_thermal[i] = mult((thermalDeform),identity());
			}

			// Compute the Renormalization Matrix / density gradient using the renormalization matrix.
			for (const auto& setName_j : setNames){

				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);

				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){

					// Index of neighbor particle in setID_j
					Uint const j = ps_i.neighbor(setID_j, i, _j);

					// This should only be computed for a solidified fluid particle.
					if( !pData[setName_j]->isSolid[j] ) continue;

					// Define all the ingedients for the particle interaction
					Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
					Real     dist = length(relpos);

					if ( dist > smoothingLength ) continue;
					Real3  reldir = divide(relpos,dist);
					Real    vol_j = pData[setName_j]->vol[j];
					Real3    gWij = gW_ij(dist, reldir, smoothingLength);
					

          			Real3x3 v_j_gWij_rij = mult(- vol_j,tensorProduct(gWij,relpos));

					// Compute First derivative Renormalization Tensor
					L_i = add(L_i, v_j_gWij_rij);

					// Compute the reference Map Gradient
					pData[setName_i]->defoGrad[i][0] = add(pData[setName_i]->defoGrad[i][0], mult( (pData[setName_j]->originPos[j][0] - pData[setName_i]->originPos[i][0]) * vol_j, gWij));
					pData[setName_i]->defoGrad[i][1] = add(pData[setName_i]->defoGrad[i][1], mult( (pData[setName_j]->originPos[j][1] - pData[setName_i]->originPos[i][1]) * vol_j, gWij)); 
					pData[setName_i]->defoGrad[i][2] = add(pData[setName_i]->defoGrad[i][2], mult( (pData[setName_j]->originPos[j][2] - pData[setName_i]->originPos[i][2]) * vol_j, gWij));

				}
			}

			// Compute the inverse of the renormalization tensor.
			setDims(L_i,dims);
			L_i = inv3x3(L_i, -1);

			pData[setName_i]->defoGrad[i][0] = mult(pData[setName_i]->L[i],pData[setName_i]->defoGrad[i][0]);
			pData[setName_i]->defoGrad[i][1] = mult(pData[setName_i]->L[i],pData[setName_i]->defoGrad[i][1]);
			pData[setName_i]->defoGrad[i][2] = mult(pData[setName_i]->L[i],pData[setName_i]->defoGrad[i][2]);
 
			// // Invert the gradient of the reference map
			setDims(pData[setName_i]->defoGrad[i], dims);
			pData[setName_i]->defoGrad[i] = inv3x3(pData[setName_i]->defoGrad[i],-1);
			pData[setName_i]->defoGrad_withoutThermal[i] = mult(1.0 / thermalDeform, pData[setName_i]->defoGrad[i]);
			setDims(pData[setName_i]->defoGrad_withoutThermal[i], dims);

			// Compute the creep component			
			MatrixXd U = MatrixXd::Zero(3,3), R = MatrixXd::Zero(3,3);
			MatrixXd F_withoutThermal = MatrixXd::Zero(3,3); toMatrix3d(pData[setName_i]->defoGrad_withoutThermal[i], F_withoutThermal);

			polarDecompose(F_withoutThermal, R, U);						
			getStretch(U, pData[setName_i]->stretch_total[i], pData[setName_i]->stretch_dirs[i], pData[setName_i]->stretch_dirs_before[i]);



			toMatrix3d(pData[setName_i]->defoGrad_withoutThermal[i], U);
			pData[setName_i]->defoGrad_elastic[i] = pData[setName_i]->defoGrad_withoutThermal[i];


			pData[setName_i]->vol[i] = pData[setName_i]->originVol[i] * det3x3(pData[setName_i]->defoGrad[i]);
		}
	}
}

void SPHSolver::computeInteractions(Uint t){
	using namespace RealOps;
	// The fluid.
	const Real dx = (Real) simData["dx"];
	const Real dt = (Real) simData["dt"];
	const Uint dims = simData["dimensions"];
	const Real smoothingLength = (Real) simData["smoothingLength"];

	VectorXd neg_delta_mn(6);
	neg_delta_mn(0) = -1.0;	neg_delta_mn(1) = -1.0;	neg_delta_mn(2) = -1.0;
	neg_delta_mn(3) =  0.0; neg_delta_mn(4) =  0.0; neg_delta_mn(5) =  0.0;

	Real3 neg_delta_mn_2D = Real3{-1.0,  // (0,0)
								  -1.0,  // (2,2)
		   							 0}; // (0,2)

	

	std::cout << "		|------ ... Performing Fluid - (Fluid / Boundary) Interactions"<< std::endl;

	#if CALC_HEAT
	for (const auto& setName_i : setNames){

		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);
		
		const Real oneOverLambda = (Real) pData[setName_i]->getOneOverLambda();
		const Real eta = (Real) pData[setName_i]->getEta();

		#pragma omp parallel for num_threads(NUMTHREADS)


		for (int i = 0; i < ps_i.n_points(); ++i){
			// std::cout << i << std::endl;
			const Real& vol_i_o = pData[setName_i]->originVol[i];
			Real    m_i  = pData[setName_i]->mass[i];
			Real  rho_i  = pData[setName_i]->dens[i];
			Real3& pos_i = pData[setName_i]->pos[i];
			Real3 u_i = sub(pos_i,pData[setName_i]->originPos[i]);

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
			// Initialize Velocity Gradient
			pData[setName_i]->velGrad[i] = zeromat;
			// Initialize Tau Gradient
			pData[setName_i]->tauGrad[i] = zeroijk;
			// Initialize Particle Density Gradient
			pData[setName_i]->particleDensityGrad[i] = zerovec;


			// Clear the time derivatives.ff
			pData[setName_i]->tauDot[i]  = zeromat;
			pData[setName_i]->acc[i] = Real3{0.0,0.0,0.0};
			pData[setName_i]->densdot[i] = 0.0;
			pData[setName_i]->enthalpydot[i] = 0.0;
			pData[setName_i]->normalVec[i] = zerovec;


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


					// Compute First derivative Renormalization Tensor
					if (pData[setName_j]->isSolid[i]){
						Real3x3 v_j_gWij_rij = mult(- vol_j,tensorProduct(gWij,relpos));
						L_i = add(L_i, v_j_gWij_rij);
					} 



					// Compute Renormalized Density Gradient
					pData[setName_i]->densGrad[i]   = add(pData[setName_i]->densGrad[i], mult((rho_j - rho_i) * vol_j,gWij));
					pData[setName_i]->tempGrad[i]   = add(pData[setName_i]->tempGrad[i], mult( (pData[setName_j]->temp[j] - pData[setName_i]->temp[i]) * vol_j, gWij));

					// Compute the velocity gradient
					pData[setName_i]->velGrad[i][0] = add(pData[setName_i]->velGrad[i][0], mult( (pData[setName_j]->vel[j][0] - pData[setName_i]->vel[i][0]) * vol_j, gWij));
					pData[setName_i]->velGrad[i][1] = add(pData[setName_i]->velGrad[i][1], mult( (pData[setName_j]->vel[j][1] - pData[setName_i]->vel[i][1]) * vol_j, gWij));
					pData[setName_i]->velGrad[i][2] = add(pData[setName_i]->velGrad[i][2], mult( (pData[setName_j]->vel[j][2] - pData[setName_i]->vel[i][2]) * vol_j, gWij));

					// Compute the gradient of the deviatoric stress				
					pData[setName_i]->tauGrad[i][0][0] = add(pData[setName_i]->tauGrad[i][0][0], mult( (pData[setName_j]->tau[j][0][0] - pData[setName_i]->tau[i][0][0]) * vol_j, gWij));
					pData[setName_i]->tauGrad[i][0][1] = add(pData[setName_i]->tauGrad[i][0][1], mult( (pData[setName_j]->tau[j][0][1] - pData[setName_i]->tau[i][0][1]) * vol_j, gWij));
					pData[setName_i]->tauGrad[i][0][2] = add(pData[setName_i]->tauGrad[i][0][2], mult( (pData[setName_j]->tau[j][0][2] - pData[setName_i]->tau[i][0][2]) * vol_j, gWij));

					pData[setName_i]->tauGrad[i][1][0] = add(pData[setName_i]->tauGrad[i][1][0], mult( (pData[setName_j]->tau[j][1][0] - pData[setName_i]->tau[i][1][0]) * vol_j, gWij));
					pData[setName_i]->tauGrad[i][1][1] = add(pData[setName_i]->tauGrad[i][1][1], mult( (pData[setName_j]->tau[j][1][1] - pData[setName_i]->tau[i][1][1]) * vol_j, gWij));
					pData[setName_i]->tauGrad[i][1][2] = add(pData[setName_i]->tauGrad[i][1][2], mult( (pData[setName_j]->tau[j][1][2] - pData[setName_i]->tau[i][1][2]) * vol_j, gWij));

					pData[setName_i]->tauGrad[i][2][0] = add(pData[setName_i]->tauGrad[i][2][0], mult( (pData[setName_j]->tau[j][2][0] - pData[setName_i]->tau[i][2][0]) * vol_j, gWij));
					pData[setName_i]->tauGrad[i][2][1] = add(pData[setName_i]->tauGrad[i][2][1], mult( (pData[setName_j]->tau[j][2][1] - pData[setName_i]->tau[i][2][1]) * vol_j, gWij));
					pData[setName_i]->tauGrad[i][2][2] = add(pData[setName_i]->tauGrad[i][2][2], mult( (pData[setName_j]->tau[j][2][2] - pData[setName_i]->tau[i][2][2]) * vol_j, gWij));



					// Compute the normal
					pData[setName_i]->normalVec[i] = add(pData[setName_i]->normalVec[i], 
										mult(Wij,relpos)
									);

				}
			}



			// Compute the inverse of the renormalization tensor.
			setDims(L_i,dims);
			toMatrix3d(L_i,_L_i);

			JacobiSVD<MatrixXd> svd_L_i(_L_i);
			// svd_L_i.singularValues()(0) / svd_L_i.singularValues()(svd_L_i.singularValues().size()-1
			pData[setName_i]->conditionNumber[i] = 
				std::max(svd_L_i.singularValues()(0) / svd_L_i.singularValues()(2),
						 svd_L_i.singularValues()(0) / svd_L_i.singularValues()(1));		
			
			_L_i = _L_i.inverse();
			toReal3x3(_L_i,L_i);

			// Assign this to be the renormalization tensor for the initial configuration for t == 0
			if(t == 0){
				pData[setName_i]->L_o[i] = L_i;

			} 

			// Compute the displacement gradient.
			// This loop goes over the neighbor map defined at the initial configuration.
			// Initialize Displacement GRadient
			// pData[setName_i]->dispGrad[i] = zeromat;
			// for(const auto& j_ : pData[setName_i]->nMap[i]){
			// 	// std::cout << j << std::endl;
			// 	const std::string setName_j = std::get<0>(j_);
			// 	const Uint j = std::get<1>(j_);
			// 	Real3  relpos_o   = sub(pData[setName_i]->originPos[i],pData[setName_j]->originPos[j]);
			// 	Real     dist_o   = length(relpos_o);
			// 	if ( dist_o > smoothingLength ) continue;
			// 	Real3   reldir_o  = divide(relpos_o,dist_o);
			// 	Real3    gWij_o = gW_ij(dist_o, reldir_o, smoothingLength);				
			// 	Real3 u_j = sub(pData[setName_j]->pos[j],pData[setName_j]->originPos[j]);
			// 	Real3 u_ji = sub(u_j, u_i);
			// 	Real vol_j_o = pData[setName_j]->originVol[j];

			// 	pData[setName_i]->dispGrad[i][0] = add(pData[setName_i]->dispGrad[i][0], mult( u_ji[0] * vol_j_o, gWij_o));
			// 	pData[setName_i]->dispGrad[i][1] = add(pData[setName_i]->dispGrad[i][1], mult( u_ji[1] * vol_j_o, gWij_o));
			// 	pData[setName_i]->dispGrad[i][2] = add(pData[setName_i]->dispGrad[i][2], mult( u_ji[2] * vol_j_o, gWij_o));

			// }
			// pData[setName_i]->dispGrad[i][0] = mult(pData[setName_i]->L_o[i],pData[setName_i]->dispGrad[i][0]);
			// pData[setName_i]->dispGrad[i][1] = mult(pData[setName_i]->L_o[i],pData[setName_i]->dispGrad[i][1]);
			// pData[setName_i]->dispGrad[i][2] = mult(pData[setName_i]->L_o[i],pData[setName_i]->dispGrad[i][2]);
			// // Now that we have the displacement gradient, compute the deformation gradient
			// pData[setName_i]->defoGrad[i] = add(pData[setName_i]->dispGrad[i],identity());

			// Normalize the normal vector.
			pData[setName_i]->normalVec[i] = normalize(pData[setName_i]->normalVec[i]);
			const Real3 pos_T = add(pData[setName_i]->pos[i], mult(dx, pData[setName_i]->normalVec[i]));

			// Set Neumann BCs for the solid on the free-surfaces
			// if( pData[setName_i]->particleDensity[i] < 0.9 ) {pData[setName_i]->secondPKStress[i] = zeromat; }

			// Compute Tau Dot
			const Real3x3 velGrad_i_T =  transpose(pData[setName_i]->velGrad[i]);

			pData[setName_i]->tauDot[i] = 
				add(
					mult(
						oneOverLambda,
						sub(
							mult(eta, add(pData[setName_i]->velGrad[i], velGrad_i_T)),
							pData[setName_i]->tau[i]
						)
					),
					add(
						mult(pData[setName_i]->velGrad[i], pData[setName_i]->tau[i]), 
						mult(pData[setName_i]->tau[i], velGrad_i_T)
					)
				);

			
			// Apply the renormalization matrix
			// pData[setName_i]->densGrad[i]  = mult(L_i,pData[setName_i]->densGrad[i]);
			pData[setName_i]->tempGrad[i]  = mult(L_i,pData[setName_i]->tempGrad[i]);

			// Normalize the normal vector.
			pData[setName_i]->velGrad[i][0] = mult(pData[setName_i]->L[i], pData[setName_i]->velGrad[i][0]);
			pData[setName_i]->velGrad[i][1] = mult(pData[setName_i]->L[i], pData[setName_i]->velGrad[i][1]);
			pData[setName_i]->velGrad[i][2] = mult(pData[setName_i]->L[i], pData[setName_i]->velGrad[i][2]);


			// For solidified fluid, compute the cauchy stress and apply some damping
			Real3x3 strainTensor_i{zeromat};
			Real3x3& defoGrad_i = pData[setName_i]->defoGrad[i];

			// if(pData[setName_i]->pos[i][0] < 1.0){
			// 	defoGrad_i = identity();
			// } else{
			// 	defoGrad_i = mult(2.0, identity());
			// }
			
			pData[setName_i]->secondPKStress[i] = zeromat;
			Real3x3& secondPKStress_i = pData[setName_i]->secondPKStress[i];

			if( pData[setName_i]->isSolid[i] ){	

				// Infinitesimal strain
				strainTensor_i = mult(0.5, sub( mult(transpose(pData[setName_i]->defoGrad_elastic[i]),pData[setName_i]->defoGrad_elastic[i]), identity()));
				Real3x3 strainTensorVol_i = mult(0.333333333333 * trace(strainTensor_i), identity());
				Real3x3 strainTensorDev_i = sub(strainTensor_i, strainTensorVol_i);

				// Second Piola-Kirchhoff Stress				
				secondPKStress_i = add(
										mult(     (Real) pData[setName_i]->getBulkModulus() , strainTensorVol_i),
										mult(2. * (Real) pData[setName_i]->getShearModulus(), strainTensorDev_i)
									);

			}

			if (pData[setName_i]->particleDensity[i] < 0.95){
				pData[setName_i]->isFS[i] = true;			
			} else{
				pData[setName_i]->isFS[i] = false;			
			}



			for (const auto& setName_j : setNames){
				
				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);

				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){

					// Index of neighbor particle in setID_j
					Uint const j = ps_i.neighbor(setID_j, i, _j);

					// Define all the ingedients for the particle interaction
					Real    rho_j = pData[setName_j]->dens[j];
					Real3&  pos_j = pData[setName_j]->pos[j];
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

					// Determine whether the particle is a real free-surface particle:
					if (pData[setName_i]->isFS[i]){

						Real3 vec_Ti = sub(pos_T,pos_i);
						Real3 vec_ji = sub(pos_j,pos_i);
						Real abcos = dot(vec_Ti,vec_ji);
						Real theta = std::acos(abcos / (dx * dist));
						// area = (0.333*h)*(dist)*sin(theta) = (0.333*h) * sdotvec_jT
						Real sdotvec_jT = dist * std::sin(theta);
						if (dist >= (Real)simData["alpha"] * 0.47140452079 * smoothingLength && length(pos_j, pos_T) < (Real)simData["alpha"] * 0.333333 * smoothingLength){
							pData[setName_i]->isFS[i] = false;
						} else if (dist < (Real)simData["alpha"] * 0.47140452079 * smoothingLength && std::abs(dot(pData[setName_i]->normalVec[i],sub(pos_j, pos_T))) + std::abs(sdotvec_jT) < (Real)simData["alpha"] * 0.333333 * smoothingLength){
							pData[setName_i]->isFS[i] = false;
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

			if(dims == 2){
				G(2,2) = 1.0; G(4,4) = 1.0; G(5,5) = 1.0;
			} else if(dims == 1){
				G(1,1) = 1.0; G(2,2) = 1.0;
				G(3,3) = 1.0; G(4,4) = 1.0;
				G(5,5) = 1.0;
			}

			JacobiSVD<MatrixXd> svd_G_i(G);
			Real cond_second_i = svd_G_i.singularValues()(0) / svd_G_i.singularValues()(svd_G_i.singularValues().size()-1);		
			
			pData[setName_i]->conditionNumber[i] = cond_second_i;			
			
			VectorXd L2_i_vec = G.fullPivLu().solve(neg_delta_mn);
			L2_i = toReal3x3From6(L2_i_vec,dims);
			setDims2(L2_i,dims);

		}

	} 
	#else
		for (const auto& setName_i : setNames){

		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);

		#pragma omp parallel for num_threads(NUMTHREADS)
		for (int i = 0; i < ps_i.n_points(); ++i){
			// Clear the time derivatives.
			pData[setName_i]->acc[i] = Real3{0.0,0.0,0.0};
			pData[setName_i]->densdot[i] = 0.0;
			pData[setName_i]->enthalpydot[i] = 0.0;			

			pData[setName_i]->tauDot = Real3x3{Real3{0,0,0},Real3{0,0,0},Real3{0,0,0}};			
		}

	} 
	#endif


	smearStress(t);

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
			const Real& vol_i_o = pData[setName_i]->originVol[i];

			const Real&   P_i = EOS(T_i, pData[setName_i]->getT0(),
									     rho_i, pData[setName_i]->getRho0(),
									     pData[setName_i]->getSoundSpeed(),
										 pData[setName_i]->getThermalExpansion());

			const Real3& vel_i   = pData[setName_i]->vel[i];
			const Real3& pos_i_o = pData[setName_i]->originPos[i];

			const Real3& tempGrad_i = pData[setName_i]->tempGrad[i];			
			const Real3& densGrad_i = pData[setName_i]->densGrad[i];
	        const Real3& normal_i =  pData[setName_i]->normalVec[i];
			const Real3x3& L_i  = pData[setName_i]->L[i];
			const Real3x3& L2_i = pData[setName_i]->L2[i];
			const Real3x3& velGrad_i = pData[setName_i]->velGrad[i];			
			const Real& cond_i = pData[setName_i]->conditionNumber[i];
			const Real& pDens_i = pData[setName_i]->particleDensity[i];

			const Real3x3 tensorViscosity{zeromat};

			// Temporary storage to keep the summation of the tensor product
			Real3x3x3 sig_ij_gW_ij_i_vol_j{zeromat,zeromat,zeromat};

			// Sum of (\sigma_i - \sigma_j) (X) \nabla W_ij

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
					
					const Real&      m_j   = pData[setName_j]->mass[j];
					const Real&    vol_j   = pData[setName_j]->vol[j];

					const Real3& vel_j = pData[setName_j]->vel[j];
					const Real3& normal_j =  pData[setName_j]->normalVec[j];
					const Real3& tempGrad_j = pData[setName_j]->tempGrad[j];
					const Real& pDens_j = pData[setName_j]->particleDensity[j];

					Real   rho_i_temp, rho_j_temp;
					Real3  relvel = sub(pData[setName_i]->vel[i],pData[setName_j]->vel[j]);
					Real3  reldir = divide(relpos,dist);

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

					
					// Compute Particle Density Gradient
					pData[setName_i]->particleDensityGrad[i] = add(pData[setName_i]->particleDensityGrad[i], mult((pDens_j - pDens_i) * vol_j,gWij));

					// Delta SPH Diffusion
					pData[setName_i]->densdot[i] += dot(mult(rho_i, relvel), gWij) * vol_j;

					pData[setName_i]->densdot[i] += diffusiveTerm_ij(rho_i,rho_j,
																	vol_j, (Real) simData["delta"],
																	pData[setName_j]->getSoundSpeed(), smoothingLength,
																	m_j,dist,
																	relpos,gWij,
																	densGrad_i,densGrad_j);

					
					#if CALC_HEAT
						// Heat Transfer between particles. Note that the boundary densities must be set to the
						// Actual density of the boundary material, to account for the correct thermal diffusivity.
						rho_i_temp = isBoundary_i ? pData[setName_i]->getMaterialDensity() : rho_i ;
						rho_j_temp = isBoundary_j ? pData[setName_j]->getMaterialDensity() : rho_j ;
	

						// If the laplacian corrector cannot be defined, use the conventional operator.						
						Real heat;

						if (cond_i < 40.0){
							heat = consistentHeatTransfer(L2_i,tempGrad_i,
														T_i, T_j, m_j, rho_i_temp, rho_j_temp, k_i, k_j,
														relpos, reldir,
														dist, gWij, vol_j);
						} else{
							heat = inconsistentHeatTransfer(tempGrad_i, tempGrad_j,
															T_i, T_j, m_j, rho_i_temp, rho_j_temp, k_i, k_j,
															relpos, reldir,
															dist, gWij, vol_j);		
						}

						// pData[setName_i]->enthalpydot[i] += heat * 15.0 / (7900.0 * 450.0);
						pData[setName_i]->enthalpydot[i] += heat;
					#endif					

					// Acceleration Due to pressure gradient. (For newtonian Fluid)
					// Real3 a = add(pressureAcc_ij(m_j,gWij,P_i,P_j,rho_i,rho_j,vol_j), viscosityAcc_ij(m_j, relpos, gWij, rho_i, rho_j,dist, relvel, mu_j));
					
					// Real3 a = add( pressureAcc_ij(m_j,gWij,P_i,P_j,rho_i,rho_j,vol_j), 
					// 			   extraStress_acc_ij(m_j,gWij,rho_i,rho_j,
					// 			   					  zeromat,
					// 							   	  zeromat)
					// 			  );

					// Interparticle forces, (IFF surface tension model).
					// a = add(a, iif_acc(simData["IIFCoeff"], reldir, dist, smoothingLength, m_i, m_j));


					// Damping for solid model
					if( pData[setName_i]->isSolid[i] && pData[setName_j]->isSolid[j]){
						pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i],
														mult(	- pData[setName_j]->mass[j] / pData[setName_i]->dens[i],
																mult(
																	scalarViscosity_ij( (Real) simData["Cl"],
																						(Real) simData["Cq"],
																							smoothingLength,
																						(Real) pData[setName_i]->getSoundSpeed(),
																							relpos,
																							relvel
																					),
																	gWij
																)
														)
							 						   );

						pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i],
								mult(	- (Real) simData["damping"] * pData[setName_j]->vol[j],
										mult(Wij, relvel))
								);

					}

					// Acceleration Due to fluid stress, Upper convected model
					// This contribution should only happens between non-solid fluid particles
					// if ( !pData[setName_i]->isSolid[i] ){
					// 	if ( !pData[setName_j]->isSolid[j] ){
					// 		pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i], 
					// 										add( pressureAcc_ij(m_j,gWij,P_i,P_j,rho_i,rho_j,vol_j), 
					// 												extraStress_acc_ij(m_j,gWij,rho_i,rho_j,
					// 																pData[setName_i]->tau[i],
					// 																pData[setName_j]->tau[j])
					// 										)
					// 									);
					// 	} else{
					// 		// Implement this part for FSI (non-solid fluid (i) <> solidified fluid (j) )
					// 	}
					// }						

				}

			}


			// If particle i is not a boundary particle, accumulate the momentum contribution.
			if (!isBoundary_i){
				// Loop over the initial neighbors.
				for(const auto& j_ : pData[setName_i]->nMap[i]){
					// If i is a solid particle and j is a solid particle, accumulate the cauchy stress contributions
					const std::string setName_j = std::get<0>(j_);
					const Uint j = std::get<1>(j_);

					if ( pData[setName_i]->isSolid[i] ){
						if( pData[setName_j]->isSolid[j] ){

							const Real&    vol_j_o = pData[setName_j]->originVol[j];
							const Real3&   pos_j_o = pData[setName_j]->originPos[j];
							Real3  relpos_o = sub(pos_i_o,pos_j_o);
							Real3  gWij_o = gW_ij(length(relpos_o), dir(relpos_o), smoothingLength);

							// pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i], 			
							// 								mult(
							// 									add(											
							// 										 mult(pData[setName_i]->defoGrad[i],pData[setName_i]->secondPKStress[i]),
							// 										 mult(pData[setName_j]->defoGrad[j],pData[setName_j]->secondPKStress[j])
							// 										),
							// 										 mult( 1.0 * vol_i_o * vol_j_o / m_i, mult(pData[setName_i]->L_o[i], gWij_o))
							// 									)
							// 							  );

							// from ganzenmuller
							pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i], 			
															mult(   1.0 * vol_i_o * vol_j_o / m_i,
																	add(											
																		mult(mult(pData[setName_i]->defoGrad[i],pData[setName_i]->secondPKStress[i]), mult(pData[setName_i]->L_o[i], gWij_o)),
																		mult(mult(pData[setName_j]->defoGrad[j],pData[setName_j]->secondPKStress[j]), mult(pData[setName_j]->L_o[j], gWij_o))
																	)
																)
														  );


							// Real J_i = det3x3(pData[setName_i]->defoGrad[i]);
							// Real J_j = det3x3(pData[setName_j]->defoGrad[j]);
							// pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i], 			
							// 								mult(
							// 									add(
							// 										mult(											
							// 											mult(mult(1./J_i,pData[setName_i]->defoGrad[i]),
							// 												pData[setName_i]->secondPKStress[i]
							// 												),
							// 											transpose(pData[setName_i]->defoGrad[i])
							// 										),
							// 										mult(
							// 											mult(mult(1./J_j,pData[setName_j]->defoGrad[j]),
							// 												pData[setName_j]->secondPKStress[j]
							// 												),
							// 											transpose(pData[setName_j]->defoGrad[j])
							// 										)
							// 										),
							// 										mult( 1.0 * (vol_i * J_i) * (vol_j * J_j) / m_i, mult(pData[setName_i]->L[i], gWij))
							// 									)
							// 							  );


						} 
						else{
							// Implement this part for FSI (non-solid fluid (j) <> solidified fluid (i) )
						}
					}

				}
			}

				// // // If the boundary particle is a sensor, accumulate the data.
				// if(isSensor_i && !isBoundary_j){
				// 	Real3 f = add(pressureAcc_ij(m_j,gWij,P_i,P_j,rho_i,rho_j,vol_j), viscosityAcc_ij(m_j, relpos, gWij, rho_i, rho_j,dist, relvel, mu_j));
				// 	f = mult(m_i,f);
				// 	pData[setName_i]->forceSensed[i] = add(pData[setName_i]->forceSensed[i], f);
				// 	#if CALC_HEAT
				// 		pData[setName_i]->heatSensed[i] += heat;
				// 	#endif
				// }
				

			// Compute shifting
			
			if(pData[setName_i]->isFS[i]){
				pData[setName_i]->shift[i] = mult(- 1.666666 * smoothingLength * pData[setName_i]->vol[i] * dt, 
													mult(
														sub(identity(),tensorProduct(pData[setName_i]->normalVec[i],pData[setName_i]->normalVec[i])),
														pData[setName_i]->particleDensityGrad[i]
												    )
												 );
			}else{
				pData[setName_i]->shift[i] = mult( -1.666666 * smoothingLength * pData[setName_i]->vol[i] * dt, pData[setName_i]->particleDensityGrad[i]);
			}
						
			if (!isBoundary_i){
				// pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i],Real3{0.0,-0.2,0.0});
				//prescribed boundary
				// if(pData[setName_i]->originPos[i][0] >= 1.9 && pData[setName_i]->originPos[i][1] >= 0.474 ){
				// 	pData[setName_i]->acc[i][1] = (-10.0 * std::min((Real)t/(100.0),(Real)1.0));
				// 	// pData[setName_i]->pos[i][1] = pData[setName_i]->originPos[i][1] - std::min(0.1 * t / (200.0)  ,0.1);
				// }
			} else{
				pData[setName_i]->acc[i] = zerovec;
				pData[setName_i]->vel[i] = zerovec;
			}
			
		}

		// std::cout << "foo" << std::endl;
	
	}
			// std::cout << "Breakpoint 2 " << std::endl;


}



void SPHSolver::smearStress(Uint t){
	using namespace RealOps;

	currentTime += (Real) simData["dt"];
	const Real dt = (Real) simData["dt"];
	const Real dx = (Real) simData["dx"];
	const Real smoothingLength = (Real) simData["smoothingLength"];

	std::cout << "   |----------------- Smearing PK2 with diffusion" << std::endl;

	for (const auto& setName : setNames){

		const int setID = ids[setName];
		const auto& ps = nsearch->point_set(setID);
		#pragma omp parallel for num_threads(NUMTHREADS) 
		
		for (int i = 0; i < ps.n_points(); ++i){	

			Real3x3 stressCorrection{zeromat};	
			pData[setName]->particleDensity[i] = 0.0;

			for (const auto& setName_j : setNames){

				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);
				const bool isBoundary_j = (setName_j == "boundary") ? true : false;

				// Loop over the neighbors of particl i in the indexed particle set.
				for (int _j = 0; _j < ps.n_neighbors(setID_j,i); _j++){

					// Index of neighbor particle in setID_j
					Uint const j = ps.neighbor(setID_j, i, _j);
					// Define all the ingedients for the particle interaction
					Real3  relpos = sub(pData[setName]->pos[i],pData[setName_j]->pos[j]);
					Real     dist = length(relpos);
					Real3  reldir = divide(relpos,dist);
					Real      Wij =  W_ij(dist, smoothingLength);
					Real3    gWij = gW_ij(dist, reldir, smoothingLength);
					Real3   v_ij = sub(pData[setName]->vel[i],pData[setName_j]->vel[j]);
					Real3x3 S_ij = sub(pData[setName]->secondPKStress[i],pData[setName_j]->secondPKStress[j]);
					
					// Computethe particle density
					pData[setName]->particleDensity[i] += pData[setName_j]->vol[j] * Wij;

					stressCorrection = add(stressCorrection, 
												mult(2.0 * dot(reldir,gWij) * pData[setName_j]->vol[j] / dist, S_ij)
											);

											
										
				}
			}		

			// Compute the particle density (add contribution from self)
			pData[setName]->particleDensity[i] += pData[setName]->vol[i] * W_ij(0, smoothingLength);
			stressCorrection = mult( (Real) simData["smearing_S"] * dt / pData[setName]->particleDensity[i] , stressCorrection);
			pData[setName]->secondPKStress[i] = add(pData[setName]->secondPKStress[i], stressCorrection);

		}
	}
	
}



void SPHSolver::XSPH(Uint t){
	using namespace RealOps;

	currentTime += (Real) simData["dt"];
	const Real dt = (Real) simData["dt"];
	const Real dx = (Real) simData["dx"];
	const Real smoothingLength = (Real) simData["smoothingLength"];

	std::cout << "   |----------------- Smoothing with XSPH" << std::endl;
	for (const auto& setName : setNames){
		const int setID = ids[setName];
		const auto& ps = nsearch->point_set(setID);
		#pragma omp parallel for num_threads(NUMTHREADS) 
		for (int i = 0; i < ps.n_points(); ++i){	
			// Real3x3 defoGradCorrection{zeromat};		
			Real3 velocityCorrection{zerovec};		
			pData[setName]->particleDensity[i] = 0.0;
			for (const auto& setName_j : setNames){

				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);
				const bool isBoundary_j = (setName_j == "boundary") ? true : false;

				// Compute the particle density while we are at it				

				// Loop over the neighbors of particl i in the indexed particle set.
				for (int _j = 0; _j < ps.n_neighbors(setID_j,i); _j++){

					// Index of neighbor particle in setID_j
					Uint const j = ps.neighbor(setID_j, i, _j);

					if (!pData[setName]->isSolid[i] && !pData[setName_j]->isSolid[j]) continue;

					// Define all the ingedients for the particle interaction
					Real3  relpos = sub(pData[setName]->pos[i],pData[setName_j]->pos[j]);
					Real3 relposj = mult(-1.0,relpos);

					Real     dist = length(relpos);
					Real      Wij = W_ij(dist, smoothingLength);

					// Real3x3 F_ji = sub(pData[setName_j]->defoGrad[j],pData[setName]->defoGrad[i]);
					Real3 v_ji = sub(pData[setName_j]->vel[j],pData[setName]->vel[i]);
					Real oneOverDensAve = 0.5 * (1./pData[setName_j]->dens[j] + 1./pData[setName]->dens[i]);
					// defoGradCorrection = add(defoGradCorrection, 
					// 						mult( Wij * (Real) simData["XSPH"] * pData[setName_j]->mass[j] * oneOverDensAve,
					// 							  F_ji)
					// 						);
					velocityCorrection = add(velocityCorrection, 
											mult( Wij * (Real) simData["XSPH"] * pData[setName_j]->mass[j] * oneOverDensAve,
												  v_ji)
											);

					// Compute the particle density
					pData[setName]->particleDensity[i] += pData[setName_j]->vol[j] * Wij;

				}
			}			

			// Compute the particle density (add contribution from self)
			pData[setName]->particleDensity[i] += pData[setName]->vol[i] * W_ij(0, smoothingLength);
			// defoGradCorrection = mult(1./pData[setName]->particleDensity[i], defoGradCorrection);
			velocityCorrection = mult(1./pData[setName]->particleDensity[i], velocityCorrection);
			// pData[setName]->defoGrad[i] = add(pData[setName]->defoGrad[i],defoGradCorrection);
			pData[setName]->vel[i] = add(pData[setName]->vel[i],velocityCorrection);
			
		}

	}
}

void SPHSolver::fixedPointIteration(Uint t){
	using namespace RealOps;

	currentTime += (Real) simData["dt"];
	const Real dt = (Real) simData["dt"];
	const Real dx = (Real) simData["dx"];
	
	std::cout << "   |----------------- Saving x(t)" << std::endl;
	for (const auto& setName : setNames){
		const int setID = ids[setName];
		const auto& ps = nsearch->point_set(setID);
		#pragma omp parallel for num_threads(NUMTHREADS) 
		for (int i = 0; i < ps.n_points(); ++i){

			pData[setName]->posBefore[i] = pData[setName]->pos[i];
			pData[setName]->velBefore[i] = pData[setName]->vel[i];			
			pData[setName]->accBefore[i] = pData[setName]->acc[i];			

			pData[setName]->stretch_plastic_dot_before[i] = pData[setName]->stretch_plastic_dot[i];
			pData[setName]->stretch_plastic_before[i]     = pData[setName]->stretch_plastic[i];
			pData[setName]->stretch_dirs_before[i]        = pData[setName]->stretch_dirs[i];

		}
	}

	for(Uint n = 0; n < (Uint) simData["fixedPointIterations"]; n ++){
		std::cout << "-------------------- Performing Fixed Point Iteration ... iteration " << n << std::endl;
		Real norm = 0;
		if( simData["computeElasticity"] == "Yes" ){
			computeDeformationGradient(t);
			// smearDefGrad(t);
		} 													
		std::cout << "	|--- Calculating Interatctions" << std::endl;
		computeInteractions(t);
		for (const auto& setName : setNames){

			const int setID = ids[setName];
			const auto& ps = nsearch->point_set(setID);

			if( setName == "fluid"){
				// #pragma omp parallel for num_threads(NUMTHREADS) 
				for (int i = 0; i < ps.n_points(); ++i){

					pData[setName]->pos[i] = add(
												pData[setName]->posBefore[i],
												 add(mult(dt, pData[setName]->velBefore[i]), 
													 mult(dt * dt / 4.0,add(pData[setName]->acc[i], pData[setName]->accBefore[i]))
													 )
												);	


					pData[setName]->vel[i] = add(
												pData[setName]->velBefore[i],
												mult(dt / 2.0, add(pData[setName]->acc[i], pData[setName]->accBefore[i]))
											);

					pData[setName]->stretch_plastic[i] = add(
												pData[setName]->stretch_plastic_before[i],
												mult(dt / 2.0, add(pData[setName]->stretch_plastic_dot[i], pData[setName]->stretch_plastic_dot_before[i]))
												);

					norm += length(pData["fluid"]->psi[i],pData[setName]->pos[i]);
					pData[setName]->psi[i] = pData[setName]->pos[i];
				}
			} else{
				// Do stuff for the boundary particles??
			}

		}
		
		std::cout << "dumb norm : " << norm << std::endl;
		// if(norm < 1.0E-9){
			// std::cout << "converged?" << std::endl;
			// break;
		// }
	}
	// XSPH(t);


	// for (const auto& setName : setNames){

	// 	const int setID = ids[setName];
	// 	const auto& ps = nsearch->point_set(setID);

	// 	if( setName == "fluid"){
	// 		#pragma omp parallel for num_threads(NUMTHREADS) 
	// 		for (int i = 0; i < ps.n_points(); ++i){

	// 			pData[setName]->vel[i] = add(
	// 									pData[setName]->vel[i],
	// 									mult(dt / 2.0, add(pData[setName]->acc[i], pData[setName]->accBefore[i]))
	// 									);

	// 		}
	// 	}
	// }

	
	
}




void SPHSolver::marchTime(Uint t){
	using namespace RealOps;

	std::cout << "-------------------- Timestep #" << t << std::endl;

	fixedPointIteration(t);
	// std::cout << "	|--- Updating Position " << std::endl;		

	// for (const auto& setName : setNames){

	// 	const int setID = ids[setName];
	// 	const auto& ps = nsearch->point_set(setID);


	// 	if( simData["computeElasticity"] == "Yes" ) computeDeformationGradient(t);
	// 	std::cout << "	|--- Calculating Interatctions" << std::endl;
	// 	computeInteractions(t);
	// 	std::cout << "	|--- Updating Vectors" << std::endl;

		

	// 	if( setName == "fluid"){
	// 	// Fluid Particles    : march the position, velocity, density.
	// 		#pragma omp parallel for num_threads(NUMTHREADS) 
	// 		for (int i = 0; i < ps.n_points(); ++i){
	// 			pData[setName]->vel[i]  = add(pData[setName]->vel[i], mult(dt,pData[setName]->acc[i]));
	// 			pData[setName]->pos[i]  = add(pData[setName]->pos[i], mult(dt,pData[setName]->vel[i]));

	// 			// ######################################
	// 			// Test For Deformation Gradient
	// 			// pData[setName]->pos[i] = add(pData[setName]-> pos[i] , Real3{(0.0 * pData[setName]->originPos[i][0] + 0.0  * pData[setName]->originPos[i][1]),
	// 																		//  (dt * 0.05 * pData[setName]->originPos[i][0] + 0.0  * pData[setName]->originPos[i][1]),0} );
	// 			// pData[setName]->pos[i] = add(pData[setName]-> pos[i] , Real3{dt * (0.5 * pData[setName]->originPos[i][0] + 0.0  * pData[setName]->originPos[i][1]),
	// 																		//  dt * (0.5 * pData[setName]->originPos[i][0] + 0.25 * pData[setName]->originPos[i][1]),0} );
	// 			// pData[setName]->pos[i] = Real3{(std::cos(currentTime) * pData[setName]->originPos[i][0] - std::sin(currentTime) * pData[setName]->originPos[i][1]),
	// 										//    (std::sin(currentTime) * pData[setName]->originPos[i][0] + std::cos(currentTime) * pData[setName]->originPos[i][1]),0};
	// 			// ######################################

	// 			// Density updates for WCSPH should only be considered for non-solid particles
	// 			if(!pData[setName]->isSolid[i]) pData[setName]->dens[i] = pData[setName]->dens[i] + dt * pData[setName]->densdot[i];
	// 			pData[setName]->temp[i] = pData[setName]->temp[i] + dt * pData[setName]->enthalpydot[i];					
	// 			pData[setName]->tau[i] = add(pData[setName]->tau[i],mult(dt, pData[setName]->tauDot[i]));				
	// 		}
	// 	} else if (setName == "boundary"){
	// 	// Boundary Particles : march the density only.
	// 		#pragma omp parallel for num_threads(NUMTHREADS)
	// 		for (int i = 0; i < ps.n_points(); ++i){
	// 			pData[setName]->dens[i] = pData[setName]->dens[i]  + dt * pData[setName]->densdot[i];
	// 			pData[setName]->temp[i] = pData[setName]->temp[i] + dt * pData[setName]->enthalpydot[i] / (pData[setName]->getSpecificHeat());
	// 		}
	// 	}
	// }


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
		