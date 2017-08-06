#ifndef __HELPERS__
#define __HELPERS__

#include <iostream>
#include <fstream>
#include <CompactNSearch>
#include <Eigen/Dense>
#include "json.hpp"

#define NUMTHREADS 36
#define MAX_NEIGHBORS 200
#define EPSL_SMALL 1E-10
// using namespace CompactNSearch;
using json = nlohmann::json;
// using Real = double;

inline void readJSON(json& jsonobj, std::string fname){
	std::ifstream input(fname);
	input >> jsonobj;
}

namespace RealOps{
	using Uint = unsigned int;
	using namespace CompactNSearch;
	using namespace Eigen;

	using Real3 = std::array<Real,3>;
	using Real3x3 = std::array<Real3,3>;

	const Uint IDXPAIR[6][2] = {{0,0},{1,1},{2,2},{0,1},{1,2},{0,2}};
	// const Real NEG_DELTA_MN[6]   = {-1.0,-1.0,-1.0,0,0,0};
	// Eigen::Matrix3d A((Eigen::Matrix3d() << 1, 2, 3, 4, 5, 6, 7, 8, 9).finished());

	inline Real3x3 toReal3x3 (VectorXd vec){
		return Real3x3{Real3{vec(0),vec(3),vec(5)},
					   Real3{vec(3),vec(1),vec(4)},
					   Real3{vec(5),vec(4),vec(2)}};
	}

	inline Real maxEig( VectorXd& a){
		Real result = a[0];
		if (a[1] > result)
		  result = a[1];

		if (a[2] > result)
		  result = a[2];

		return result;
	}

	inline void assign( MatrixXd& a, const Real3x3& b ){
		for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			a(i,j)=b[i][j];
	}

	inline std::array<Real,3> add( const std::array<Real,3>& a,  const std::array<Real,3>& b) {
	    return std::array<Real,3>{a[0]+b[0],a[1]+b[1],a[2]+b[2]};
	}

	inline std::array<Real,3> sub( const std::array<Real,3>& a,  const std::array<Real,3>& b) {
	    return std::array<Real,3>{a[0]-b[0],a[1]-b[1],a[2]-b[2]};
	}

	inline std::array<Real,3> mult( Real c, const std::array<Real,3>& a) {
	    return std::array<Real,3>{c * a[0], c * a[1], c * a[2]};
	}

	inline std::array<Real,3> div(const std::array<Real,3>& a, Real c) {
	    return std::array<Real,3>{a[0] / c, a[1] / c, a[2] / c};
	}

	inline std::array<Real,3> neg(const std::array<Real,3>& a) {
	    return std::array<Real,3>{ -a[0], -a[1], -a[2]};
	}

	inline std::array<Real,3> dir( const std::array<Real,3>& a ) {
	    return div(a,a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	}

	inline std::array<Real,3> cross( const std::array<Real,3>& a,  const std::array<Real,3>& b) {
		return std::array<Real,3>{ a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
	}

	inline Real dot( const std::array<Real,3>& a,  const std::array<Real,3>& b) {
	    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
	}

	inline Real length( const std::array<Real,3>& a ) {
	    return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	}

	inline Real length2( const std::array<Real,3>& a ) {
	    return (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	}

	inline bool iseq( const std::array<Real,3>& a, const std::array<Real,3>& b){
		return (a[0]==b[0] && a[1]==b[1] && a[2]==b[2]) ? true : false;
	}
	inline void print( const Real3x3& a){
		std::cout << " Real3x3 --------------------------------------------------" << std::endl;
		std::cout << " Real3 : { " << a[0][0] << ", " << a[0][1] << ", " << a[0][2] << " }" << std::endl;
		std::cout << " Real3 : { " << a[1][0] << ", " << a[1][1] << ", " << a[1][2] << " }" << std::endl;
		std::cout << " Real3 : { " << a[2][0] << ", " << a[2][1] << ", " << a[2][2] << " }" << std::endl;

	}

	inline void print(  const std::string& s , const std::array<Real,3>& a ){
		std::cout << s << " Real3 : { " << a[0] << ", " << a[1] << ", " << a[2] << " }" << std::endl;
 	}

	inline void print(  const std::string& s , const Real& a){
		std::cout << s << " Real : { " << a << " }"  << std::endl;
 	}


 	inline Real3x3 mult(Real a, Real3x3 b){
 		return Real3x3{mult(a,b[0]),mult(a,b[1]),mult(a,b[2])};
 	}

	inline Real3x3 add( const Real3x3& a,  const Real3x3& b) {
	    return Real3x3 { add(a[0],b[0]), add(a[1],b[1]), add(a[2],b[2])};
	}

	inline Real doubleDot(const Real3x3& a, const Real3x3& b){
		Real bar=0;
		for(int i=0;i<3;i++)for(int j=0;j<3;j++)
			bar+=(a[i][j]*b[i][j]);
		return bar;
	}
	// inline Real3x3 add( const Real3x3& a,  const Real3x3& b) {
	//     return Real3x3 { Real3{a[0][0] + b[0][0],a[0][1] + b[0][1],a[0][2] + b[0][2]},
	// 				     Real3{a[1][0] + b[1][0],a[1][1] + b[1][1],a[1][2] + b[1][2]},
	// 				     Real3{a[2][0] + b[2][0],a[2][1] + b[2][1],a[2][2] + b[2][2]} };
	// }

	inline Real3x3 inv( const Real3x3& a){
		Real3x3 b;
		Real a11 = a[0][0], a12 = a[0][1], a13 = a[0][2];
		Real a21 = a[1][0], a22 = a[1][1], a23 = a[1][2];
		Real a31 = a[2][0], a32 = a[2][1], a33 = a[2][2];
		Real d = a11 * (a33 * a22 - a32 * a23) - a21 * (a33 * a12 - a32 * a13) + a31 * (a23 * a12 - a22 * a13);
		b[0][0] =  (a33 * a22 - a32 * a23)/d,  b[0][1] = - (a33 * a12 - a32 * a13)/d, b[0][2] =  (a23 * a12 - a22 * a13)/d;
		b[1][0] = -(a33 * a21 - a31 * a23)/d,  b[1][1] =   (a33 * a11 - a31 * a13)/d, b[1][2] = -(a23 * a11 - a21 * a13)/d;
		b[2][0] =  (a32 * a21 - a31 * a22)/d,  b[2][1] = - (a32 * a11 - a31 * a12)/d, b[2][2] =  (a22 * a11 - a21 * a12)/d;
		return b;
	}

 	inline Real3x3 tensorProduct(Real3& a, Real3& b){
 		Real3x3 c;

 		c[0][0] = a[0]*b[0], c[0][1] = a[0]*b[1], c[0][2] = a[0]*b[2];
 		c[1][0] = a[1]*b[0], c[1][1] = a[1]*b[1], c[1][2] = a[1]*b[2];
 		c[2][0] = a[2]*b[0], c[2][1] = a[2]*b[1], c[2][2] = a[2]*b[2];

		// if (c[0][0] <= 1.0E-10 && c[0][0] >= -1.0E-10 ) c[0][0] = 1.0;
		// if (c[1][1] <= 1.0E-10 && c[1][1] >= -1.0E-10 ) c[1][1] = 1.0;
		// if (c[2][2] <= 1.0E-10 && c[2][2] >= -1.0E-10 ) c[2][2] = 1.0;

 		return c;
 	}

	inline void checkSingularity(Real3x3& L_i){

		if (abs(L_i[0][0]) < 1.0E-10 &&
				abs(L_i[1][0]) < 1.0E-10 && abs(L_i[0][1]) < 1.0E-10 &&
				abs(L_i[1][2]) < 1.0E-10 && abs(L_i[2][1]) < 1.0E-10){
			L_i[0][0] = 1.0;
		}
		if (abs(L_i[1][1]) < 1.0E-10 &&
				abs(L_i[1][2]) < 1.0E-10 && abs(L_i[2][1]) < 1.0E-10 &&
				abs(L_i[0][1]) < 1.0E-10 && abs(L_i[1][0]) < 1.0E-10){
			L_i[1][1] = 1.0;
		}
		if (abs(L_i[2][2]) < 1.0E-10 &&
				abs(L_i[2][0]) < 1.0E-10 && abs(L_i[2][1]) < 1.0E-10 &&
				abs(L_i[0][2]) < 1.0E-10 && abs(L_i[1][2]) < 1.0E-10){
			L_i[2][2] = 1.0;
		}

	}

 	inline Real3 mult(Real3x3& a, Real3& b){
 		return Real3{a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2],
 					 a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2],
 					 a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2]};
 	}

 	inline Real W_Wendland_2D_2h(Real dist, Real smoothingLength){
 		Real h = smoothingLength * 0.5;
 		Real q = (dist / h);
 		if (q > 2.0){
 			return 0.0;
 		} else{
 			Real foo = (1.0 - 0.5 * q);
			Real bar = foo * foo;
 			return (0.55704230082163367519109317180380 / (h * h)) * bar * bar * (2.0 * q + 1.0);
		}
 	}

 	inline Real W_Wendland_3D_2h(Real dist, Real smoothingLength){
 		Real h = smoothingLength * 0.5;
 		Real q = (dist / h);
 		if (q > 2.0){
 			return 0.0;
 		} else{
 			Real foo = (1.0 - 0.5 * q);
			Real bar = foo * foo;
 			return (0.41778172561622525639 / (h * h)) * bar * bar * (2.0 * q + 1.0);
		}
 	}


 	inline Real W_Wendland_2D(Real dist, Real smoothingLength){
 		Real q = dist / smoothingLength;
 		if (q > 1.0){
 			return 0.0;
 		} else{
 			Real foo = (1.0 - q);
			Real bar = foo * foo;
			// return ((4.456338406572776) / (smoothingLength*smoothingLength)) * 0.5 * ((1.0-q)*(1.0-q)*(1.0-q)*(1.0-q)) * (1.0 + 4.0 * q);
	 		return (2.228169203286388) / ( smoothingLength * smoothingLength ) * bar * bar * (1.0 + 4.0 * q);
		}
 	}

 	inline Real3 gW_Wendland_3D_2h(Real dist, std::array<Real,3> dir, Real smoothingLength){
 		Real h = smoothingLength * 0.5;
 		Real q = (dist / h);
 		if (q > 2.0){
 			return std::array<Real,3>{0,0,0};
 		} else{
 			Real foo = 1.0 - 0.5 * q;
 			// This is the "2h" version.
 			// (Referred from the DualSphysics Guidbook, \frac{alpha_D}{h} * (-5 * q) * (1-\frac{1}{2}q)^3)
	 		return mult(((0.41778172561622525639)/(h*h*h))*(-5.0 * q)*foo*foo*foo, dir);
		}
 	}

 	inline Real3 gW_Wendland_2D_2h(Real dist, std::array<Real,3> dir, Real smoothingLength){
 		Real h = smoothingLength * 0.5;
 		Real q = (dist / h);
 		if (q > 2.0){
 			return std::array<Real,3>{0,0,0};
 		} else{
 			Real foo = 1.0 - 0.5 * q;
 			// This is the "2h" version.
 			// (Referred from the DualSphysics Guidbook, \frac{alpha_D}{h} * (-5 * q) * (1-\frac{1}{2}q)^3)
	 		return mult(((0.557042300821633675191)/(h*h*h))*(-5.0 * q)*foo*foo*foo, dir);
		}
 	}


 	inline Real3 gW_Wendland_2D(Real dist, std::array<Real,3> dir, Real smoothingLength){
 		Real q = dist / smoothingLength;
 		if (q > 1.0){
 			return std::array<Real,3>{0,0,0};
 		} else{
	 		return mult( ( (4.456338406572776) / ( smoothingLength * smoothingLength * smoothingLength ) ) * (-10.0 * q) * (1.0 - q) * (1.0 - q) * (1.0 - q), dir);
		}
 	}


 	inline Real delta_SPH(Real  rho_i,      Real rho_j,
 						  Real  vol_j,        Real delta,
 						  Real  soundSpeed, Real smoothingLength,
 						  Real  dist,
 						  Real3 relpos,    	Real3 gWij,
						  Real3 densGrad_i, Real3 densGrad_j){
 		const Real psi = 2.0 * ( - rho_j + rho_i) - dot(add(densGrad_j,densGrad_i),relpos);
 		const Real qoo = dot(relpos, gWij) * vol_j * soundSpeed * smoothingLength * delta / (dist * dist);
 		return psi * qoo;
 	}

 	// inline Real temperatureDependentViscosity(Real mu, Real T){return 0;}
 	inline Real equal (Real T){
 		return T;
 	}

 	inline Real fixedViscosity(Real mu, Real T){return mu;}

 	inline Real linearEOS(Real T, Real T0, Real rho, Real rho0, Real c, Real alpha){
 		return c * c * ( (rho - rho0) + alpha * (T - T0) );
 	}

 	inline Real taitEOS(Real T, Real T0, Real rho, Real rho0, Real c, Real alpha){
 		return (1./7.) * rho0 * c * c * (std::pow((rho/rho0),7) - 1.0);
 	}

 	inline Real3 viscosity_acc_ij_Shao(Real mj, Real3 relpos, Real3 gWij, Real rho_i, Real rho_j, Real dist, Real3 relvel, Real mu){

 		return mult(( 4. * mj * mu * dot(relpos, gWij) / ((rho_i + rho_j) * (dist * dist + EPSL_SMALL )) ), relvel );
 	}

 	inline Real3 pressure_acc_ij_1(Real mj, Real3 gWij, Real P_i, Real P_j, Real rho_i, Real rho_j, Real vol_j){
 		return mult( - mj * (P_i + P_j) / (rho_i * rho_j), gWij);
 	}

 	inline Real3 pressure_acc_ij_2(Real mj, Real3 gWij, Real P_i, Real P_j, Real rho_i, Real rho_j, Real vol_j){
 		return mult( - (1./rho_i) * (P_i + P_j) * vol_j, gWij);
 	}

 	inline Real3 pressure_acc_ij_3(Real mj, Real3 gWij, Real P_i, Real P_j, Real rho_i, Real rho_j, Real vol_j){
 		return mult( - mj * ( P_i / (rho_i * rho_i) + P_j / (rho_j * rho_j) ), gWij);
 	}

 	inline Real constantConductivity(std::string type, Real T){
 		return 1.0;
 	}

 	// inline Real cleary(Real Ti, Real Tj, Real mj, Real rho_i, Real rho_j, Real ki, Real kj, Real3 relpos, Real dist, Real3 gWij){
 	// 	return 4.0 * mj * ki * kj * (Ti - Tj) * dot(relpos, gWij) / (rho_i * rho_j * (ki + kj) * (dist * dist));
 	// }
	inline Real consistentHeatTransfer(Real3x3 L2_i, Real3 gradT_i,
						   Real Ti, Real Tj, Real mj, Real rho_i, Real rho_j, Real ki, Real kj,
						   Real3 relpos, Real3 reldir,
						   Real dist, Real3 gWij, Real vol_j){
		Real3x3 eij_gWij = tensorProduct(reldir,gWij);

		// return 4.0 * mj * ki * kj * (Ti - Tj) * dot(relpos, gWij) / (rho_i * rho_j * (ki + kj) * (dist * dist));

		// return 4.0 * (ki * kj / (ki + kj)) / (rho_i) * doubleDot(L2_i,eij_gWij) * ((Ti - Tj) / (dist) - (ki * kj / (ki + kj)) * dot(reldir,gradT_i)) * vol_j;
		return (4.0 / (rho_i)) * doubleDot(L2_i,eij_gWij) * (ki * kj / (ki + kj)) * ((Ti - Tj) / (dist) - dot(reldir,gradT_i)) * vol_j;
		return (2.0 / (rho_i)) * doubleDot(L2_i,eij_gWij) * ((Ti - Tj) / (dist) - dot(reldir,gradT_i)) * vol_j;
		// return 2.0 * ((Ti-Tj) / dist) * dot(reldir,gWij) * vol_j;
	}

 	inline Real3 gravity_acc(){
 		return Real3{0,0,-10.0};
 	}

 	struct SymTensor3{
 		Real data[3][3][3];
 		SymTensor3(){
 			for(int i=0;i<3;i++)
 				for(int j=0;j<3;j++)
 					for(int k=0;k<3;k++)
 						data[i][j][k] = 0;
 		}
 		void add(Uint k, Uint m, Uint n, Real value){
 			data[k][m][n] += value;
 		}
 		void mult(Uint k, Uint m, Uint n, Real value){
 			data[k][m][n] = data[k][m][n] * value;
 		}
 		inline Real operator()(Uint k, Uint m, Uint n){
 			return data[k][m][n];
 		}
 	};


}

#endif
