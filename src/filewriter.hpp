#include <H5Part.h>
#include <iostream>
#include "helpers.hpp"
#include "particleAttrib.hpp"
#include "json.hpp"

using namespace RealOps;

class FileWriter{
	H5PartFile *fileWriter;
	bool isOpen;
	unsigned int fileNumber;
	unsigned int timeStep;
	std::string outputSuffix;
	std::string currentFileName;
	std::vector<Real> mass,x,y,z,vx,vy,vz,dens,temp,heatSensed,fxSensed,fySensed,fzSensed,isSensor;

public:
	FileWriter(std::string outputFile){
		isOpen = false;
		outputFile = outputFile.substr(1, outputFile.size()); outputFile = outputFile.substr(0, outputFile.size()-1);
		outputSuffix = outputFile;
		fileNumber = 0;
		timeStep = 0;
		currentFileName = 	"../outputs/" + outputSuffix + "_" + std::to_string(fileNumber) + ".h5part";
		fileWriter = H5PartOpenFile(currentFileName.c_str(),H5PART_WRITE);
	}
	void openNextFile(){
		H5PartCloseFile(fileWriter);
		fileNumber += 1;
		currentFileName = 	"../outputs/" + outputSuffix + "_" + std::to_string(fileNumber) + ".h5part";
		fileWriter = H5PartOpenFile(currentFileName.c_str(),H5PART_WRITE);
		timeStep = 0;
	}
	void close(){
		H5PartCloseFile(fileWriter);
	}
	void write(std::map<std::string,ParticleAttributes*> &pData, Real smoothingLength){
		int totalParticles = 0;
		std::vector<Real> h,x,y,z,vx,vy,vz,dens,temp,heatSensed,fxSensed,fySensed,fzSensed,isSensor;
		std::vector<Real> nx,ny,nz,isFS,tgx,tgy,tgz,curvature,hdot;

		for (const auto& dataType : pData) totalParticles += dataType.second->numParticles;
		std::cout << "... Writing Results to " << currentFileName << std::endl;
		std::cout << "... Total Particles to Write : " << totalParticles << std::endl;

		x.resize(totalParticles); vx.resize(totalParticles);
		y.resize(totalParticles); vy.resize(totalParticles);
		z.resize(totalParticles); vz.resize(totalParticles);

		isFS.resize(totalParticles);
		curvature.resize(totalParticles);
		nx.resize(totalParticles);
		ny.resize(totalParticles);
		nz.resize(totalParticles);

		tgx.resize(totalParticles);
		tgy.resize(totalParticles);
		tgz.resize(totalParticles);

		hdot.resize(totalParticles);

		mass.resize(totalParticles);
		h.resize(totalParticles);
		dens.resize(totalParticles); temp.resize(totalParticles);
		std::cout << "foo" << std::endl;
		isSensor.resize(totalParticles);
		heatSensed.resize(totalParticles);
		fxSensed.resize(totalParticles);
		fySensed.resize(totalParticles);
		fzSensed.resize(totalParticles);
		int idx = 0;

		for (const auto& dataType : pData){

			#pragma omp parallel for num_threads(NUMTHREADS)
			for (int i = idx; i < idx + dataType.second->numParticles; i++ ){

				x[i] = dataType.second->pos[i - idx][0]; vx[i] = dataType.second->vel[i - idx][0];
				y[i] = dataType.second->pos[i - idx][1]; vy[i] = dataType.second->vel[i - idx][1];
				z[i] = dataType.second->pos[i - idx][2]; vz[i] = dataType.second->vel[i - idx][2];

				nx[i] = dataType.second->normalVec[i - idx][0];
				ny[i] = dataType.second->normalVec[i - idx][1];
				nz[i] = dataType.second->normalVec[i - idx][2];

				tgx[i] = dataType.second->tempGrad[i - idx][0];
				tgy[i] = dataType.second->tempGrad[i - idx][1];
				tgz[i] = dataType.second->tempGrad[i - idx][2];

				hdot[i] = dataType.second->enthalpydot[i - idx];

				isFS[i] = dataType.second->isFS[i - idx];
				curvature[i] = dataType.second->curvature[i - idx];
				dens[i] = dataType.second->dens[i - idx];
				mass[i] = dataType.second->mass[i - idx];
				h[i] = smoothingLength;
				// dens[i] = dataType.second->tempGrad[i-idx][0];
				// dens[i] = dataType.second->densdot[i-idx];
				temp[i] = dataType.second->temp[i - idx];
				isSensor[i] = (dataType.second->isSensor[i-idx] ? 1 : 0);

				//Sensor values. These are accumulated valued over the output interval.
				heatSensed[i] = dataType.second->heatSensed[i - idx];
				fxSensed[i] = dataType.second->forceSensed[i - idx][0];
				fySensed[i] = dataType.second->forceSensed[i - idx][1];
				fzSensed[i] = dataType.second->forceSensed[i - idx][2];

				//Wipe the sensor values.
				dataType.second->heatSensed[i - idx] = 0;
				dataType.second->forceSensed[i - idx] = Real3{0,0,0};
			}

			idx += dataType.second->numParticles;
		}
		std::cout << "Copied into memory. Writing with H5Part ..." << std::endl;

		H5PartSetStep(fileWriter,timeStep);
		H5PartSetNumParticles(fileWriter,totalParticles);

		H5PartWriteDataFloat64(fileWriter,"x",&x[0]);
		H5PartWriteDataFloat64(fileWriter,"y",&y[0]);
		H5PartWriteDataFloat64(fileWriter,"z",&z[0]);
		H5PartWriteDataFloat64(fileWriter,"vx",&vx[0]);
		H5PartWriteDataFloat64(fileWriter,"vy",&vy[0]);
		H5PartWriteDataFloat64(fileWriter,"vz",&vz[0]);

		H5PartWriteDataFloat64(fileWriter,"nx",&nx[0]);
		H5PartWriteDataFloat64(fileWriter,"ny",&ny[0]);
		H5PartWriteDataFloat64(fileWriter,"nz",&nz[0]);
		H5PartWriteDataFloat64(fileWriter,"tgx",&tgx[0]);
		H5PartWriteDataFloat64(fileWriter,"tgy",&tgy[0]);
		H5PartWriteDataFloat64(fileWriter,"tgz",&tgz[0]);
		H5PartWriteDataFloat64(fileWriter,"hdot",&hdot[0]);
		H5PartWriteDataFloat64(fileWriter,"freeSurface",&isFS[0]);
		H5PartWriteDataFloat64(fileWriter,"curvature",&curvature[0]);

		H5PartWriteDataFloat64(fileWriter,"dens",&dens[0]);
		H5PartWriteDataFloat64(fileWriter,"mass",&mass[0]);
		H5PartWriteDataFloat64(fileWriter,"h",&h[0]);
		H5PartWriteDataFloat64(fileWriter,"temp",&temp[0]);
		H5PartWriteDataFloat64(fileWriter,"heatSensed",&heatSensed[0]);
		H5PartWriteDataFloat64(fileWriter,"fxSensed",&fxSensed[0]);
		H5PartWriteDataFloat64(fileWriter,"fySensed",&fySensed[0]);
		H5PartWriteDataFloat64(fileWriter,"fzSensed",&fzSensed[0]);

		std::cout << "Done." << std::endl;
		timeStep ++;
		}

};
