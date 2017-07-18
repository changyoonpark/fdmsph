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
	std::vector<Real> x,y,z,vx,vy,vz,dens,temp,heatSensed,fxSensed,fySensed,fzSensed,isSensor; 
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
	void write(std::map<std::string,ParticleAttributes*> &pData){
		int totalParticles = 0;
		std::vector<Real> x,y,z,vx,vy,vz,dens,temp,heatSensed,fxSensed,fySensed,fzSensed,isSensor;
		for (const auto& dataType : pData) totalParticles += dataType.second->numParticles;
		std::cout << "... Writing Results to " << currentFileName << std::endl;
		x.resize(totalParticles); vx.resize(totalParticles);
		y.resize(totalParticles); vy.resize(totalParticles);
		z.resize(totalParticles); vz.resize(totalParticles);
		dens.resize(totalParticles); temp.resize(totalParticles);
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
				// dens[i] = dataType.second->dens[i - idx];
				// dens[i] = dataType.second->tempGrad[i-idx][0];
				dens[i] = dataType.second->densdot[i-idx];
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

		H5PartSetStep(fileWriter,timeStep);
		H5PartSetNumParticles(fileWriter,totalParticles);

		H5PartWriteDataFloat64(fileWriter,"x",&x[0]);
		H5PartWriteDataFloat64(fileWriter,"y",&y[0]);
		H5PartWriteDataFloat64(fileWriter,"z",&z[0]);
		H5PartWriteDataFloat64(fileWriter,"vx",&vx[0]);
		H5PartWriteDataFloat64(fileWriter,"vy",&vy[0]);
		H5PartWriteDataFloat64(fileWriter,"vz",&vz[0]);
		H5PartWriteDataFloat64(fileWriter,"dens",&dens[0]);
		H5PartWriteDataFloat64(fileWriter,"temp",&temp[0]);
		H5PartWriteDataFloat64(fileWriter,"heatSensed",&heatSensed[0]);
		H5PartWriteDataFloat64(fileWriter,"fxSensed",&fxSensed[0]);
		H5PartWriteDataFloat64(fileWriter,"fySensed",&fySensed[0]);
		H5PartWriteDataFloat64(fileWriter,"fzSensed",&fzSensed[0]);
				
		timeStep ++;
		}

};

