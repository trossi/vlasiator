/*
 This file is part of Vlasiator.
 
 Copyright 2011 Finnish Meteorological Institute
 
 Vlasiator is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License version 3
 as published by the Free Software Foundation.
 
 Vlasiator is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "cell_spatial.h"
#include "common.h"
#include "project.h"
#include "parameters.h"

typedef dispersionParameters DispP;
Real DispP::BX0 = NAN;
Real DispP::BY0 = NAN;
Real DispP::BZ0 = NAN;
Real DispP::DENSITY = NAN;
Real DispP::TEMPERATURE = NAN;
//Real DispP::sigma1 = NAN;
//Real DispP::X1 = NAN;
//Real DispP::sigma2 = NAN;
//Real DispP::X2 = NAN;
//Real DispP::magPertAmp = NAN;
Real DispP::densityPertAmp = NAN;
uint DispP::nSpaceSamples = 0;
uint DispP::nVelocitySamples = 0;

bool initializeProject(void) {
   uint seed;
   typedef Readparameters RP;
   RP::add("Dispersion.BX0", "Background field value (T)", 1.0e-9);
   RP::add("Dispersion.BY0", "Background field value (T)", 2.0e-9);
   RP::add("Dispersion.BZ0", "Background field value (T)", 3.0e-9);
   RP::add("Dispersion.rho", "Number density (m^-3)", 1.0e7);
   RP::add("Dispersion.Temperature", "Temperature (K)", 2.0e6);
//   RP::add("Dispersion.sigma1", "Width of 1st magnetic field Gaussian (m)", 10000.0);
//   RP::add("Dispersion.X1", "Center abscissa of 1st magnetic field Gaussian (m)",-20000.0);
//   RP::add("Dispersion.sigma2", "Width of 2nd magnetic field Gaussian (m)", 10000.0);
//   RP::add("Dispersion.X2", "Center abscissa of 2nd magnetic field Gaussian (m)",20000.0);
//   RP::add("Dispersion.magPertAmp", "Amplitude factor of the magnetic perturbation", 1.0);
   RP::add("Dispersion.densityPertAmp", "Amplitude factor of the density perturbation", 0.1);
   RP::add("Dispersion.seed","Seed integer for the srand() function", 42);
   RP::add("Dispersion.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
   RP::add("Dispersion.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   RP::parse();
   RP::get("Dispersion.BX0", DispP::BX0);
   RP::get("Dispersion.BY0", DispP::BY0);
   RP::get("Dispersion.BZ0", DispP::BZ0);
   RP::get("Dispersion.rho", DispP::DENSITY);
   RP::get("Dispersion.Temperature", DispP::TEMPERATURE);
//   RP::get("Dispersion.sigma1", DispP::sigma1);
//   RP::get("Dispersion.X1", DispP::X1);
//   RP::get("Dispersion.sigma2", DispP::sigma2);
//   RP::get("Dispersion.X2", DispP::X2);
//   RP::get("Dispersion.magPertAmp", DispP::magPertAmp);
   RP::get("Dispersion.densityPertAmp", DispP::densityPertAmp);
   RP::get("Dispersion.seed", seed);
   RP::get("Dispersion.nSpaceSamples", DispP::nSpaceSamples);
   RP::get("Dispersion.nVelocitySamples", DispP::nVelocitySamples);
   
   srand(seed);
   
   return true;
}

bool cellParametersChanged(creal& t) {return false;}


Real getDistribValue(creal& vx,creal& vy, creal& vz) {
   creal k = 1.3806505e-23; // Boltzmann
   creal mass = 1.67262171e-27; // m_p in kg
   return exp(- mass * (vx*vx + vy*vy + vz*vz) / (2.0 * k * DispP::TEMPERATURE));
}

Real calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
   creal mass = 1.67262171e-27; // m_p in kg
   creal k = 1.3806505e-23; // Boltzmann
   creal mu0 = 1.25663706144e-6; // mu_0
   creal q = 1.60217653e-19; // q_i
   
   static int spaceIndexOld[3] = {0};
   static int spaceIndex[3] = {0};
   static int rnd = 0;
   static int cpt = 0;
   spaceIndex[0] = (int) (x / dx);
   spaceIndex[1] = (int) (y / dy);
   spaceIndex[2] = (int) (z / dz);
   if(spaceIndex[0] != spaceIndexOld[0] ||
      spaceIndex[1] != spaceIndexOld[1] ||
      spaceIndex[2] != spaceIndexOld[2]) {
//      if(++cpt%(rand()%10+1) == 0)
      {
	 rnd = rand();
      }
   }
   spaceIndexOld[0] = spaceIndex[0];
   spaceIndexOld[1] = spaceIndex[1];
   spaceIndexOld[2] = spaceIndex[2];
   
   creal d_vx = dvx / (DispP::nVelocitySamples-1);
   creal d_vy = dvy / (DispP::nVelocitySamples-1);
   creal d_vz = dvz / (DispP::nVelocitySamples-1);
   Real avg = 0.0;
   #pragma omp parallel for collapse(3) reduction(+:avg)
   for (uint vi=0; vi<DispP::nVelocitySamples; ++vi)
      for (uint vj=0; vj<DispP::nVelocitySamples; ++vj)
	 for (uint vk=0; vk<DispP::nVelocitySamples; ++vk)
         {
	    avg += getDistribValue(vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz);
         }
   
   return avg *
   DispP::DENSITY * (1.0 + DispP::densityPertAmp * (0.5 - (double)rnd / (double)RAND_MAX)) *
   pow(mass / (2.0 * M_PI * k * DispP::TEMPERATURE), 1.5) /
//   (Parameters::vzmax - Parameters::vzmin) / 
(DispP::nVelocitySamples*DispP::nVelocitySamples*DispP::nVelocitySamples);
}

void calcBlockParameters(Real* blockParams) {
   //blockParams[BlockParams::Q_PER_M] = 1.0;
}

void calcCellParameters(Real* cellParams,creal& t) {
   creal x = cellParams[CellParams::XCRD];
   creal dx = cellParams[CellParams::DX];
   cellParams[CellParams::EX   ] = 0.0;
   cellParams[CellParams::EY   ] = 0.0;
   cellParams[CellParams::EZ   ] = 0.0;
   cellParams[CellParams::BX   ] = DispP::BX0;
   cellParams[CellParams::BY   ] = DispP::BY0;
   cellParams[CellParams::BZ   ] = DispP::BZ0;
/*   cellParams[CellParams::BZ   ] = DispP::B0 * (1.0 + 
   (-DispP::magPertAmp*(x+0.5*dx-DispP::X1) * exp(-0.5*(x+0.5*dx-DispP::X1)*(x+0.5*dx-DispP::X1) / (DispP::sigma1*DispP::sigma1))) - 
   (-DispP::magPertAmp*(x+0.5*dx-DispP::X2) * exp(-0.5*(x+0.5*dx-DispP::X2)*(x+0.5*dx-DispP::X2) / (DispP::sigma2*DispP::sigma2))));*/
}

// TODO use this instead: template <class Grid, class CellData> void calcSimParameters(Grid<CellData>& mpiGrid...
void calcSimParameters(dccrg::Dccrg<SpatialCell>& mpiGrid, creal& t, Real& /*dt*/) {
   std::vector<uint64_t> cells = mpiGrid.get_cells();
   for (uint i = 0; i < cells.size(); ++i) {
      calcCellParameters(mpiGrid[cells[i]]->cpu_cellParams, t);
   }
}
