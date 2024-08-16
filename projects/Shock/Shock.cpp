/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../object_wrapper.h"
#include "../../velocity_mesh_parameters.h"
#include "../../backgroundfield/backgroundfield.h"

#include "Shock.h"


namespace projects {
   Shock::Shock(): Project() { }
   Shock::~Shock() { }

   bool Shock::initialize(void) {return Project::initialize();}

   void Shock::addParameters() {
      typedef Readparameters RP;
      RP::add("Shock.BX0", "Background field value (T)", 1.0e-9);
      RP::add("Shock.BY0", "Background field value (T)", 2.0e-9);
      RP::add("Shock.BZ0", "Background field value (T)", 3.0e-9);
      RP::add("Shock.EX0", "Background electric field", 0.0);
      RP::add("Shock.VX0", "Bulk velocity in x", 0.0);
      RP::add("Shock.VY0", "Bulk velocity in y", 0.0);
      RP::add("Shock.VZ0", "Bulk velocuty in z", 0.0);
      RP::add("Shock.rho", "Number density (m^-3)", 1.0e7);
      RP::add("Shock.Temperature", "Temperature (K)", 2.0e6);
      RP::add("Shock.magPertAmp", "Amplitude of the magnetic perturbation", 1.0e-9);
      RP::add("Shock.densityPertAmp", "Amplitude factor of the density perturbation", 0.1);
      RP::add("Shock.velocityPertAmp", "Amplitude of the velocity perturbation", 1.0e6);
      RP::add("Shock.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
      RP::add("Shock.Scale_x", "Scale length in x (m)", 2.0e6);
      RP::add("Shock.Scale_y", "Scale length in y (m)", 2.0e6);
      RP::add("Shock.Sharp_Y", "Sharpness of tannh", 0.1);
   }

   void Shock::getParameters() {
      Project::getParameters();
      typedef Readparameters RP;

      if(getObjectWrapper().particleSpecies.size() > 1) {
         std::cerr << "The selected project does not support multiple particle populations! Aborting in " << __FILE__ << " line " << __LINE__ << std::endl;
         abort();
      }
      RP::get("Shock.BX0", this->BX0);
      RP::get("Shock.BY0", this->BY0);
      RP::get("Shock.BZ0", this->BZ0);
      RP::get("Shock.EX0", this->EX0);
      RP::get("Shock.VX0", this->VX0);
      RP::get("Shock.VY0", this->VY0);
      RP::get("Shock.VZ0", this->VZ0);
      RP::get("Shock.rho", this->DENSITY);
      RP::get("Shock.Temperature", this->TEMPERATURE);
      RP::get("Shock.magPertAmp", this->magPertAmp);
      RP::get("Shock.densityPertAmp", this->densityPertAmp);
      RP::get("Shock.velocityPertAmp", this->velocityPertAmp);
      RP::get("Shock.maxwCutoff", this->maxwCutoff);
      RP::get("Shock.Scale_x", this->SCA_X);
      RP::get("Shock.Scale_y", this->SCA_Y);
      RP::get("Shock.Sharp_Y", this->Sharp_Y);
   }

   inline Real Shock::getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, const uint popID) const {
      creal kb = physicalconstants::K_B;
      creal mass = physicalconstants::MASS_PROTON;
      return exp(- mass * ((vx-this->VX0)*(vx-this->VX0) + (vy-this->VY0)*(vy-this->VY0)+ (vz-this->VZ0)*(vz-this->VZ0)) / (2.0 * kb * this->TEMPERATURE));
      //*exp(-pow(x-Parameters::xmax/2.0, 2.0)/pow(this->SCA_X, 2.0))*exp(-pow(y-Parameters::ymax/4.0, 2.0)/pow(this->SCA_Y, 2.0));
   }

   Real Shock::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz,
           creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,const uint popID) const {
      const size_t meshID = getObjectWrapper().particleSpecies[popID].velocityMesh;
      vmesh::MeshParameters& velocityMeshParams = vmesh::getMeshWrapper()->velocityMeshes->at(meshID);
      if (vx < velocityMeshParams.meshMinLimits[0] + 0.5*dvx ||
          vy < velocityMeshParams.meshMinLimits[1] + 0.5*dvy ||
          vz < velocityMeshParams.meshMinLimits[2] + 0.5*dvz ||
          vx > velocityMeshParams.meshMaxLimits[0] - 1.5*dvx ||
          vy > velocityMeshParams.meshMaxLimits[1] - 1.5*dvy ||
          vz > velocityMeshParams.meshMaxLimits[2] - 1.5*dvz) {
         return 0.0;
      }
      
      creal mass = physicalconstants::MASS_PROTON;
      creal kb = physicalconstants::K_B;
      
      creal result = getDistribValue(x+0.5*dx, y+0.5*dy, z+0.5*dz, vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz, popID)
         * this->DENSITY * pow(mass / (2.0 * M_PI * kb * this->TEMPERATURE), 1.5);               

      if(result < this->maxwCutoff) {
         return 0.0;
      } else {
         return result;
      }
   }

   void Shock::calcCellParameters(spatial_cell::SpatialCell* cell,creal& t) { }
   
   void Shock::setProjectBField(
      BFieldFsGrid & perBGrid,
      BgBFsGrid & BgBGrid,
      TechnicalFsGrid & technicalGrid
   ) {
      setBackgroundFieldToZero(BgBGrid);
      
      if(!P::isRestart) {
         auto localSize = perBGrid.getLocalSize();
         
#pragma omp parallel for collapse(3)
         for (FsGridTools::FsIndex_t x = 0; x < localSize[0]; ++x) {
            for (FsGridTools::FsIndex_t y = 0; y < localSize[1]; ++y) {
               for (FsGridTools::FsIndex_t z = 0; z < localSize[2]; ++z) {
                  const auto xyz = perBGrid.getPhysicalCoords(x, y, z);
                  auto cell = perBGrid.get(x, y, z);
                  
                  cell[fsgrids::bfield::PERBX] = 0.0;
                  cell[fsgrids::bfield::PERBY] = 0.0;
                  cell[fsgrids::bfield::PERBZ] = this->BZ0*(3.0 + 2.0*tanh((xyz[1] - FSParams.ymax/2.0)/(this->Sharp_Y*FSParams.ymax)));
               }
            }
         }
      }
   }
}//namespace projects
