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

#include "fs_common.h"
#include "derivatives.hpp"
#include "fs_limiters.h"
#include "../arch/arch_sysboundary_api.h"

/*! \brief Low-level spatial derivatives calculation.
 *
 * For the cell with ID cellID calculate the spatial derivatives or apply the derivative boundary conditions defined in project.h.
 * Uses RHO, V[XYZ] and B[XYZ] in the first-order time accuracy method and in the second step of the second-order method,
 * and RHO_DT2, V[XYZ]1 and B[XYZ]1 in the first step of the second-order method.
 *
 * For sysBoundaryLayer 1 or 2, we are near a boundary, and we wish to use regular centered differences instead of slope limiter-adjusted values.
 * This is to minimize oscillations as a smooth behaviour is required near artificial boundaries,
 * unlike at boundaries and shocks inside the simulation domain.
 *
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param momentsGrid fsGrid holding the moment quantities
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * 
 * \sa calculateDerivativesSimple calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivatives(
   cint i,
   cint j,
   cint k,
   const arch::buf<BFieldFsGrid> & perBGrid,
   const arch::buf<MomentsFsGrid> & momentsGrid,
   const arch::buf<DPerBFsGrid> & dPerBGrid,
   const arch::buf<DMomentsFsGrid> & dMomentsGrid,
   const arch::buf<TechnicalFsGrid> & technicalGrid,
   const arch::buf<SysBoundary>& sysBoundaries,
   cint& RKCase
) {
   auto dPerB = dPerBGrid.get(i,j,k);
   auto dMoments = dMomentsGrid.get(i,j,k);

   // Get boundary flag for the cell:
   cuint sysBoundaryFlag  = technicalGrid.get(i,j,k)->sysBoundaryFlag;
   cuint sysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;

   // Constants for electron pressure derivatives
   // Upstream pressure
   Real Peupstream = Parameters::electronTemperature * Parameters::electronDensity * physicalconstants::K_B;
   Real Peconst = Peupstream * pow(Parameters::electronDensity, -Parameters::electronPTindex);

   auto centMoments = momentsGrid.get(i,j,k);
   auto centPerB = perBGrid.get(i,j,k);
   #ifdef DEBUG_SOLVERS
   if (centMoments[fsgrids::moments::RHOM] <= 0) {
      std::cerr << __FILE__ << ":" << __LINE__
         << (centMoments[fsgrids::moments::RHOM] < 0 ? " Negative" : " Zero") << " density in spatial cell at (" << i << " " << j << " " << k << ")"
         << std::endl;
      abort();
   }
   #endif
   
   // Calculate x-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      auto leftPerB = perBGrid.get(i-1,j,k);
      auto rghtPerB = perBGrid.get(i+1,j,k);
      auto leftMoments = momentsGrid.get(i-1,j,k);
      auto rghtMoments = momentsGrid.get(i+1,j,k);
      #ifdef DEBUG_SOLVERS
      if (leftMoments[fsgrids::moments::RHOM] <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (leftMoments[fsgrids::moments::RHOM] < 0 ? " Negative" : " Zero") << " density in spatial cell " //<< leftNbrID
            << std::endl;
         abort();
      }
      if (rghtMoments[fsgrids::moments::RHOM] <= 0) {
         std::cerr << __FILE__ << ":" << __LINE__
            << (rghtMoments[fsgrids::moments::RHOM] < 0 ? " Negative" : " Zero") << " density in spatial cell " //<< rightNbrID
            << std::endl;
         abort();
      }
      #endif
      
      if(sysBoundaryLayer == 1 || sysBoundaryLayer == 2) {
         dMoments[fsgrids::dmoments::drhomdx] = (rghtMoments[fsgrids::moments::RHOM]-leftMoments[fsgrids::moments::RHOM])/2;
         dMoments[fsgrids::dmoments::drhoqdx] = (rghtMoments[fsgrids::moments::RHOQ]-leftMoments[fsgrids::moments::RHOQ])/2;
         dMoments[fsgrids::dmoments::dp11dx] = (rghtMoments[fsgrids::moments::P_11]-leftMoments[fsgrids::moments::P_11])/2;
         dMoments[fsgrids::dmoments::dp22dx] = (rghtMoments[fsgrids::moments::P_22]-leftMoments[fsgrids::moments::P_22])/2;
         dMoments[fsgrids::dmoments::dp33dx] = (rghtMoments[fsgrids::moments::P_33]-leftMoments[fsgrids::moments::P_33])/2;

         dMoments[fsgrids::dmoments::dVxdx]  = (rghtMoments[fsgrids::moments::VX]-leftMoments[fsgrids::moments::VX])/2;
         dMoments[fsgrids::dmoments::dVydx]  = (rghtMoments[fsgrids::moments::VY]-leftMoments[fsgrids::moments::VY])/2;
         dMoments[fsgrids::dmoments::dVzdx]  = (rghtMoments[fsgrids::moments::VZ]-leftMoments[fsgrids::moments::VZ])/2;
         dPerB[fsgrids::dperb::dPERBydx]  = (rghtPerB[fsgrids::bfield::PERBY]-leftPerB[fsgrids::bfield::PERBY])/2;
         dPerB[fsgrids::dperb::dPERBzdx]  = (rghtPerB[fsgrids::bfield::PERBZ]-leftPerB[fsgrids::bfield::PERBZ])/2;
      } else {
         dMoments[fsgrids::dmoments::drhomdx] = limiter(leftMoments[fsgrids::moments::RHOM],centMoments[fsgrids::moments::RHOM],rghtMoments[fsgrids::moments::RHOM]);
         dMoments[fsgrids::dmoments::drhoqdx] = limiter(leftMoments[fsgrids::moments::RHOQ],centMoments[fsgrids::moments::RHOQ],rghtMoments[fsgrids::moments::RHOQ]);
         dMoments[fsgrids::dmoments::dp11dx] = limiter(leftMoments[fsgrids::moments::P_11],centMoments[fsgrids::moments::P_11],rghtMoments[fsgrids::moments::P_11]);
         dMoments[fsgrids::dmoments::dp22dx] = limiter(leftMoments[fsgrids::moments::P_22],centMoments[fsgrids::moments::P_22],rghtMoments[fsgrids::moments::P_22]);
         dMoments[fsgrids::dmoments::dp33dx] = limiter(leftMoments[fsgrids::moments::P_33],centMoments[fsgrids::moments::P_33],rghtMoments[fsgrids::moments::P_33]);

         dMoments[fsgrids::dmoments::dVxdx]  = limiter(leftMoments[fsgrids::moments::VX], centMoments[fsgrids::moments::VX], rghtMoments[fsgrids::moments::VX]);
         dMoments[fsgrids::dmoments::dVydx]  = limiter(leftMoments[fsgrids::moments::VY], centMoments[fsgrids::moments::VY], rghtMoments[fsgrids::moments::VY]);
         dMoments[fsgrids::dmoments::dVzdx]  = limiter(leftMoments[fsgrids::moments::VZ], centMoments[fsgrids::moments::VZ], rghtMoments[fsgrids::moments::VZ]);
         dPerB[fsgrids::dperb::dPERBydx]  = limiter(leftPerB[fsgrids::bfield::PERBY],centPerB[fsgrids::bfield::PERBY],rghtPerB[fsgrids::bfield::PERBY]);
         dPerB[fsgrids::dperb::dPERBzdx]  = limiter(leftPerB[fsgrids::bfield::PERBZ],centPerB[fsgrids::bfield::PERBZ],rghtPerB[fsgrids::bfield::PERBZ]);
      }

      // pres_e = const * np.power(rho_e, index)
      dMoments[fsgrids::dmoments::dPedx] = Peconst * limiter(pow(leftMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex));      
      
      if (FSParams.ohmHallTerm < 2 || sysBoundaryLayer == 1) {
        dPerB[fsgrids::dperb::dPERBydxx] = 0.0;
        dPerB[fsgrids::dperb::dPERBzdxx] = 0.0;
      } else {
        dPerB[fsgrids::dperb::dPERBydxx] = leftPerB[fsgrids::bfield::PERBY] + rghtPerB[fsgrids::bfield::PERBY] - 2.0*centPerB[fsgrids::bfield::PERBY];
        dPerB[fsgrids::dperb::dPERBzdxx] = leftPerB[fsgrids::bfield::PERBZ] + rghtPerB[fsgrids::bfield::PERBZ] - 2.0*centPerB[fsgrids::bfield::PERBZ];
      }
   } else {
      SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 0);
   }

   // Calculate y-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      auto leftPerB = perBGrid.get(i,j-1,k);
      auto rghtPerB = perBGrid.get(i,j+1,k);
      auto leftMoments = momentsGrid.get(i,j-1,k);
      auto rghtMoments = momentsGrid.get(i,j+1,k);

      if(sysBoundaryLayer == 1 || sysBoundaryLayer == 2) {
         dMoments[fsgrids::dmoments::drhomdy] = (rghtMoments[fsgrids::moments::RHOM]-leftMoments[fsgrids::moments::RHOM])/2;
         dMoments[fsgrids::dmoments::drhoqdy] = (rghtMoments[fsgrids::moments::RHOQ]-leftMoments[fsgrids::moments::RHOQ])/2;
         dMoments[fsgrids::dmoments::dp11dy] = (rghtMoments[fsgrids::moments::P_11]-leftMoments[fsgrids::moments::P_11])/2;
         dMoments[fsgrids::dmoments::dp22dy] = (rghtMoments[fsgrids::moments::P_22]-leftMoments[fsgrids::moments::P_22])/2;
         dMoments[fsgrids::dmoments::dp33dy] = (rghtMoments[fsgrids::moments::P_33]-leftMoments[fsgrids::moments::P_33])/2;
         dMoments[fsgrids::dmoments::dVxdy]  = (rghtMoments[fsgrids::moments::VX]-leftMoments[fsgrids::moments::VX])/2;
         dMoments[fsgrids::dmoments::dVydy]  = (rghtMoments[fsgrids::moments::VY]-leftMoments[fsgrids::moments::VY])/2;
         dMoments[fsgrids::dmoments::dVzdy]  = (rghtMoments[fsgrids::moments::VZ]-leftMoments[fsgrids::moments::VZ])/2;

         dPerB[fsgrids::dperb::dPERBxdy]  = (rghtPerB[fsgrids::bfield::PERBX]-leftPerB[fsgrids::bfield::PERBX])/2;
         dPerB[fsgrids::dperb::dPERBzdy]  = (rghtPerB[fsgrids::bfield::PERBZ]-leftPerB[fsgrids::bfield::PERBZ])/2;
      } else {
         dMoments[fsgrids::dmoments::drhomdy] = limiter(leftMoments[fsgrids::moments::RHOM],centMoments[fsgrids::moments::RHOM],rghtMoments[fsgrids::moments::RHOM]);
         dMoments[fsgrids::dmoments::drhoqdy] = limiter(leftMoments[fsgrids::moments::RHOQ],centMoments[fsgrids::moments::RHOQ],rghtMoments[fsgrids::moments::RHOQ]);
         dMoments[fsgrids::dmoments::dp11dy] = limiter(leftMoments[fsgrids::moments::P_11],centMoments[fsgrids::moments::P_11],rghtMoments[fsgrids::moments::P_11]);
         dMoments[fsgrids::dmoments::dp22dy] = limiter(leftMoments[fsgrids::moments::P_22],centMoments[fsgrids::moments::P_22],rghtMoments[fsgrids::moments::P_22]);
         dMoments[fsgrids::dmoments::dp33dy] = limiter(leftMoments[fsgrids::moments::P_33],centMoments[fsgrids::moments::P_33],rghtMoments[fsgrids::moments::P_33]);
         dMoments[fsgrids::dmoments::dVxdy]  = limiter(leftMoments[fsgrids::moments::VX], centMoments[fsgrids::moments::VX], rghtMoments[fsgrids::moments::VX]);
         dMoments[fsgrids::dmoments::dVydy]  = limiter(leftMoments[fsgrids::moments::VY], centMoments[fsgrids::moments::VY], rghtMoments[fsgrids::moments::VY]);
         dMoments[fsgrids::dmoments::dVzdy]  = limiter(leftMoments[fsgrids::moments::VZ], centMoments[fsgrids::moments::VZ], rghtMoments[fsgrids::moments::VZ]);

         dPerB[fsgrids::dperb::dPERBxdy]  = limiter(leftPerB[fsgrids::bfield::PERBX],centPerB[fsgrids::bfield::PERBX],rghtPerB[fsgrids::bfield::PERBX]);
         dPerB[fsgrids::dperb::dPERBzdy]  = limiter(leftPerB[fsgrids::bfield::PERBZ],centPerB[fsgrids::bfield::PERBZ],rghtPerB[fsgrids::bfield::PERBZ]);
      }

      // pres_e = const * np.power(rho_e, index)
      dMoments[fsgrids::dmoments::dPedy] = Peconst * limiter(pow(leftMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex));      

      if (FSParams.ohmHallTerm < 2 || sysBoundaryLayer == 1) {
         dPerB[fsgrids::dperb::dPERBxdyy] = 0.0;
         dPerB[fsgrids::dperb::dPERBzdyy] = 0.0;
      } else {
         dPerB[fsgrids::dperb::dPERBxdyy] = leftPerB[fsgrids::bfield::PERBX] + rghtPerB[fsgrids::bfield::PERBX] - 2.0*centPerB[fsgrids::bfield::PERBX];
         dPerB[fsgrids::dperb::dPERBzdyy] = leftPerB[fsgrids::bfield::PERBZ] + rghtPerB[fsgrids::bfield::PERBZ] - 2.0*centPerB[fsgrids::bfield::PERBZ];
      }
      
   } else {
      SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 1);
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if ((sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) || (sysBoundaryLayer == 1)) {
      auto leftPerB = perBGrid.get(i,j,k-1);
      auto rghtPerB = perBGrid.get(i,j,k+1);
      auto leftMoments = momentsGrid.get(i,j,k-1);
      auto rghtMoments = momentsGrid.get(i,j,k+1);
      if(sysBoundaryLayer == 1 || sysBoundaryLayer == 2) {
         dMoments[fsgrids::dmoments::drhomdz] = (rghtMoments[fsgrids::moments::RHOM]-leftMoments[fsgrids::moments::RHOM])/2;
         dMoments[fsgrids::dmoments::drhoqdz] = (rghtMoments[fsgrids::moments::RHOQ]-leftMoments[fsgrids::moments::RHOQ])/2;
         dMoments[fsgrids::dmoments::dp11dz] = (rghtMoments[fsgrids::moments::P_11]-leftMoments[fsgrids::moments::P_11])/2;
         dMoments[fsgrids::dmoments::dp22dz] = (rghtMoments[fsgrids::moments::P_22]-leftMoments[fsgrids::moments::P_22])/2;
         dMoments[fsgrids::dmoments::dp33dz] = (rghtMoments[fsgrids::moments::P_33]-leftMoments[fsgrids::moments::P_33])/2;
         dMoments[fsgrids::dmoments::dVxdz]  = (rghtMoments[fsgrids::moments::VX]-leftMoments[fsgrids::moments::VX])/2;
         dMoments[fsgrids::dmoments::dVydz]  = (rghtMoments[fsgrids::moments::VY]-leftMoments[fsgrids::moments::VY])/2;
         dMoments[fsgrids::dmoments::dVzdz]  = (rghtMoments[fsgrids::moments::VZ]-leftMoments[fsgrids::moments::VZ])/2;
         
         dPerB[fsgrids::dperb::dPERBxdz]  = (rghtPerB[fsgrids::bfield::PERBX]-leftPerB[fsgrids::bfield::PERBX])/2;
         dPerB[fsgrids::dperb::dPERBydz]  = (rghtPerB[fsgrids::bfield::PERBY]-leftPerB[fsgrids::bfield::PERBY])/2;
      } else {
         dMoments[fsgrids::dmoments::drhomdz] = limiter(leftMoments[fsgrids::moments::RHOM],centMoments[fsgrids::moments::RHOM],rghtMoments[fsgrids::moments::RHOM]);
         dMoments[fsgrids::dmoments::drhoqdz] = limiter(leftMoments[fsgrids::moments::RHOQ],centMoments[fsgrids::moments::RHOQ],rghtMoments[fsgrids::moments::RHOQ]);
         dMoments[fsgrids::dmoments::dp11dz] = limiter(leftMoments[fsgrids::moments::P_11],centMoments[fsgrids::moments::P_11],rghtMoments[fsgrids::moments::P_11]);
         dMoments[fsgrids::dmoments::dp22dz] = limiter(leftMoments[fsgrids::moments::P_22],centMoments[fsgrids::moments::P_22],rghtMoments[fsgrids::moments::P_22]);
         dMoments[fsgrids::dmoments::dp33dz] = limiter(leftMoments[fsgrids::moments::P_33],centMoments[fsgrids::moments::P_33],rghtMoments[fsgrids::moments::P_33]);
         dMoments[fsgrids::dmoments::dVxdz]  = limiter(leftMoments[fsgrids::moments::VX], centMoments[fsgrids::moments::VX], rghtMoments[fsgrids::moments::VX]);
         dMoments[fsgrids::dmoments::dVydz]  = limiter(leftMoments[fsgrids::moments::VY], centMoments[fsgrids::moments::VY], rghtMoments[fsgrids::moments::VY]);
         dMoments[fsgrids::dmoments::dVzdz]  = limiter(leftMoments[fsgrids::moments::VZ], centMoments[fsgrids::moments::VZ], rghtMoments[fsgrids::moments::VZ]);
         
         dPerB[fsgrids::dperb::dPERBxdz]  = limiter(leftPerB[fsgrids::bfield::PERBX],centPerB[fsgrids::bfield::PERBX],rghtPerB[fsgrids::bfield::PERBX]);
         dPerB[fsgrids::dperb::dPERBydz]  = limiter(leftPerB[fsgrids::bfield::PERBY],centPerB[fsgrids::bfield::PERBY],rghtPerB[fsgrids::bfield::PERBY]);
      }

      // pres_e = const * np.power(rho_e, index)
      dMoments[fsgrids::dmoments::dPedz] = Peconst * limiter(pow(leftMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex),pow(centMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex),pow(rghtMoments[fsgrids::moments::RHOQ]/physicalconstants::CHARGE,Parameters::electronPTindex));      

      if (FSParams.ohmHallTerm < 2 || sysBoundaryLayer == 1) {
        dPerB[fsgrids::dperb::dPERBxdzz] = 0.0;
        dPerB[fsgrids::dperb::dPERBydzz] = 0.0;
      } else {
        dPerB[fsgrids::dperb::dPERBxdzz] = leftPerB[fsgrids::bfield::PERBX] + rghtPerB[fsgrids::bfield::PERBX] - 2.0*centPerB[fsgrids::bfield::PERBX];
        dPerB[fsgrids::dperb::dPERBydzz] = leftPerB[fsgrids::bfield::PERBY] + rghtPerB[fsgrids::bfield::PERBY] - 2.0*centPerB[fsgrids::bfield::PERBY];
      }
      
   } else {
      SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 2);
   }
   
   if (FSParams.ohmHallTerm < 2 || sysBoundaryLayer == 1) {
      dPerB[fsgrids::dperb::dPERBxdyz] = 0.0;
      dPerB[fsgrids::dperb::dPERBydxz] = 0.0;
      dPerB[fsgrids::dperb::dPERBzdxy] = 0.0;
   } else {
      // Calculate xy mixed derivatives:
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         auto botLeft = perBGrid.get(i-1,j-1,k);
         auto botRght = perBGrid.get(i+1,j-1,k);
         auto topLeft = perBGrid.get(i-1,j+1,k);
         auto topRght = perBGrid.get(i+1,j+1,k);
         dPerB[fsgrids::dperb::dPERBzdxy] = FOURTH * (botLeft[fsgrids::bfield::PERBZ] + topRght[fsgrids::bfield::PERBZ] - botRght[fsgrids::bfield::PERBZ] - topLeft[fsgrids::bfield::PERBZ]);
      } else {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 3);
      }
      
      // Calculate xz mixed derivatives:
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         auto botLeft = perBGrid.get(i-1,j,k-1);
         auto botRght = perBGrid.get(i+1,j,k-1);
         auto topLeft = perBGrid.get(i-1,j,k+1);
         auto topRght = perBGrid.get(i+1,j,k+1);
         dPerB[fsgrids::dperb::dPERBydxz] = FOURTH * (botLeft[fsgrids::bfield::PERBY] + topRght[fsgrids::bfield::PERBY] - botRght[fsgrids::bfield::PERBY] - topLeft[fsgrids::bfield::PERBY]);
      } else {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 4);
      }
      
      // Calculate yz mixed derivatives:
      if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
         auto botLeft = perBGrid.get(i,j-1,k-1);
         auto botRght = perBGrid.get(i,j+1,k-1);
         auto topLeft = perBGrid.get(i,j-1,k+1);
         auto topRght = perBGrid.get(i,j+1,k+1);
         dPerB[fsgrids::dperb::dPERBxdyz] = FOURTH * (botLeft[fsgrids::bfield::PERBX] + topRght[fsgrids::bfield::PERBX] - botRght[fsgrids::bfield::PERBX] - topLeft[fsgrids::bfield::PERBX]);
      } else {
         SBC::SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, 5);
      }
   }
}


/*! \brief High-level derivative calculation wrapper function.
 * 

 * B has to be updated because after the system boundary update in propagateMagneticFieldSimple there is no consistent state of B yet everywhere.
 * 
 * Then the derivatives are calculated.
 * 
 * \param perBGrid fsGrid holding the perturbed B quantities
 * \param perBDt2Grid fsGrid holding the perturbed B quantities at runge-kutta t=0.5
 * \param momentsGrid fsGrid holding the moment quantities
 * \param momentsDt2Grid fsGrid holding the moment quantities at runge-kutta t=0.5
 * \param dPerBGrid fsGrid holding the derivatives of perturbed B
 * \param dMomentsGrid fsGrid holding the derviatives of moments
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 * \param RKCase Element in the enum defining the Runge-Kutta method steps
 * \param communicateMoments If true, the derivatives of moments (rho, V, P) are communicated to neighbours.
 
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateBVOLDerivatives
 */
void calculateDerivativesSimple(
   arch::buf<BFieldFsGrid> & perBGrid,
   arch::buf<BFieldFsGrid> & perBDt2Grid,
   arch::buf<MomentsFsGrid> & momentsGrid,
   arch::buf<MomentsFsGrid> & momentsDt2Grid,
   arch::buf<DPerBFsGrid> & dPerBGrid,
   arch::buf<DMomentsFsGrid> & dMomentsGrid,
   arch::buf<TechnicalFsGrid> & technicalGrid,
   arch::buf<SysBoundary>& sysBoundaries,
   cint& RKCase,
   const bool communicateMoments) {
   int timer;
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const int* gridDims = &technicalGrid.grid()->getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   
   phiprof::start("Calculate face derivatives");
   
   timer=phiprof::initializeTimer("MPI","MPI");
   phiprof::start(timer);
   
   switch (RKCase) {
      case RK_ORDER1:
         // Means initialising the solver as well as RK_ORDER1
         // standard case Exchange PERB* with neighbours
         // The update of PERB[XYZ] is needed after the system
         // boundary update of propagateMagneticFieldSimple.
         perBGrid.syncHostData();
         perBGrid.grid()->updateGhostCells();
         perBGrid.syncDeviceData();
         if(communicateMoments) {
            momentsGrid.syncHostData();
            momentsGrid.grid()->updateGhostCells();
            momentsGrid.syncDeviceData();
         }
         break;
      case RK_ORDER2_STEP1:
         // Exchange PERB*_DT2,RHO_DT2,V*_DT2 with neighbours The
         // update of PERB[XYZ]_DT2 is needed after the system
         // boundary update of propagateMagneticFieldSimple.
         perBDt2Grid.syncHostData();
         perBDt2Grid.grid()->updateGhostCells();
         perBDt2Grid.syncDeviceData();
         if(communicateMoments) {
            momentsDt2Grid.syncHostData();
            momentsDt2Grid.grid()->updateGhostCells();
            momentsDt2Grid.syncDeviceData();
         }
         break;
      case RK_ORDER2_STEP2:
         // Exchange PERB*,RHO,V* with neighbours The update of B
         // is needed after the system boundary update of
         // propagateMagneticFieldSimple.
         perBGrid.syncHostData();
         perBGrid.grid()->updateGhostCells();
         perBGrid.syncDeviceData();
         if(communicateMoments) {
            momentsGrid.syncHostData();
            momentsGrid.grid()->updateGhostCells();
            momentsGrid.syncDeviceData();
         }
         break;
      default:
         cerr << __FILE__ << ":" << __LINE__ << " Went through switch, this should not happen." << endl;
         abort();
   }
   
   phiprof::stop(timer);

   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);

   // Calculate derivatives
   arch::parallel_for({(uint)gridDims[0], (uint)gridDims[1], (uint)gridDims[2]}, ARCH_LOOP_LAMBDA(int i, int j, int k) { 
      if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) return;
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         calculateDerivatives(i,j,k, perBGrid, momentsGrid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, RKCase);
      } else {
         calculateDerivatives(i,j,k, perBDt2Grid, momentsDt2Grid, dPerBGrid, dMomentsGrid, technicalGrid, sysBoundaries, RKCase);
      }
   });

   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   phiprof::stop("Calculate face derivatives",N_cells,"Spatial Cells");   
}

/*! \brief Low-level spatial derivatives calculation.
 *
 * Calculate the spatial derivatives of BVOL or set them to zero.
 *
 * For sysBoundaryLayer 1 or 2, we are near a boundary, and we wish to use regular centered differences instead of slope limiter-adjusted values.
 * This is to minimize oscillations as a smooth behaviour is required near artificial boundaries,
 * unlike at boundaries and shocks inside the simulation domain.
 *
 * \param volGrid fsGrid holding the volume averaged fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param sysBoundaries System boundary conditions existing
 *
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateDerivativesSimple
 */

void calculateBVOLDerivatives(
   const arch::buf<VolFsGrid> & volGrid,
   const arch::buf<TechnicalFsGrid> & technicalGrid,
   cint i,
   cint j,
   cint k,
   const arch::buf<SysBoundary>& sysBoundaries
) {
   auto array = volGrid.get(i,j,k);
   
   cuint sysBoundaryFlag = technicalGrid.get(i,j,k)->sysBoundaryFlag;
   cuint sysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;

   // Calculate x-derivatives (is not TVD for AMR mesh):
   if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY || sysBoundaryLayer == 1) {

      auto left = volGrid.get(i-1,j,k);
      auto rght = volGrid.get(i+1,j,k);
      
      if (sysBoundaryLayer == 1 || sysBoundaryLayer == 2) {
         array[fsgrids::volfields::dPERBXVOLdx] = (rght[fsgrids::volfields::PERBXVOL]-left[fsgrids::volfields::PERBXVOL])/2;
         array[fsgrids::volfields::dPERBYVOLdx] = (rght[fsgrids::volfields::PERBYVOL]-left[fsgrids::volfields::PERBYVOL])/2;
         array[fsgrids::volfields::dPERBZVOLdx] = (rght[fsgrids::volfields::PERBZVOL]-left[fsgrids::volfields::PERBZVOL])/2;
      } else {
         array[fsgrids::volfields::dPERBXVOLdx] = limiter(left[fsgrids::volfields::PERBXVOL],array[fsgrids::volfields::PERBXVOL],rght[fsgrids::volfields::PERBXVOL]);
         array[fsgrids::volfields::dPERBYVOLdx] = limiter(left[fsgrids::volfields::PERBYVOL],array[fsgrids::volfields::PERBYVOL],rght[fsgrids::volfields::PERBYVOL]);
         array[fsgrids::volfields::dPERBZVOLdx] = limiter(left[fsgrids::volfields::PERBZVOL],array[fsgrids::volfields::PERBZVOL],rght[fsgrids::volfields::PERBZVOL]);
      }
   } else {
      SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 0);
   }
   
   // Calculate y-derivatives (is not TVD for AMR mesh):
   if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY || sysBoundaryLayer == 1) {
      auto left = volGrid.get(i,j-1,k);
      auto rght = volGrid.get(i,j+1,k);
      
      if (sysBoundaryLayer == 1 || sysBoundaryLayer == 2) {
         array[fsgrids::volfields::dPERBXVOLdy] = (rght[fsgrids::volfields::PERBXVOL]-left[fsgrids::volfields::PERBXVOL])/2;
         array[fsgrids::volfields::dPERBYVOLdy] = (rght[fsgrids::volfields::PERBYVOL]-left[fsgrids::volfields::PERBYVOL])/2;
         array[fsgrids::volfields::dPERBZVOLdy] = (rght[fsgrids::volfields::PERBZVOL]-left[fsgrids::volfields::PERBZVOL])/2;
      } else {
         array[fsgrids::volfields::dPERBXVOLdy] = limiter(left[fsgrids::volfields::PERBXVOL],array[fsgrids::volfields::PERBXVOL],rght[fsgrids::volfields::PERBXVOL]);
         array[fsgrids::volfields::dPERBYVOLdy] = limiter(left[fsgrids::volfields::PERBYVOL],array[fsgrids::volfields::PERBYVOL],rght[fsgrids::volfields::PERBYVOL]);
         array[fsgrids::volfields::dPERBZVOLdy] = limiter(left[fsgrids::volfields::PERBZVOL],array[fsgrids::volfields::PERBZVOL],rght[fsgrids::volfields::PERBZVOL]);
      }
   } else {
      SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 1);
   }
   
   // Calculate z-derivatives (is not TVD for AMR mesh):
   if (sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY || sysBoundaryLayer == 1) {
      auto left = volGrid.get(i,j,k-1);
      auto rght = volGrid.get(i,j,k+1);
      
      if (sysBoundaryLayer == 1 || sysBoundaryLayer == 2) {
         array[fsgrids::volfields::dPERBXVOLdz] = (rght[fsgrids::volfields::PERBXVOL]-left[fsgrids::volfields::PERBXVOL])/2;
         array[fsgrids::volfields::dPERBYVOLdz] = (rght[fsgrids::volfields::PERBYVOL]-left[fsgrids::volfields::PERBYVOL])/2;
         array[fsgrids::volfields::dPERBZVOLdz] = (rght[fsgrids::volfields::PERBZVOL]-left[fsgrids::volfields::PERBZVOL])/2;
      } else {
         array[fsgrids::volfields::dPERBXVOLdz] = limiter(left[fsgrids::volfields::PERBXVOL],array[fsgrids::volfields::PERBXVOL],rght[fsgrids::volfields::PERBXVOL]);
         array[fsgrids::volfields::dPERBYVOLdz] = limiter(left[fsgrids::volfields::PERBYVOL],array[fsgrids::volfields::PERBYVOL],rght[fsgrids::volfields::PERBYVOL]);
         array[fsgrids::volfields::dPERBZVOLdz] = limiter(left[fsgrids::volfields::PERBZVOL],array[fsgrids::volfields::PERBZVOL],rght[fsgrids::volfields::PERBZVOL]);
      }
   } else {
      SBC::SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, 2);
   }
}

/*! \brief High-level derivative calculation wrapper function.
 * 
 * BVOL has been calculated locally by calculateVolumeAveragedFields but not communicated.
 * For the acceleration step one needs the cross-derivatives of BVOL
 * 
 * \param volGrid fsGrid holding the volume averaged fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 * 
 * \sa calculateDerivatives calculateBVOLDerivatives calculateDerivativesSimple
 */
void calculateBVOLDerivativesSimple(
   arch::buf<VolFsGrid> & volGrid,
   arch::buf<TechnicalFsGrid> & technicalGrid,
   arch::buf<SysBoundary>& sysBoundaries
) {
   int timer;
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const int* gridDims = &technicalGrid.grid()->getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   
   phiprof::start("Calculate volume derivatives");
   
   timer=phiprof::initializeTimer("Start comm","MPI");
   phiprof::start(timer);
   volGrid.syncHostData();
   volGrid.grid()->updateGhostCells();
   // volGrid.syncDeviceData();
   
   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   
   // Calculate derivatives
   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);
   
   //ARCH_TODO
   // arch::parallel_for({(uint)gridDims[0], (uint)gridDims[1], (uint)gridDims[2]}, ARCH_LOOP_LAMBDA(int i, int j, int k) {
   //       if (technicalGrid.get(i,j,k)->sysBoundaryFlag != sysboundarytype::DO_NOT_COMPUTE) {
   //          calculateBVOLDerivatives(volGrid,technicalGrid,i,j,k,sysBoundaries);
   //       }
   //    });
   technicalGrid.syncHostData();
   #pragma omp parallel for collapse(3)
   for (int k=0; k<gridDims[2]; k++) {
      for (int j=0; j<gridDims[1]; j++) {
         for (int i=0; i<gridDims[0]; i++) {
            if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
               continue;
            }
            calculateBVOLDerivatives(volGrid,technicalGrid,i,j,k,sysBoundaries);
         }
      }
   }
   volGrid.syncDeviceData();
   technicalGrid.syncDeviceData();

   phiprof::stop(timer,N_cells,"Spatial Cells");

   phiprof::stop("Calculate volume derivatives",N_cells,"Spatial Cells");
}

/*! \brief Low-level curvature calculation.
 * 
 * 
 * \param volGrid fsGrid holding the volume averaged fields
 * \param bgbGrid fsGrid holding the background fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param i,j,k fsGrid cell coordinates for the current cell
 * \param sysBoundaries System boundary conditions existing
 *
 * http://fusionwiki.ciemat.es/wiki/Magnetic_curvature
 * 
 * \sa calculateDerivatives calculateBVOLDerivativesSimple calculateDerivativesSimple
 */

void calculateCurvature(
   const arch::buf<VolFsGrid> & volGrid,
   const arch::buf<BgBFsGrid> & bgbGrid,
   const arch::buf<TechnicalFsGrid> & technicalGrid,
   cint i,
   cint j,
   cint k,
   const arch::buf<SysBoundary>& sysBoundaries
) {
   if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY && technicalGrid.get(i,j,k)->sysBoundaryLayer != 1 && technicalGrid.get(i,j,k)->sysBoundaryLayer != 2) {
      auto vol = volGrid.get(i,j,k);
      auto bg = bgbGrid.get(i,j,k);
      auto vol_left_x = volGrid.get(i-1,j,k);
      auto vol_rght_x = volGrid.get(i+1,j,k);
      auto vol_left_y = volGrid.get(i,j-1,k);
      auto vol_rght_y = volGrid.get(i,j+1,k);
      auto vol_left_z = volGrid.get(i,j,k-1);
      auto vol_rght_z = volGrid.get(i,j,k+1);
      auto bg_left_x = bgbGrid.get(i-1,j,k);
      auto bg_rght_x = bgbGrid.get(i+1,j,k);
      auto bg_left_y = bgbGrid.get(i,j-1,k);
      auto bg_rght_y = bgbGrid.get(i,j+1,k);
      auto bg_left_z = bgbGrid.get(i,j,k-1);
      auto bg_rght_z = bgbGrid.get(i,j,k+1);
      
      Real bx = bg[fsgrids::bgbfield::BGBXVOL] + vol[fsgrids::volfields::PERBXVOL];
      Real by = bg[fsgrids::bgbfield::BGBYVOL] + vol[fsgrids::volfields::PERBYVOL];
      Real bz = bg[fsgrids::bgbfield::BGBZVOL] + vol[fsgrids::volfields::PERBZVOL];
      creal bnorm = sqrt(bx*bx + by*by + bz*bz);
      bx /= bnorm;
      by /= bnorm;
      bz /= bnorm;
      Real left_x_bx = bg_left_x[fsgrids::bgbfield::BGBXVOL] + vol_left_x[fsgrids::volfields::PERBXVOL];
      Real left_x_by = bg_left_x[fsgrids::bgbfield::BGBYVOL] + vol_left_x[fsgrids::volfields::PERBYVOL];
      Real left_x_bz = bg_left_x[fsgrids::bgbfield::BGBZVOL] + vol_left_x[fsgrids::volfields::PERBZVOL];
      creal left_x_bnorm = sqrt(left_x_bx*left_x_bx + left_x_by*left_x_by + left_x_bz*left_x_bz);
      left_x_bx /= left_x_bnorm;
      left_x_by /= left_x_bnorm;
      left_x_bz /= left_x_bnorm;
      
      Real rght_x_bx = bg_rght_x[fsgrids::bgbfield::BGBXVOL] + vol_rght_x[fsgrids::volfields::PERBXVOL];
      Real rght_x_by = bg_rght_x[fsgrids::bgbfield::BGBYVOL] + vol_rght_x[fsgrids::volfields::PERBYVOL];
      Real rght_x_bz = bg_rght_x[fsgrids::bgbfield::BGBZVOL] + vol_rght_x[fsgrids::volfields::PERBZVOL];
      creal rght_x_bnorm = sqrt(rght_x_bx*rght_x_bx + rght_x_by*rght_x_by + rght_x_bz*rght_x_bz);
      rght_x_bx /= rght_x_bnorm;
      rght_x_by /= rght_x_bnorm;
      rght_x_bz /= rght_x_bnorm;
      
      Real left_y_bx = bg_left_y[fsgrids::bgbfield::BGBXVOL] + vol_left_y[fsgrids::volfields::PERBXVOL];
      Real left_y_by = bg_left_y[fsgrids::bgbfield::BGBYVOL] + vol_left_y[fsgrids::volfields::PERBYVOL];
      Real left_y_bz = bg_left_y[fsgrids::bgbfield::BGBZVOL] + vol_left_y[fsgrids::volfields::PERBZVOL];
      creal left_y_bnorm = sqrt(left_y_bx*left_y_bx + left_y_by*left_y_by + left_y_bz*left_y_bz);
      left_y_bx /= left_y_bnorm;
      left_y_by /= left_y_bnorm;
      left_y_bz /= left_y_bnorm;
      
      Real rght_y_bx = bg_rght_y[fsgrids::bgbfield::BGBXVOL] + vol_rght_y[fsgrids::volfields::PERBXVOL];
      Real rght_y_by = bg_rght_y[fsgrids::bgbfield::BGBYVOL] + vol_rght_y[fsgrids::volfields::PERBYVOL];
      Real rght_y_bz = bg_rght_y[fsgrids::bgbfield::BGBZVOL] + vol_rght_y[fsgrids::volfields::PERBZVOL];
      creal rght_y_bnorm = sqrt(rght_y_bx*rght_y_bx + rght_y_by*rght_y_by + rght_y_bz*rght_y_bz);
      rght_y_bx /= rght_y_bnorm;
      rght_y_by /= rght_y_bnorm;
      rght_y_bz /= rght_y_bnorm;
      
      Real left_z_bx = bg_left_z[fsgrids::bgbfield::BGBXVOL] + vol_left_z[fsgrids::volfields::PERBXVOL];
      Real left_z_by = bg_left_z[fsgrids::bgbfield::BGBYVOL] + vol_left_z[fsgrids::volfields::PERBYVOL];
      Real left_z_bz = bg_left_z[fsgrids::bgbfield::BGBZVOL] + vol_left_z[fsgrids::volfields::PERBZVOL];
      creal left_z_bnorm = sqrt(left_z_bx*left_z_bx + left_z_by*left_z_by + left_z_bz*left_z_bz);
      left_z_bx /= left_z_bnorm;
      left_z_by /= left_z_bnorm;
      left_z_bz /= left_z_bnorm;
      
      Real rght_z_bx = bg_rght_z[fsgrids::bgbfield::BGBXVOL] + vol_rght_z[fsgrids::volfields::PERBXVOL];
      Real rght_z_by = bg_rght_z[fsgrids::bgbfield::BGBYVOL] + vol_rght_z[fsgrids::volfields::PERBYVOL];
      Real rght_z_bz = bg_rght_z[fsgrids::bgbfield::BGBZVOL] + vol_rght_z[fsgrids::volfields::PERBZVOL];
      creal rght_z_bnorm = sqrt(rght_z_bx*rght_z_bx + rght_z_by*rght_z_by + rght_z_bz*rght_z_bz);
      rght_z_bx /= rght_z_bnorm;
      rght_z_by /= rght_z_bnorm;
      rght_z_bz /= rght_z_bnorm;
      
      vol[fsgrids::volfields::CURVATUREX] = bx * 0.5*(left_x_bx-rght_x_bx) / technicalGrid.grid()->DX + by * 0.5*(left_y_bx-rght_y_bx) / technicalGrid.grid()->DY + bz * 0.5*(left_z_bx-rght_z_bx) / technicalGrid.grid()->DZ;
      vol[fsgrids::volfields::CURVATUREY] = bx * 0.5*(left_x_by-rght_x_by) / technicalGrid.grid()->DX + by * 0.5*(left_y_by-rght_y_by) / technicalGrid.grid()->DY + bz * 0.5*(left_z_by-rght_z_by) / technicalGrid.grid()->DZ;
      vol[fsgrids::volfields::CURVATUREZ] = bx * 0.5*(left_x_bz-rght_x_bz) / technicalGrid.grid()->DX + by * 0.5*(left_y_bz-rght_y_bz) / technicalGrid.grid()->DY + bz * 0.5*(left_z_bz-rght_z_bz) / technicalGrid.grid()->DZ;
   }
}

/*! \brief High-level curvature calculation wrapper function.
 * 
 * \param volGrid fsGrid holding the volume averaged fields
 * \param bgbGrid fsGrid holding the background fields
 * \param technicalGrid fsGrid holding technical information (such as boundary types)
 * \param sysBoundaries System boundary conditions existing
 * 
 * \sa calculateDerivatives calculateBVOLDerivatives calculateDerivativesSimple
 */
void calculateCurvatureSimple(
   arch::buf<VolFsGrid> & volGrid,
   arch::buf<BgBFsGrid> & bgbGrid,
   arch::buf<TechnicalFsGrid> & technicalGrid,
   arch::buf<SysBoundary>& sysBoundaries
) {
   int timer;
   //const std::array<int, 3> gridDims = technicalGrid.getLocalSize();
   const int* gridDims = &technicalGrid.grid()->getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   
   phiprof::start("Calculate curvature");
   
   timer=phiprof::initializeTimer("Start comm","MPI");
   phiprof::start(timer);
   volGrid.syncHostData();
   volGrid.grid()->updateGhostCells();
   // volGrid.syncDeviceData();
   phiprof::stop(timer,N_cells,"Spatial Cells");

   //ARCH_TODO
   // arch::parallel_for({(uint)gridDims[0], (uint)gridDims[1], (uint)gridDims[2]}, ARCH_LOOP_LAMBDA(int i, int j, int k) {
   //       if (technicalGrid.get(i,j,k)->sysBoundaryFlag != sysboundarytype::DO_NOT_COMPUTE) {
   //          calculateCurvature(volGrid,bgbGrid,technicalGrid,i,j,k,sysBoundaries);
   //       }
   //    });
   bgbGrid.syncHostData();
   technicalGrid.syncHostData();
   #pragma omp parallel for collapse(3)
   for (int k=0; k<gridDims[2]; k++) {
      for (int j=0; j<gridDims[1]; j++) {
         for (int i=0; i<gridDims[0]; i++) {
            if (technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
               continue;
            }
            calculateCurvature(volGrid,bgbGrid,technicalGrid,i,j,k,sysBoundaries);
         }
      }
   }
   volGrid.syncDeviceData();
   bgbGrid.syncDeviceData();
   technicalGrid.syncDeviceData();

   phiprof::stop("Calculate curvature",N_cells,"Spatial Cells");
}

/*! \brief Returns perturbed volumetric B of cell
 *
 */
static std::array<Real, 3> getPerB(SpatialCell* cell)
{
   return std::array<Real, 3> { {cell->parameters[CellParams::PERBXVOL], cell->parameters[CellParams::PERBYVOL], cell->parameters[CellParams::PERBZVOL]} };
}

// /*! \brief Returns volumetric E of cell
//  *
//  */
// static std::array<Real, 3> getE(SpatialCell* cell)
// {
//    return std::array<Real, 3> { {cell->parameters[CellParams::EXVOL], cell->parameters[CellParams::EYVOL], cell->parameters[CellParams::EZVOL]} };
// }

// /*! \brief Returns volumetric B of cell
//  *
//  */
// static std::array<Real, 3> getB(SpatialCell* cell)
// {
//    return std::array<Real, 3> {
//       {
//          cell->parameters[CellParams::BGBXVOL] + cell->parameters[CellParams::PERBXVOL],
//          cell->parameters[CellParams::BGBYVOL] + cell->parameters[CellParams::PERBYVOL],
//          cell->parameters[CellParams::BGBZVOL] + cell->parameters[CellParams::PERBZVOL]
//       }
//    };
// }

/*! \brief Calculates momentum density of cell
 *
 */
static std::array<Real, 3> getMomentumDensity(SpatialCell* cell)
{
   Real rho = cell->parameters[CellParams::RHOM];
   return std::array<Real, 3> { {rho * cell->parameters[CellParams::VX], rho * cell->parameters[CellParams::VY], rho * cell->parameters[CellParams::VZ]} };
}

/*! \brief Calculates energy density for spatial cell with only perturbated magnetic field
 *
 */
static Real calculateU1(SpatialCell* cell)
{
   std::array<Real, 3> p = getMomentumDensity(cell);
   std::array<Real, 3> B = getPerB(cell);
   return (pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2)) / (2.0 * cell->parameters[CellParams::RHOM]) + (pow(B[0], 2) + pow(B[1], 2) + pow(B[2], 2)) / (2.0 * physicalconstants::MU_0);
}

/*! \brief Low-level scaled gradients calculation
 * 
 * For the SpatialCell* cell and its neighbors, calculate scaled gradients and their maximum alpha
 * The gradients are the same as in the GUMICS simulation, see
 * Janhunen, P., Palmroth, M., Laitinen, T., Honkonen, I., Juusola, L., Facsko, G., & Pulkkinen, T. I. (2012). The GUMICS-4 global MHD magnetosphere-ionosphere coupling simulation. Journal of Atmospheric and Solar - Terrestrial Physics, 80, 48-59. https://doi.org/10.1016/j.jastp.2012.03.006
 *
 */
void calculateScaledDeltas(
   SpatialCell* cell,
   std::vector<SpatialCell*>& neighbors)
{
   Real dRho {0};
   Real dU {0};
   Real dPsq {0};
   Real dBsq {0};
   Real dB {0};

   Real myRho {cell->parameters[CellParams::RHOM]};
   Real myU {calculateU1(cell)};
   std::array<Real, 3> myP = getMomentumDensity(cell);
   std::array<Real, 3> myB = getPerB(cell);
   for (SpatialCell* neighbor : neighbors) {
      Real otherRho = neighbor->parameters[CellParams::RHOM];
      Real otherU = calculateU1(neighbor);
      std::array<Real, 3> otherP = getMomentumDensity(neighbor);
      std::array<Real, 3> otherB = getPerB(neighbor);
      Real deltaBsq = pow(myB[0] - otherB[0], 2) + pow(myB[1] - otherB[1], 2) + pow(myB[2] - otherB[2], 2);

      // Assignment intentional
      if (Real maxRho = std::max(myRho, otherRho)) {
         dRho = std::max(fabs(myRho - otherRho) / maxRho, dRho);
      }
      if (Real maxU = std::max(myU, otherU)) {
         dU = std::max(fabs(myU - otherU) / maxU, dU);
         dPsq = std::max((pow(myP[0] - otherP[0], 2) + pow(myP[1] - otherP[1], 2) + pow(myP[2] - otherP[2], 2)) / (2 * myRho * maxU), dPsq) / 4.0;
         dBsq = std::max(deltaBsq / (2 * physicalconstants::MU_0 * maxU), dBsq) / 4.0;
      }
      if(Real maxB = sqrt(std::max(pow(myB[0], 2) + pow(myB[1], 2) + pow(myB[2], 2), pow(otherB[0], 2) + pow(otherB[1], 2) + pow(otherB[2], 2)))) {
         dB = std::max(sqrt(deltaBsq) / maxB, dB) / 2.0;
      }
   }
   
   Real alpha = dRho;
   if (dU > alpha) {
      alpha = dU;
   }
   if (dPsq > alpha) {
      alpha = dPsq;
   }
   if (dBsq > alpha) {
      alpha = dBsq;
   }
   if (dB > alpha) {
      alpha = dB;
   }

   Real dBXdy {cell->derivativesBVOL[bvolderivatives::dPERBXVOLdy]};
   Real dBXdz {cell->derivativesBVOL[bvolderivatives::dPERBXVOLdz]};
   Real dBYdx {cell->derivativesBVOL[bvolderivatives::dPERBYVOLdx]};
   Real dBYdz {cell->derivativesBVOL[bvolderivatives::dPERBYVOLdz]};
   Real dBZdx {cell->derivativesBVOL[bvolderivatives::dPERBZVOLdx]};
   Real dBZdy {cell->derivativesBVOL[bvolderivatives::dPERBZVOLdy]};

   // Note missing factor of mu_0, since we want B and J in same units later
   std::array<Real, 3> myJ = {dBZdy - dBYdz, dBXdz - dBZdx, dBYdx - dBXdy};
   Real BdotJ {0.0};
   Real Bsq {0.0};
   for (int i = 0; i < 3; ++i) {
      BdotJ += myB[i] * myJ[i];
      Bsq += myB[i] * myB[i];
   }

   Real Bperp {0.0};
   Real J {0.0};
   for (int i = 0; i < 3; ++i) {
      Bperp += std::pow(myB[i] * (1 - BdotJ / Bsq), 2);
      J += myJ[i] * myJ[i];
   }
   Bperp = std::sqrt(Bperp);
   J = std::sqrt(J);

   cell->parameters[CellParams::AMR_DRHO] = dRho;
   cell->parameters[CellParams::AMR_DU] = dU;
   cell->parameters[CellParams::AMR_DPSQ] = dPsq;
   cell->parameters[CellParams::AMR_DBSQ] = dBsq;
   cell->parameters[CellParams::AMR_DB] = dB;
   cell->parameters[CellParams::AMR_ALPHA] = alpha;
   cell->parameters[CellParams::AMR_JPERB] = J / Bperp;
}

/*! \brief High-level scaled gradient calculation wrapper function.
 * 
 * Calculates gradients needed for alpha everywhere in the grid
 * 
 */

void calculateScaledDeltasSimple(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid)
{
   const vector<CellID>& cells = getLocalCells();
   int N_cells = cells.size();
   int timer;
   phiprof::start("Calculate volume gradients");
   
   timer=phiprof::initializeTimer("Start comm","MPI");
   phiprof::start(timer);

   // We only need nearest neighbourhood and spatial data here
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_SPATIAL_DATA);
   mpiGrid.update_copies_of_remote_neighbors(NEAREST_NEIGHBORHOOD_ID);
   
   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   // Calculate derivatives
   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);

   #pragma omp parallel for
   for (uint i = 0; i < cells.size(); ++i) {
   //for (CellID id : cells) {
      CellID id = cells[i];
      SpatialCell* cell = mpiGrid[id];
      std::vector<SpatialCell*> neighbors;
      for (auto neighPair : mpiGrid.get_face_neighbors_of(id)) {
         neighbors.push_back(mpiGrid[neighPair.first]);
      }
      calculateScaledDeltas(cell, neighbors);
   }

   phiprof::stop(timer,N_cells,"Spatial Cells");

   phiprof::stop("Calculate volume gradients",N_cells,"Spatial Cells");
}
