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

#include "fs_common.h"
#include "ldz_gradpe.hpp"
#include "../arch/arch_sysboundary_api.h"

#ifdef DEBUG_VLASIATOR
   #define DEBUG_FSOLVER
#endif

using namespace std;

void calculateEdgeGradPeTermXComponents(
   const arch::buf<EGradPeFsGrid> & EGradPeGrid,
   const arch::buf<MomentsFsGrid> & momentsGrid,
   const arch::buf<DMomentsFsGrid> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
   Real hallRhoq = 0.0;
   Real rhoq = 0.0;
   switch (FSParams.ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if FSParams.ohmGradPeTerm == 0." << endl;
         break;
         
      case 1:
         rhoq = momentsGrid.get(i,j,k)[fsgrids::moments::RHOQ];
         hallRhoq = (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;
         //EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EXGRADPE] = -physicalconstants::K_B*Parameters::electronTemperature*dMomentsGrid.get(i,j,k)[fsgrids::dmoments::drhoqdx] / (hallRhoq*EGradPeGrid.DX);
         EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EXGRADPE] = - dMomentsGrid.get(i,j,k)[fsgrids::dmoments::dPedx] / (hallRhoq*EGradPeGrid.grid()->DX);
	 break;
         
      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

void calculateEdgeGradPeTermYComponents(
   const arch::buf<EGradPeFsGrid> & EGradPeGrid,
   const arch::buf<MomentsFsGrid> & momentsGrid,
   const arch::buf<DMomentsFsGrid> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
   Real hallRhoq = 0.0;
   Real rhoq = 0.0;
   switch (FSParams.ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if FSParams.ohmGradPeTerm == 0." << endl;
         break;
         
      case 1:
         rhoq = momentsGrid.get(i,j,k)[fsgrids::moments::RHOQ];
         hallRhoq = (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;
         //EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EYGRADPE] = -physicalconstants::K_B*Parameters::electronTemperature*dMomentsGrid.get(i,j,k)[fsgrids::dmoments::drhoqdy] / (hallRhoq*EGradPeGrid.DY);
         EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EYGRADPE] = - dMomentsGrid.get(i,j,k)[fsgrids::dmoments::dPedy] / (hallRhoq*EGradPeGrid.grid()->DY);
         break;
         
      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

void calculateEdgeGradPeTermZComponents(
   const arch::buf<EGradPeFsGrid> & EGradPeGrid,
   const arch::buf<MomentsFsGrid> & momentsGrid,
   const arch::buf<DMomentsFsGrid> & dMomentsGrid,
   cint i,
   cint j,
   cint k
) {
   Real hallRhoq = 0.0;
   Real rhoq = 0.0;
   switch (FSParams.ohmGradPeTerm) {
      case 0:
         cerr << __FILE__ << __LINE__ << "You shouldn't be in a electron pressure gradient term function if FSParams.ohmGradPeTerm == 0." << endl;
         break;
         
      case 1:
         rhoq = momentsGrid.get(i,j,k)[fsgrids::moments::RHOQ];
         hallRhoq = (rhoq <= Parameters::hallMinimumRhoq ) ? Parameters::hallMinimumRhoq : rhoq ;
         //EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EZGRADPE] = -physicalconstants::K_B*Parameters::electronTemperature*dMomentsGrid.get(i,j,k)[fsgrids::dmoments::drhoqdz] / (hallRhoq*EGradPeGrid.DZ);
         EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EZGRADPE] = - dMomentsGrid.get(i,j,k)[fsgrids::dmoments::dPedz] / (hallRhoq*EGradPeGrid.grid()->DZ);
         break;
         
      default:
         cerr << __FILE__ << ":" << __LINE__ << "You are welcome to code higher-order Hall term correction terms." << endl;
         break;
   }
}

/** Calculate the electron pressure gradient term on all given cells.
 * @param sysBoundaries System boundary condition functions.
 */
void calculateGradPeTerm(
   const arch::buf<EGradPeFsGrid> & EGradPeGrid,
   const arch::buf<MomentsFsGrid> & momentsGrid,
   const arch::buf<DMomentsFsGrid> & dMomentsGrid,
   const arch::buf<TechnicalFsGrid> & technicalGrid,
   cint i,
   cint j,
   cint k,
   const arch::buf<SysBoundary>& sysBoundaries
) {
   #ifdef DEBUG_FSOLVER
   if (technicalGrid.get(i,j,k) == NULL) {
      cerr << "NULL pointer in " << __FILE__ << ":" << __LINE__ << endl;
      exit(1);
   }
   #endif
   
   cuint cellSysBoundaryFlag = technicalGrid.get(i,j,k)->sysBoundaryFlag;
   
   if (cellSysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) return;
   
   cuint cellSysBoundaryLayer = technicalGrid.get(i,j,k)->sysBoundaryLayer;
   
   if ((cellSysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY) && (cellSysBoundaryLayer != 1)) {
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag).fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid,i,j,k,0);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag).fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid,i,j,k,1);
      sysBoundaries.getSysBoundary(cellSysBoundaryFlag).fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid,i,j,k,2);
   } else {
      calculateEdgeGradPeTermXComponents(EGradPeGrid,momentsGrid,dMomentsGrid,i,j,k);
      calculateEdgeGradPeTermYComponents(EGradPeGrid,momentsGrid,dMomentsGrid,i,j,k);
      calculateEdgeGradPeTermZComponents(EGradPeGrid,momentsGrid,dMomentsGrid,i,j,k);
   }
}

void calculateGradPeTermSimple(
   arch::buf<EGradPeFsGrid> & EGradPeGrid,
   arch::buf<MomentsFsGrid> & momentsGrid,
   arch::buf<MomentsFsGrid> & momentsDt2Grid,
   arch::buf<DMomentsFsGrid> & dMomentsGrid,
   arch::buf<TechnicalFsGrid> & technicalGrid,
   arch::buf<SysBoundary>& sysBoundaries,
   cint& RKCase
) {
   int timer;
   const auto gridDims = &technicalGrid.grid()->getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];
   phiprof::start("Calculate GradPe term");

   timer=phiprof::initializeTimer("MPI","MPI");
   phiprof::start(timer);
   dMomentsGrid.syncHostData();
   dMomentsGrid.grid()->updateGhostCells();
   dMomentsGrid.syncDeviceData();
   phiprof::stop(timer);

   // Calculate GradPe term
   timer=phiprof::initializeTimer("Compute cells");
   phiprof::start(timer);
   arch::parallel_for({(uint)gridDims[0], (uint)gridDims[1], (uint)gridDims[2]}, ARCH_LOOP_LAMBDA(int i, int j, int k) { 
      if (RKCase == RK_ORDER1 || RKCase == RK_ORDER2_STEP2) {
         calculateGradPeTerm(EGradPeGrid, momentsGrid, dMomentsGrid, technicalGrid, i, j, k, sysBoundaries);
      } else {
         calculateGradPeTerm(EGradPeGrid, momentsDt2Grid, dMomentsGrid, technicalGrid, i, j, k, sysBoundaries);
      }
   });
   phiprof::stop(timer,N_cells,"Spatial Cells");
   
   phiprof::stop("Calculate GradPe term",N_cells,"Spatial Cells");
}
