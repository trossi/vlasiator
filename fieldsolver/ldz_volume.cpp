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
#include "ldz_volume.hpp"

#ifdef DEBUG_VLASIATOR
   #define DEBUG_FSOLVER
#endif

using namespace std;
using index_t = FsGridTools::FsIndex_t;


template<typename SomeFsGrid>
inline static auto get(
   index_t i, index_t j, index_t k,
   SomeFsGrid & grid,
   const TechnicalFsGrid & technicalGrid
) {
   #ifdef DEBUG_FSOLVER
   if (technicalGrid.get(i,j,k) == NULL) {
      stringstream ss;
      ss << "ERROR, got NULL neighbor in " << __FILE__ << ":" << __LINE__ << endl;
      cerr << ss.str(); exit(1);
   }
   #endif
   return grid.get(i,j,k);
}


template<long unsigned int N, typename T>
inline static auto calculateAverage(
   const T & grid_0, const T & grid_1, const T & grid_2, const T & grid_3
) {
   CHECK_FLOAT(grid_0->at(N));
   CHECK_FLOAT(grid_1->at(N));
   CHECK_FLOAT(grid_2->at(N));
   CHECK_FLOAT(grid_3->at(N));
   auto val = FOURTH*(grid_0->at(N) + grid_1->at(N) + grid_2->at(N) + grid_3->at(N));
   CHECK_FLOAT(val);
   return val;
}


void calculateVolumeAveragedFields(
   BFieldFsGrid & perBGrid,
   EFieldFsGrid & EGrid,
   DPerBFsGrid & dPerBGrid,
   VolFsGrid & volGrid,
   TechnicalFsGrid & technicalGrid
) {
   const auto gridDims = &technicalGrid.getLocalSize()[0];
   const size_t N_cells = gridDims[0]*gridDims[1]*gridDims[2];

   phiprof::Timer timer {"Calculate volume averaged fields"};
   int parallelTimerId {phiprof::initializeTimer("volume averaged fields compute cells")};
   #pragma omp parallel
   {
      phiprof::Timer parallelTimer {parallelTimerId};
      #pragma omp for collapse(2)
      for (index_t k=0; k<gridDims[2]; k++) {
         for (index_t j=0; j<gridDims[1]; j++) {
            for (index_t i=0; i<gridDims[0]; i++) {
               std::array<Real, Rec::N_REC_COEFFICIENTS> perturbedCoefficients;
               auto volGrid0 = volGrid.get(i,j,k);

               // Calculate reconstruction coefficients for this cell:
               // This handles domain edges so no need to skip DO_NOT_COMPUTE or OUTER_BOUNDARY_PADDING cells.
               reconstructionCoefficients(
                  perBGrid,
                  dPerBGrid,
                  perturbedCoefficients,
                  i,
                  j,
                  k,
                  2
               );

               // Calculate volume average of B:
               volGrid0->at(fsgrids::volfields::PERBXVOL) = perturbedCoefficients[Rec::a_0];
               volGrid0->at(fsgrids::volfields::PERBYVOL) = perturbedCoefficients[Rec::b_0];
               volGrid0->at(fsgrids::volfields::PERBZVOL) = perturbedCoefficients[Rec::c_0];

               // This avoids out of domain accesses below.
               if(technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE || technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::OUTER_BOUNDARY_PADDING) continue;

               // Calculate volume average of E
               if ( technicalGrid.get(i,j,k)->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY ||
                    (technicalGrid.get(i,j,k)->sysBoundaryFlag != sysboundarytype::NOT_SYSBOUNDARY && technicalGrid.get(i,j,k)->sysBoundaryLayer == 1)
                  ) {
                  auto EGrid_i1j1k1 = get(i  ,j  ,k  , EGrid, technicalGrid);
                  auto EGrid_i1j1k2 = get(i  ,j  ,k+1, EGrid, technicalGrid);
                  auto EGrid_i1j2k1 = get(i  ,j+1,k  , EGrid, technicalGrid);
                  auto EGrid_i1j2k2 = get(i  ,j+1,k+1, EGrid, technicalGrid);
                  auto EGrid_i2j1k1 = get(i+1,j  ,k  , EGrid, technicalGrid);
                  auto EGrid_i2j1k2 = get(i+1,j  ,k+1, EGrid, technicalGrid);
                  auto EGrid_i2j2k1 = get(i+1,j+1,k  , EGrid, technicalGrid);

                  volGrid0->at(fsgrids::volfields::EXVOL) =
                     calculateAverage<fsgrids::efield::EX>(
                        EGrid_i1j1k1, EGrid_i1j2k1, EGrid_i1j1k2, EGrid_i1j2k2
                     );
                  volGrid0->at(fsgrids::volfields::EYVOL) =
                     calculateAverage<fsgrids::efield::EY>(
                        EGrid_i1j1k1, EGrid_i2j1k1, EGrid_i1j1k2, EGrid_i2j1k2
                     );
                  volGrid0->at(fsgrids::volfields::EZVOL) =
                     calculateAverage<fsgrids::efield::EZ>(
                        EGrid_i1j1k1, EGrid_i2j1k1, EGrid_i1j2k1, EGrid_i2j2k1
                     );
               } else {
                  volGrid0->at(fsgrids::volfields::EXVOL) = 0.0;
                  volGrid0->at(fsgrids::volfields::EYVOL) = 0.0;
                  volGrid0->at(fsgrids::volfields::EZVOL) = 0.0;
               }
            }
         }
      }
   }

   timer.stop(N_cells, "Spatial Cells");
}
