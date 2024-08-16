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

#ifndef OUTFLOW_H
#define OUTFLOW_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"
#include "outflowFieldBoundary.h"

namespace SBC {

	struct OutflowSpeciesParameters {
      /*! Array of bool telling which faces are going to be skipped by the Vlasov system boundary condition.*/
			std::array<bool, 6> facesToSkipVlasov;
      /*! List of schemes to use for the Vlasov outflow boundary conditions on each face ([xyz][+-]). */
      std::array<uint, 6> faceVlasovScheme;
      /*! List of faces on which outflow boundary conditions are to be reapplied upon restart ([xyz][+-]). */
      std::vector<std::string> faceToReapplyUponRestartList;

      /*! Factor by which to quench the inflowing parts of the velocity distribution function.*/
      Real quenchFactor;
	};

   /*!\brief Outflow is a class applying copy/outflow boundary conditions.
    * 
    * Outflow is a class handling cells tagged as sysboundarytype::OUTFLOW by this system boundary condition. It applies copy/outflow boundary conditions.
    * 
    * These consist in:
    * - Copy the distribution and moments from the nearest NOT_SYSBOUNDARY cell;
    * - Copy the perturbed B components from the nearest NOT_SYSBOUNDARY cell. EXCEPTION: the face components adjacent to the simulation domain at the +x/+y/+z faces are propagated still.
    */
   class Outflow: public OuterBoundaryCondition {
   public:
      Outflow();
      virtual ~Outflow();
      
      static void addParameters();
      virtual void getParameters();
      
      virtual bool initSysBoundary(
         creal& t,
         Project &project
      );
      bool initFieldBoundary();
      OutflowFieldBoundary* getFieldBoundary() {return fieldBoundary;} 
      virtual bool applyInitialState(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         TechnicalFsGrid & technicalGrid,
         BFieldFsGrid & perBGrid,
         Project &project
      );
      ARCH_HOSTDEV Real fieldSolverBoundaryCondMagneticField(
         const arch::buf<BFieldFsGrid> & bGrid,
         const arch::buf<TechnicalFsGrid> & technicalGrid,
         cint i,
         cint j,
         cint k,
         creal& dt,
         cuint& component
      ) {
         return fieldBoundary->fieldSolverBoundaryCondMagneticField(bGrid, technicalGrid, i, j, k, dt, component);
      }
      ARCH_HOSTDEV void fieldSolverBoundaryCondMagneticFieldProjection(
         const arch::buf<BFieldFsGrid> & bGrid,
         const arch::buf<TechnicalFsGrid> & technicalGrid,
         cint i,
         cint j,
         cint k
      ) {
         fieldBoundary->fieldSolverBoundaryCondMagneticFieldProjection(bGrid, technicalGrid, i, j, k);
      }
      ARCH_HOSTDEV void fieldSolverBoundaryCondElectricField(
         const arch::buf<EFieldFsGrid> & EGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) {
         fieldBoundary->fieldSolverBoundaryCondElectricField(EGrid, i, j, k, component);
      }
      ARCH_HOSTDEV void fieldSolverBoundaryCondHallElectricField(
         const arch::buf<EHallFsGrid> & EHallGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) {
         fieldBoundary->fieldSolverBoundaryCondHallElectricField(EHallGrid, i, j, k, component);
      }
      ARCH_HOSTDEV void fieldSolverBoundaryCondGradPeElectricField(
         const arch::buf<EGradPeFsGrid> & EGradPeGrid,
         cint i,
         cint j,
         cint k,
         cuint component
      ) {
         fieldBoundary->fieldSolverBoundaryCondGradPeElectricField(EGradPeGrid, i, j, k, component);
      }
      ARCH_HOSTDEV void fieldSolverBoundaryCondDerivatives(
         const arch::buf<DPerBFsGrid> & dPerBGrid,
         const arch::buf<DMomentsFsGrid> & dMomentsGrid,
         cint i,
         cint j,
         cint k,
         cuint& RKCase,
         cuint& component
      ) {
         fieldBoundary->fieldSolverBoundaryCondDerivatives(dPerBGrid, dMomentsGrid, i, j, k, RKCase, component);
      }
      ARCH_HOSTDEV void fieldSolverBoundaryCondBVOLDerivatives(
         const arch::buf<VolFsGrid> & volGrid,
         cint i,
         cint j,
         cint k,
         cuint& component
      ) {
         fieldBoundary->fieldSolverBoundaryCondBVOLDerivatives(volGrid, i, j, k, component);
      }
      virtual void vlasovBoundaryCondition(
         const dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
         const CellID& cellID,
         const uint popID,
         const bool calculate_V_moments
      );
      
      virtual void getFaces(bool* faces);
      virtual std::string getName() const;
      virtual uint getIndex() const;
      
   protected:
      /*! Array of bool telling which faces are going to be processed by the fields system boundary condition.*/
      bool facesToSkipFields[6];
      /*! Array of bool telling which faces are going to be reapplied upon restart.*/
      bool facesToReapply[6];
      /*! List of faces on which outflow boundary conditions are to be applied ([xyz][+-]). */
      std::vector<std::string> faceList;
      /*! List of faces on which no fields outflow boundary conditions are to be applied ([xyz][+-]). */
      std::vector<std::string> faceNoFieldsList;
      std::vector<OutflowSpeciesParameters> speciesParams;
      
      /*! Factor by which to quench the inflowing parts of the velocity distribution function.*/
      Real quenchFactor;
      
      enum vlasovscheme {
         NONE,
         COPY,
         LIMIT,
         N_SCHEMES
      };

      OutflowFieldBoundary* fieldBoundary;
      
   }; // class Outflow
} // namespace SBC

#endif
