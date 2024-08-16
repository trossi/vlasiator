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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#ifndef GRID_H
#define GRID_H

#include "definitions.h"
#include "spatial_cell_wrapper.hpp"
#include <dccrg.hpp>
#include <dccrg_cartesian_geometry.hpp>
#include "sysboundary/sysboundary.h"
#include "projects/project.h"
#include <string>

/*!
  \brief Initialize DCCRG and fsgrids
*/
void initializeGrids(
   int argn,
   char **argc,
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   BFieldFsGrid & perBGrid,
   BgBFsGrid & BgBGrid,
   MomentsFsGrid & momentsGrid,
   MomentsFsGrid & momentsDt2Grid,
   EFieldFsGrid & EGrid,
   EGradPeFsGrid & EGradPeGrid,
   VolFsGrid & volGrid,
   TechnicalFsGrid & technicalGrid,
   SysBoundary& sysBoundaries,
   Project& project
);

/*!
  \brief Balance load

    \param[in,out] mpiGrid The DCCRG grid with spatial cells
*/
void balanceLoad(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, SysBoundary& sysBoundaries);

/*!

Updates velocity block lists between remote neighbors and
prepares local copies of remote neighbors to receive velocity block
data. This is needed if one has locally adjusted velocity blocks

\param mpiGrid   The DCCRG grid with spatial cells
*/
void updateRemoteVelocityBlockLists(
   dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
   const uint popID,
   const uint neighborhood=DIST_FUNC_NEIGHBORHOOD_ID
);

/*! Deallocates all blocks in remote cells in order to save
 *  memory. 
 * \param mpiGrid Spatial grid
 */
void deallocateRemoteCellBlocks(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

/*! Adjust sparse velocity space to make it consistent in all 6 dimensions

 1) Compute which blocks have content (done for all cells in mpiGrid)
 2) Adjust local velocity blocks. That is, make sure blocks exist which have content, or have
 neighbors with content in all 6-dimensions. This is done for cells in cellsToAdjust list.
 3) Make sure remote cells are up-to-date and ready to receive data, if doPrepareToReceiveBlocks is true.

 Note that block existence does not use vlasov stencil as it is important to also include diagonals to avoid massloss

 \param mpiGrid  Parallel grid with spatial cells
 \param cellsToAdjust  List of cells that are adjusted, that is cells which blocks are added or removed. 
 \param doPrepareToReceiveBlocks If true, then remote cells are set up so that velocity space data can be received. Global operation, value has to be the same for all processes.
 
*/
bool adjustVelocityBlocks(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid,
                          const std::vector<CellID>& cellsToAdjust,
                          bool doPrepareToReceiveBlocks,
                            const uint popID);

/*! Estimates memory consumption and writes it into logfile. Collective operation on MPI_COMM_WORLD
 * \param mpiGrid Spatial grid
 */
void report_grid_memory_consumption(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

/*! Shrink to fit velocity space data to save memory.
 * \param mpiGrid Spatial grid
 */
void shrink_to_fit_grid_data(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

void setFaceNeighborRanks( dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid );

/*! Map grid refinement to FsGrid
 */
void mapRefinement(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, TechnicalFsGrid & technicalGrid);

/*! Refine spatial cells and update necessary information
 * \param mpiGrid Spatial grid
 * \param technicalGrid Technical grid
 * \param sysBoundaries System boundaries
 * \param project Project used
 * \param useStatic Used for forcing static refinement on restart. Negative values use adaptive refinement, non-negative values correspond to static refinement pass in Project::forceRefinement
 */
bool adaptRefinement(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid, TechnicalFsGrid & technicalGrid, SysBoundary& sysBoundaries, Project& project, int useStatic = -1);

#endif
