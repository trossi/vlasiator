/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
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

#ifndef DERIVATIVES_HPP
#define DERIVATIVES_HPP

#include <vector>
#include <array>

#include "../definitions.h"
#include "../spatial_cell.hpp"
#include "../sysboundary/sysboundary.h"

#include "fs_limiters.h"

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
   const bool communicateMoments);

void calculateBVOLDerivativesSimple(
   arch::buf<VolFsGrid> & volGrid,
   arch::buf<TechnicalFsGrid> & technicalGrid,
   arch::buf<SysBoundary>& sysBoundaries
);

void calculateCurvatureSimple(
   arch::buf<VolFsGrid> & volGrid,
   arch::buf<BgBFsGrid> & bgbGrid,
   arch::buf<TechnicalFsGrid> & technicalGrid,
   arch::buf<SysBoundary>& sysBoundaries
);


void calculateScaledDeltasSimple(dccrg::Dccrg<SpatialCell,dccrg::Cartesian_Geometry>& mpiGrid);

#endif
