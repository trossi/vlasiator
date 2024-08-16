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

#ifndef LDZ_MAGNETIC_FIELD_HPP
#define LDZ_MAGNETIC_FIELD_HPP

#include <vector>

#include "../definitions.h"
#include "../common.h"
#include "../spatial_cell.hpp"
#include "../arch/arch_sysboundary_api.h"

#include "fs_common.h"

void propagateMagneticField(
   arch::buf<BFieldFsGrid> & perBGrid,
   arch::buf<BFieldFsGrid> & perBDt2Grid,
   arch::buf<EFieldFsGrid> & EGrid,
   arch::buf<EFieldFsGrid> & EDt2Grid,
   cint i,
   cint j,
   cint k,
   creal& dt,
   cint& RKCase,
   const bool doX=true,
   const bool doY=true,
   const bool doZ=true
);

void propagateMagneticFieldSimple(
   arch::buf<BFieldFsGrid> & perBGrid,
   arch::buf<BFieldFsGrid> & perBDt2Grid,
   arch::buf<EFieldFsGrid> & EGrid,
   arch::buf<EFieldFsGrid> & EDt2Grid,
   arch::buf<TechnicalFsGrid> & technicalGrid,
   arch::buf<SysBoundary> & sysBoundaries, 
   creal& dt,
   cint& RKCase
);

#endif
