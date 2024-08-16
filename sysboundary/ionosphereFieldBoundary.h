#include "../common.h"
#include "../fieldsolver/fs_common.h"
#include "../parameters.h"
#include "fieldBoundary.h"
#include <iostream>

using namespace std;

namespace SBC {
    class IonosphereFieldBoundary : FieldBoundary {

    public:
    IonosphereFieldBoundary(Real center[3], Real radius, uint geometry) {
        this->center[0] = center[0];
        this->center[1] = center[1];
        this->center[2] = center[2];
        this->radius = radius;
        this->geometry = geometry;
    }
 
    /*! We want here to
    * 
    * -- Average perturbed face B from the nearest neighbours
    */
    ARCH_HOSTDEV Real fieldSolverBoundaryCondMagneticField(
        const arch::buf<BFieldFsGrid> & bGrid,
        const arch::buf<TechnicalFsGrid> & technicalGrid,
        cint i,
        cint j,
        cint k,
        creal& dt,
        cuint& component
    ) {
        if (technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) {
            switch(component) {
                case 0:
                if (  ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BX) == compute::BX)
                    && ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BX) == compute::BX)
                ) {
                    return 0.5 * (bGrid.get(i-1,j,k)[fsgrids::bfield::PERBX] + bGrid.get(i+1,j,k)[fsgrids::bfield::PERBX]);
                } else if ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BX) == compute::BX) {
                    return bGrid.get(i-1,j,k)[fsgrids::bfield::PERBX];
                } else if ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BX) == compute::BX) {
                    return bGrid.get(i+1,j,k)[fsgrids::bfield::PERBX];
                } else {
                    Real retval = 0.0;
                    uint nCells = 0;
                    if ((technicalGrid.get(i,j-1,k)->SOLVE & compute::BX) == compute::BX) {
                        retval += bGrid.get(i,j-1,k)[fsgrids::bfield::PERBX];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j+1,k)->SOLVE & compute::BX) == compute::BX) {
                        retval += bGrid.get(i,j+1,k)[fsgrids::bfield::PERBX];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j,k-1)->SOLVE & compute::BX) == compute::BX) {
                        retval += bGrid.get(i,j,k-1)[fsgrids::bfield::PERBX];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j,k+1)->SOLVE & compute::BX) == compute::BX) {
                        retval += bGrid.get(i,j,k+1)[fsgrids::bfield::PERBX];
                        nCells++;
                    }
                    if (nCells == 0) {
                        for (int a=i-1; a<i+2; a++) {
                            for (int b=j-1; b<j+2; b++) {
                            for (int c=k-1; c<k+2; c++) {
                                if ((technicalGrid.get(a,b,c)->SOLVE & compute::BX) == compute::BX) {
                                    retval += bGrid.get(a,b,c)[fsgrids::bfield::PERBX];
                                    nCells++;
                                }
                            }
                            }
                        }
                    }
                    if (nCells == 0) {
                        #ifndef __CUDA_ARCH__
                        cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                        #endif
                        return 0.0;
                    }
                    return retval / nCells;
                }
                case 1:
                if (  (technicalGrid.get(i,j-1,k)->SOLVE & compute::BY) == compute::BY
                    && (technicalGrid.get(i,j+1,k)->SOLVE & compute::BY) == compute::BY
                ) {
                    return 0.5 * (bGrid.get(i,j-1,k)[fsgrids::bfield::PERBY] + bGrid.get(i,j+1,k)[fsgrids::bfield::PERBY]);
                } else if ((technicalGrid.get(i,j-1,k)->SOLVE & compute::BY) == compute::BY) {
                    return bGrid.get(i,j-1,k)[fsgrids::bfield::PERBY];
                } else if ((technicalGrid.get(i,j+1,k)->SOLVE & compute::BY) == compute::BY) {
                    return bGrid.get(i,j+1,k)[fsgrids::bfield::PERBY];
                } else {
                    Real retval = 0.0;
                    uint nCells = 0;
                    if ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BY) == compute::BY) {
                        retval += bGrid.get(i-1,j,k)[fsgrids::bfield::PERBY];
                        nCells++;
                    }
                    if ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BY) == compute::BY) {
                        retval += bGrid.get(i+1,j,k)[fsgrids::bfield::PERBY];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j,k-1)->SOLVE & compute::BY) == compute::BY) {
                        retval += bGrid.get(i,j,k-1)[fsgrids::bfield::PERBY];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j,k+1)->SOLVE & compute::BY) == compute::BY) {
                        retval += bGrid.get(i,j,k+1)[fsgrids::bfield::PERBY];
                        nCells++;
                    }
                    if (nCells == 0) {
                        for (int a=i-1; a<i+2; a++) {
                            for (int b=j-1; b<j+2; b++) {
                            for (int c=k-1; c<k+2; c++) {
                                if ((technicalGrid.get(a,b,c)->SOLVE & compute::BY) == compute::BY) {
                                    retval += bGrid.get(a,b,c)[fsgrids::bfield::PERBY];
                                    nCells++;
                                }
                            }
                            }
                        }
                    }
                    if (nCells == 0) {
                        #ifndef __CUDA_ARCH__
                        cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                        #endif
                        return 0.0;
                    }
                    return retval / nCells;
                }
                case 2:
                if (  (technicalGrid.get(i,j,k-1)->SOLVE & compute::BZ) == compute::BZ
                    && (technicalGrid.get(i,j,k+1)->SOLVE & compute::BZ) == compute::BZ
                ) {
                    return 0.5 * (bGrid.get(i,j,k-1)[fsgrids::bfield::PERBZ] + bGrid.get(i,j,k+1)[fsgrids::bfield::PERBZ]);
                } else if ((technicalGrid.get(i,j,k-1)->SOLVE & compute::BZ) == compute::BZ) {
                    return bGrid.get(i,j,k-1)[fsgrids::bfield::PERBZ];
                } else if ((technicalGrid.get(i,j,k+1)->SOLVE & compute::BZ) == compute::BZ) {
                    return bGrid.get(i,j,k+1)[fsgrids::bfield::PERBZ];
                } else {
                    Real retval = 0.0;
                    uint nCells = 0;
                    if ((technicalGrid.get(i-1,j,k)->SOLVE & compute::BZ) == compute::BZ) {
                        retval += bGrid.get(i-1,j,k)[fsgrids::bfield::PERBZ];
                        nCells++;
                    }
                    if ((technicalGrid.get(i+1,j,k)->SOLVE & compute::BZ) == compute::BZ) {
                        retval += bGrid.get(i+1,j,k)[fsgrids::bfield::PERBZ];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j-1,k)->SOLVE & compute::BZ) == compute::BZ) {
                        retval += bGrid.get(i,j-1,k)[fsgrids::bfield::PERBZ];
                        nCells++;
                    }
                    if ((technicalGrid.get(i,j+1,k)->SOLVE & compute::BZ) == compute::BZ) {
                        retval += bGrid.get(i,j+1,k)[fsgrids::bfield::PERBZ];
                        nCells++;
                    }
                    if (nCells == 0) {
                        for (int a=i-1; a<i+2; a++) {
                            for (int b=j-1; b<j+2; b++) {
                            for (int c=k-1; c<k+2; c++) {
                                if ((technicalGrid.get(a,b,c)->SOLVE & compute::BZ) == compute::BZ) {
                                    retval += bGrid.get(a,b,c)[fsgrids::bfield::PERBZ];
                                    nCells++;
                                }
                            }
                            }
                        }
                    }
                    if (nCells == 0) {
                        #ifndef __CUDA_ARCH__ 
                        cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                        #endif 
                        return 0.0;
                    }
                    return retval / nCells;
                }
                default:
                #ifndef __CUDA_ARCH__
                cerr << "ERROR: ionosphere boundary tried to copy nonsensical magnetic field component " << component << endl;
                #endif 
                return 0.0;
            }
        } else { // L2 cells
            Real retval = 0.0;
            uint nCells = 0;
            for (int a=i-1; a<i+2; a++) {
                for (int b=j-1; b<j+2; b++) {
                for (int c=k-1; c<k+2; c++) {
                    if (technicalGrid.get(a,b,c)->sysBoundaryLayer == 1) {
                        retval += bGrid.get(a,b,c)[fsgrids::bfield::PERBX + component];
                        nCells++;
                    }
                }
                }
            }
            if (nCells == 0) {
                #ifndef __CUDA_ARCH__
                cerr << __FILE__ << ":" << __LINE__ << ": ERROR: this should not have fallen through." << endl;
                #endif 
                return 0.0;
            }
            return retval / nCells;
        }
    }

    /*! We want here to
        *
        * -- Retain only the boundary-normal projection of perturbed face B
        */
    ARCH_HOSTDEV void fieldSolverBoundaryCondMagneticFieldProjection(
        const arch::buf<BFieldFsGrid> & bGrid,
        const arch::buf<TechnicalFsGrid> & technicalGrid,
        cint i,
        cint j,
        cint k
    ) {
        // Projection of B-field to normal direction
        Real BdotN = 0;
        Real normalDirection[3]; 
        fieldSolverGetNormalDirection(technicalGrid, i, j, k, normalDirection);
        for(uint component=0; component<3; component++) {
            BdotN += bGrid.get(i,j,k)[fsgrids::bfield::PERBX+component] * normalDirection[component];
        }
        // Apply to any components that were not solved
        if ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 2) ||
            ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k)->SOLVE & compute::BX) != compute::BX))
            ) {
            bGrid.get(i,j,k)[fsgrids::bfield::PERBX] = BdotN*normalDirection[0];
        }
        if ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 2) ||
            ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k)->SOLVE & compute::BY) != compute::BY))
            ) {
            bGrid.get(i,j,k)[fsgrids::bfield::PERBY] = BdotN*normalDirection[1];
        }
        if ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 2) ||
            ((technicalGrid.get(i,j,k)->sysBoundaryLayer == 1) && ((technicalGrid.get(i,j,k)->SOLVE & compute::BZ) != compute::BZ))
            ) {
            bGrid.get(i,j,k)[fsgrids::bfield::PERBZ] = BdotN*normalDirection[2];
        }
    }

    ARCH_HOSTDEV void fieldSolverBoundaryCondElectricField(
      const arch::buf<EFieldFsGrid> & EGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGrid.get(i,j,k)[fsgrids::efield::EX+component] = 0.0;
   }
   
   ARCH_HOSTDEV void fieldSolverBoundaryCondHallElectricField(
      const arch::buf<EHallFsGrid> & EHallGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      auto cp = EHallGrid.get(i,j,k);
      switch (component) {
         case 0:
            cp[fsgrids::ehall::EXHALL_000_100] = 0.0;
            cp[fsgrids::ehall::EXHALL_010_110] = 0.0;
            cp[fsgrids::ehall::EXHALL_001_101] = 0.0;
            cp[fsgrids::ehall::EXHALL_011_111] = 0.0;
            break;
         case 1:
            cp[fsgrids::ehall::EYHALL_000_010] = 0.0;
            cp[fsgrids::ehall::EYHALL_100_110] = 0.0;
            cp[fsgrids::ehall::EYHALL_001_011] = 0.0;
            cp[fsgrids::ehall::EYHALL_101_111] = 0.0;
            break;
         case 2:
            cp[fsgrids::ehall::EZHALL_000_001] = 0.0;
            cp[fsgrids::ehall::EZHALL_100_101] = 0.0;
            cp[fsgrids::ehall::EZHALL_010_011] = 0.0;
            cp[fsgrids::ehall::EZHALL_110_111] = 0.0;
            break;
         default:
            #ifndef __CUDA_ARCH__
            cerr << __FILE__ << ":" << __LINE__ << ":" << " Invalid component" << endl;
            #endif
            return;
      }
   }
   
   ARCH_HOSTDEV void fieldSolverBoundaryCondGradPeElectricField(
      const arch::buf<EGradPeFsGrid> & EGradPeGrid,
      cint i,
      cint j,
      cint k,
      cuint component
   ) {
      EGradPeGrid.get(i,j,k)[fsgrids::egradpe::EXGRADPE+component] = 0.0;
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
      SysBoundaryCondition::setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
      return;
   }
   
   ARCH_HOSTDEV void fieldSolverBoundaryCondBVOLDerivatives(
      const arch::buf<VolFsGrid> & volGrid,
      cint i,
      cint j,
      cint k,
      cuint& component
   ) {
      // FIXME This should be OK as the BVOL derivatives are only used for Lorentz force JXB, which is not applied on the ionosphere cells.
      SysBoundaryCondition::setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
   } 

  }; 

}
