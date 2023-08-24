#ifndef CUDAGRID_CUH
#define CUDAGRID_CUH

#include "pyGrid.h"
#include "pyMol2.h"

void invoke_compute_grid_softcore_HB_CUDA(Grid* grid, Mol2* rec); 

void invoke_compute_grid_hardcore_HB_CUDA(Grid* grid, Mol2* rec); 

#endif