#ifndef CUDAGRID_CUH
#define CUDAGRID_CUH

#include "../../src/pyGrid.h"
#include "../../src/pyMol2.h"

void invoke_compute_grid_softcore_HB_omp(Grid* grid, Mol2* rec); 
void print_difference(const std::vector<std::vector<std::vector<double>>>& vec, double* arr);

#endif