#ifndef CUDA_ENERGY2_CUH
#define CUDA_ENERGY2_CUH

typedef struct st_DeviceGridInterpol {

  double vdwA;
  double vdwB;
  double elec;
  double solv_gauss;
  double rec_solv_gauss;
  double hb_donor;
  double hb_acceptor;

} DeviceGridInterpol;

double invoke_compute_ene_from_grids_softcore_solvation(Grid* grid, Mol2 *lig, 
                                                 const std::vector<std::vector<double>>& xyz);

#endif 
