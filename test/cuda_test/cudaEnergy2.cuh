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

#endif CUDA_ENERGY2_CUH
