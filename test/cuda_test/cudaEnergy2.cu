#include "../../src/pyEnergy2.h"
#include "cudaEnergy2.cuh"
#include <__clang_cuda_builtin_vars.h>

__device__ double trilinear_interpolation(double *grid, int x0, int y0, int z0,
                                          int x1, int y1, int z1, double xd,
                                          double yd, double zd, int npointsx,
                                          int npointsy) {

  double c00 =
      (grid[(x0 * npointsx * npointsy) + (npointsx * y0) + z0] * (1 - xd)) +
      (grid[(x1 * npointsx * npointsy) + (npointsx * y0) + z0] * xd);
  double c10 =
      (grid[(x0 * npointsx * npointsy) + (npointsx * y1) + z0] * (1 - xd)) +
      (grid[(x1 * npointsx * npointsy) + (npointsx * y1) + z0] * xd);
  double c01 =
      (grid[(x0 * npointsx * npointsy) + (npointsx * y0) + z1] * (1 - xd)) +
      (grid[(x1 * npointsx * npointsy) + (npointsx * y0) + z1] * xd);
  double c11 =
      (grid[(x0 * npointsx * npointsy) + (npointsx * y1) + z1] * (1 - xd)) +
      (grid[(x1 * npointsx * npointsy) + (npointsx * y1) + z1] * xd);

  double c0 = (c00 * (1 - yd)) + (c10 * yd);
  double c1 = (c01 * (1 - yd)) + (c11 * yd);

  return (c0 * (1 - zd)) + (c1 * zd);
}

__device__ DeviceGridInterpol grids_trilinear_interpolation(
    double *d_elec_grid, double *d_vdwA, double *d_vdwB,
    double *d_hb_donor_grid, double *d_acceptor_grid, double *d_rec_solv_gauss,
    double *d_solv_gauss, double x, double y, double z, int x0, int y0, int z0,
    int x1, int y1, int z1, int npointsx, int npointsy) {
  double xd, yd, zd;
  double c00, c10, c01, c11, c0, c1;

  xd = double((x - x0) / (x1 - x0));
  yd = double((y - y0) / (y1 - y0));
  zd = double((z - z0) / (z1 - z0));

  DeviceGridInterpol GI;

  GI.elec = trilinear_interpolation(d_elec_grid, x0, y0, z0, x1, y1, z1, xd, yd,
                                    zd, npointsx, npointsy);
  GI.vdwA = trilinear_interpolation(d_vdwA, x0, y0, z0, x1, y1, z1, xd, yd, zd,
                                    npointsx, npointsy);
  GI.vdwB = trilinear_interpolation(d_vdwB, x0, y0, z0, x1, y1, z1, xd, yd, zd,
                                    npointsx, npointsy);
  GI.hb_donor = trilinear_interpolation(d_hb_donor_grid, x0, y0, z0, x1, y1, z1,
                                        xd, yd, zd, npointsx, npointsy);
  GI.hb_acceptor = trilinear_interpolation(d_acceptor_grid, x0, y0, z0, x1, y1,
                                           z1, xd, yd, zd, npointsx, npointsy);
  GI.rec_solv_gauss = trilinear_interpolation(
      d_rec_solv_gauss, x0, y0, z0, x1, y1, z1, xd, yd, zd, npointsx, npointsy);
  GI.solv_gauss = trilinear_interpolation(d_solv_gauss, x0, y0, z0, x1, y1, z1,
                                          xd, yd, zd, npointsx, npointsy);

  return GI;
}

__global__ void compute_ene_from_grids_softcore_solvation() {

  double elec = 0.0, vdwA = 0.0, vdwB = 0.00, rec_solv = 0.00, lig_solv = 0.00,
         lig_affinity, hb_donor = 0.0, hb_acceptor = 0.0;
  int a1 = 0, b1 = 0, c1 = 0, a2 = 0, b2 = 0, c2 = 0, donor_index;
  double af, bf, cf, angle_term;

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if (i < N) {

    if (a1 > 0 and b1 > 0 and c1 > 0 and a2 < npointsx and b2 < npointsy and
        c2 < npointsz) {
      GridInterpol *GI = new GridInterpol;
      trilinear_interpolation(Grids, af, bf, cf, a1, b1, c1, a2, b2, c2, GI);

      if (use_pbsa and pbsa_loaded) {
        elec += charges[i] * GI->pbsa;
      } else if (use_delphi and delphi_loaded) {
        elec += charges[i] * GI->delphi;
      } else {
        elec += charges[i] * GI->elec;
      }
      vdwA += epsilons_sqrt[i] * pow(radii[i], 6) * GI->vdwA;
      vdwB += epsilons_sqrt[i] * pow(radii[i], 3) * GI->vdwB;

      lig_affinity =
          (solvation_alpha * charges[i] * charges[i]) + solvation_beta;
      rec_solv += GI->solv_gauss * (4.0 / 3.0) * PI * pow(radii[i], 3);
      lig_solv += lig_affinity * GI->rec_solv_gauss;

      if (atom_is_acceptor(i, Lig)) {
        hb_donor += GI->hb_donor;
      } else {
        if (HBdonors.size() > 0) {
          donor_index = H_is_donor(i, Lig);
          if (donor_index >= 0) { // atom is a H donor
            double x = (af * grid_spacing) + xbegin;
            double y = (bf * grid_spacing) + ybegin;
            double z = (cf * grid_spacing) + zbegin;
            angle_term = angle(xyz[donor_index][0], xyz[donor_index][1],
                               xyz[donor_index][2], xyz[i][0], xyz[i][1],
                               xyz[i][2], x, y, z);
            angle_term =
                cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) *
                cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0);
            hb_acceptor += GI->hb_acceptor * angle_term;
          }
        }
      }
      delete GI;
    } else {
      elec += 999999.9;
      vdwA += 999999.9;
      vdwB += 999999.9;
    }
  }
  energy_result->elec = scale_elec_energy * elec;
  energy_result->vdw = scale_vdw_energy * (vdwA - vdwB);
  energy_result->rec_solv = rec_solv;
  energy_result->lig_solv = lig_solv;
  energy_result->hb_donor = hb_donor;
  energy_result->hb_acceptor = hb_acceptor;
  energy_result->total = (scale_elec_energy * elec) +
                         (scale_vdw_energy * (vdwA - vdwB)) + rec_solv +
                         lig_solv + hb_donor + hb_acceptor;

  out = ((scale_elec_energy * elec) + (scale_vdw_energy * (vdwA - vdwB)) +
         rec_solv + lig_solv + hb_donor + hb_acceptor);
}
