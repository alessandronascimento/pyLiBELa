#include "../../src/pyEnergy2.h"
#include "../../src/pyGrid.h"
#include "../../src/pyMol2.h"
#include <iostream> 
#include <vector>
#include "cudaEnergy2.cuh"

__device__ double distance_squared(double x1, double x2, double y1, double y2,
                                   double z1, double z2) {

  return pow(x2 - x1, 2.0) + pow(y2 - y1, 2.0) + pow(z2 - z1, 2.0);
}

__device__ double distance(double x1, double x2, double y1, double y2,
                           double z1, double z2) {

  return sqrt(pow(x2 - x1, 2.0) + pow(y2 - y1, 2.0) + pow(z2 - z1, 2.0));
}

__device__ double angle(double x1, double y1, double z1, double x2, double y2,
                        double z2, double x3, double y3, double z3) {

  double ab = sqrt(distance_squared(x1, x2, y1, y2, z1, z2));
  if (ab == 0.0)
    ab = 0.0001;
  double ac = sqrt(distance_squared(x1, x3, y1, y3, z1, z3));
  if (ac == 0.0)
    ac = 0.0001;
  double bc = sqrt(distance_squared(x2, x3, y2, y3, z2, z3));
  if (bc == 0.0)
    bc = 0.0001;

  double angle = acos(((ab * ab) + (bc * bc) - (ac * ac)) / (2 * ab * bc));
  angle = angle * 180.0 / PI;
  return (angle);
}

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

__device__ void grids_trilinear_interpolation(
    double *d_elec_grid, double *d_vdwA, double *d_vdwB,
    double *d_hb_donor_grid, double *d_acceptor_grid, double *d_rec_solv_gauss,
    double *d_solv_gauss, double x, double y, double z, int x0, int y0, int z0,
    int x1, int y1, int z1, int npointsx, int npointsy,
    DeviceGridInterpol *GI) {

  double xd, yd, zd;
  double c00, c10, c01, c11, c0, c1;

  xd = double((x - x0) / (x1 - x0));
  yd = double((y - y0) / (y1 - y0));
  zd = double((z - z0) / (z1 - z0));

  GI->elec = trilinear_interpolation(d_elec_grid, x0, y0, z0, x1, y1, z1, xd,
                                     yd, zd, npointsx, npointsy);
  GI->vdwA = trilinear_interpolation(d_vdwA, x0, y0, z0, x1, y1, z1, xd, yd, zd,
                                     npointsx, npointsy);
  GI->vdwB = trilinear_interpolation(d_vdwB, x0, y0, z0, x1, y1, z1, xd, yd, zd,
                                     npointsx, npointsy);
  GI->hb_donor = trilinear_interpolation(d_hb_donor_grid, x0, y0, z0, x1, y1,
                                         z1, xd, yd, zd, npointsx, npointsy);
  GI->hb_acceptor = trilinear_interpolation(d_acceptor_grid, x0, y0, z0, x1, y1,
                                            z1, xd, yd, zd, npointsx, npointsy);
  GI->rec_solv_gauss = trilinear_interpolation(
      d_rec_solv_gauss, x0, y0, z0, x1, y1, z1, xd, yd, zd, npointsx, npointsy);
  GI->solv_gauss = trilinear_interpolation(d_solv_gauss, x0, y0, z0, x1, y1, z1,
                                           xd, yd, zd, npointsx, npointsy);
}

__device__ bool atom_is_acceptor(int a, int *d_HBacceptors, int HBacceptors_size) {

  for (int i = 0; i < HBacceptors_size; i++) {
    if (d_HBacceptors[i] == a)
      return true;
  }
  return false;
}

__device__ bool H_is_donor(int a, int* d_HBdonors, int HBdonors_size) {

  for (int i = 0 ; i < HBdonors_size ; i++) {
    if (d_HBdonors[i] == a)
      return true;
  }
  return false;
}

__global__ void compute_ene_from_grids_softcore_solvation(
    int N, int npointsx, int npointsy, int npointsz, double solvation_alpha, double solvation_beta, 
    int HBacceptors_size, int HBdonors_size, double scale_elec_energy, double scale_vdw_energy, 
    double xbegin, double ybegin, double zbegin, double grid_spacing, 
    double *d_elec_grid, double *d_vdwA, double *d_vdwB,
    double *d_hb_donor_grid, double *d_acceptor_grid, double *d_rec_solv_gauss,
    double *d_solv_gauss, double *d_charges, double *d_epsilons_sqrt, double *d_xyz,
    double *d_radii, int *d_HBacceptors, int *d_HBdonors, double *output) {

  double elec = 0.0, vdwA = 0.0, vdwB = 0.00, rec_solv = 0.00, lig_solv = 0.00,
         lig_affinity, hb_donor = 0.0, hb_acceptor = 0.0;
  int a1 = 0, b1 = 0, c1 = 0, a2 = 0, b2 = 0, c2 = 0, donor_index;
  double af, bf, cf, angle_term;

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if (i < N) {

    if (a1 > 0 and b1 > 0 and c1 > 0 and a2 < npointsx and b2 < npointsy and
        c2 < npointsz) {

      DeviceGridInterpol GI;
      grids_trilinear_interpolation(
          d_elec_grid, d_vdwA, d_vdwB, d_hb_donor_grid, d_acceptor_grid, d_rec_solv_gauss, d_solv_gauss,
          af, bf, cf, a1, b1, c1, a2, b2, c2, npointsx, npointsy, &GI);

      /*
            if (use_pbsa and pbsa_loaded) {
              elec += charges[i] * GI->pbsa;
            } else if (use_delphi and delphi_loaded) {
              elec += charges[i] * GI->delphi;
            } else {
              elec += charges[i] * GI->elec;
            }
      */
      elec += d_charges[i] * GI.elec;

      vdwA += d_epsilons_sqrt[i] * pow(d_radii[i], 6) * GI.vdwA;
      vdwB += d_epsilons_sqrt[i] * pow(d_radii[i], 3) * GI.vdwB;

      lig_affinity =
          (solvation_alpha * d_charges[i] * d_charges[i]) + solvation_beta;
      rec_solv += GI.solv_gauss * (4.0 / 3.0) * PI * pow(d_radii[i], 3);
      lig_solv += lig_affinity * GI.rec_solv_gauss;

      if (atom_is_acceptor(i, d_HBacceptors, HBacceptors_size )) {
        hb_donor += GI.hb_donor;
      } 
      else {
        if (HBdonors_size > 0) {
          donor_index = H_is_donor(i, d_HBdonors, HBdonors_size);
          if (donor_index >= 0) { // atom is a H donor
            double x = (af * grid_spacing) + xbegin;
            double y = (bf * grid_spacing) + ybegin;
            double z = (cf * grid_spacing) + zbegin;
            angle_term = angle(d_xyz[donor_index * 3 + 0], d_xyz[donor_index * 3 + 1],
                               d_xyz[donor_index * 3 + 2], d_xyz[i * 3 + 0], d_xyz[i * 3 + 1],
                               d_xyz[i * 3 + 2], x, y, z);
            angle_term =
                cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) *
                cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0);
            hb_acceptor += GI.hb_acceptor * angle_term;
          }
        }
      }
    } else {
      elec += 999999.9;
      vdwA += 999999.9;
      vdwB += 999999.9;
    }
  }

  *output = ((scale_elec_energy * elec) + (scale_vdw_energy * (vdwA - vdwB)) +
         rec_solv + lig_solv + hb_donor + hb_acceptor);
}

double invoke_compute_ene_from_grids_softcore_solvation(Grid* Grids, Mol2 *Lig, 
                                                 std::vector<std::vector<double>> xyz) {

  double *d_elec_grid, *d_vdwA, *d_vdwB, *d_hb_donor_grid, *d_acceptor_grid, *d_rec_solv_gauss, *d_solv_gauss;
  double *d_charges, *d_epsilons_sqrt, *d_xyz, *d_radii, *output;
  int *d_HBdonors, *d_HBacceptors; 

  
}
