#include "pyEnergy2.h"
#include "pyGrid.h"
#include "pyMol2.h"
#include <cuda_runtime_api.h>
#include <device_atomic_functions.h>
#include <driver_types.h>
#include <iostream> 
#include <vector>
#include "cudaEnergy2.cuh"

__device__ double ene_distance_squared(double x1, double x2, double y1, double y2,
                                   double z1, double z2) {

  return pow(x2 - x1, 2.0) + pow(y2 - y1, 2.0) + pow(z2 - z1, 2.0);
}

__device__ double ene_angle(double x1, double y1, double z1, double x2, double y2,
                        double z2, double x3, double y3, double z3) {

  double ab = sqrt(ene_distance_squared(x1, x2, y1, y2, z1, z2));
  if (ab == 0.0)
    ab = 0.0001;
  double ac = sqrt(ene_distance_squared(x1, x3, y1, y3, z1, z3));
  if (ac == 0.0)
    ac = 0.0001;
  double bc = sqrt(ene_distance_squared(x2, x3, y2, y3, z2, z3));
  if (bc == 0.0)
    bc = 0.0001;

  double angle = acos(((ab * ab) + (bc * bc) - (ac * ac)) / (2 * ab * bc));
  angle = angle * 180.0 / PI;
  return (angle);
}

__device__ bool atom_is_acceptor(int a, int *d_HBacceptors, int HBacceptors_size) {

  for (int i = 0; i < HBacceptors_size; i++) {
    if (d_HBacceptors[i] == a)
      return true;
  }
  return false;
}

__device__ int H_is_donor(int a, int* d_HBdonors0, int* d_HBdonors1, int HBdonors_size) {

    for (int i = 0; i < HBdonors_size; i++) {
      if (d_HBdonors1[i] == a) return d_HBdonors0[i];
    }

    return -1;
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
    int x1, int y1, int z1, int npointsx, int npointsy, int scoring_function,
    DeviceGridInterpol *GI) {

  double xd, yd, zd;

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

  if ((scoring_function == 0) or (scoring_function == 2) or (scoring_function == 4)){

    GI->rec_solv_gauss = trilinear_interpolation(
        d_rec_solv_gauss, x0, y0, z0, x1, y1, z1, xd, yd, zd, npointsx, npointsy);
    GI->solv_gauss = trilinear_interpolation(d_solv_gauss, x0, y0, z0, x1, y1, z1,
                                            xd, yd, zd, npointsx, npointsy);
  }

  else {
    GI->rec_solv_gauss = 0.0;
    GI->solv_gauss = 0.0;
  }
}


__global__ void compute_ene_from_grids_softcore_solvation(
    int N, int npointsx, int npointsy, int npointsz, int scoring_function, 
    double solvation_alpha, double solvation_beta, int HBacceptors_size, int HBdonors_size,  
    double xbegin, double ybegin, double zbegin, double grid_spacing, 
    double *d_elec_grid, double *d_vdwA, double *d_vdwB,
    double *d_hb_donor_grid, double *d_acceptor_grid, double *d_rec_solv_gauss,
    double *d_solv_gauss, double *d_charges, double *d_epsilons_sqrt, double *d_xyz,
    double *d_radii, int *d_HBacceptors, int *d_HBdonors0, int *d_HBdonors1,  
    double *elec, double *vdwA, double *vdwB, double *rec_solv, 
    double *lig_solv, double *hb_donor, double *hb_acceptor) {

  int a1 = 0, b1 = 0, c1 = 0, a2 = 0, b2 = 0, c2 = 0, donor_index;
  double af, bf, cf, lig_affinity, angle_term;

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if (i >= N) return;

  af = (d_xyz[i * 3 + 0] - xbegin)/grid_spacing;
  bf = (d_xyz[i * 3 + 1] - ybegin)/grid_spacing;
  cf = (d_xyz[i * 3 + 2] - zbegin)/grid_spacing;
  a1 = floor(af);
  b1 = floor(bf);
  c1 = floor(cf);
  a2 = ceil(af);
  if (a2 <= a1){
      a2++;
  }
  b2 = ceil(bf);
  if (b2 <= b1){
      b2++;
  }
  c2 = ceil(cf);
  if (c2 <= c1){
      c2++;
  }

  // If it is outside the Grid boundaries, penalize with a high potential 
  if (a1 <= 0 or b1 <= 0 or c1 <= 0 or a2 >= npointsx or b2 >= npointsy or c2 >= npointsz) {
    atomicAdd(elec, 999999.9); 
    atomicAdd(vdwA, 999999.9); 
    atomicAdd(vdwB, 999999.9); 
    return;
  }

  DeviceGridInterpol GI;
  grids_trilinear_interpolation(
      d_elec_grid, d_vdwA, d_vdwB, d_hb_donor_grid, d_acceptor_grid, d_rec_solv_gauss, d_solv_gauss,
      af, bf, cf, a1, b1, c1, a2, b2, c2, npointsx, npointsy, scoring_function, &GI);

  /*
        if (use_pbsa and pbsa_loaded) {
          elec += charges[i] * GI->pbsa;
        } else if (use_delphi and delphi_loaded) {
          elec += charges[i] * GI->delphi;
        } else {
          elec += charges[i] * GI->elec;
        }
  */
  atomicAdd(elec, d_charges[i] * GI.elec); 

  atomicAdd(vdwA, d_epsilons_sqrt[i] * pow(d_radii[i], 6) * GI.vdwA); 
  atomicAdd(vdwB, d_epsilons_sqrt[i] * pow(d_radii[i], 3) * GI.vdwB); 

  lig_affinity = (solvation_alpha * d_charges[i] * d_charges[i]) + solvation_beta;
  atomicAdd(rec_solv, GI.solv_gauss * (4.0 / 3.0) * PI * pow(d_radii[i], 3));
  atomicAdd(lig_solv, lig_affinity * GI.rec_solv_gauss); 

  if (atom_is_acceptor(i, d_HBacceptors, HBacceptors_size )) {
    atomicAdd(hb_donor, GI.hb_donor); 
  } 

  else if (HBdonors_size > 0) {

    donor_index = H_is_donor(i, d_HBdonors0, d_HBdonors1, HBdonors_size);

    if (donor_index >= 0) { // atom is a H donor
      double x = (af * grid_spacing) + xbegin;
      double y = (bf * grid_spacing) + ybegin;
      double z = (cf * grid_spacing) + zbegin;
      angle_term = ene_angle(d_xyz[donor_index * 3 + 0], d_xyz[donor_index * 3 + 1],
                          d_xyz[donor_index * 3 + 2], d_xyz[i * 3 + 0], d_xyz[i * 3 + 1],
                          d_xyz[i * 3 + 2], x, y, z);
      angle_term =
          cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) *
          cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0);
      atomicAdd(hb_acceptor, GI.hb_acceptor * angle_term);
    }
  }
}


double invoke_compute_ene_from_grids_softcore_solvation(Grid* grid, Mol2 *lig, 
                                                 const std::vector<std::vector<double>>& xyz) {

  double *d_elec_grid, *d_vdwA, *d_vdwB, *d_hb_donor_grid, *d_acceptor_grid, *d_rec_solv_gauss, *d_solv_gauss;
  double *d_charges, *d_epsilons_sqrt, *d_xyz, *d_radii;
  double *elec, *vdwA, *vdwB, *rec_solv, *lig_solv, *hb_donor, *hb_acceptor; 
  int *d_HBdonors0, *d_HBdonors1, *d_HBacceptors; 

  size_t size_bytes{grid->npointsx * grid->npointsy * grid->npointsz * sizeof(double)};

  cudaMalloc(&d_elec_grid, size_bytes);
  cudaMalloc(&d_vdwA, size_bytes);
  cudaMalloc(&d_vdwB, size_bytes);
  cudaMalloc(&d_solv_gauss, size_bytes);
  cudaMalloc(&d_rec_solv_gauss, size_bytes);
  cudaMalloc(&d_hb_donor_grid, size_bytes);
  cudaMalloc(&d_acceptor_grid, size_bytes);

  cudaMemcpy(d_elec_grid, grid->r_elec_grid, size_bytes, cudaMemcpyHostToDevice); 
  cudaMemcpy(d_vdwA, grid->r_vdwA_grid, size_bytes, cudaMemcpyHostToDevice); 
  cudaMemcpy(d_vdwB, grid->r_vdwB_grid, size_bytes, cudaMemcpyHostToDevice); 
  cudaMemcpy(d_solv_gauss, grid->r_solv_gauss, size_bytes, cudaMemcpyHostToDevice); 
  cudaMemcpy(d_rec_solv_gauss, grid->r_rec_solv_gauss, size_bytes, cudaMemcpyHostToDevice); 
  cudaMemcpy(d_hb_donor_grid, grid->r_hb_donor_grid, size_bytes, cudaMemcpyHostToDevice); 
  cudaMemcpy(d_acceptor_grid, grid->r_hb_acceptor_grid, size_bytes, cudaMemcpyHostToDevice); 

  cudaMalloc(&elec, sizeof(double));
  cudaMalloc(&vdwA, sizeof(double));
  cudaMalloc(&vdwB, sizeof(double));
  cudaMalloc(&rec_solv, sizeof(double));
  cudaMalloc(&lig_solv, sizeof(double));
  cudaMalloc(&hb_donor, sizeof(double));
  cudaMalloc(&hb_acceptor, sizeof(double));

  double *out_elec = (double*) malloc(sizeof(double)); 
  double *out_vdwA = (double*) malloc(sizeof(double)); 
  double *out_vdwB = (double*) malloc(sizeof(double)); 
  double *out_rec_solv = (double*) malloc(sizeof(double)); 
  double *out_lig_solv = (double*) malloc(sizeof(double)); 
  double *out_hb_donor = (double*) malloc(sizeof(double)); 
  double *out_hb_acceptor = (double*) malloc(sizeof(double)); 

  cudaMemset(elec, 0, sizeof(double));
  cudaMemset(vdwA, 0, sizeof(double));
  cudaMemset(vdwB, 0, sizeof(double));
  cudaMemset(rec_solv, 0, sizeof(double));
  cudaMemset(lig_solv, 0, sizeof(double));
  cudaMemset(hb_donor, 0, sizeof(double));
  cudaMemset(hb_acceptor, 0, sizeof(double));

  cudaMalloc(&d_xyz, xyz.size() * xyz[0].size() * sizeof(double));
  cudaMalloc(&d_charges, lig->charges.size() * sizeof(double));
  cudaMalloc(&d_radii, lig->radii.size() * sizeof(double));
  cudaMalloc(&d_epsilons_sqrt, lig->epsilons_sqrt.size() * sizeof(double));
  cudaMalloc(&d_HBdonors0, lig->HBdonors.size() * sizeof(int));
  cudaMalloc(&d_HBdonors1, lig->HBdonors.size() * sizeof(int));
  cudaMalloc(&d_HBacceptors, lig->HBacceptors.size() * sizeof(int));

  double *xyz_temp = new double[xyz.size() * 3];

  for (int i = 0; i < xyz.size(); i++) {
    xyz_temp[3 * i] = xyz[i][0];
    xyz_temp[(3 * i) + 2] = xyz[i][2];
    xyz_temp[(3 * i) + 1] = xyz[i][1];
  }

  int *HBdonor0_temp = new int[lig->HBdonors.size()];
  int *HBdonor1_temp = new int[lig->HBdonors.size()];

  for (int i = 0; i < lig->HBdonors.size(); i++) {
    HBdonor0_temp[i] = lig->HBdonors[i][0];
    HBdonor1_temp[i] = lig->HBdonors[i][1];
  }
  
  cudaMemcpy(d_xyz, xyz_temp, xyz.size() * 3 * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_charges, lig->charges.data(), lig->charges.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_radii, lig->radii.data(), lig->radii.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_epsilons_sqrt, lig->epsilons_sqrt.data(), lig->epsilons_sqrt.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_HBdonors0, HBdonor0_temp, lig->HBdonors.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_HBdonors1, HBdonor1_temp, lig->HBdonors.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_HBacceptors, lig->HBacceptors.data(), lig->HBacceptors.size() * sizeof(int), cudaMemcpyHostToDevice);

  // ------------ KERNEL START -----------------------------------------

  compute_ene_from_grids_softcore_solvation<<<(lig->N + 255) /256 , 256>>>(
      lig->N, grid->npointsx, grid->npointsy, grid->npointsz, grid->Input->scoring_function,
      grid->Input->solvation_alpha, grid->Input->solvation_beta, lig->HBacceptors.size(), 
      lig->HBdonors.size(), grid->xbegin, grid->ybegin, grid->zbegin, grid->grid_spacing,
      d_elec_grid, d_vdwA, d_vdwB, d_hb_donor_grid, d_acceptor_grid, d_rec_solv_gauss,
      d_solv_gauss, d_charges, d_epsilons_sqrt, d_xyz, d_radii, d_HBacceptors,
      d_HBdonors0, d_HBdonors1, elec, vdwA, vdwB, rec_solv, lig_solv, hb_donor, hb_acceptor
      );

  cudaMemcpy(out_elec, elec, sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(out_vdwA, vdwA, sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(out_vdwB, vdwB, sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(out_rec_solv, rec_solv, sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(out_lig_solv, lig_solv, sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(out_hb_donor, hb_donor, sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(out_hb_acceptor, hb_acceptor, sizeof(double), cudaMemcpyDeviceToHost);

  double energy = ((grid->Input->scale_elec_energy * *out_elec) + 
                  (grid->Input->scale_vdw_energy*(*out_vdwA - *out_vdwB)) + 
                  *out_rec_solv + *out_lig_solv + *out_hb_donor + *out_hb_acceptor);


  cudaFree(d_elec_grid);
  cudaFree(d_vdwA);
  cudaFree(d_vdwB);
  cudaFree(d_solv_gauss);
  cudaFree(d_rec_solv_gauss);
  cudaFree(d_hb_donor_grid);
  cudaFree(d_acceptor_grid);

  cudaFree(elec);
  cudaFree(vdwA);
  cudaFree(vdwB);
  cudaFree(rec_solv);
  cudaFree(lig_solv);
  cudaFree(hb_donor);
  cudaFree(hb_acceptor);

  free(out_elec);
  free(out_vdwA);
  free(out_vdwB);
  free(out_rec_solv);
  free(out_lig_solv);
  free(out_hb_donor);
  free(out_hb_acceptor);

  cudaFree(d_xyz);
  cudaFree(d_charges);
  cudaFree(d_radii);
  cudaFree(d_epsilons_sqrt);
  cudaFree(d_HBdonors0);
  cudaFree(d_HBdonors1);
  cudaFree(d_HBacceptors);

  delete[] (xyz_temp);
  delete[] (HBdonor0_temp);
  delete[] (HBdonor1_temp);

  return energy;
}

// TODO Error handling
// TODO Extract distance, distance_squared and angle functions
