
#include "pyGrid.h"
#include "pyMol2.h"
#include <cassert>
#include <iostream>
#include <vector>

#define HB_C12 55332.873
#define HB_C10 18393.199

typedef enum class e_DielectricModel {
  constant,
  r4,
  none,
} DielectricModel;

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
  double ac = sqrt(distance_squared(x1, x3, y1, y3, z1, z3));
  double bc = sqrt(distance_squared(x2, x3, y2, y3, z2, z3));

  double angle = acos(((ab * ab) + (bc * bc) - (ac * ac)) / (2 * ab * bc));
  angle = angle * 180.0 / PI;
  return (angle);
}

__global__ void compute_grid_softcore_HB_CUDA(
    int npointsx, int npointsy, int npointsz, double grid_spacing,
    double xbegin, double ybegin, double zbegin, double *xyz, int N, //
    double deltaij_es6, double deltaij6, double diel, double sigma,
    double solvation_alpha, double solvation_beta, double *charges,
    double *radii, double *epsilons_sqrt, int NHBdonors, int NHBacceptors,
    int *HBdonor1, int *HBdonor2, int *HBacceptors, double *out_elec_grid,
    double *out_vdwA_grid, double *out_vdwB_grid, double *out_rec_solv_grid,
    double *out_solv_grid, double *out_hbd_grid, double *out_hba_grid) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;
  int k = threadIdx.z + blockIdx.z * blockDim.z;

  if (i < npointsx && j < npointsy && k < npointsz) {
    int index = k + (j * npointsx) + (npointsx * npointsy * i);
    double x = i * grid_spacing + xbegin;
    double y = j * grid_spacing + ybegin;
    double z = k * grid_spacing + zbegin;

    double elec = 0.0;
    double solv = 0.0;
    double rec_solv = 0.0;
    double vdwA = 0.0;
    double vdwB = 0.0;
    double hb_donor = 0.0;
    double hb_acceptor = 0.0;

    for (int a = 0; a < N; a++) {
      double d2 = distance_squared(x, xyz[a * 3 + 0], y, xyz[a * 3 + 1], z,
                                   xyz[a * 3 + 2]);
      double d6 = d2 * d2 * d2;
      double denom_vdw = (d6 + deltaij6);
      double denom = cbrt(d6 + deltaij_es6);

      elec += (332.0 * charges[a]) / (diel * denom);
      solv += ((solvation_alpha * charges[a] * charges[a]) + solvation_beta) *
              exp((-denom / (2 * sigma * sigma))) / (sigma * sigma * sigma);
      rec_solv += (4.0 / 3.0) * PI * pow(radii[a], 3) *
                  exp((-denom / (2 * sigma * sigma))) / (sigma * sigma * sigma);

      vdwA += (4096.0 * epsilons_sqrt[a] * pow(radii[a], 6)) /
              (denom_vdw * denom_vdw);
      vdwB += (128.0 * epsilons_sqrt[a] * pow(radii[a], 3)) / denom_vdw;
    }

    for (int a = 0; a < NHBdonors; a++) {
      double d2 = distance_squared(xyz[HBdonor2[a] * 3 + 0], x,
                                   xyz[HBdonor2[a] * 3 + 1], y,
                                   xyz[HBdonor2[a] * 3 + 2], z);
      double d10 = d2 * d2 * d2 * d2 * d2;
      double ang =
          angle(xyz[HBdonor1[a] * 3 + 0], xyz[HBdonor1[a] * 3 + 1],
                xyz[HBdonor1[a] * 3 + 2], xyz[HBdonor2[a] * 3 + 0],
                xyz[HBdonor2[a] * 3 + 1], xyz[HBdonor2[a] * 3 + 2], x, y, z);
      double angle_term = cos(ang * PI / 180.0) * cos(ang * PI / 180.0) *
                          cos(ang * PI / 180.0) * cos(ang * PI / 180.0);
      hb_donor += ((HB_C12 / (d10 * d2)) - (HB_C10 / d10)) * angle_term;
    }

    for (unsigned a = 0; a < NHBacceptors; a++) {
      double d2 = distance_squared(xyz[HBacceptors[a] * 3 + 0], x,
                                   xyz[HBacceptors[a] * 3 + 1], y,
                                   xyz[HBacceptors[a] * 3 + 2], z);
      double d10 = d2 * d2 * d2 * d2 * d2;
      hb_acceptor += ((HB_C12 / (d10 * d2)) - (HB_C10 / d10));
    }

    out_elec_grid[index] = elec;
    out_vdwA_grid[index] = vdwA;
    out_vdwB_grid[index] = vdwB;
    out_rec_solv_grid[index] = rec_solv;
    out_solv_grid[index] = solv;
    out_hbd_grid[index] = hb_donor;
    out_hba_grid[index] = hb_acceptor;
  }
}

__global__ void compute_grid_hardcore_HB_CUDA(
    int npointsx, int npointsy, int npointsz, double grid_spacing,
    double xbegin, double ybegin, double zbegin,
    DielectricModel dielectric_model, double *xyz, int N, //
    double deltaij_es6, double deltaij6, double diel, double sigma,
    double solvation_alpha, double solvation_beta, double *charges,
    double *radii, double *epsilons_sqrt, int NHBdonors, int NHBacceptors,
    int *HBdonor1, int *HBdonor2, int *HBacceptors, double *out_elec_grid,
    double *out_vdwA_grid, double *out_vdwB_grid, double *out_rec_solv_grid,
    double *out_solv_grid, double *out_hbd_grid, double *out_hba_grid) {
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int j = threadIdx.y + blockIdx.y * blockDim.y;
  int k = threadIdx.z + blockIdx.z * blockDim.z;

  if (i < npointsx && j < npointsy && k < npointsz) {
    int index = k + (j * npointsx) + (npointsx * npointsy * i);
    double x = i * grid_spacing + xbegin;
    double y = j * grid_spacing + ybegin;
    double z = k * grid_spacing + zbegin;

    double elec = 0.0;
    double solv = 0.0;
    double rec_solv = 0.0;
    double vdwA = 0.0;
    double vdwB = 0.0;
    double hb_donor = 0.0;
    double hb_acceptor = 0.0;

    for (int a = 0; a < N; a++) {
      double d2 = distance_squared(x, xyz[a * 3 + 0], y, xyz[a * 3 + 1], z,
                                   xyz[a * 3 + 2]);
      double d6 = d2 * d2 * d2;

      if (dielectric_model == DielectricModel::constant) {
        double d = sqrt(d2);
        elec += 332.0 * ((charges[a]) / (d * diel));
      } else if (dielectric_model == DielectricModel::r4) {
        elec += 332.0 * (charges[a] / (4 * d2));
      } else {
        elec += 332.0 * (charges[a] / d2);
      }

      vdwA += (64.0 * epsilons_sqrt[a] * pow(radii[a], 6)) / (d6 * d6);
      vdwB += (8.0 * M_SQRT2 * epsilons_sqrt[a] * pow(radii[a], 3)) / d6;

      solv += ((solvation_alpha * charges[a] * charges[a]) + solvation_beta) *
              exp((-d2 / (2 * sigma * sigma))) / (sigma * sigma * sigma);
      rec_solv += (4.0 / 3.0) * PI * pow(radii[a], 3) *
                  exp((-d2 / (2 * sigma * sigma))) / (sigma * sigma * sigma);
    }

    for (int a = 0; a < NHBdonors; a++) {
      double d2 = distance_squared(xyz[HBdonor2[a] * 3 + 0], x,
                                   xyz[HBdonor2[a] * 3 + 1], y,
                                   xyz[HBdonor2[a] * 3 + 2], z);
      double d10 = d2 * d2 * d2 * d2 * d2;
      double ang =
          angle(xyz[HBdonor1[a] * 3 + 0], xyz[HBdonor1[a] * 3 + 1],
                xyz[HBdonor1[a] * 3 + 2], xyz[HBdonor2[a] * 3 + 0],
                xyz[HBdonor2[a] * 3 + 1], xyz[HBdonor2[a] * 3 + 2], x, y, z);
      double angle_term = cos(ang * PI / 180.0) * cos(ang * PI / 180.0) *
                          cos(ang * PI / 180.0) * cos(ang * PI / 180.0);
      hb_donor += ((HB_C12 / (d10 * d2)) - (HB_C10 / d10)) * angle_term;
    }

    for (unsigned a = 0; a < NHBacceptors; a++) {
      double d2 = distance_squared(xyz[HBacceptors[a] * 3 + 0], x,
                                   xyz[HBacceptors[a] * 3 + 1], y,
                                   xyz[HBacceptors[a] * 3 + 2], z);
      double d10 = d2 * d2 * d2 * d2 * d2;
      hb_acceptor += ((HB_C12 / (d10 * d2)) - (HB_C10 / d10));
    }

    out_elec_grid[index] = elec;
    out_vdwA_grid[index] = vdwA;
    out_vdwB_grid[index] = vdwB;
    out_rec_solv_grid[index] = rec_solv;
    out_solv_grid[index] = solv;
    out_hbd_grid[index] = hb_donor;
    out_hba_grid[index] = hb_acceptor;
  }
}

void invoke_compute_grid_softcore_HB_CUDA(Grid *grid, Mol2 *rec) {

  double *d_xyz, *d_charges, *d_radii, *d_epsilon_sqrt;
  double *d_out_elec_grid, *d_out_vdwA_grid, *d_out_vdwB_grid,
      *d_out_rec_solv_grid, *d_out_solv_grid, *d_out_hbd_grid, *d_out_hba_grid;

  cudaMalloc(&d_xyz, rec->xyz.size() * rec->xyz[0].size() * sizeof(double));
  cudaMalloc(&d_charges, rec->charges.size() * sizeof(double));
  cudaMalloc(&d_radii, rec->radii.size() * sizeof(double));
  cudaMalloc(&d_epsilon_sqrt, rec->epsilons_sqrt.size() * sizeof(double));

  size_t size_bytes{grid->npointsx * grid->npointsy * grid->npointsz *
                    sizeof(double)};

  cudaMalloc(&d_out_elec_grid, size_bytes);
  cudaMalloc(&d_out_vdwA_grid, size_bytes);
  cudaMalloc(&d_out_vdwB_grid, size_bytes);
  cudaMalloc(&d_out_solv_grid, size_bytes);
  cudaMalloc(&d_out_rec_solv_grid, size_bytes);
  cudaMalloc(&d_out_hbd_grid, size_bytes);
  cudaMalloc(&d_out_hba_grid, size_bytes);

  double *xyz_arr = new double[rec->xyz.size() * 3];

  for (int i = 0; i < rec->xyz.size(); i++) {
    xyz_arr[3 * i] = rec->xyz[i][0];
    xyz_arr[(3 * i) + 1] = rec->xyz[i][1];
    xyz_arr[(3 * i) + 2] = rec->xyz[i][2];
  }

  int *d_HBdonnor1;
  int *d_HBdonnor2;
  int *HBdonnor1 = new int[rec->HBdonors.size()];
  int *HBdonnor2 = new int[rec->HBdonors.size()];
  int *d_HBacceptors = new int[rec->HBacceptors.size()];

  for (int i = 0; i < rec->HBdonors.size(); i++) {
    HBdonnor1[i] = rec->HBdonors[i][0];
    HBdonnor2[i] = rec->HBdonors[i][1];
  }

  cudaMalloc(&d_HBdonnor1, rec->HBdonors.size() * sizeof(int));
  cudaMalloc(&d_HBdonnor2, rec->HBdonors.size() * sizeof(int));
  cudaMalloc(&d_HBacceptors, rec->HBacceptors.size() * sizeof(int));

  cudaMemcpy(d_xyz, xyz_arr,
             rec->xyz.size() * rec->xyz[0].size() * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_charges, rec->charges.data(),
             rec->charges.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_radii, rec->radii.data(), rec->radii.size() * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_epsilon_sqrt, rec->epsilons_sqrt.data(),
             rec->epsilons_sqrt.size() * sizeof(double),
             cudaMemcpyHostToDevice);

  cudaMemcpy(d_HBdonnor1, HBdonnor1, rec->HBdonors.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_HBdonnor2, HBdonnor2, rec->HBdonors.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_HBacceptors, rec->HBacceptors.data(),
             rec->HBacceptors.size() * sizeof(int), cudaMemcpyHostToDevice);

  grid->r_elec_grid = (double *)malloc(size_bytes);
  grid->r_vdwA_grid = (double *)malloc(size_bytes);
  grid->r_vdwB_grid = (double *)malloc(size_bytes);
  grid->r_solv_gauss = (double *)malloc(size_bytes);
  grid->r_rec_solv_gauss = (double *)malloc(size_bytes);
  grid->r_hb_donor_grid = (double *)malloc(size_bytes);
  grid->r_hb_acceptor_grid = (double *)malloc(size_bytes);

  dim3 blockdims(8, 8, 8);
  dim3 griddims(8, 8, 8);

  printf("Entering kernel\n");
  compute_grid_softcore_HB_CUDA<<<griddims, blockdims>>>(
      grid->npointsx, grid->npointsy, grid->npointsz, grid->grid_spacing,
      grid->xbegin, grid->ybegin, grid->zbegin, d_xyz, rec->N,
      grid->Input->deltaij_es6, grid->Input->deltaij6, grid->Input->diel,
      grid->Input->sigma, grid->Input->solvation_alpha,
      grid->Input->solvation_beta, d_charges, d_radii, d_epsilon_sqrt,
      rec->HBdonors.size(), rec->HBacceptors.size(), d_HBdonnor1, d_HBdonnor2,
      d_HBacceptors, d_out_elec_grid, d_out_vdwA_grid, d_out_vdwB_grid,
      d_out_rec_solv_grid, d_out_solv_grid, d_out_hbd_grid, d_out_hba_grid);

  printf("Kernel has ended\n");

  cudaMemcpy(grid->r_elec_grid, d_out_elec_grid, size_bytes,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(grid->r_vdwA_grid, d_out_vdwA_grid, size_bytes,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(grid->r_vdwB_grid, d_out_vdwB_grid, size_bytes,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(grid->r_solv_gauss, d_out_solv_grid, size_bytes,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(grid->r_rec_solv_gauss, d_out_rec_solv_grid, size_bytes,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(grid->r_hb_donor_grid, d_out_hbd_grid, size_bytes,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(grid->r_hb_acceptor_grid, d_out_hba_grid, size_bytes,
             cudaMemcpyDeviceToHost);

  int nx{grid->npointsx};
  int ny{grid->npointsy};
  int nz{grid->npointsz};

  grid->elec_grid = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));
  grid->vdwA_grid = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));
  grid->vdwB_grid = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));
  grid->rec_solv_gauss = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));
  grid->solv_gauss = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));
  grid->hb_acceptor_grid = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));
  grid->hb_donor_grid = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));

  for (int i = 0; i < nx; i++) {

    for (int j = 0; j < ny; j++) {
      int idx = (nx * ny * i + nx * j);
      grid->elec_grid[i][j].insert(grid->elec_grid[i][j].end(),
                                   grid->r_elec_grid + idx,
                                   grid->r_elec_grid + idx + nz);
      grid->vdwA_grid[i][j].insert(grid->vdwA_grid[i][j].end(),
                                   grid->r_vdwA_grid + idx,
                                   grid->r_vdwA_grid + idx + nz);
      grid->vdwB_grid[i][j].insert(grid->vdwB_grid[i][j].end(),
                                   grid->r_vdwB_grid + idx,
                                   grid->r_vdwB_grid + idx + nz);
      grid->rec_solv_gauss[i][j].insert(grid->rec_solv_gauss[i][j].end(),
                                        grid->r_rec_solv_gauss + idx,
                                        grid->r_rec_solv_gauss + idx + nz);
      grid->solv_gauss[i][j].insert(grid->solv_gauss[i][j].end(),
                                    grid->r_solv_gauss + idx,
                                    grid->r_solv_gauss + idx + nz);
      grid->hb_acceptor_grid[i][j].insert(grid->hb_acceptor_grid[i][j].end(),
                                          grid->r_hb_acceptor_grid + idx,
                                          grid->r_hb_acceptor_grid + idx + nz);
      grid->hb_donor_grid[i][j].insert(grid->hb_donor_grid[i][j].end(),
                                       grid->r_hb_donor_grid + idx,
                                       grid->r_hb_donor_grid + idx + nz);
    }
  }

  delete[] (xyz_arr);
  delete (HBdonnor1);
  delete (HBdonnor2);
  cudaFree(d_xyz);

  cudaFree(d_charges);
  cudaFree(d_radii);
  cudaFree(d_epsilon_sqrt);
  cudaFree(d_HBacceptors);
  cudaFree(d_HBdonnor1);
  cudaFree(d_HBdonnor2);
  cudaFree(d_out_elec_grid);
  cudaFree(d_out_vdwA_grid);
  cudaFree(d_out_vdwB_grid);
  cudaFree(d_out_solv_grid);
  cudaFree(d_out_rec_solv_grid);
  cudaFree(d_out_solv_grid);
  cudaFree(d_out_hba_grid);
  cudaFree(d_out_hbd_grid);

  printf("Invoking finished\n");
}

void invoke_compute_grid_hardcore_HB_CUDA(Grid *grid, Mol2 *rec) {

  double *d_xyz, *d_charges, *d_radii, *d_epsilon_sqrt;
  double *d_out_elec_grid, *d_out_vdwA_grid, *d_out_vdwB_grid,
      *d_out_rec_solv_grid, *d_out_solv_grid, *d_out_hbd_grid, *d_out_hba_grid;

  cudaMalloc(&d_xyz, rec->xyz.size() * rec->xyz[0].size() * sizeof(double));
  cudaMalloc(&d_charges, rec->charges.size() * sizeof(double));
  cudaMalloc(&d_radii, rec->radii.size() * sizeof(double));
  cudaMalloc(&d_epsilon_sqrt, rec->epsilons_sqrt.size() * sizeof(double));

  size_t size_bytes{grid->npointsx * grid->npointsy * grid->npointsz *
                    sizeof(double)};

  cudaMalloc(&d_out_elec_grid, size_bytes);
  cudaMalloc(&d_out_vdwA_grid, size_bytes);
  cudaMalloc(&d_out_vdwB_grid, size_bytes);
  cudaMalloc(&d_out_solv_grid, size_bytes);
  cudaMalloc(&d_out_rec_solv_grid, size_bytes);
  cudaMalloc(&d_out_hbd_grid, size_bytes);
  cudaMalloc(&d_out_hba_grid, size_bytes);

  double *xyz_arr = new double[rec->xyz.size() * 3];

  for (int i = 0; i < rec->xyz.size(); i++) {
    xyz_arr[3 * i] = rec->xyz[i][0];
    xyz_arr[(3 * i) + 1] = rec->xyz[i][1];
    xyz_arr[(3 * i) + 2] = rec->xyz[i][2];
  }

  int *d_HBdonnor1;
  int *d_HBdonnor2;
  int *HBdonnor1 = new int[rec->HBdonors.size()];
  int *HBdonnor2 = new int[rec->HBdonors.size()];
  int *d_HBacceptors = new int[rec->HBacceptors.size()];

  for (int i = 0; i < rec->HBdonors.size(); i++) {
    HBdonnor1[i] = rec->HBdonors[i][0];
    HBdonnor2[i] = rec->HBdonors[i][1];
  }

  cudaMalloc(&d_HBdonnor1, rec->HBdonors.size() * sizeof(int));
  cudaMalloc(&d_HBdonnor2, rec->HBdonors.size() * sizeof(int));
  cudaMalloc(&d_HBacceptors, rec->HBacceptors.size() * sizeof(int));

  cudaMemcpy(d_xyz, xyz_arr,
             rec->xyz.size() * rec->xyz[0].size() * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_charges, rec->charges.data(),
             rec->charges.size() * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_radii, rec->radii.data(), rec->radii.size() * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_epsilon_sqrt, rec->epsilons_sqrt.data(),
             rec->epsilons_sqrt.size() * sizeof(double),
             cudaMemcpyHostToDevice);

  cudaMemcpy(d_HBdonnor1, HBdonnor1, rec->HBdonors.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_HBdonnor2, HBdonnor2, rec->HBdonors.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(d_HBacceptors, rec->HBacceptors.data(),
             rec->HBacceptors.size() * sizeof(int), cudaMemcpyHostToDevice);

  grid->r_elec_grid = (double *)malloc(size_bytes);
  grid->r_vdwA_grid = (double *)malloc(size_bytes);
  grid->r_vdwB_grid = (double *)malloc(size_bytes);
  grid->r_solv_gauss = (double *)malloc(size_bytes);
  grid->r_rec_solv_gauss = (double *)malloc(size_bytes);
  grid->r_hb_donor_grid = (double *)malloc(size_bytes);
  grid->r_hb_acceptor_grid = (double *)malloc(size_bytes);

  dim3 blockdims(8, 8, 8);
  dim3 griddims(8, 8, 8);

  DielectricModel dielectric_model;

  if (grid->Input->dielectric_model == "constant")
    dielectric_model = DielectricModel::constant;
  else if (grid->Input->dielectric_model == "4r")
    dielectric_model = DielectricModel::r4;
  else
    dielectric_model = DielectricModel::none;

  printf("Entering kernel\n");
  compute_grid_hardcore_HB_CUDA<<<griddims, blockdims>>>(
      grid->npointsx, grid->npointsy, grid->npointsz, grid->grid_spacing,
      grid->xbegin, grid->ybegin, grid->zbegin, dielectric_model, d_xyz, rec->N,
      grid->Input->deltaij_es6, grid->Input->deltaij6, grid->Input->diel,
      grid->Input->sigma, grid->Input->solvation_alpha,
      grid->Input->solvation_beta, d_charges, d_radii, d_epsilon_sqrt,
      rec->HBdonors.size(), rec->HBacceptors.size(), d_HBdonnor1, d_HBdonnor2,
      d_HBacceptors, d_out_elec_grid, d_out_vdwA_grid, d_out_vdwB_grid,
      d_out_rec_solv_grid, d_out_solv_grid, d_out_hbd_grid, d_out_hba_grid);

  printf("Kernel has ended\n");

  cudaMemcpy(grid->r_elec_grid, d_out_elec_grid, size_bytes,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(grid->r_vdwA_grid, d_out_vdwA_grid, size_bytes,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(grid->r_vdwB_grid, d_out_vdwB_grid, size_bytes,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(grid->r_solv_gauss, d_out_solv_grid, size_bytes,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(grid->r_rec_solv_gauss, d_out_rec_solv_grid, size_bytes,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(grid->r_hb_donor_grid, d_out_hbd_grid, size_bytes,
             cudaMemcpyDeviceToHost);
  cudaMemcpy(grid->r_hb_acceptor_grid, d_out_hba_grid, size_bytes,
             cudaMemcpyDeviceToHost);

  int nx{grid->npointsx};
  int ny{grid->npointsy};
  int nz{grid->npointsz};

  grid->elec_grid = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));
  grid->vdwA_grid = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));
  grid->vdwB_grid = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));
  grid->rec_solv_gauss = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));
  grid->solv_gauss = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));
  grid->hb_acceptor_grid = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));
  grid->hb_donor_grid = std::vector<std::vector<std::vector<double>>>(
      nx, std::vector<std::vector<double>>(ny));

  for (int i = 0; i < nx; i++) {

    for (int j = 0; j < ny; j++) {
      int idx = (nx * ny * i + nx * j);
      grid->elec_grid[i][j].insert(grid->elec_grid[i][j].end(),
                                   grid->r_elec_grid + idx,
                                   grid->r_elec_grid + idx + nz);
      grid->vdwA_grid[i][j].insert(grid->vdwA_grid[i][j].end(),
                                   grid->r_vdwA_grid + idx,
                                   grid->r_vdwA_grid + idx + nz);
      grid->vdwB_grid[i][j].insert(grid->vdwB_grid[i][j].end(),
                                   grid->r_vdwB_grid + idx,
                                   grid->r_vdwB_grid + idx + nz);
      grid->rec_solv_gauss[i][j].insert(grid->rec_solv_gauss[i][j].end(),
                                        grid->r_rec_solv_gauss + idx,
                                        grid->r_rec_solv_gauss + idx + nz);
      grid->solv_gauss[i][j].insert(grid->solv_gauss[i][j].end(),
                                    grid->r_solv_gauss + idx,
                                    grid->r_solv_gauss + idx + nz);
      grid->hb_acceptor_grid[i][j].insert(grid->hb_acceptor_grid[i][j].end(),
                                          grid->r_hb_acceptor_grid + idx,
                                          grid->r_hb_acceptor_grid + idx + nz);
      grid->hb_donor_grid[i][j].insert(grid->hb_donor_grid[i][j].end(),
                                       grid->r_hb_donor_grid + idx,
                                       grid->r_hb_donor_grid + idx + nz);
    }
  }

  delete[] (xyz_arr);
  delete[] (HBdonnor1);
  delete[] (HBdonnor2);

  cudaFree(d_xyz);
  cudaFree(d_charges);
  cudaFree(d_radii);
  cudaFree(d_epsilon_sqrt);
  cudaFree(d_HBacceptors);
  cudaFree(d_HBdonnor1);
  cudaFree(d_HBdonnor2);
  cudaFree(d_out_elec_grid);
  cudaFree(d_out_vdwA_grid);
  cudaFree(d_out_vdwB_grid);
  cudaFree(d_out_solv_grid);
  cudaFree(d_out_rec_solv_grid);
  cudaFree(d_out_solv_grid);
  cudaFree(d_out_hba_grid);
  cudaFree(d_out_hbd_grid);

  printf("Invoking finished\n");
}
