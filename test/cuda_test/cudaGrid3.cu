#include <iostream>
#include <vector>
#include <cassert>
#include "../../src/pyGrid.h"
#include "../../src/pyMol2.h"

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

#define HB_C12 55332.873
#define HB_C10 18393.199

typedef enum class e_DieletricModel {
    constant,
    four_r,
    none,
} DieletricModel;

__device__ double distance_squared(double x1, double x2,
                                   double y1, double y2,
                                   double z1, double z2) {

    return pow(x2 - x1, 2.0) + pow(y2 - y1, 2.0) + pow(z2 - z1, 2.0);
}

__device__ double distance(double x1, double x2,
                                   double y1, double y2,
                                   double z1, double z2) {

    return sqrt(pow(x2 - x1, 2.0) + pow(y2 - y1, 2.0) + pow(z2 - z1, 2.0));
}

__device__ double angle(double x1, double x2, double x3,
                        double y1, double y2, double y3,
                        double z1, double z2, double z3 ) {

    double ab = sqrt(distance_squared(x1, x2, y1, y2, z1, z2));
    double ac = sqrt(distance_squared(x1, x3, y1, y3, z1, z3));
    double bc = sqrt(distance_squared(x2, x3, y2, y3, z2, z3));

    return 180.0/M_PI * acos((pow(ab, 2.0) + pow(bc, 2.0) - pow(ac, 2.0))/ (2*ab*bc));
}

__global__
void compute_grid(int npointsx, int npointsy, int npointsz,
                                  double grid_spacing,
                                  double xbegin, double ybegin, double zbegin,
                                  double* xyz, int N,//
                                  double deltaij_es6,
                                  double deltaij6,
                                  double diel,
                                  double sigma,
                                  double solvation_alpha,
                                  double solvation_beta,
                                  double* charges,
                                  double* radii,
                                  double* epsilons_sqrt,
                                  int NHBdonors,
                                  int NHBacceptors,
                                  int* HBdonor1,
                                  int* HBdonor2,
                                  int* HBacceptors,
                                  double* out_elec_grid,
                                  double* out_vdwA_grid,
                                  double* out_vdwB_grid,
                                  double* out_rec_solv_grid,
                                  double* out_solv_grid,
                                  double* out_hbd_grid,
                                  double* out_hba_grid) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    if (i < npointsx && j < npointsy && k < npointsz){
        int index = k + (j*npointsx) + (npointsx*npointsy*i);
        double x = i*grid_spacing + xbegin;
        double y = j*grid_spacing + ybegin;
        double z = k*grid_spacing + zbegin;

        double elec = 0.0;
        double solv = 0.0;
        double rec_solv=0.0;
        double vdwA = 0.0;
        double vdwB = 0.0;
        double hb_donor = 0.0;
        double hb_acceptor = 0.0;

        for (int a = 0; a < N; a++){
            double d2 = distance_squared(x, xyz[a*3 + 0], y, xyz[a*3 + 1], z, xyz[a*3 + 2]);
            double d6 = d2*d2*d2;
            double denom_vdw = (d6 + deltaij6);
            double denom = cbrt(d6 + deltaij_es6);

            elec += (332.0*charges[a])/(diel*denom);
            solv += ((solvation_alpha * charges[a] * charges[a])+ solvation_beta) *  exp((-denom/(2*sigma*sigma))) / (sigma*sigma*sigma);
            rec_solv += (4.0/3.0) * PI * pow(radii[a], 3) * exp((-denom/(2*sigma*sigma))) / (sigma*sigma*sigma);

            
            vdwA += (4096.0 * epsilons_sqrt[a] * pow(radii[a], 6)) / (denom_vdw*denom_vdw);
            vdwB += ( 128.0 * epsilons_sqrt[a] * pow(radii[a], 3)) / denom_vdw;
        }

        for (int a=0; a<NHBdonors; a++){
            double d2 = distance_squared(xyz[HBdonor2[a]*3+0], x, xyz[HBdonor2[a]*3+1], y, xyz[HBdonor2[a]*3+2], z);
            double d10 = d2*d2*d2*d2*d2;
            double ang = angle(xyz[HBdonor1[a]*3+0], xyz[HBdonor1[a]*3+1], xyz[HBdonor1[a]*3+2], xyz[HBdonor2[a]*3+0],
            xyz[HBdonor2[a]*3+1], xyz[HBdonor2[a]*3+2], x, y, z);
            double angle_term = cos(ang * PI / 180.0) * cos(ang * PI / 180.0) * cos(ang * PI / 180.0) * cos(ang * PI / 180.0);
            hb_donor += ((HB_C12/(d10*d2)) - (HB_C10/d10)) * angle_term;
        }

        for (unsigned a=0; a<NHBacceptors; a++){
            double d2 = distance_squared(xyz[HBacceptors[a]*3+0], x, xyz[HBacceptors[a]*3+1], y, xyz[HBacceptors[a]*3+2], z);
            double d10 = d2*d2*d2*d2*d2;
            hb_acceptor += ((HB_C12/(d10*d2)) - (HB_C10/d10));
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
 
void invoke_compute_grid_softcore_HB_omp(Grid* grid, Mol2* rec) {

    double *d_xyz, *d_charges, *d_radii, *d_epsilon_sqrt ;
    double *d_out_elec_grid, *d_out_vdwA_grid, *d_out_vdwB_grid, *d_out_rec_solv_grid, *d_out_solv_grid, *d_out_hbd_grid, *d_out_hba_grid;

    cudaMalloc(&d_xyz, rec->xyz.size() * rec->xyz[0].size() * sizeof(double));
    cudaMalloc(&d_charges, rec->charges.size() * sizeof(double));
    cudaMalloc(&d_radii, rec->radii.size()* sizeof(double));
    cudaMalloc(&d_epsilon_sqrt, rec->epsilons_sqrt.size()*sizeof(double));
    

    size_t size_bytes {grid->npointsx * grid->npointsy * grid->npointsz * sizeof(double)};

    cudaMalloc(&d_out_elec_grid, size_bytes);
    cudaMalloc(&d_out_vdwA_grid, size_bytes);
    cudaMalloc(&d_out_vdwB_grid, size_bytes);
    cudaMalloc(&d_out_solv_grid, size_bytes);
    cudaMalloc(&d_out_rec_solv_grid, size_bytes);
    cudaMalloc(&d_out_hbd_grid, size_bytes);
    cudaMalloc(&d_out_hba_grid, size_bytes);

    double* xyz_arr = new double[rec->xyz.size()*3];

    for (int i = 0; i < rec->xyz.size(); i++){
        xyz_arr[3*i] = rec->xyz[i][0];
        xyz_arr[(3*i)+1] = rec->xyz[i][1];
        xyz_arr[(3*i)+2] = rec->xyz[i][2];
    }

    int* d_HBdonnor1; 
    int* d_HBdonnor2; 
    int* HBdonnor1 = new int[rec->HBdonors.size()];
    int* HBdonnor2 = new int[rec->HBdonors.size()];
    int* d_HBacceptors = new int[rec->HBacceptors.size()];

    for (int i=0; i< rec->HBdonors.size(); i++){
        HBdonnor1[i] = rec->HBdonors[i][0];
        HBdonnor2[i] = rec->HBdonors[i][1];
    }

    cudaMalloc(&d_HBdonnor1, rec->HBdonors.size() * sizeof(int));
    cudaMalloc(&d_HBdonnor2, rec->HBdonors.size() * sizeof(int));
    cudaMalloc(&d_HBacceptors, rec->HBacceptors.size() * sizeof(int));

    cudaMemcpy(d_xyz, xyz_arr, rec->xyz.size() * rec->xyz[0].size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_charges, rec->charges.data(), rec->charges.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_radii, rec->radii.data(), rec->radii.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_epsilon_sqrt, rec->epsilons_sqrt.data(), rec->epsilons_sqrt.size() * sizeof(double), cudaMemcpyHostToDevice);

    cudaMemcpy(d_HBdonnor1, HBdonnor1, rec->HBdonors.size()*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_HBdonnor2, HBdonnor2, rec->HBdonors.size()*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_HBacceptors, rec->HBacceptors.data(), rec->HBacceptors.size()*sizeof(int), cudaMemcpyHostToDevice);

    //Só para teste, preferencialmente irá escrever diretamente para o Grid
    double* out_elec_grid = (double*) malloc(size_bytes);
    double* out_vdwA_grid = (double*) malloc(size_bytes);
    double* out_vdwB_grid = (double*) malloc(size_bytes);
    double* out_solv_grid = (double*) malloc(size_bytes);
    double* out_rec_solv_grid = (double*) malloc(size_bytes);
    double* out_hbd_grid = (double*) malloc(size_bytes);
    double* out_hba_grid = (double*) malloc(size_bytes);

    dim3 blockdims(8,8,8);
    dim3 griddims(8,8,8);

    printf("Entering kernel\n");
    printf("Last error: %s\n", cudaGetErrorString(cudaGetLastError()));
    //TODO: ver melhor estes parâmetros de launch
    compute_grid<<<griddims, blockdims>>>(grid->npointsx, grid->npointsy, grid->npointsz,
                                                           grid->grid_spacing,
                                                           grid->xbegin, grid->ybegin, grid->zbegin,
                                                           d_xyz, rec->N,
                                                           grid->Input->deltaij_es6,
                                                           grid->Input->deltaij6,
                                                           grid->Input->diel,
                                                           grid->Input->sigma,
                                                           grid->Input->solvation_alpha,
                                                           grid->Input->solvation_beta,
                                                           d_charges,
                                                           d_radii,
                                                           d_epsilon_sqrt,
                                                           rec->HBdonors.size(),
                                                           rec->HBacceptors.size(),
                                                           d_HBdonnor1,
                                                           d_HBdonnor2,
                                                           d_HBacceptors,
                                                           d_out_elec_grid,
                                                           d_out_vdwA_grid,
                                                           d_out_vdwB_grid,
                                                           d_out_rec_solv_grid,
                                                           d_out_solv_grid,
                                                           d_out_hbd_grid,
                                                           d_out_hba_grid);

    printf("Kernel has ended\n");
    printf("Last error: %s\n", cudaGetErrorString(cudaGetLastError()));

    cudaMemcpy(out_elec_grid, d_out_elec_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_vdwA_grid, d_out_vdwA_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_vdwB_grid, d_out_vdwB_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_solv_grid, d_out_solv_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_rec_solv_grid, d_out_rec_solv_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_hbd_grid, d_out_hbd_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_hba_grid, d_out_hba_grid, size_bytes, cudaMemcpyDeviceToHost);

    for (int i=0; i< grid->npointsx/10; i++){
        for (int j=0; j< grid->npointsy/10; j++){
            for (int k=0; k< grid->npointsz/10; k++){
                int index = k + (j*grid->npointsx) + (grid->npointsx*grid->npointsy*i);
                printf("ElecGrid[%3d][%3d][%3d] = %8.3f \t CudaGrid[%8d] = %8.3f \t Diff = %8.3f\n", i, j, k, grid->elec_grid[i][j][k], index, out_elec_grid[index], abs(grid->elec_grid[i][j][k]-abs(out_elec_grid[index])));
                printf("VdwAGrid[%3d][%3d][%3d] = %8.3f \t CudaGrid[%8d] = %8.3f \t Diff = %8.3f\n", i, j, k, grid->vdwA_grid[i][j][k], index, out_vdwA_grid[index], abs(grid->vdwA_grid[i][j][k]-abs(out_vdwA_grid[index])));
                printf("VdwBGrid[%3d][%3d][%3d] = %8.3f \t CudaGrid[%8d] = %8.3f \t Diff = %8.3f\n", i, j, k, grid->vdwB_grid[i][j][k], index, out_vdwB_grid[index], abs(grid->vdwB_grid[i][j][k]-abs(out_vdwB_grid[index])));
                printf("RecSGrid[%3d][%3d][%3d] = %8.3f \t CudaGrid[%8d] = %8.3f \t Diff = %8.3f\n", i, j, k, grid->rec_solv_gauss[i][j][k], index, out_rec_solv_grid[index], abs(grid->rec_solv_gauss[i][j][k]-abs(out_rec_solv_grid[index])));
                printf("SolvGrid[%3d][%3d][%3d] = %8.3f \t CudaGrid[%8d] = %8.3f \t Diff = %8.3f\n", i, j, k, grid->solv_gauss[i][j][k], index, out_solv_grid[index], abs(grid->solv_gauss[i][j][k]-abs(out_solv_grid[index])));
                printf("HBAcGrid[%3d][%3d][%3d] = %8.3f \t CudaGrid[%8d] = %8.3f \t Diff = %8.3f\n", i, j, k, grid->hb_acceptor_grid[i][j][k], index, out_hba_grid[index], abs(grid->hb_acceptor_grid[i][j][k]-abs(out_hba_grid[index])));
                printf("HBDoGrid[%3d][%3d][%3d] = %8.3f \t CudaGrid[%8d] = %8.3f \t Diff = %8.3f\n", i, j, k, grid->hb_donor_grid[i][j][k], index, out_hbd_grid[index], abs(grid->hb_donor_grid[i][j][k]-abs(out_hbd_grid[index])));
                printf("\n\n");
            }
        }
    }

    delete [] (xyz_arr);
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
