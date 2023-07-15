#include <iostream>
#include <vector>
#include <cassert>
#include "../../src/pyGrid.h"
#include "../../src/pyMol2.h"


// static void HandleError( cudaError_t err,
//                          const char *file,
//                          int line ) {
//     if (err != cudaSuccess) {
//         printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
//                 file, line );
//         exit( EXIT_FAILURE );
//     }
// 

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#define DISTANCE_SQUARED(x1, x2, y1, y2, z1, z2) (pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0) + pow(z1 - z2, 2.0))

#define HB_C12 55332.873
#define HB_C10 18393.199

// Will interpret the ditric_model attribute as int because
// string handling through the kernel may not be straight forward
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

__device__ double angle(double x1, double x2, double x3,
                        double y1, double y2, double y3,
                        double z1, double z2, double z3 ) {

    double ab = sqrt(distance_squared(x1, x2, y1, y2, z1, z2));
    double ac = sqrt(distance_squared(x1, x3, y1, y3, z1, z3));
    double bc = sqrt(distance_squared(x2, x3, y2, y3, z2, z3));

    return 180.0/PI * acos((pow(ab, 2.0) + pow(bc, 2.0) - pow(ac, 2.0))/ (2*ab*bc));
}

template <typename T>
void print_values_3D(const T* array, size_t size_x, size_t size_y, size_t size_z) {

    for (int i = 0 ; i < size_x ; i++)
    {
        for (int j = 0 ; j < size_y ; j++)
        {
            for (int k = 0 ; k < size_z ; k++)
            {
                printf("%.4lf\t", static_cast<float>(array[(i * size_x + j) * size_y + k]));
            }
        }
        printf("\n");
    }
    printf("\n");
}

void print_difference(const std::vector<std::vector<std::vector<double>>>& vec, double* arr) {

    size_t x {vec.size()};
    size_t y {vec[0].size()};
    size_t z {vec[0][0].size()};

    printf("x %ld, y %ld, z %ld\n", x, y, z);

    double diff{};
    double max{0};
    double sum{0};
    int count{0};
    int original_abnormal{0};
    int new_abnormal{0};
    int lin_index{};
    int zeroes_original{0};
    int zeroes_new{0};

    double tolerance{0};

    for (int i = 0 ; i < x ; i++)
    {
        for (int j = 0 ; j < y ; j++)
        {
            for (int k = 0 ; k < z ; k++)
            {
                lin_index = k + (j*y) + (x*y*i);
                //printf("index %d\n", lin_index);
                if (!std::isfinite(vec[i][j][k]))
                {
                    original_abnormal++;
                    continue;
                }
                if (!std::isfinite(arr[lin_index]))
                {
                    new_abnormal++;
                    continue;
                }
                if (abs(vec[i][j][k]) <= tolerance)
                {
                    zeroes_original++;
                    continue;
                }
                if (abs(arr[lin_index]) <= tolerance) 
                {
                    zeroes_new++;
                    continue;
                }

                diff = abs(vec[i][j][k] - arr[lin_index]) / std::max(vec[i][j][k], arr[lin_index]);
                sum += diff;

                // if (diff > 1) printf("Diff > 1 --> original: %lf\t new: %lf\n", vec[i][j][k], arr[lin_index]);
                //  printf("Diff > 1 --> original: %lf\t new: %lf\n", vec[i][j][k], arr[lin_index]);

                if (diff > max) max = diff;
                count++;
            }
        }
    }
    printf("Avg: %lf\nMax diff: %lf\nInfs or NaNs in original: %d\nInfs or Nans in new: %d\nzeroes_original: %d\nzeroes_new: %d\n\n",
           sum/count, max, original_abnormal, new_abnormal, zeroes_original, zeroes_new);
}

__global__
void compute_grid_softcore_HB_omp(int npointsx, int npointsy, int npointsz,
                                   double grid_spacing,
                                   double xbegin, double ybegin, double zbegin,
                                   DieletricModel dielectric_model,
                                   double deltaij_es6, double deltaij6,
                                   double solvation_alpha, double solvation_beta,
                                   double sigma, double diel,
                                   int N,
                                   int xyz_w, double* xyz, //
                                   double* charges,
                                   double* radii,
                                   double* epsilons_sqrt,
                                   int HBacceptors_size, int* HBacceptors,
                                   int HBdonors_w, int HBdonors_h, int* HBdonors, //
                                   double* out_elec_grid,
                                   double* out_vdwA_grid,
                                   double* out_vdwB_grid,
                                   double* out_solv_gauss,
                                   double* out_rec_solv_gauss,
                                   double* out_hb_donor_grid,
                                   double* out_hb_acceptor_grid) {

    int idxi=threadIdx.x + blockIdx.x * blockDim.x;
    int idxj=threadIdx.y + blockIdx.y * blockDim.y;
    int idxk=threadIdx.z + blockIdx.z * blockDim.z;
    int xstride = blockDim.x * gridDim.x; // Assuming they are the same for x,y,z
    int ystride = blockDim.y * gridDim.y; // Assuming they are the same for x,y,z
    int zstride = blockDim.z * gridDim.z; // Assuming they are the same for x,y,z
    double deltaij_es3 = sqrt(deltaij_es6);

    for (int i=idxi; i<npointsx; i+= xstride){
        double x = i*grid_spacing + xbegin;

        for(int j=idxj; j<npointsy; j+= ystride){
            double y = j*grid_spacing + ybegin;

            for (int k=idxk; k<npointsz; k+= zstride){
                double z = k*grid_spacing + zbegin;

                double elec = 0.0;
                double vdwA = 0.0;
                double vdwB = 0.0;
                double solv = 0.0;
                double rec_solv = 0.0;

                for (int a = 0; a < N; a++) {
                    double d2 = distance_squared(x, xyz[a*xyz_w + 0], y, xyz[a*xyz_w + 1], z, xyz[a*xyz_w + 2]);
                    // double d2 = DISTANCE_SQUARED(x, xyz[a*N + 0], y, xyz[a*N + 1], z, xyz[a*N + 2]);
                    // double d2 = DISTANCE_SQUARED(x, xyz[a*N + 0], y, xyz[a*N + 1], z, xyz[a*N + 2]);
                    double d6 = d2*d2*d2;
                    double denom = 0.0;

                    if (dielectric_model == DieletricModel::constant) {
                        denom = cbrt(d6 + deltaij_es6);
                        elec += (332.0*charges[a])/(diel*denom);
                        solv += ((solvation_alpha * charges[a] * charges[a]) + solvation_beta)
                                * exp((-denom/(2*pow(sigma, 2.0)))) / (pow(sigma, 3.0));
                        rec_solv += (4.0/3.0) * M_PI * pow(radii[a], 3.0)
                                * exp((-denom/(2*pow(sigma, 2.0)))) / (pow(sigma, 3.0));
                    }
                    else {
                        // if (i == 0 && j == 0 && k == 0) printf("iter: %d, d2: %lf\n", a, d2);
                        denom = cbrt(d6 + deltaij_es6);
                        // if (a == 0) printf("%lf\n", denom);
                        elec += (332.0*charges[a])/(diel*denom);
                        solv += ((solvation_alpha * charges[a] * charges[a]) + solvation_beta)
                                * exp((-denom/(2*pow(sigma, 2.0)))) / (pow(sigma, 3.0));
                        rec_solv += (4.0/3.0) * M_PI * pow(radii[a], 3.0)
                                * exp((-denom/(2*pow(sigma, 2.0)))) / (pow(sigma, 3.0));
                    }
                    
                    denom = (d6 + deltaij6);
                    vdwA += (4096*epsilons_sqrt[a] * pow(radii[a], 6.0)) / pow(denom, 2.0);
                    vdwB += (128*epsilons_sqrt[a] * pow(radii[a], 3.0)) / denom;
                }
                // printf("ended inner for loop\n");
                double hb_donor = 0;
                double hb_acceptor = 0;

                for (int m = 0; m < HBdonors_h; m++) {
                    int HB0 = HBdonors[m * HBdonors_w + 0];
                    int HB1 = HBdonors[m * HBdonors_w + 1];

                    double d2 = distance_squared(xyz[HB1*xyz_w + 0], x, xyz[HB1*xyz_w + 1], y, xyz[HB1*xyz_w + 2], z);
                    double d10 = d2*d2*d2*d2*d2;
                    double ang = angle(xyz[HB0*xyz_w + 0], xyz[HB0*xyz_w + 1], xyz[HB0*xyz_w + 2],
                            xyz[HB1*xyz_w + 0], xyz[HB1*xyz_w + 1], xyz[HB1*xyz_w + 2], x, y, z);
                    double angle_term = (pow(cos(ang * PI / 180.0), 4.0));
                    hb_donor += (HB_C12/(d10*d2)) - (HB_C10/d10);
                }


                for (int n = 0 ; n < HBacceptors_size ; n++) {
                    double d2 = distance_squared(xyz[HBacceptors[n]*xyz_w + 0], x, xyz[HBacceptors[n]*xyz_w + 1], y, xyz[HBacceptors[n]*xyz_w + 2], z);
                    double d10 = pow(d2, 5.0);
                    hb_acceptor += (HB_C12/(d10*d2)) - (HB_C10/d10);
                }
                // printf("%lf\n", elec);
                out_elec_grid[(i * npointsz + j) * npointsy + k] = elec;
                out_vdwA_grid[(i * npointsz + j) * npointsy + k] = vdwA;
                out_vdwB_grid[(i * npointsz + j) * npointsy + k] = vdwB;
                out_solv_gauss[(i * npointsz + j) * npointsy + k] = solv;
                out_rec_solv_gauss[(i * npointsz + j) * npointsy + k] = rec_solv;
                out_hb_donor_grid[(i * npointsz + j) * npointsy + k] = hb_donor;
                out_hb_acceptor_grid[(i * npointsz + j) * npointsy + k] = hb_acceptor;
            }
        }
    }
}

void invoke_compute_grid_softcore_HB_omp(Grid* grid, Mol2* rec) {

    double *d_xyz, *d_charges, *d_radii, *d_epsilons_sqrt;
    int *d_hbacceptors, *d_hbdonors;

    double *d_out_elec_grid, *d_out_vdwa_grid, *d_out_vdwb_grid, *d_out_solv_gauss, *d_out_rec_solv_gauss, *d_out_hb_donor_grid, *d_out_hb_acceptor_grid;

    cudaMalloc(&d_xyz, rec->xyz.size() * rec->xyz[0].size() * sizeof(double));
    cudaMalloc(&d_charges, rec->charges.size() * sizeof(double));
    cudaMalloc(&d_radii, rec->radii.size() * sizeof(double));
    cudaMalloc(&d_epsilons_sqrt, rec->epsilons_sqrt.size() * sizeof(double));

    cudaMalloc(&d_hbdonors, rec->HBdonors.size() * rec->HBdonors[0].size() * sizeof(int));
    cudaMalloc(&d_hbacceptors, rec->HBacceptors.size() * sizeof(int));

    size_t size_bytes {grid->npointsx * grid->npointsy * grid->npointsz * sizeof(double)};

    cudaMalloc(&d_out_elec_grid, size_bytes);
    cudaMalloc(&d_out_vdwa_grid, size_bytes);
    cudaMalloc(&d_out_vdwb_grid, size_bytes);
    cudaMalloc(&d_out_solv_gauss, size_bytes);
    cudaMalloc(&d_out_rec_solv_gauss, size_bytes);
    cudaMalloc(&d_out_hb_donor_grid, size_bytes);
    cudaMalloc(&d_out_hb_acceptor_grid, size_bytes);

    // double* xyz_arr = (double*) malloc(rec->xyz.size() * 3 * sizeof(double));
    
    double* xyz_arr = new double[rec->xyz.size()*3];

    for (int i = 0; i < rec->xyz.size(); i++){
        xyz_arr[3*i] = rec->xyz[i][0];
        xyz_arr[(3*i)+1] = rec->xyz[i][1];
        xyz_arr[(3*i)+2] = rec->xyz[i][2];
    }

    double* HBdonors_arr = new double[rec->HBdonors.size()*2];
    for (int i = 0; i < rec->HBdonors.size(); i++){
        for (int j = 0; j < rec->HBdonors[0].size(); j++){
            HBdonors_arr[rec->HBdonors[0].size()*i + j] = rec->HBdonors[i][j];
        }
    }
   // for (int i = 0 ; i < rec->HBdonors.size(); i++)
   // {
   //          HBdonors_arr[2*i + 0] = rec->HBdonors[i][0];
   //          HBdonors_arr[2*i + 1] = rec->HBdonors[i][1];
   // }

    rec->charges.shrink_to_fit();
    rec->radii.shrink_to_fit();
    rec->epsilons_sqrt.shrink_to_fit();
    rec->HBacceptors.shrink_to_fit();

    cudaMemcpy(d_xyz, xyz_arr, rec->xyz.size() * rec->xyz[0].size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_charges, rec->charges.data(), rec->charges.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_radii, rec->radii.data(), rec->radii.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_epsilons_sqrt, rec->epsilons_sqrt.data(), rec->epsilons_sqrt.size() * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_hbacceptors, rec->HBacceptors.data(), rec->HBacceptors.size() * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_hbdonors, HBdonors_arr, rec->HBdonors.size() * rec->HBdonors[0].size() * sizeof(int), cudaMemcpyHostToDevice);

    DieletricModel dieletric_model{};

    if (grid->Input->dielectric_model == "constant")
        dieletric_model = DieletricModel::constant;
    else
        dieletric_model = DieletricModel::none;


    //Só para teste, preferencialmente irá escrever diretamente para o Grid
    double* out_elec_grid = (double*) malloc(size_bytes);
    double* out_vdwa_grid = (double*) malloc(size_bytes);
    double* out_vdwb_grid = (double*) malloc(size_bytes);
    double* out_solv_gauss = (double*) malloc(size_bytes);
    double* out_rec_solv_gauss = (double*) malloc(size_bytes);
    double* out_hb_donor_grid = (double*) malloc(size_bytes);
    double* out_hb_acceptor_grid = (double*) malloc(size_bytes);

    dim3 blockdims(4,4,4);
    dim3 griddims(8,8,8);

    printf("Entering kernel\n");
    printf("Last error: %s\n", cudaGetErrorString(cudaGetLastError()));
    //TODO: ver melhor estes parâmetros de launch
    compute_grid_softcore_HB_omp<<<griddims, blockdims>>>(grid->npointsx, grid->npointsy, grid->npointsz,
                                                           grid->grid_spacing,
                                                           grid->xbegin, grid->ybegin, grid->zbegin,
                                                           dieletric_model,
                                                           grid->Input->deltaij_es6, grid->Input->deltaij6,
                                                           grid->Input->solvation_alpha, grid->Input->solvation_beta,
                                                           grid->Input->sigma, grid->Input->diel,
                                                           rec->N,
                                                           rec->xyz[0].size(), d_xyz,
                                                           d_charges,
                                                           d_radii,
                                                           d_epsilons_sqrt,
                                                           rec->HBacceptors.size(),d_hbacceptors,
                                                           rec->HBdonors[0].size(), rec->HBdonors.size(), d_hbdonors,
                                                           d_out_elec_grid,
                                                           d_out_vdwa_grid,
                                                           d_out_vdwb_grid,
                                                           d_out_solv_gauss,
                                                           d_out_rec_solv_gauss,
                                                           d_out_hb_donor_grid,
                                                           d_out_hb_acceptor_grid);

    double out_rec_si = 0.0;
    for (int a = 0; a < rec->N ; a++)
    {
        out_rec_si += (grid->Input->solvation_alpha * pow(rec->charges[a], 2.0)) + grid->Input->solvation_beta;
    }

    printf("Kernel has ended\n");
    printf("Last error: %s\n", cudaGetErrorString(cudaGetLastError()));

    cudaMemcpy(out_elec_grid, d_out_elec_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_vdwa_grid, d_out_vdwa_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_vdwb_grid, d_out_vdwb_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_solv_gauss, d_out_solv_gauss, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_rec_solv_gauss, d_out_rec_solv_gauss, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_hb_donor_grid, d_out_hb_donor_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_hb_acceptor_grid, d_out_hb_acceptor_grid, size_bytes, cudaMemcpyDeviceToHost);

   //  print_values_3D(out_elec_grid, grid->npointsx, grid->npointsy, grid->npointsz);
    // print_values_3D(out_vdwa_grid, grid->npointsx, grid->npointsy, grid->npointsz);
    // print_values_3D(out_vdwb_grid, grid->npointsx, grid->npointsy, grid->npointsz);
    // print_values_3D(out_solv_gauss, grid->npointsx, grid->npointsy, grid->npointsz);
    // print_values_3D(out_rec_solv_gauss, grid->npointsx, grid->npointsy, grid->npointsz);
    // print_values_3D(out_hb_donor_grid, grid->npointsx, grid->npointsy, grid->npointsz);
    // print_values_3D(out_hb_acceptor_grid, grid->npointsx, grid->npointsy, grid->npointsz);
    printf("\nrec_si from invoke function: %lf\n\n", out_rec_si);

    print_difference(grid->elec_grid, out_elec_grid);
//    print_difference(grid->vdwA_grid, out_vdwa_grid);
//    print_difference(grid->vdwB_grid, out_vdwb_grid);
    print_difference(grid->solv_gauss, out_solv_gauss);
    print_difference(grid->rec_solv_gauss, out_rec_solv_gauss);
//    print_difference(grid->hb_donor_grid, out_hb_donor_grid);
//    print_difference(grid->hb_acceptor_grid, out_hb_acceptor_grid);

    delete [] (xyz_arr);
    delete [] (HBdonors_arr);

    free(out_vdwa_grid);
    free(out_vdwb_grid);
    free(out_solv_gauss);
    free(out_hb_donor_grid);
    free(out_hb_acceptor_grid);

    cudaFree(d_xyz);
    cudaFree(d_charges);
    cudaFree(d_radii);
    cudaFree(d_epsilons_sqrt);
    cudaFree(d_hbdonors);
    cudaFree(d_hbacceptors);
    cudaFree(d_out_vdwa_grid);
    cudaFree(d_out_vdwb_grid);
    cudaFree(d_out_solv_gauss);
    cudaFree(d_out_hb_donor_grid);
    cudaFree(d_out_hb_acceptor_grid);

    printf("Invoking finished\n");

}
