#include "math.h"
#include <iostream>
#include <vector>
#include <cmath>
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
// }

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

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

    return 180.0/M_PI * acos((pow(ab, 2.0) + pow(bc, 2.0) - pow(ac, 2.0))/ (2*ab*bc));
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
    int original_anormal{0};
    int new_anormal{0};
    int lin_index{};
    int zeroes_original{0};
    int zeroes_new{0};

    for (int i = 0 ; i < x ; i++)
    {
        for (int j = 0 ; j < y ; j++)
        {
            for (int k = 0 ; k < z ; k++)
            {
                lin_index = (i * z + j) * y + k;
                //printf("index %d\n", lin_index);
                if (!std::isfinite(vec[i][j][k]))
                {
                    original_anormal++;
                    continue;
                }
                if (!std::isfinite(arr[lin_index]))
                {
                    new_anormal++;
                    continue;
                }
                if (vec[i][j][k] == 0.0)
                {
                    zeroes_original++;
                    continue;
                }
                if (arr[lin_index] == 0.0) zeroes_new++;

                diff = abs(vec[i][j][k] - arr[lin_index]) / vec[i][j][k];
                sum += diff;

                if (diff > max) max = diff;
                count++;
            }
        }
    }
    printf("Avg: %lf\nMax diff: %lf\nInfs or NaNs in original: %d\nInfs or Nans in new: %d\nzeroes_original: %d\nzeroes_new: %d\n\n",
           sum/count, max, original_anormal, new_anormal, zeroes_original, zeroes_new);
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
                                  int HBdonors_w, int* HBdonors, //
                                  double* out_elec_grid,
                                  double* out_vdwA_grid,
                                  double* out_vdwB_grid,
                                  double* out_solv_gauss,
                                  double* out_rec_solv_gauss,
                                  double* out_hb_donor_grid,
                                  double* out_hb_acceptor_grid,
                                  int* out_rec_si) {

    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int I = 0;
    int J = 0;
    int K = 0;

    if (i < npointsx && j < npointsy && k < npointsz)
    {
        double x = i*grid_spacing + xbegin;
        double y = j*grid_spacing + ybegin;
        double z = k*grid_spacing + zbegin;

        double elec = 1.0;
        double vdwA = 0.0;
        double vdwB = 0.0;
        double solv = 0.0;
        double rec_solv = 0.0;

        for (int a = 0; a < N; a++)
        {
            double d2 = distance_squared(x, xyz[a*xyz_w + 0], y, xyz[a*xyz_w + 1], z, xyz[a*xyz_w + 2]);
            double d6 = d2*d2*d2;
            double denom = 0.0;

            if (dielectric_model == DieletricModel::constant)
            {
                denom = pow(d6 + deltaij_es6, 1/3);
                elec += (332.0*charges[a])/(diel*denom);
                solv += ((solvation_alpha * charges[a] * charges[a]) + solvation_beta)
                        * exp((-denom/(2*pow(sigma, 2.0)))) / (pow(sigma, 3.0));
                rec_solv += (4.0/3.0) * M_PI * pow(radii[a], 3.0)
                        * exp((-denom/(2*pow(sigma, 2.0)))) / (pow(sigma, 3.0));
            }
            //TODO: condicional desnecessário?
            else
            {
                denom = pow(d6 + deltaij_es6, 1/3);
                elec += (332.0*charges[a])/(diel*denom);
                solv += ((solvation_alpha * charges[a] * charges[a]) + solvation_beta)
                        * exp((-denom/(2*pow(sigma, 2.0)))) / (pow(sigma, 3.0));
                rec_solv += (4.0/3.0) * M_PI * pow(radii[a], 3.0)
                        * exp((-denom/(2*pow(sigma, 2.0)))) / (pow(sigma, 3.0));
            }

            denom = (d6 + deltaij6);
            vdwA = (4096*epsilons_sqrt[a] * pow(radii[a], 6.0)) / pow(denom, 2.0);
            vdwB = (128*epsilons_sqrt[a] * pow(radii[a], 3.0)) / denom;
        }

        double hb_donor = 0;

        for (int m = 0; m < sizeof(HBdonors)/sizeof(HBdonors[0]); m++)
        {
            int HB0 = HBdonors[m * HBdonors_w + 0];
            int HB1 = HBdonors[m * HBdonors_w + 1];

            double d2 = distance_squared(xyz[HB1*xyz_w + 0], x, xyz[HB1*xyz_w + 0], y, xyz[HB1*xyz_w + 0], z);
            double d10 = d2*d2*d2*d2*d2;
            double ang = angle(xyz[HB0*xyz_w + 0], xyz[HB0*xyz_w + 1], xyz[HB0*xyz_w + 2],
                    xyz[HB1*xyz_w + 0], xyz[HB1*xyz_w + 1], xyz[HB1*xyz_w + 2], x, y, z);
            double angle_term = (pow(cos(ang * M_PI / 180.0), 4.0));
            hb_donor += (HB_C12/(d10*d2)) - (HB_C10/d10);
        }


        double hb_acceptor = 0;
        for (int n = 0 ; n < HBacceptors_size ; n++)
        {
            double d2 = distance_squared(xyz[HBacceptors[n]*xyz_w + 0], x, xyz[HBacceptors[n]*xyz_w + 1], y, xyz[HBacceptors[n]*xyz_w + 2], z);
            double d10 = pow(d2, 5.0);
            hb_acceptor += (HB_C12/(d10*d2)) - (HB_C10/d10);
        }

        printf("%lf", elec);

        out_elec_grid[(i * npointsx + j) * npointsy + k] = elec;
        out_vdwA_grid[(i * npointsx + j) * npointsy + k] = vdwA;
        out_vdwB_grid[(i * npointsx + j) * npointsy + k] = vdwB;
        out_solv_gauss[(i * npointsx + j) * npointsy + k] = solv;
        out_rec_solv_gauss[(i * npointsx + j) * npointsy + k] = rec_solv;
        out_hb_donor_grid[(i * npointsx + j) * npointsy + k] = hb_donor;
        out_hb_acceptor_grid[(i * npointsx + j) * npointsy + k] = hb_acceptor;
    }

    out_rec_si[0] = 0.0;
    for (int a = 0; a < N ; a++)
    {
        out_rec_si[0] += (solvation_alpha * pow(charges[a], 2.0)) + solvation_beta;
    }
}

__global__
void compute_grid_softcore_HB_omp2(int npointsx, int npointsy, int npointsz,
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
                                   int HBdonors_w, int* HBdonors, //
                                   double* out_elec_grid,
                                   double* out_vdwA_grid,
                                   double* out_vdwB_grid,
                                   double* out_solv_gauss,
                                   double* out_rec_solv_gauss,
                                   double* out_hb_donor_grid,
                                   double* out_hb_acceptor_grid,
                                   double* out_rec_si) {

    int idxi=threadIdx.x, idxj=threadIdx.y, idxk=threadIdx.x;
    double deltaij_es3 = sqrt(deltaij_es6);

    for (int i=idxi; i<npointsx; i++){
        double x = i*grid_spacing + xbegin;

        for(int j=idxj; j<npointsy; j++){
            double y = j*grid_spacing + ybegin;

            for (int k=idxk; k<npointsz; k++){
                double z = k*grid_spacing + zbegin;

                double elec = 0.0;
                double vdwA = 0.0;
                double vdwB = 0.0;
                double solv = 0.0;
                double rec_solv = 0.0;

                for (int a = 0; a < N; a++){
                    double d2 = distance_squared(x, xyz[a*xyz_w + 0], y, xyz[a*xyz_w + 1], z, xyz[a*xyz_w + 2]);
                    double d6 = d2*d2*d2;
                    double denom = 0.0;

                    if (dielectric_model == DieletricModel::constant){
                        double d3 = sqrt(d6);
                        denom = pow(d3 + deltaij_es3, 1/3);
                        elec += (332.0*charges[a])/(diel*denom);
                        solv += ((solvation_alpha * charges[a] * charges[a]) + solvation_beta)
                                * exp((-denom/(2*pow(sigma, 2.0)))) / (pow(sigma, 3.0));
                        rec_solv += (4.0/3.0) * M_PI * pow(radii[a], 3.0)
                                * exp((-denom/(2*pow(sigma, 2.0)))) / (pow(sigma, 3.0));
                    }
                    else{
                        denom = pow(d6 + deltaij_es6, 1/3);
                        elec += (332.0*charges[a])/(diel*denom);
                        solv += ((solvation_alpha * charges[a] * charges[a]) + solvation_beta)
                                * exp((-denom/(2*pow(sigma, 2.0)))) / (pow(sigma, 3.0));
                        rec_solv += (4.0/3.0) * M_PI * pow(radii[a], 3.0)
                                * exp((-denom/(2*pow(sigma, 2.0)))) / (pow(sigma, 3.0));
                    }

                    denom = (d6 + deltaij6);
                    vdwA = (4096*epsilons_sqrt[a] * pow(radii[a], 6.0)) / pow(denom, 2.0);
                    vdwB = (128*epsilons_sqrt[a] * pow(radii[a], 3.0)) / denom;
                }

                double hb_donor = 0;

                for (int m = 0; m < HBdonors_w; m++)
                {
                    int HB0 = HBdonors[m * HBdonors_w + 0];
                    int HB1 = HBdonors[m * HBdonors_w + 1];

                    double d2 = distance_squared(xyz[HB1*xyz_w + 0], x, xyz[HB1*xyz_w + 1], y, xyz[HB1*xyz_w + 2], z);
                    double d10 = d2*d2*d2*d2*d2;
                    double ang = angle(xyz[HB0*xyz_w + 0], xyz[HB0*xyz_w + 1], xyz[HB0*xyz_w + 2],
                            xyz[HB1*xyz_w + 0], xyz[HB1*xyz_w + 1], xyz[HB1*xyz_w + 2], x, y, z);
                    double angle_term = (pow(cos(ang * M_PI / 180.0), 4.0));
                    hb_donor += (HB_C12/(d10*d2)) - (HB_C10/d10);
                }

                double hb_acceptor = 0;
                for (int n = 0 ; n < HBacceptors_size ; n++)
                {
                    double d2 = distance_squared(xyz[HBacceptors[n]*xyz_w + 0], x, xyz[HBacceptors[n]*xyz_w + 1], y, xyz[HBacceptors[n]*xyz_w + 2], z);
                    double d10 = pow(d2, 5.0);
                    hb_acceptor += (HB_C12/(d10*d2)) - (HB_C10/d10);
                }

                out_elec_grid[k + (npointsy*j) + (npointsx*npointsy*i)] = elec;
                out_vdwA_grid[k + (npointsy*j) + (npointsx*npointsy*i)] = vdwA;
                out_vdwB_grid[k + (npointsy*j) + (npointsx*npointsy*i)] = vdwB;
                out_solv_gauss[k + (npointsy*j) + (npointsx*npointsy*i)] = solv;
                out_rec_solv_gauss[k + (npointsy*j) + (npointsx*npointsy*i)] = rec_solv;
                out_hb_donor_grid[k + (npointsy*j) + (npointsx*npointsy*i)] = hb_donor;
                out_hb_acceptor_grid[k + (npointsy*j) + (npointsx*npointsy*i)] = hb_acceptor;
            }
        }
    }
    out_rec_si[0] = 0.0;
    for (int a = 0; a < N ; a++)
    {
        out_rec_si[0] += (solvation_alpha * pow(charges[a], 2.0)) + solvation_beta;
    }
}

void invoke_compute_grid_softcore_HB_omp(Grid* grid, Mol2* rec) {

    double *d_xyz, *d_charges, *d_radii, *d_epsilons_sqrt;
    int *d_hbacceptors, *d_hbdonors;

    double *d_out_elec_grid, *d_out_vdwa_grid, *d_out_vdwb_grid, *d_out_solv_gauss, *d_out_rec_solv_gauss, *d_out_hb_donor_grid, *d_out_hb_acceptor_grid;
    double *d_out_rec_si;

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

    cudaMalloc(&d_out_rec_si, 1 * sizeof(double));

    //Deep copying C++ vectors into linearized arrays
    //TODO: xyz.size() vem antes de xyz[0].size() mesmo?
    double* xyz_arr = (double*) malloc(rec->xyz.size() * rec->xyz[0].size() * sizeof(double));
    for (int i = 0; i < rec->xyz.size(); i++)
    {
        for (int j = 0; j < rec->xyz[0].size(); j++)
        {
            xyz_arr[rec->xyz[0].size()*i + j] = rec->xyz[i][j];
        }
    }

    double* HBdonors_arr = (double*) malloc(rec->HBdonors.size() * rec->HBdonors[0].size() * sizeof(double));
    for (int i = 0; i < rec->HBdonors.size(); i++)
    {
        for (int j = 0; j < rec->HBdonors[0].size(); j++)
        {
            HBdonors_arr[rec->HBdonors[0].size()*i + j] = rec->HBdonors[i][j];
        }
    }


    //Guarantee that size matches the vector' capacity
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

    dim3 blockdims(1,1,1);
    dim3 griddims(1,1,1);

    printf("Entering kernel\n");
    printf("Last error: %s\n", cudaGetErrorString(cudaGetLastError()));
    //TODO: ver melhor estes parâmetros de launch
    compute_grid_softcore_HB_omp2<<<griddims, blockdims>>>(grid->npointsx, grid->npointsy, grid->npointsz,
                                                          grid->grid_spacing,
                                                          grid->xbegin, grid->ybegin, grid->zbegin,
                                                          dieletric_model,
                                                          grid->Input->deltaij_es6, grid->Input->deltaij6,
                                                          grid->Input->solvation_alpha, grid->Input->solvation_beta,
                                                          grid->Input->sigma, grid->Input->diel,
                                                          rec->N,
                                                          rec->xyz[0].size(), d_xyz,//fixme: pode ser o shape errado
            d_charges,
            d_radii,
            d_epsilons_sqrt,
            rec->HBacceptors.size(),d_hbacceptors,
            rec->HBdonors[0].size(), d_hbdonors,
            d_out_elec_grid,
            d_out_vdwa_grid,
            d_out_vdwb_grid,
            d_out_solv_gauss,
            d_out_rec_solv_gauss,
            d_out_hb_donor_grid,
            d_out_hb_acceptor_grid,
            d_out_rec_si);

    printf("Kernel has ended\n");
    printf("Last error: %s\n", cudaGetErrorString(cudaGetLastError()));

    cudaMemcpy(out_elec_grid, d_out_elec_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_vdwa_grid, d_out_vdwa_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_vdwb_grid, d_out_vdwb_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_solv_gauss, d_out_solv_gauss, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_rec_solv_gauss, d_out_rec_solv_gauss, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_hb_donor_grid, d_out_hb_donor_grid, size_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(out_hb_acceptor_grid, d_out_hb_acceptor_grid, size_bytes, cudaMemcpyDeviceToHost);

    //TODO: verificar como funciona. se for move semantics, checar se é memory safe

    // print_values_3D(out_elec_grid, grid->npointsx, grid->npointsy, grid->npointsz);
    // print_values_3D(out_vdwa_grid, grid->npointsx, grid->npointsy, grid->npointsz);
    // print_values_3D(out_vdwb_grid, grid->npointsx, grid->npointsy, grid->npointsz);
    // print_values_3D(out_solv_gauss, grid->npointsx, grid->npointsy, grid->npointsz);
    // print_values_3D(out_rec_solv_gauss, grid->npointsx, grid->npointsy, grid->npointsz);
    // print_values_3D(out_hb_donor_grid, grid->npointsx, grid->npointsy, grid->npointsz);
    // print_values_3D(out_hb_acceptor_grid, grid->npointsx, grid->npointsy, grid->npointsz);

    print_difference(grid->elec_grid, out_elec_grid);
    print_difference(grid->vdwA_grid, out_vdwa_grid);
    print_difference(grid->vdwB_grid, out_vdwb_grid);
    print_difference(grid->solv_gauss, out_solv_gauss);
    print_difference(grid->rec_solv_gauss, out_rec_solv_gauss);
    print_difference(grid->hb_donor_grid, out_hb_donor_grid);
    print_difference(grid->hb_acceptor_grid, out_hb_acceptor_grid);

    free(xyz_arr);
    free(HBdonors_arr);

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
    cudaFree(d_out_rec_si);

    printf("Invoking finished\n");

}
