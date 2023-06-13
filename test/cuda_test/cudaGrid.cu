#include "math.h"
#include <stdio.h>
#include <vector>
#include "../../src/pyGrid.h"
#include "../../src/pyMol2.h"

#define HB_C12 55332.873 
#define HB_C10 18393.199

// Will interpret the dielectric_model attribute as int because 
// string handling through the kernel may not be straight forward
typedef enum class e_DieletricModel {
    Constant,
    Four_r,
    None
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
                                int* HBacceptors,
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

    if (i < npointsx && j < npointsy && k < npointsz)
    {
        double x = i*grid_spacing + xbegin;
        double y = j*grid_spacing + ybegin;
        double z = k*grid_spacing + zbegin;

        double elec = 0.0;
        double vdwA = 0.0;
        double vdwB = 0.0;
        double solv = 0.0;
        double rec_solv = 0.0;
        
        for (int a = 0; a < N; a++) 
        {
            double d2 = distance_squared(x, xyz[a*xyz_w + 0], y, xyz[a*xyz_w + 1], z, xyz[a*xyz_w + 2]);
            double d6 = d2*d2*d2;
            double denom = 0.0;

            if (dielectric_model == Constant)
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
            double HB0 = HBdonors[m * HBdonors_w + 0];
            double HB1 = HBdonors[m * HBdonors_w + 1];

            double d2 = distance_squared(xyz[HB1*xyz_w + 0], x, xyz[HB1*xyz_w + 0], y, xyz[HB1*xyz_w + 0], z);
            double d10 = d2*d2*d2*d2*d2;
            double ang = angle(xyz[HB0*xyz_w + 0], xyz[HB0*xyz_w + 1], xyz[HB0*xyz_w + 2],
                      xyz[HB1*xyz_w + 0], xyz[HB1*xyz_w + 1], xyz[HB1*xyz_w + 2], x, y, z);
            double angle_term = (pow(cos(ang * M_PI / 180.0), 4.0));
            hb_donor += (HB_C12/(d10*d2)) - (HB_C10/d10);
        }

        double hb_acceptor = 0;
        {
           double d2 = distance_squared(xyz[HBacceptors[n]*xyz_w + 0], x, xyz[HBacceptors[n]*xyz_w + 1], y, xyz[HBacceptors[n]*xyz_w + 2], z);
           double d10 = pow(d2, 5.0);
           hb_acceptor += (HB_C12/(d10*d2)) - (HB_C10/d10);
        }

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

void invoke_compute_grid_softcore_HB_omp(const Grid& grid, const Mol2& rec) {

   double* d_xyz, d_charges, d_radii, d_epsilons_sqrt; 
   int* d_HBacceptors, d_HBdonors;

   double* out_elec_grid, out_vdwA_grid, out_vdwB_grid, out_solv_gauss, out_rec_solv_gauss, out_hb_donor_grid, out_hb_acceptor_grid;
   int* out_rec_si;

   cudaMalloc(&d_xyz, rec.xyz.size() * rec.xyz[0].size() * sizeof(double));
   cudaMalloc(&d_charges, rec.charges.size() * sizeof(double));
   cudaMalloc(&d_radii, rec.radii.size() * sizeof(double));
   cudaMalloc(&d_epsilons_sqrt, rec.epsilons_sqrt.size() * sizeof(double));

   cudaMalloc(&d_HBdonors, rec.HBdonors.size() * rec.HBdonors[0]*size() * sizeof(int));
   cudaMalloc(&d_HBacceptors, rec.HBacceptors.size() * sizeof(int));

   size_t size {(grid.npointsx * grid.npointsy * grid.npointsz) * sizeof(double)};
   cudaMalloc(&out_rec_si, size);
   cudaMalloc(&out_vdwA_grid, size);
   cudaMalloc(&out_vdwB_grid, size);
   cudaMalloc(&out_solv_gauss, size);
   cudaMalloc(&out_rec_solv_gauss, size);
   cudaMalloc(&out_hb_donor_grid, size);
   cudaMalloc(&out_hb_acceptor_grid, size);

   cudaMalloc(&out_rec_si, 1 * sizeof(int));

   //FIXME: possível memory error
   cudaMemcpy(rec.xyz.data(), d_xyz, rec.xyz.size() * rec.xyz[0].size() * sizeof(double), cudaMemcpyHostToDevice);
   cudaMemcpy(rec.charges.data(), d_charges, rec.charges.size() * sizeof(double), cudaMemcpyHostToDevice);
   cudaMemcpy(rec.radii.data(), d_radii, rec.radii.size() * sizeof(double), cudaMemcpyHostToDevice);
   cudaMemcpy(rec.epsilons_sqrt.data(), d_epsilons_sqrt, rec.epsilons_sqrt.size() * sizeof(double, cudaMemcpyHostToDevice));

   cudaMemcpy(rec.HBacceptors.data(), d_HBacceptors, rec.HBacceptors.size() * sizeof(int), cudaMemcpyHostToDevice);
   cudaMemcpy(rec.HBdonors.data(), d_HBdonors, rec.HBdonors.size() * rec.HBdonors[0].size() * sizeof(int), cudaMemcpyHostToDevice);

   DieletricModel dieletric_model{};

   switch (grid.Input->dielectric_model)
   {
   case "constant":
      dieletric_model = DieletricModel::Constant;
      break;
   
   default:
      dieletric_model = DieletricModel::None;
   }

   //TODO: ver melhor estes parâmetros de launch
   compute_grid_softcore_HB_omp<<<(size+255)/256, 256>>>(grid.npointsx, grid.npointsy, grid.npointsz,
                                                         grid.grid_spacing,
                                                         grid.xbegin, grid.ybegin, grid.zbegin,
                                                         dieletric_model,
                                                         grid.Input->deltaij_es6, grid.Input->deltaij6,
                                                         grid.Input->solvation_alpha, grid.Input->solvation_beta,
                                                         grid.Input->sigma, grid.Input->diel,
                                                         rec.N,
                                                         rec.xyz[0].size(), d_xyz,//FIXME: pode ser o shape errado 
                                                         d_charges,
                                                         d_radii,
                                                         d_epsilons_sqrt,
                                                         d_HBacceptors,
                                                         rec.HBdonors[0].size(), d_HBdonors,
                                                         out_elec_grid,
                                                         out_vdwA_grid,
                                                         out_vdwB_grid,
                                                         out_solv_gauss,
                                                         out_rec_solv_gauss,
                                                         out_hb_donor_grid,
                                                         out_hb_acceptor_grid,
                                                         out_rec_si)   

   //TODO: verificar como funciona. Se for move semantics, checar se é memory safe
   


   cudaFree();
   cudaFree();
   cudaFree();
   cudaFree();
   cudaFree();
   cudaFree();
   cudaFree();
   cudaFree();
   cudaFree();
   cudaFree();
   cudaFree();
   cudaFree();
   
}

class Bar{
public:
    std::vector<int> val(3, 10);
    Bar() = default;
};

int main() {

    Bar test{};
    Mol2 mol{};
    printf("%d, %d\n", test.val[1], mol.Nbonds);

    return 0;
}
