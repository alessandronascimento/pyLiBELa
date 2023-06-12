#include "math.h"
#include <stdio.h>
#include <vector>
#include "../../src/pyMol2.h"

#define HB_C12 55332.873 
#define HB_C10 18393.199

// Will interpret the dielectric_model attribute as int because 
// string handling through the kernel may not be straight forward
typedef enum e_DieletricModel {
    Constant,
    Four_r
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
                                double** xyz,
                                double* charges,
                                double* radii,
                                double* epsilons_sqrt,
                                int* HBacceptors,
                                int** HBdonors,
                                double*** out_elec_grid,
                                double*** out_vdwA_grid,
                                double*** out_vdwB_grid,
                                double*** out_solv_gauss,
                                double*** out_rec_solv_gauss,
                                double*** out_hb_donor_grid,
                                double*** out_hb_acceptor_grid,
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
            double d2 = distance_squared(x, xyz[a][0], y, xyz[a][1], z, xyz[a][2]);
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
        double hb_acceptor = 0;

        for (int m = 0; m < sizeof(HBdonors)/sizeof(HBdonors[0]); m++)
        {
            double d2 = distance_squared(xyz[HBdonors[m][1]][0], x, xyz[HBdonors[m][1]][0], y, xyz[HBdonors[m][1]][0], z);
            double d10 = d2*d2*d2*d2*d2;
            double ang = angle(xyz[HBdonors[m][0]][0], xyz[HBdonors[m][0]][1], xyz[HBdonors[m][0]][2],
                      xyz[HBdonors[m][1]][0], xyz[HBdonors[m][1]][1], xyz[HBdonors[m][1]][2], x, y, z);
            double angle_term = (pow(cos(ang * M_PI / 180.0), 4.0));
            hb_donor += (HB_C12/(d10*d2)) - (HB_C10/d10);
        }

        for (int n = 0 ; n < sizeof(HBacceptors)/sizeof(HBacceptors[0]) ; n++)
        {
           double d2 = distance_squared(xyz[HBacceptors[n]][0], x, xyz[HBacceptors[n]][1], y, xyz[HBacceptors[n]][2], z);
           double d10 = pow(d2, 5.0);
           hb_acceptor += (HB_C12/(d10*d2)) - (HB_C10/d10);
        }

        out_elec_grid[i][j][k] = elec;
        out_vdwA_grid[i][j][k] = vdwA;
        out_vdwB_grid[i][j][k] = vdwB;
        out_solv_gauss[i][j][k] = solv;
        out_rec_solv_gauss[i][j][k] = rec_solv;
        out_hb_donor_grid[i][j][k] = hb_donor;
        out_hb_acceptor_grid[i][j][k] = hb_acceptor;
    }

    out_rec_si[0] = 0.0;
    for (int a = 0; a < N ; a++)
    {
        out_rec_si[0] += (solvation_alpha * pow(charges[a], 2.0)) + solvation_beta; 
    }
}

//TODO: talvez tenha uma forma direta de fazer isso com C++ 
void invoke_compute_grid_softcore_HB_omp() {

}

class Bar{
public:
    std::vector<int> val{3, 10};
    Bar() = default;
};

int main() {

    Bar test{};
    Mol2 mol{};
    printf("%d, %d\n", test.val[1], mol.Nbonds);

    return 0;
}
