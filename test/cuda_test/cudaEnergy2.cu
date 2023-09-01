#include "../../src/pyEnergy2.h"

__global__ void compute_energy_softcore_solvation(
    int N, double solvation_alpha, double solvation_beta, double sigma,
    double deltaij_es6, double deltaij6, double diel, double out,
    double *charges, double *xyz, double *lig_xyz, double *radii,
    double *epsilons_sqrt) {

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if (i < N) {

    double lig_solv = 0.00;
    double vdw = 0.0;
    double elec = 0.0;
    double acoef, bcoef, lig_solv_affinity, lig_solv_distf, rec_solv_distf,
        rec_solv_affinity, dij, dij2, deff, dij6; //, rij, eij;
    double sqrt2 = sqrt(2);

    rec_solv_affinity =
        (solvation_alpha * ((charges[i]) * (charges[i]))) + solvation_beta;
    for (int j = 0; j < N; j++) {

      //! distances

      dij2 = distance_squared(xyz[i][0], lig_xyz[j][0], xyz[i][1],
                              lig_xyz[j][1], xyz[i][2], lig_xyz[j][2]);
      dij6 = dij2 * dij2 * dij2;
      dij = sqrt(dij2);
      deff = pow(((dij * dij * dij) + deltaij_es3), (1.0 / 3.0));

      //! electrostactic energy (softcore)
      if (dielectric_model == "constant") {
        elec += (332.0 * charges[i] * charges[j]) / (diel * deff);

        lig_solv_affinity =
            (solvation_alpha * ((charges[j]) * (charges[j]))) + solvation_beta;
        lig_solv_distf = ((4.0 / 3.0) * PI * (radii[i] * radii[i] * radii[i]));
        lig_solv_distf =
            (lig_solv_distf * exp(-(deff * deff) / (2 * (sigma * sigma)))) /
            (sigma * sigma * sigma);
        rec_solv_distf = ((4.0 / 3.0) * PI * (radii[j] * radii[j] * radii[j]));
        rec_solv_distf =
            (rec_solv_distf * exp(-(deff * deff) / (2 * (sigma * sigma)))) /
            (sigma * sigma * sigma);
      } else {
        deff = pow(((dij6) + deltaij_es6), (1.0 / 3.0));
        elec += (332.0 * charges[i] * charges[j]) /
                (diel * (pow((dij6 + deltaij_es6), (1.0 / 3.0))));

        lig_solv_affinity =
            (solvation_alpha * ((charges[j]) * (charges[j]))) + solvation_beta;
        lig_solv_distf = ((4.0 / 3.0) * PI * (radii[i] * radii[i] * radii[i]));
        lig_solv_distf = (lig_solv_distf * exp(-deff / (2 * (sigma * sigma)))) /
                         (sigma * sigma * sigma);
        rec_solv_distf = ((4.0 / 3.0) * PI * (radii[j] * radii[j] * radii[j]));
        rec_solv_distf = (rec_solv_distf * exp(-deff / (2 * (sigma * sigma)))) /
                         (sigma * sigma * sigma);
      }

      //! VDW energy (softcore)

      acoef = (epsilons_sqrt[i] * pow(2 * radii[i], 6)) *
              (epsilons_sqrt[j] * pow(2 * radii[j], 6));
      bcoef = (sqrt2 * epsilons_sqrt[i] * pow(2 * radii[i], 3)) *
              (sqrt2 * epsilons_sqrt[j] * pow(2 * radii[j], 3));
      vdw +=
          ((acoef / pow((dij6 + deltaij6), 2)) - (bcoef / (dij6 + deltaij6)));

      //! Solvation energy

      rec_solv += rec_solv_affinity * rec_solv_distf;
      lig_solv += lig_solv_affinity * lig_solv_distf;
    }
  }
}
