#ifndef ENERGY2_H
#define ENERGY2_H

#include<vector>
#include<cmath>
#include "pyPARSER.h"
#include "pyMol2.h"
#include "iMcLiBELa.h"
#include "pyGrid.h"

/**
 * @brief The Energy2 class is used for energy calculations during docking and MC simulations.
 *
 */

class Energy2 {
public:
    struct GridInterpol{
        double vdwA;
        double vdwB;
        double elec;
        double pbsa;
        double delphi;
        double solv_gauss;
        double rec_solv_gauss;
        double hb_donor;
        double hb_acceptor;
    };

    //! class PARSER provided to this class;
    PARSER* Input;

    Energy2(PARSER* _Input);
    //!
    //! \brief distance Computes the distance between two atoms
    //! \param x1 X coordinate for atom 1
    //! \param x2 X coordinate for atom 2
    //! \param y1 Y coordinate for atom 1
    //! \param y2 Y coordinate for atom 2
    //! \param z1 Z coordinate for atom 1
    //! \param z2 Z coordinate for atom 2
    //! \return the distance between the two atoms
    //!

    bool atom_is_acceptor(int a, Mol2* Lig);
    int H_is_donor(int a, Mol2* Lig);
    double angle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);
    double distance(double x1, double x2, double y1, double y2, double z1, double z2);

    //!
    //! \brief distance_squared Computes the square of the distance between two atoms
    //! \param x1 X coordinate for atom 1
    //! \param x2 X coordinate for atom 2
    //! \param y1 Y coordinate for atom 1
    //! \param y2 Y coordinate for atom 2
    //! \param z1 Z coordinate for atom 1
    //! \param z2 Z coordinate for atom 2
    //! \return the square o the distance between two atoms
    //!
    double distance_squared(double x1, double x2, double y1, double y2, double z1, double z2);

    //!
    //! \brief compute_energy_softcore_solvation Computes the interaction energy between ligand and
    //! receptor using a softcore model with solvation (SV model)
    //! \param Rec Object from class MOL2 with receptor description
    //! \param Lig Object from class MOL2 with ligand description
    //! \param lig_xyz C++ vector of vector with ligand coordinates
    //! \return Rhe interaction energy in kcal/mol.
    //!
    double compute_energy_softcore_solvation(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz);
    double compute_energy_softcore(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz);
    double compute_energy_hardcore_solvation(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz);
    double compute_energy_hardcore(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz);

    double compute_ene_from_grids_softcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);
    double compute_ene_from_grids_softcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);
    double compute_ene_from_grids_hardcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);
    double compute_ene_from_grids_hardcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);


    double compute_energy_softcore_solvation(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result);
    double compute_energy_softcore(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result);
    double compute_energy_hardcore_solvation(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result);
    double compute_energy_hardcore(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result);

    double compute_ene(Mol2* Rec,   Mol2* Lig, vector<vector<double> >lig_xyz);
    double compute_ene(Grid *Grids, Mol2* Lig, vector<vector<double> >lig_xyz);
    double compute_ene(Grid* Grids, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result);
    double compute_ene(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result);

    double compute_ene_from_grids_hardcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result);
    double compute_ene_from_grids_hardcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result);
    double compute_ene_from_grids_softcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result);
    double compute_ene_from_grids_softcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result);

    double trilinear_interpolation(vector<vector<vector<double> > > grid, double x, double y, double z, int x0, int y0, int z0, int x1, int y1, int z1);
    void trilinear_interpolation(Grid* Grids, double x, double y, double z, int x0, int y0, int z0, int x1, int y1, int z1, GridInterpol* GI);

    double evaluate_forces_hardcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);
};

#endif // ENERGY2_H
