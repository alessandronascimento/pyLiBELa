/*
 * PARSER.h
 *
 *  Created on: 28/05/2010
 *      Author: Nascimento
 */

#ifndef PARSER_H_
#define PARSER_H_

#include<stdio.h>
#include<stdlib.h>
#include<fstream>
#include<iostream>
#include<string>
#include"string.h"
#include<math.h>
#include<vector>
#ifdef HAS_GUI
#include <QStringList>
#include <QThread>
#endif

using namespace std;

/**
 * @brief The PARSER class is used to get the input parameters from the input file. These
 * parameters are used in almost every class of the program.
 */


class PARSER {

public:

//!	Used to parse the input file.
	string param;
//!	Number of steps of the equilibrium simulation.
	int number_steps;
//! Temperature of the equilibrium simulation. Used as final temperature of the simulated annealing also.
	double temp;
//! Number of steps for each sweep
    int sweep_steps;
//! Number of Equilibration Steps
    int eq_steps;
//! Defines the maximum amplitude of the shift during the simulation.
	double cushion;
//! Soft-core value for the dispersion-repulsion potential.
	double deltaij6;
//! Soft-core value for the electrostatic potential**6.
	double deltaij_es6; // = 1.75*1.75*1.75*1.75*1.75*1.75;
//! Soft-core value for the electrostatic potential ** 3.
    double deltaij_es3;
//! Dielectric constant
	double diel;
//! Weighting scheme in rec/lig_solv_distf. From Verkhivker (1999), it has a value of 3.5 Angstrons.
	double sigma;
//! Temporary information.
	string rprmtop, rcrd, lprmtop, lcrd;
//!	PRMTOP file with receptor information (AMBER format).
	const char* rec_prmtop;
//! Receptor coordinates in AMBER crd format.
	const char* rec_crd;
//! PRMTOP file with ligand information (AMBER format).
	const char* lig_prmtop;
//! Ligand coordinates in AMBER crd format.
	const char* lig_crd;
//! Used to parse and yes or no answer to the next variable.
	string sa;
//! Will SA optimization procedure be used?
	bool sa_scheme;
//! Start temperature for a Simulated Annealing simulation.
	double sa_start_temp;
//! Number of simulation steps between each cycle of temperature decrease.
	int sa_steps;
//! Defines the prefix for output files.
	string output;
//! Defines the rotation step during the MC sorting.
	double rotation_step;
//! Defines if torsion angles should be sampled in MC.
    bool sample_torsions;
//! Defines the step for torsion sampling in MC sorting.
    double torsion_step;
	//! Maximum dimension in X direction to the sampling box
	double x_dim;
	//! Maximum dimension in Y direction to the sampling box
	double y_dim;
	//! Maximum dimension in Z direction to the sampling box
	double z_dim;
	/*! Filename with ligand "trajectory" to use it in a semi-flexible approach.
	 * The file format is the same as in AMBER trajectory.
	 */
	string lig_traj;
	//! Filename with receptor "trajectory" to use it in a semi-flexible approach.
	//! The file format is the same as in AMBER trajectory.
	string rec_traj;
	/*! Flag to define if the ligand will be defined as flexible. The flag is a
	* automatically turned on when the input command lig_traj is used.
	*/
	bool flex_lig;
	/*! Flag to define if the receptor will be defined as flexible. The flag is a
	 * automatically turned on when the input command rec_traj is used.
	 */
	bool flex_rec;
	/*! Defines a timeout function (in seconds) in Equilibrium simulations. Above this
	 * limit, the program will abort the execution.
	 */
	int timeout;
    //! String to parse temporary information
	string tmp;
	/*!
	 * Defines the scoring function to be used throughout the code.
	 *  0 = Softcore + Solvation
	 *  1 = Softcore
	 *  2 = Amber FF + Solvation
	 *  3 = Amber FF
	 */
	int scoring_function;
	//! Defines wether or not the equilibrium MC will use a dynamic temperature scaling dependent on the acceptance rate.
	bool dynamic_rate;
	//! Lower limit for the dynamic rate
	double dyn_rate_inf_limit;
	//! Higher limit for the dynamic rate
	double dyn_rate_sup_limit;
	//! Number of steps in the dynamic rate scheme.
	int dyn_rate_steps;
	//! C++ String used to parse the receptor MOL2 filename.
	string rec_mol2;
	//! C++ string used to parse the ligand MOL2 filename.
	string lig_mol2;
	//! Defines if mol2 file (for both ligand and receptor) uses Amber (true) or Sybyl (false) atom types.
	bool mol2_aa;
    //! Optimization initial move
	double min_delta;
    //! Optimization tolerance used for overlay mathcing
	double min_tol;
    //! Minimization timeout (in seconds)
	int min_timeout;
    //! Mode to use (eq, sa or dock);
	string mode;
    //! Electrostatic scale in the overlay matching.
	double elec_scale;
    //! VDW scale in the overlay mathcing
	double vdw_scale;
    //! Mol2 file with the reference ligand (already posed in the binding pocket).
	string reflig_mol2;
    //! File with the list of molecules to be docked.
	string multifile;
    //! flag to use smiles
    bool use_smiles;
    //! File with the list of smiles molecules to be docked.
    string smiles_multifile;
    /*! Algorithm for overlay optimization. Current options are:
     * lbfgs
     * ln_auglag
     * cobyla
     * mma
     * subplex
     * ld_auglag
     * crs
     * direct
     */
	string overlay_optimizer;
    /*!
     * \brief energy_optimizer
     * Algorithm for binding energy optimization. Current options are:
     * lbfgs
     * ln_auglag
     * cobyla
     * mma
     * subplex
     * ld_auglag
     * crs
     * direct
     * stogo
     * SA
     * none
     */
    string energy_optimizer;
    //! Ignore hydrogen atoms during overlay matching
	bool dock_no_h;
	//! Size of population for density estimation algorithm
	int deal_pop_size;
	//! Number of generations in DEAl
	int deal_generations;
	//! Number of top-ranked poses used to estimate density
	int deal_top_ranked;
    //! Use binding energy optimization through Density Estimation Algorithm
	bool deal;
    //! Generate conformers using Openbabel API
	bool generate_conformers;
    //! Number of conformers to generate
	int lig_conformers;
    //! Number of conformers that should be optimized
	int conformers_to_evaluate;
    //! Number of minimization steps before conformer generation, when using WRS.
    int conformer_min_steps;
    //! Number of search trials for conformer generation in confab
    int conf_search_trials;
    //! Tolerance for binding energy optimization
	double dock_min_tol;
    //! Turns on the equilibrium mode
	bool eq_mode;
    //! Turns on the simulated annealing mode
	bool sa_mode;
    //! Turns on the docking mode.
	bool dock_mode;
    //! True if parallel execution is active (through OpenMP)
	bool dock_parallel;
    //! Number of parallel threads to use
	unsigned parallel_jobs;
    //! Write docked molecules as mol2 files ?
	bool write_mol2;
    //! SA mu
	double sa_mu_t;
    //! Grid spacing when using grid computation
	double grid_spacing;
    //! Use energy grids
    bool use_grids;
    //! Maximal box to search during optimization
	double search_box_x, search_box_y, search_box_z;
    //! Load pre-computed grids from file
    bool load_grid_from_file;
    //! Grid file prefix
	string grid_prefix;
    //! Write computed grids to file
	bool write_grids;
    //! Show RMSD. Used for debuggin purposes only.
	bool show_rmsd;
    //! Sort conformers using initial binding energy. If set to false, conformers
    //! will be ranked by Overlay objective function.
	bool sort_by_energy;
    //! Only scores ligands, without optimization.
    bool only_score;
    //! Alpha parameter in solvent affinity.
    double solvation_alpha;
    //! Beta parameter in solvent affinity.
    double solvation_beta;
    //! seed for random number generation
    int seed;
    //! ligand simulation
    bool ligsim;
    //! An "stride" parameter for writing trajectory file. Coordinates are written at mc_stride intervals.
    int mc_stride;
    //! step for translation in full search
    double translation_step;
    //!
    bool full_search_mode;
    //!
    double scale_vdw_energy;
    //!
    double scale_elec_energy;
    //!
    bool use_score_optimization;
    //! If true, will not take into account the ligand internal energy during MC simulations.
    bool use_only_binding_energy;


    /*!
     * \brief dielectric_model This keyword brings the model for dieletric constant usage
     * in electrostatic evaluation within the code. There are three available models:
     * \"constant\" : will use a single constant value for dieletric. This value can be
     * defined in PARSER::diel;
     * \"r\" : will use a dielectric function where D is equal to r, i.e., the distance between
     * the two atoms.
     * \"4r\" : uses a 4 times r function for D. This is the model used in UCSF DOC6.
     */
    string dielectric_model;

    //! Size for Monte Carlo Recursion (MCR) series
    int mcr_size;

    //! Vector with MCR coefficients (bi values from Li and Scheraga's paper);
    vector<double> mcr_coefficients;

    //! Activates MCR mode
    bool mcr_mode;

    //! b coefficient for MCR in the ith iteration;
    double bi;
    //! Defines whether extensive output should be given
    bool verbose;
    //! Number of bins for histograming the MC data for ligand rotation
    int entropy_rotation_bins;
    //! Number of bins for histograming the MC data for ligand translation
    int entropy_translation_bins;
    //! There are two currently implemented functions for evaluating the ligand internal energy.
    //! The possible choices are:
    //! \"gaff\" or
    //! \"MMFF94\"
    string ligand_energy_model;
    //! Defines whether atoms should be defined as GAFF types or AMBER types
    string atomic_model_ff;
    //! Defines an electrostatic potential grid as computed by PBSA (AMBERTools)
    string pbsa_grid;
    //! Defines an electrostatic potential grid as computed by DelPhi
    string delphi_grid;
    //! Defines whether PBSA potential from AmberTools is to be used
    bool use_pbsa;
    //! Defines whether potential from Delphi (binary file) is to be used
    bool use_delphi;
    //! Number of points in delphi grid in each direction. Must be an odd number.
    int delphi_gsize;
    //! Defines an electrostatic potential in DelPhi CUBE format
    string delphi_cube_grid;
    //! Activates the full flexible approach in MC simulations
    bool mc_full_flex;
    //!
    bool compute_rotation_entropy;
    //! Maximal atomic displacement in MC full flex
    double max_atom_displacement;
    //!
    bool use_writeMol2_score_cutoff;
    //!
    bool use_writeMol2_energy_cutoff;
    //!
    double writeMol2_score_cutoff;
    //!
    double writeMol2_energy_cutoff;
    //! Sigma parameter for Gaussian-Weight in Coulomb potential scaling
    double coulomb_sigma;
    //! Sigma parameter for Gaussian-Weight in Lennard-Jones scaling
    double LJ_sigma;
    //! Use Gaussian-weighting for attractive LJ potential
    bool use_GW_LJ6;
    //! Use Gaussian-weighting for repulsive LJ potential
    bool use_GW_LJ12;
    //! Use Gaussian-weighting for Coulomb potential
    bool use_GW_Coulomb;
    //!
    bool use_Erestraints;
    //!
    double restraints_weight;
    //! Defines whether there should be used a cutoff in overlay to energy optimization
    bool use_overlay_cutoff;
    //! Cutoff value in similatiy index for energy optimization
    double overlay_cutoff;
    //! Maximal number of atoms allowed for smiles molecules
    int atom_limit;

#ifdef HAS_GUI
	QStringList docking_molecules;
#endif


	// FUNCTIONS



	PARSER();
	/*!
	 * This function sets some default parameters and then calls the next function
	 * while parses the input files to grab the parameters of the simulation;
	 */
	void set_parameters(char* arg);

	/*!
	 * This function actually processes the parameter parsed in the input file.
	 */
	void comparing (string param, ifstream &input);

	/*!
     * This function check the consistency of some of the given parameters.
	 */
	void check_parameters_sanity(void);
};

#endif /* PARSER_H_ */
