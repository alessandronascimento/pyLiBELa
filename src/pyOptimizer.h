/*
 * Optimizer.h
 *
 *  Created on: 21/03/2012
 *      Author: Nascimento
 */

#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_

#include"nlopt.h"
#include"nlopt.hpp"
#include"pyMol2.h"
//#include"ENERGY.h"
#include "iMcLiBELa.h"
#include "pyEnergy2.h"
#include"pyPARSER.h"
#include"pyCOORD_MC.h"
#include"pyGaussian.h"
#include "pyWRITER.h"
#include<gsl/gsl_rng.h>
#include<time.h>
#include <memory>
#include<vector>
#include "pySA.h"
#include "pyMC.h"

using namespace std;

class Optimizer {
public:

	static Mol2* Rec;
	static Mol2* RefLig;
	static PARSER* Parser;
	static Grid* Grids;

    struct opt_result_t{
		int optimization_status;
		double f_min;
        energy_result_t* energy_result;
		vector<vector<double> > optimized_xyz;
	};

    struct align_result_t{
        double rmsd;
        vector<double> translation;
        vector<double> rotation;
        int opt_status;
    };

    struct align_t{
        vector<vector<double> > ref_xyz;
        vector<vector<double> > current_xyz;
    };


	Optimizer(Mol2* _Rec, Mol2* _RefLig, PARSER* _Parser);
	Optimizer(Mol2* _Rec, Mol2* _RefLig, PARSER* _Parser, Grid* _Grids);
	virtual ~Optimizer();

	static double evaluate_rmsd(Mol2* Lig1, Mol2* Lig2);
	static void Simulated_Annealing(Mol2* Lig, opt_result_t* opt_result);

	static vector<vector<double> > update_coords(const std::vector<double> &x, Mol2* Lig2);
    static double distance(vector<double> atom1, vector<double> atom2);
    static double distance_squared(vector<double> atom1, vector<double> atom2);

    static double evaluate_energy(Mol2* Lig2, vector<vector<double> > new_xyz);
    static void evaluate_energy2(Mol2* Lig2, vector<vector<double> > new_xyz, energy_result_t* energy_result);

	static double objective_energy_function(const vector<double> &x, vector<double> &grad, void *data);
	static double objective_overlay_function(const vector<double> &x, vector<double> &grad, void *data);
    static double objective_prealign_function(const vector<double> &x, vector<double> &grad, void *data);
    static double superpose_function(const vector<double> &x, vector<double> &grad, void *data);

	static void minimize_energy_nlopt_cobyla(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_energy_nlopt_lbfgs(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_energy_nlopt_ln_auglag(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_energy_nlopt_ld_auglag(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_energy_nlopt_mma(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_energy_nlopt_subplex(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_energy_nlopt_simplex(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_energy_nlopt_crs(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_energy_nlopt_direct(Mol2* Lig2, opt_result_t* opt_result);
    static void minimize_energy_nlopt_direct_only(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_energy_nlopt_stogo(Mol2* Lig2, opt_result_t* opt_result);
    static void minimize_energy_nlopt_isres(Mol2* Lig2, opt_result_t* opt_result);
    static void minimize_energy_nlopt_esch(Mol2* Lig2, opt_result_t* opt_result);
    static void minimize_energy_adaptative(Mol2* Lig2, opt_result_t* opt_result);

	static void minimize_overlay_nlopt_cobyla(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_overlay_nlopt_subplex(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_overlay_nlopt_lbfgs(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_overlay_nlopt_ln_auglag(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_overlay_nlopt_ld_auglag(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_overlay_nlopt_mma(Mol2* Lig2, opt_result_t* opt_result);
	static void minimize_overlay_nlopt_crs(Mol2* Lig2, opt_result_t* opt_result);
    static void minimize_overlay_nlopt_direct(Mol2* Lig2, opt_result_t* opt_result);

    static void minimize_alignment_nlopt_simplex(align_t* align_data, align_result_t* opt_result, vector<double> current_com);

    static void pre_align(Mol2* Lig2, opt_result_t* opt_result);

/*	
    void set_rec(Mol2* rec) {
            Rec = rec;
        }

    void set_ref_lig(Mol2* ref_lig) {
            RefLig = ref_lig;
        }

    void set_parser(PARSER* parser) {
            Parser = parser;
        }

    void set_grids(Grid* grids) {
            Grids = grids;
        }
    void set_writer(WRITER* writer) {
            WRITER = writer;
        }


     void run() {
            // some code to run optimization
            cout << "Running optimization..." << endl;
        }
*/	
};
/*
Mol2* Optimizer::Rec = nullptr;
Mol2* Optimizer::RefLig = nullptr;
PARSER* Optimizer::Parser = nullptr;
Grid* Optimizer::Grids = nullptr;
WRITER* Optimizer::Writer = nullptr;
*/
#endif /* OPTIMIZER_H_ */
