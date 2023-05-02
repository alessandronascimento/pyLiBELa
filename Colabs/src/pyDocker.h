/*
 * Docker.h
 *
 *  Created on: 24/03/2012
 *      Author: Nascimento
 */

#ifndef DOCKER_H_
#define DOCKER_H_

#include "pyMol2.h"
#include "pyOptimizer.h"
#include "pyPARSER.h"
#include "pyWRITER.h"
#include "pyCOORD_MC.h"
#include <vector>
#include <string>
#include <string.h>
#include "pyGrid.h"

#ifdef HAS_GUI
#include "GUI/QtWriter.h"
#include "Deal.h"
#endif

class Docker{
public:
	//! Objective function value
	double f;
	//! Objective function for molecular overlay
	double fo;
	//! Info for printing
    char* info;//[98];
	//!
	double best_ene;
	//!
	int best_conf;
    //!
    WRITER* Writer;
#ifdef HAS_GUI
    QtWriter* QWriter;
    Docker(QtWriter* _QWriter);
#endif

    Docker(WRITER* _Writer);
    void print_line();
    void print_info(char info[98]);
    void write_mol2(Mol2* Lig, vector<vector<double> > new_xyz, double ene, double si);
    void write_mol2(Mol2* Lig, vector<vector<double> > new_xyz, energy_result_t* result, double si);
	//! Optimizer class
	//! Class constructor. Automatically starts the docking procedure.
    void run(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, unsigned counter);
    //! Class constructor. Automatically starts the docking procedure.
    void run(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, Grid* Grids, unsigned counter);

	//! Method for docking a set of conformers of a molecule. All the conformers are docked and a selected number are optimized in the active site.
    void Dock_conformers(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, unsigned counter);
    //! Method for docking a set of conformers of a molecule. All the conformers are docked and a selected number are optimized in the active site.
    void Dock_conformers(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, Grid* Grids, unsigned counter);

#ifdef HAS_GUI
	//! Constructor for docking. Automatically starts the docking procedure. Used in the Qt GUI.
    Docker(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, QtWriter* Writer, unsigned counter);
//    Docker(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, QtWriter* Writer);
    //! Constructor for docking. Automatically starts the docking procedure. Used in the Qt GUI.
    Docker(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, QtWriter* Writer, Grid* Grids, unsigned counter);
    //! Constructor class for ligand docking with automatic conformer generation
    void Dock_conformers(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, QtWriter* Writer, Grid* Grids, unsigned counter);
    //! Constructor class for ligand docking with automatic conformer generation
    void Dock_conformers(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, QtWriter* Writer, unsigned counter);
#endif

	virtual ~Docker();

	//! C++ vector with the docked (oriented) coordinates.
	vector<vector<double> > new_coord;

	/*! Method to sort a set of conformers and rank the N best conformers using the binding energy criteria.
	 * N is defined in the PARSER class as PARSER::conformers_to_evaluate.
	 */
	vector<unsigned> sort_vector(vector<double> vec);
	vector<unsigned> sort_vector_inv(vector<double> vec);

    /*!
     * \brief minimize_overlay This method calls the methods in the Optimizer class to optimize ligand
     * overlay with respect to the reference ligand.
     * \param Input PARSER object with execution parameters
     * \param Opt Optimizer object
     * \param Lig Mol2 object with docking ligand.
     * \param opt_result Type defined in Optimizer class to get all the optimization results.
     * \return false if optimization has any problem.
     */
    bool minimize_overlay(PARSER* Input, Optimizer* Opt, Mol2* Lig, Optimizer::opt_result_t* opt_result);

    /*!
     * \brief minimize_energy This method calls the methods in the Optimizer class to optimize ligand
     * binding within receptor active site.
     * \param Input PARSER object with execution parameters
     * \param Opt Optimizer object
     * \param Rec Mol2 object with biological receptor.
     * \param Lig Mol2 object with docking ligand.
     * \param opt_result Type defined in Optimizer class to get all the optimization results.
     * \return false if optimization has any problem
     */
    bool minimize_energy(PARSER* Input, Optimizer* Opt, Mol2* Rec, Mol2* Lig, Optimizer::opt_result_t* opt_result);
    bool minimize_energy(PARSER* Input, Optimizer* Opt, Mol2* Rec, Mol2* Lig, Mol2* Reflig, Optimizer::opt_result_t* opt_result);

};

#endif /* DOCKER_H_ */
