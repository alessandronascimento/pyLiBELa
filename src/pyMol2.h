/*
 * Mol2.h
 *
 *  Created on: 10/10/2011
 *      Author: Nascimento
 */

#ifndef MOL2_H_
#define MOL2_H_

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<cmath>
#include"pyPARSER.h"
#include <zlib.h>

using namespace std;

class Mol2 {
public:
    struct atom_param{
        string type;
        double radius;
        double epsilon;
        double mass;
    };


// variaveis

	//! Number of atoms
	int N ;
	//! Number of residues
	int Nres;
	//! Number of atomtypes
	int Natomtypes;
	//! Number of bonds;
	int Nbonds;
	//!
	string molname;
	//! Charges of the atoms. Must be divided by 18.2223 to get electron charges
	vector<double>charges;
	//! Atomic masses for the system.
	vector<double>masses;
	//! Atom names for the system according to AMBER FF.
	vector<string>amberatoms;
	//! Atomic parameters
	vector<string>atomtypes_prm;
	//! Atomic names
	vector<string> atomnames;
	//! Atomic coordinates
	vector<vector<double> > xyz;
	//! Atomic coordinates after rotation/translation
	vector<vector<double> > new_xyz;
	//! Atomic coordinates from the last step accepted.
	vector<vector<double> > old_xyz;
	//! Atomic coordinates from the optimized pose.
	vector<vector<double> > opt_overlay_xyz;
	//! Atomic coordinates from the optimized pose.
	vector<vector<double> > opt_energy_xyz;
	//! Trajectory vector
	vector<vector<vector<double> > > mcoords;
	//! Atomic parameter for LJ computation
	vector<double>epsilons;
	//! Squared root of atoms epsilons
	vector<double> epsilons_sqrt;
	//! Atomic radii
	vector<double>radii;
	//! C++ vector with the names of the residues
	vector<string> resnames;
	//! Pointers to the number of residues / number of atoms.
	vector<int> residue_pointer;
    //! C++ vector containing the index for atoms that are HB acceptors
    vector<int> HBacceptors;
    //! C++ vector containing the index for atoms that are HB acceptors. For a X-H donor, mapping is X in element [0] and H in element [1];
    vector<vector<int> > HBdonors;
	//! Temporary string to parse prmtop file.
	string line;
	//!
    char* str ;
	//! Keeps Vaa for RefMol/CompMol
	double self_obj_function;
	//!
	vector<vector<string> >bonds;
	//!
	vector<string> sybyl_atoms;

	vector<vector<vector<double> > > new_mcoords;
	string atm;
    //! Energies evaluated for the conformers generated. Uses GAFF.
    vector<double> conformer_energies;
    //! energy used when parsing a trajectory ensemble
    double energy;
    //! RMSD used when parsing a trajectory ensemble
    double rmsd;
    //! Longest distance among any pair of atoms in the molecule;
    vector<int> longest_axis;
    //! Molecule radius;
    double radius;

	/*!
	 * Initializer. This class has, as arguments, a pointer to the class PARSER.
	 * The class uses some information given by the user there. The filename is also
	 * given as argument. This makes possible to use the same object to reference and
	 * comparing molecules.
	 */
	Mol2();
    Mol2(PARSER Input, string molfile);
    bool parse_gzipped_file(PARSER Input, string molfile);
    bool parse_mol2file(PARSER Input, string molfile);
    bool parse_gzipped_ensemble(PARSER Input, string molfile, int skipper);
    vector<vector<double> > get_next_xyz(gzFile mol2file);

	/*!
	 *
	 */
	~Mol2();

    void initialize_gaff();
    vector<double> ensemble_energies;
    vector<double> ensemble_rmsd;
    vector<atom_param> gaff_force_field;
    void get_gaff_atomic_parameters(string gaff_atom, atom_param* ap);
    string sybyl_2_gaff(string atom);
    string sybyl_2_amber(string atom);

    //!
    void find_longest_axis();
    double distance(vector<double> atom1, vector<double> atom2);
};

#endif /* MOL2_H_ */
