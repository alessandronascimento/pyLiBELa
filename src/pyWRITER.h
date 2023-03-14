/*
 * WRITER.h
 *
 *  Created on: 23/11/2010
 *      Author: Nascimento
 */

#ifndef WRITER_H_
#define WRITER_H_

#include<string>
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<iostream>
#include<zlib.h>
#include"pyMol2.h"
#include "pyPARSER.h"
#include "iMcLiBELa.h"

using namespace std;


class WRITER {
public:
//! The output file is defined as a log file of the simulation.
	FILE *output;
//! String used to generate the output file name.
	string outputfile;
//! String used to parse the prefix for all output files.
	string output_prefix;
    //! Mol2 molecule used to get coordinates, atom names, etc.
	Mol2 *Cmol;
    //! Gzipped file used to write MC trajectory or docked molecules
    gzFile outmol2;
    //! Input parameters from class PARSER
    PARSER* Input;






/*!
 * The constructor is used only to parse the output_prefix as defined in
 * the input file and store it in the class.
 */

    WRITER(string prefix, PARSER* _Input);
    WRITER(string prefix);
    WRITER(PARSER* _Input);
    ~WRITER();

/*!
 *
 * @param center C++ vector with center of mass for box computation, i.e., comX, comY and comZ.
 * @param min_x
 * @param min_y
 * @param min_z
 * @param max_x
 * @param max_y
 * @param max_z
 */
	void write_box(vector<double>center, double min_x, double min_y, double min_z, double max_x, double max_y, double max_z);

/*!
 * The write pdb method is used to write gzipped file everytime a configuration
 * is accepted during a SA run or during an equilibrium run.
 */

	void write_pdb(vector<vector<double> >xyz, int N, vector<string> atomnames, double energy, double temp, string outname);

	/*!
	 *
	 */
	void write_pdb(Mol2 *Cmol, vector<vector<double> >xyz, double energy, double rmsd, string outname);
	void writeMol2(Mol2* Cmol, vector<vector<double> >xyz, double energy, double rmsd, string outname);
	void writeMol2(Mol2* Cmol, vector<vector<double> >xyz, double energy, double rmsd);
    void writeMol2(Mol2* Cmol, vector<vector<double> >xyz, energy_result_t* result, double rmsd);
    void writeMol2(Mol2* Cmol, vector<vector<double> >xyz, energy_result_t* result, double rmsd, string outname);
	void writeMol2_Mol_new_xyz(Mol2* Cmol, double energy, double rmsd);

    /**
     * @brief write_pqr Function to write a PQR file from a MOL2 file. Atomic radii and charges are taken ]
     * LIBELA conversion (AMBER/GAFF) and from mol2 files, respectively.
     * @param Cmol MOL2 object
     * @param outname prefix for file writting.
     */
    void write_pqr(Mol2 *Cmol, string outname);

/*!
 * The print_info method is used to, at the same time, print some useful(?)
 * information about the simulation in the screen and write it to the log
 * file.
 */

	void print_info(char info[98]);

	void print_info_no_newline(char info[98]);

	void print_line(void);

/*!
 *
 */
	void print_welcome(void);

    void print_params();

	bool FileExists(string strFilename);
    void print_dock_params();
};

#endif /* WRITER_H_ */
