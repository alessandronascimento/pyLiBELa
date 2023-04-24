/*
 * Grid.h
 *
 *  Created on: Jan 7, 2013
 *      Author: asn
 */

#ifndef GRID_H_
#define GRID_H_

#include <vector>
#include <cmath>
#include <sstream>
#include "pyMol2.h"
#include "pyPARSER.h"
#include "pyWRITER.h"
#ifdef HAS_GUI
#include "GUI/QtWriter.h"
#endif
#include <zlib.h>

class Grid {
public:
    //! Spacing between grid points
	double grid_spacing;
    //! Number of points in each grid direction.
	int npointsx, npointsy, npointsz;
    //! Coordinates for the edges of the grid cube
	double xbegin, ybegin, zbegin, xend, yend, zend;
    //! Coordinates
	vector<vector<vector<double> > > elec_grid;
    vector<vector<vector<double> > > pbsa_grid;
	vector<vector<vector<double> > > vdwA_grid;
	vector<vector<vector<double> > > vdwB_grid;
	vector<vector<vector<double> > > solv_gauss;
	vector<vector<vector<double> > > rec_solv_gauss;
    vector<vector<vector<double> > > delphi_grid;
    vector<vector<vector<double> > > hb_donor_grid;
    vector<vector<vector<double> > > hb_acceptor_grid;

	double rec_si;
    //! Pointer to the PARSER class
	PARSER* Input;
    //! Pointer to the WRITER class
    WRITER* Writer;
#ifdef HAS_GUI
    QtWriter* QWriter;
#endif
    //! Defines whether or not PBSA electrostatic grid was loaded
    bool pbsa_loaded; // = false;
    //! Defines whether or not DELPHI electrostatic grid was loaded
    bool delphi_loaded; // = false;
    //! Char array to use with Writer class
    char* info;//[98];

    /*!
     * \brief Grid This function initializes the class. It passes a copy of the PARSER
     * class to the Grid class and defines the grid spacing
     * \param _Input Pointer to the PARSER class.
     * \param _Writer Pointer to the WRITER class.
     */
    Grid(PARSER* _Input, WRITER* _Writer);
    /*!
     * \brief Grid Initializer of the Grid class used through the code. This method copies
     * a pointer to the class PARSER to this class, defines the grid spacing and calls the
     * appropriate methods to start grid computation. Also, if the keyword "write_grids" is
     * set to "yes" in the PARSER, it calls the method to write the grids to a file.
     * \param _Input Pointer to the PARSER class.
     * \param _Writer Pointer to the WRITER class.
     * \param Rec POINTER to a Mol2 class with Receptor information.
     * \param com C++ vector with the coordinates of the center of mass of the reference
     * ligand. It is used to define the center of the computation box.
     */

#ifdef HAS_GUI

    Grid(PARSER* _Input, QtWriter* _QWriter, Mol2* Rec, vector<double> com);

#endif

    Grid(PARSER* _Input, WRITER* _Writer, Mol2* Rec, vector<double> com);

	double distance(double x1, double x2, double y1, double y2, double z1, double z2) ;

    double distance_squared(double x1, double x2, double y1, double y2, double z1, double z2) ;

    double angle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);

    void generate_points(vector<double> ref_com);

    void compute_grid_softcore(Mol2* Lig);
    void compute_grid_softcore_omp(Mol2* Lig);
    void compute_grid_softcore_HB_omp(Mol2* Lig);

	void compute_grid_hardcore(Mol2* Lig);
    void compute_grid_hardcore_Gaussian(Mol2* Lig);
    /**
     * @brief write_grids_to_file This function writes a computed potential grid to a file
     */
	void write_grids_to_file(void);
    /**
     * @brief load_grids_from_file This function loads a pre-computed Grid file, computed with
     * LiBELa or McGrid
     */
	void load_grids_from_file(void);
    /**
     * @brief load_Ambergrids_from_file This function loads a PBSA potential file as computed
     * by pbsa program in AmberTools
     */
    void load_Ambergrids_from_file(void);

    /**
     * @brief load_gzAmbergrids_from_file This function loads a gzipped PBSA potential file as computed
     * by pbsa program in AmberTools.
     */
    void load_gzAmbergrids_from_file(void);

    /**
     * @brief load_Delphi_Grid_from_file This fuction loads a Delphi potential file in
     * binary BIOSYM format file
     */
    void load_Delphi_Grid_from_file(void);
    /**
     * @brief load_phimap_from_file This function loads a Delphi potential file in
     * binary (unformatted) Delphi file
     * @param gsize number of points in the grid for each direction
     */
    void load_phimap_from_file(int gsize);
    /**
     * @brief compute_grid_hardcore_omp This function computes hardcore (Amber FF) potential using multicore
     * implementation with OpenMP.
     * @param Rec class Mol2 object with receptor description
     */
    void compute_grid_hardcore_omp(Mol2* Rec);
    void compute_grid_hardcore_HB_omp(Mol2* Rec);
    void compute_grid_hardcore_omp_HB_Gaussian(Mol2* Rec);

    /**
     * @brief load_delphi_cube This function loads a DelPhi potential in
     * a CUBE (ASCII) format.
     */
    void load_delphi_cube(void);
    void load_delphi_gzcube(void);

    void print_line();
    void print_info(char info[98]);

	virtual ~Grid();

};

#endif /* GRID_H_ */
