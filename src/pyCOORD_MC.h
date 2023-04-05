/*
 * COORD_MC.h
 *
 *  Created on: 16/04/2010
 *      Author: Nascimento
 */

#ifndef COORD_MC_H_
#define COORD_MC_H_

#include<vector>
#include "pyRAND.cpp"
#include<iostream>
#include<cmath>
#include"pyMol2.h"
#include "iMcLiBELa.h"

using namespace std;

/**
 * The C++ class \a COORD_MC handles the coordinates changes necessary throughout the code.
 * It includes methods for translating and rotating cartesian coordinates around Euler angles,
 * Center of Mass (COM) computation and RMSD computation.
 */

class COORD_MC {
public:

	// VARIABLES

	//! Atomic coordinate X
	double x,
	//! Atomic coordinate Y
	y,
	//! Atomic coordinate Z
	z,
	//! X component of the center of mass
	centerx,
	//! Y component of the center of mass
	centery,
	//! Z component of the center of mass
	centerz,
	//! Sum over the total atoms in the system (mainly, the ligand). Used to compute the center of mass.
	totalmass;
	//! C++ vector that stores the atomic coordinates
	vector<double> coordinates;
	//! C++ vector that stores the modified atomic coordinates (rotated/translated)
//	vector<vector<double> >new_coordinates;
	//! C++ vector that stores the center of mass for the system (ligand)
	vector<double> com;
	//! Root mean square deviation. Used to compute the deviation among original and rotated/translated ligand coordinates.
	double rmsd;
	//! Variation in binding energy used to compute Metropolis criterium.
	double delta_energy;
	//! Alpha angle to rotate around this Euler angle.
	double alpha;
	//! Beta angle to rotate around this Euler angle.
	double beta;
	//! Gamma angle to rotate around this Euler angle.
	double gamma;
	//! X-shift for molecule translation
	double transx;
	//! Y-shift for molecule translation
	double transy;
	//! Z-shift for molecule translation
	double transz;

	//! Pointer to a PRMTOP object
	Mol2 *Cmol;

	// METHODS

	COORD_MC();
	/*! This method computed the center of mass of a ligand given by the object Cmol (Class PRMTOP)
	 * @param Cmol Pointer to the PRMTOP class with molecule properties.
	 * @return A C++ vector with three elements: X, Y and Z components of the center of mass.
	 */
	vector<double> compute_com(Mol2 *Cmol);

	/*! This overloaded method computed the center of mass of a ligand given by the object
	 *  Cmol (Class PRMTOP)
	 *  @param coords C++ vector with ligand coordinates
	 *  @param Cmol Pointer to a PRMTOP object with ligand description (atributes)
	 *  @return A C++ vector with three elements: the X, Y and Z components of the center of mass.
	 */
	vector<double> compute_com(vector<vector<double> >coords, Mol2 *Cmol);

	/*! This method shifts (i.e. translates) a molecule to any desired position.
	 * @param coordinates C++ vector with molecular (cartesian) coordinates.
	 * @param N Number of atoms in the molecule.
	 * @param tx Shift in X direction.
	 * @param ty Shift in Y direction.
	 * @param tz Shift in Z direction.
	 * @return A new C++ vector with new (shifted) coordinates.
	 */
	vector<vector<double> >translate(vector<vector<double> >coordinates, int N, double tx, double ty, double tz);

	/*!
	 * Method to rotate a 3D system in Euler angles (alpha, beta and gamma). To keep it simpler,
	 * the rotation is carried out around the system center of mass. That's the reason why COM
	 * is required in this method.
	 * @param coordinates C++ vector with molecular (cartesian) coordinates.
	 * @param N Number of atoms in the molecule.
	 * @param alpha Alpha angle for rotation.
	 * @param beta Beta angle for rotation.
	 * @param gamma Gamma angle for rotation.
	 * @return A new C++ vector with rotated coordinates.
	 */
	vector<vector<double> >rotate(vector<vector<double> >coordinates, int N, double alpha, double beta, double gamma);

	/*! Method to rotate (in Euler angles) and translate a set of coordinates at once.
	 * @param coordinates C++ vector with original ligand coordinates
	 * @param N Number of atoms
	 * @param alpha First Euler angle to define a rotation (see http://www.gregslabaugh.name/publications/euler.pdf for a reference)
	 * @param beta Second Euler angle.
	 * @param gamma Third Euler angle.
	 * @param transx Shift in X direction.
	 * @param transy Shift in Y direction.
	 * @param transz Shift in Z direction.
	 * @return C++ vector with modified coordinates
	 */
	vector<vector<double> >rototranslate(vector<vector<double> >coordinates, int N, double alpha, double beta, double gamma, double transx, double transy, double transz);
	vector<vector<double> >rototranslate(vector<vector<double> >coordinates, Mol2* Lig, double alpha, double beta, double gamma, double transx, double transy, double transz);
	vector<vector<double> >rototranslate(vector<vector<double> >coordinates, Mol2* LIG, RAND* Rand);
	/*! This method rototranslates a set of coordinates from a MD simulation, for example.
	 *  The new coordinates are written in the PRMTOP attribute 'new_mcoord'.
	 *
	 * @param Cmol Pointer to a PRMTOP object defining a ligand.
	 * @param Rand Ponter to a RAND object with sorted angles and shifts.
	 */

	void rototranslate_all(Mol2 *Cmol, RAND *Rand);

	/*! This overloaded method rototranslates a set of coordinates from a MD simulation, for example.
	 *  The new coordinates are written in the PRMTOP attribute 'new_mcoord'.
	 *
	 * @param Cmol Pointer to a PRMTOP object defining a ligand.
	 * @param alpha First Euler angle to define a rotation (see http://www.gregslabaugh.name/publications/euler.pdf for a reference)
	 * @param beta Second Euler angle.
	 * @param gamma Third Euler angle.
	 * @param transx Shift in X direction.
	 * @param transy Shift in Y direction.
	 * @param transz Shift in Z direction.
	 */

	void rototranslate_all(Mol2 *Cmol, double alpha, double beta, double gamma, double transx, double transy, double transz);

	/*! This overloaded method rototranslates a ligand applying a rotations around Euler angles and shifts
	 *  in X, Y and Z direction.
	 *  @param coordinates C++ vector with ligand coordinates.
	 *  @param N Number of atoms in the ligand.
	 *  @param rand Pointer to the RAND object with sorted (random) numbers for Euler rotation and for translation.
	 *  @return C++ vector with modified coordinates
	 */

	vector<vector<double> >rototranslate(vector<vector<double> >coordinates, int N, RAND* rand);

	/*! Method to compute the root mean square deviation among the original ligand coordinates and
	 *  the rotated/translated coordinates.
	 * @param coordinates C++ vector with ligand coordinates
	 * @param new_coord C++ vector with modified coordinates
	 * @param N Number of atoms
	 * @return The computed RMSD (double precision)
	 */
	double compute_rmsd(vector<vector<double> >coordinates, vector<vector<double> >new_coord, int N);

	/*! This method computes the probability for acceptance for a move usign Metropolis criterium.
	 * @param old_energy Energy of the system prior to the attempt to move.
	 * @param new_energy New energy of the system after the movement.
	 * @param temp Temperature of the system.
	 * @return Probability (double precision)
	 * The probability is computed as:
	 *
	 * \f$ prob = \exp(\frac{-\Delta E}{kT}) = \exp(\frac{-\Delta E}{(0.00198587752*T)})\f$.
	 *
	 * where
	 *
	 * \f$ \Delta E = new  energy - old  energy \f$
	 */
	double compute_prob(double old_energy, double new_energy, double temp);
};

#endif /* COORD_MC_H_ */
