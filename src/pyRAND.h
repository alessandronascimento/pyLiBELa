/*
 * RAND.h
 *
 *  Created on: 07/10/2010
 *      Author: Nascimento
 */

#ifndef RAND_H_
#define RAND_H_

#include "pyPARSER.h"
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<cstdlib>
#include"pyMol2.cpp"
#include<gsl/gsl_rng.h>

class RAND {
public:
//! Alpha angle in Euler notation. Used for rotation.
	double a,
//! Beta angle in Euler notation. Used for rotation.
	b,
//! Gamma angle in Euler notation. Used for rotation.
	g,
//! Shift in X axis (X translation)
	transx,
//! Shift in Y axis (Y translation)
	transy,
//!Shift in Z axis (Z translation)
	transz;
//!
	int recn;
//!
	int lign;





//! Random number. Sorted using STDLIB.
	double rnumber;

//	const gsl_rng_type *T0, *T;
	gsl_rng * r;


//	int conformation;
//! The Parser class is used as a pointer in this class

//	PARSER* Input;

	Mol2 *CRec;
	Mol2 *CLig;

/*!
 * The random method is used to sort 6 random numbers that will
 * be used to define a movement within 3 Euler angles (alpha, beta
 * and gamma) and some translation in coordinates (x, y and z shift).
 */
	RAND();
    void random(double cushion, double rotation_step, Mol2 *CRec, Mol2 *CLig);
	void random( double cushion, double rotation_step);

/*!
 * The print method prints the 6 sorted random numbers.
 */
	void print();
};

#endif /* RAND_H_ */
