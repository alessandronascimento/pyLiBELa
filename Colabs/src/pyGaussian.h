/*
 * Gaussian.h
 *
 *  Created on: 18/11/2011
 *      Author: Nascimento
 */

#ifndef GAUSSIAN_H_
#define GAUSSIAN_H_

#include<stdio.h>
#include<iostream>
#include<cmath>
#include <ctime>
#include <vector>
#include "iMcLiBELa.h"
#include"pyMol2.h"
#include"pyPARSER.h"

using namespace std;

class Gaussian {
public:
	clock_t time0, time1;

	Gaussian(void);
	Gaussian(Mol2* RefMol, Mol2* CompMol, double* Vab);
	Gaussian(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz, double* Vab);
	double compute_si(Mol2* RefMol, Mol2* CompMol);
	double compute_si(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz);
	double compute_shape_density(Mol2* RefMol, Mol2* CompMol);
	double compute_shape_density(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz);
	double compute_shape_density(Mol2* CompMol, vector<vector<double> > xyz);
	double compute_shape_and_charge_density(PARSER *Input, Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz);
	double compute_shape_and_charge_density(PARSER *Input, Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz, vector<vector<double> > cxyz);
	double dist_squared(double x1, double x2, double y1, double y2, double z1, double z2);
	virtual ~Gaussian();
};

#endif /* GAUSSIAN_H_ */
