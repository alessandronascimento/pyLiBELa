/*
 * Gaussian.cpp
 *
 *  Created on: 18/11/2011
 *      Author: Nascimento
 */

#include "pyGaussian.h"

Gaussian::Gaussian(Mol2* RefMol, Mol2* CompMol, double* Vab) {
	time0=clock();
	printf("Refmol N: %d    CompMol N: %d\n", RefMol->N, CompMol->N);
	*Vab = this->compute_si(RefMol, CompMol);
	time1=clock();

#ifdef DEBUG
	printf("Gaussian computation took %7.3f seconds.\n", ((double(time1)-double(time0))/CLOCKS_PER_SEC));
#endif

}

Gaussian::Gaussian(void){
}

Gaussian::Gaussian(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz, double* Vab) {
	time0=clock();
	printf("Refmol N: %d    CompMol N: %d\n", RefMol->N, CompMol->N);
	*Vab = this->compute_si(RefMol, CompMol, xyz);
	time1=clock();

#ifdef DEBUG
	printf("Gaussian computation took %7.3f seconds.\n", ((double(time1)-double(time0))/CLOCKS_PER_SEC));
#endif

}

Gaussian::~Gaussian(){

}

double Gaussian::compute_shape_density(Mol2* RefMol, Mol2* CompMol){
	double V=0.0;
	double pi=2.0*sqrt(2);
	double pj=2.0*sqrt(2);
	double dij2;
	double alpha_i, alpha_j;
	for (int i=0; i<RefMol->N; i++){
		alpha_i = PI*( pow( ( (3*pi)/(4*PI*pow(RefMol->radii[i], 3))), (2.0/3.0)));
		for (int j=0; j<CompMol->N; j++){
			alpha_j = PI*pow( ((3*pi)/(4*PI*pow(CompMol->radii[j], 3))), (2.0/3.0));
			dij2 = dist_squared(RefMol->xyz[i][0], CompMol->xyz[j][0],RefMol->xyz[i][1], CompMol->xyz[j][1],RefMol->xyz[i][2], CompMol->xyz[j][2]);
			V+= pi*pj*(exp((-alpha_i*alpha_j*dij2)/(alpha_i+alpha_j)))*(pow((PI/(alpha_i+alpha_j)), (3.0/2.0)));
		}
	}
	return(V);
}

double Gaussian::compute_shape_density(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz){
	double V=0.0;
	double pi=2.0*sqrt(2);
	double pj=2.0*sqrt(2);
	double dij2;
	double alpha_i, alpha_j;
	for (int i=0; i<RefMol->N; i++){
		alpha_i = PI*( pow( ( (3*pi)/(4*PI*pow(RefMol->radii[i], 3))), (2.0/3.0)));
		for (int j=0; j<CompMol->N; j++){
			alpha_j = PI*pow(  ( (3*pi)/(4*PI*pow(CompMol->radii[j], 3))), (2.0/3.0));
			dij2 = dist_squared(RefMol->xyz[i][0], xyz[j][0],RefMol->xyz[i][1], xyz[j][1],RefMol->xyz[i][2], xyz[j][2]);
			V+= pi*pj*(exp((-alpha_i*alpha_j*dij2)/(alpha_i+alpha_j)))*(pow((PI/(alpha_i+alpha_j)), (3.0/2.0)));
		}
	}
	return(V);
}

double Gaussian::compute_shape_density(Mol2* CompMol, vector<vector<double> > xyz){
	double V=0.0;
	double pi=2.0*sqrt(2);
	double pj=2.0*sqrt(2);
	double dij2;
	double alpha_i, alpha_j;
	for (int i=0; i<CompMol->N; i++){
		alpha_i = PI*( pow( ( (3*pi)/(4*PI*pow(CompMol->radii[i], 3))), (2.0/3.0)));
		for (int j=0; j<CompMol->N; j++){
			alpha_j = PI*pow(  ( (3*pi)/(4*PI*pow(CompMol->radii[j], 3))), (2.0/3.0));
			dij2 = dist_squared(xyz[i][0], xyz[j][0],xyz[i][1], xyz[j][1],xyz[i][2], xyz[j][2]);
			V+= pi*pj*(exp((-alpha_i*alpha_j*dij2)/(alpha_i+alpha_j)))*(pow((PI/(alpha_i+alpha_j)), (3.0/2.0)));
		}
	}
	return(V);
}

double Gaussian::compute_si(Mol2* RefMol, Mol2* CompMol){
	double si = ((2*this->compute_shape_density(RefMol, CompMol))/(this->compute_shape_density(RefMol, RefMol) + this->compute_shape_density(CompMol, CompMol)));
	printf("Shape density correlation: %.3f\n", si);
	return (si);
}

double Gaussian::compute_si(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz){
	double t1, t2, t3, si;
	t1= this->compute_shape_density(RefMol, CompMol, xyz);
	t2 = this->compute_shape_density(RefMol, RefMol);
	t3 = this->compute_shape_density(CompMol, xyz);
	si = 2*t1/(t2+t3);

#ifdef DEBUG
	printf("Shape density correlation: %.3f\n", si);
#endif
	return (si);
}

double Gaussian::dist_squared(double x1, double x2, double y1, double y2, double z1, double z2){
	double r2 = ((x2-x1)*(x2-x1)) + ((y2-y1)*(y2-y1)) + ((z2-z1)*(z2-z1));
	return (r2);
}

double Gaussian::compute_shape_and_charge_density(PARSER *Input, Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz){
	double pi=2.0*sqrt(2);
    double pj=2.0*sqrt(2);
    double dij2;
    double Vshape=0.0, Vpos=0.0, Vneg=0.0;
    double alphai_shape=0.0, alphai_pos=0.0, alphai_neg=0.0;
    double alphaj_shape=0.0, alphaj_pos=0.0, alphaj_neg=0.0;
    for (int i=0; i<RefMol->N; i++){
        if (RefMol->radii[i] != 0.0){
            alphai_shape = PI*( pow( ( (3*pi)/(4*PI*pow(RefMol->radii[i], 3))), (2.0/3.0)));
        }
        else {
            alphai_shape=0.0;
        }

        if (RefMol->charges[i] > 0.0){
            alphai_pos = PI*( pow( ( (3*pi)/(4*PI*pow((1.00 + RefMol->charges[i]), 3))), (2.0/3.0)));
            alphai_neg=0.0;
        }

        else if (RefMol->charges[i] < 0.0) {
            alphai_neg = PI*( pow( ( (3*pi)/(4*PI*pow((1.00 + abs(RefMol->charges[i])), 3))), (2.0/3.0)));
            alphai_pos=0.0;
        }

        for (int j=0; j<CompMol->N; j++){
            if (CompMol->radii[j] != 0.0000){
                alphaj_shape = PI*pow(  ( (3*pi)/(4*PI*pow(CompMol->radii[j], 3))), (2.0/3.0));
            }
            else {
                alphaj_shape=0.0000;
            }

            if (CompMol->charges[j] > 0.00){
                alphaj_pos = PI*pow(  ( (3*pi)/(4*PI*pow((1.00 + CompMol->charges[j]), 3))), (2.0/3.0));
                alphaj_neg=0.0;
            }

            else if (CompMol->charges[j] < 0.00){
                alphaj_neg = PI*pow(  ( (3*pi)/(4*PI*pow((1.00 + abs(CompMol->charges[j])), 3))), (2.0/3.0));
                alphaj_pos = 0.0;
            }

            dij2 = dist_squared(RefMol->xyz[i][0], xyz[j][0],RefMol->xyz[i][1], xyz[j][1],RefMol->xyz[i][2], xyz[j][2]);

            if (alphai_shape != 0.0 and alphaj_shape != 0.0 ){
                Vshape += pi*pj*(exp((-alphai_shape*alphaj_shape*dij2)/(alphai_shape+alphaj_shape)))*(pow((PI/(alphai_shape+alphaj_shape)), (3.0/2.0)));
            }

            if (alphai_pos != 0.0 and alphaj_pos != 0.0 ){
                Vpos += pi*pj*(exp((-alphai_pos*alphaj_pos*dij2)/(alphai_pos+alphaj_pos)))*(pow((PI/(alphai_pos+alphaj_pos)), (3.0/2.0)));
            }

            if (alphai_neg != 0.0 and alphaj_neg != 0.0 ){
                Vneg += pi*pj*(exp((-alphai_neg*alphaj_neg*dij2)/(alphai_neg+alphaj_neg)))*(pow((PI/(alphai_neg+alphaj_neg)), (3.0/2.0)));
            }
        }
    }
    return((Input->vdw_scale*Vshape)+(Input->elec_scale*Vpos)+(Input->elec_scale*Vneg));
}

double Gaussian::compute_shape_and_charge_density(PARSER *Input, Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz, vector<vector<double> > cxyz){
	double pi=2.0*sqrt(2);
	double pj=2.0*sqrt(2);
	double dij2;
    double Vshape=0.0, Vpos=0.0, Vneg=0.0;
    double alphai_shape=0.0, alphai_pos=0.0, alphai_neg=0.0;
    double alphaj_shape=0.0, alphaj_pos=0.0, alphaj_neg=0.0;
    for (int i=0; i<RefMol->N; i++){
        if (RefMol->radii[i] != 0.0){
            alphai_shape = PI*( pow( ( (3*pi)/(4*PI*pow(RefMol->radii[i], 3))), (2.0/3.0)));
        }
        else {
            alphai_shape=0.0;
        }

        if (RefMol->charges[i] > 0.0){
            alphai_pos = PI*( pow( ( (3*pi)/(4*PI*pow((1.00 + RefMol->charges[i]), 3))), (2.0/3.0)));
            alphai_neg=0.0;
        }

        else if (RefMol->charges[i] < 0.0) {
            alphai_neg = PI*( pow( ( (3*pi)/(4*PI*pow((1.00 + abs(RefMol->charges[i])), 3))), (2.0/3.0)));
            alphai_pos=0.0;
        }

        for (int j=0; j<CompMol->N; j++){
            if (CompMol->radii[j] != 0.0000){
                alphaj_shape = PI*pow(  ( (3*pi)/(4*PI*pow(CompMol->radii[j], 3))), (2.0/3.0));
            }
            else {
                alphaj_shape=0.0000;
            }

            if (CompMol->charges[j] > 0.00){
                alphaj_pos = PI*pow(  ( (3*pi)/(4*PI*pow((1.00 + CompMol->charges[j]), 3))), (2.0/3.0));
                alphaj_neg=0.0;
            }

            else if (CompMol->charges[j] < 0.00){
                alphaj_neg = PI*pow(  ( (3*pi)/(4*PI*pow((1.00 + abs(CompMol->charges[j])), 3))), (2.0/3.0));
                alphaj_pos = 0.0;
            }

            dij2 = this->dist_squared(xyz[i][0], cxyz[j][0],xyz[i][1], cxyz[j][1],xyz[i][2], cxyz[j][2]);

            if (alphai_shape != 0.0 and alphaj_shape != 0.0 ){
                Vshape += pi*pj*(exp((-alphai_shape*alphaj_shape*dij2)/(alphai_shape+alphaj_shape)))*(pow((PI/(alphai_shape+alphaj_shape)), (3.0/2.0)));
            }

            if (alphai_pos != 0.0 and alphaj_pos != 0.0 ){
                Vpos += pi*pj*(exp((-alphai_pos*alphaj_pos*dij2)/(alphai_pos+alphaj_pos)))*(pow((PI/(alphai_pos+alphaj_pos)), (3.0/2.0)));
            }

            if (alphai_neg != 0.0 and alphaj_neg != 0.0 ){
                Vneg += pi*pj*(exp((-alphai_neg*alphaj_neg*dij2)/(alphai_neg+alphaj_neg)))*(pow((PI/(alphai_neg+alphaj_neg)), (3.0/2.0)));
            }
        }
    }
    return((Input->vdw_scale*Vshape)+(Input->elec_scale*Vpos)+(Input->elec_scale*Vneg));
}


#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;


BOOST_PYTHON_MODULE(pyGaussian)
{
    double (Gaussian::*cs1)(Mol2* RefMol, Mol2* CompMol) = &Gaussian::compute_si;
    double (Gaussian::*cs2)(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz)= &Gaussian::compute_si;

    double (Gaussian::*csd1)(Mol2* RefMol, Mol2* CompMol) = &Gaussian::compute_shape_density;
    double (Gaussian::*csd2)(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz) = &Gaussian::compute_shape_density;
    double (Gaussian::*csd3)(Mol2* CompMol, vector<vector<double> > xyz) = &Gaussian::compute_shape_density;

    double (Gaussian::*csacd1)(PARSER *Input, Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz) = &Gaussian::compute_shape_and_charge_density;
    double (Gaussian::*csacd2)(PARSER *Input, Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz, vector<vector<double> > cxyz) = &Gaussian::compute_shape_and_charge_density;

    class_<Gaussian>("Gaussian", init<>())
        .def(init<Mol2*, Mol2*, double*>())
        .def(init<Mol2*, Mol2*, vector<vector<double> >, double*>())
        .def_readwrite("time0", & Gaussian::time0)
        .def_readwrite("time1", & Gaussian::time1)

        .def("compute_si",cs1)
        .def("compute_si",cs2)

        .def("compute_shape_density",csd1)
        .def("compute_shape_density",csd2)
        .def("compute_shape_density",csd3)

        .def("compute_shape_and_charge_density",cs1)
        .def("compute_shape_and_charge_density",cs2)

        .def("dist_squared", & Gaussian::dist_squared)
    ;
}
