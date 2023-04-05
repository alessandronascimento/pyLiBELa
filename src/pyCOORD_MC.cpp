/*
 * COORD_MC.cpp
 *
 *  Created on: 16/04/2010
 *      Author: Nascimento
 */

#include "pyCOORD_MC.h"

using namespace std;


COORD_MC::COORD_MC(){
}

vector<double> COORD_MC::compute_com(Mol2 *Cmol){
	double centerx=0.0;
	double centery=0.0;
	double centerz=0.0;
	double totalmass=0.0;
	vector<double>com (3);
	for(int i=0; i<Cmol->N; i++){
		centerx+= (Cmol->masses[i]*Cmol->xyz[i][0]);
		centery+= (Cmol->masses[i]*Cmol->xyz[i][1]);
		centerz+= (Cmol->masses[i]*Cmol->xyz[i][2]);
		totalmass+= Cmol->masses[i];
	}
	com[0] = (centerx/totalmass);
	com[1] = (centery/totalmass);
	com[2] = (centerz/totalmass);
	return(com);
}

vector<double> COORD_MC::compute_com(vector<vector<double> >coords, Mol2 *Cmol){
	double centerx=0.0;
	double centery=0.0;
	double centerz=0.0;
	double totalmass=0.0;
	vector<double> com(3);
	for(int i=0; i<Cmol->N; i++){
		centerx+= (Cmol->masses[i]*coords[i][0]);
		centery+= (Cmol->masses[i]*coords[i][1]);
		centerz+= (Cmol->masses[i]*coords[i][2]);
		totalmass+= Cmol->masses[i];
	}
	com[0] = (centerx/totalmass);
	com[1] = (centery/totalmass);
	com[2] = (centerz/totalmass);
	return(com);
}

vector<vector<double> >COORD_MC::rotate(vector<vector<double> >coordinates, int N, double alpha, double beta, double gamma){
	vector<vector<double> >new_coordinates;
	vector<double> txyz;
	for(int i=0; i < N ; i++){
		x=coordinates[i][0];
		y=coordinates[i][1];
		z=coordinates[i][2];
		txyz.push_back(((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))));
		txyz.push_back(((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))));
		txyz.push_back(((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)));
		new_coordinates.push_back(txyz);
		txyz.clear();
	}
	return(new_coordinates);
}

vector<vector<double> >COORD_MC::translate(vector<vector<double> >coordinates, int N, double tx, double ty, double tz){
	vector<vector<double> >new_coordinates;
	vector<double> txyz;
	for (int i=0; i < N; i++){
		txyz.push_back(coordinates[i][0] + tx);
		txyz.push_back(coordinates[i][1] + ty);
		txyz.push_back(coordinates[i][2] + tz);
		new_coordinates.push_back(txyz);
		txyz.clear();
	}
	return(new_coordinates);
}

vector<vector<double> >COORD_MC::rototranslate(vector<vector<double> >coordinates, Mol2* LIG, RAND* Rand){
	vector<vector<double> >new_coordinates;
	vector<double> COM;
	COM = this->compute_com(coordinates, LIG);
	vector<double> txyz;
	double alpha = Rand->a;
	double beta = Rand->b;
	double gamma = Rand->g;
	for(int i=0; i < LIG->N ; i++){
		x=coordinates[i][0]-COM[0];
		y=coordinates[i][1]-COM[1];
		z=coordinates[i][2]-COM[2];
		txyz.push_back((((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))))+Rand->transx+COM[0]);
		txyz.push_back((((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))))+Rand->transy + COM[1]);
		txyz.push_back((((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)))+Rand->transz + COM[2]);
		new_coordinates.push_back(txyz);
		txyz.clear();
	}
	return(new_coordinates);
}

vector<vector<double> >COORD_MC::rototranslate(vector<vector<double> >coordinates, int N, double alpha, double beta, double gamma, double transx, double transy, double transz){
    vector<vector<double> >new_coordinates(N);
    vector<double> txyz(3);
	double x, y, z;
    for(unsigned i=0; i < unsigned(N); i++){
		x=coordinates[i][0];
		y=coordinates[i][1];
		z=coordinates[i][2];
        txyz[0] = ((((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))))+transx);
        txyz[1] = ((((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))))+transy);
        txyz[2] = ((((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)))+transz);
        new_coordinates[i] = txyz;
	}
	return(new_coordinates);
}

vector<vector<double> >COORD_MC::rototranslate(vector<vector<double> >coordinates, Mol2* Lig, double alpha, double beta, double gamma, double transx, double transy, double transz){
    vector<vector<double> >new_coordinates(Lig->N);
	vector<double> txyz(3);
    vector<double> COM(3);
    COM = this->compute_com(coordinates, Lig);
	double x, y, z;
	for(int i=0; i < Lig->N ; i++){
		x=coordinates[i][0]-COM[0];
		y=coordinates[i][1]-COM[1];
		z=coordinates[i][2]-COM[2];
		txyz[0] = ((((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))))+transx+COM[0]);
		txyz[1] = ((((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))))+transy + COM[1]);
		txyz[2] = ((((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)))+transz + COM[2]);
        new_coordinates[i] = txyz;
	}
    return(new_coordinates);
}

void COORD_MC::rototranslate_all(Mol2 *Cmol, RAND *Rand){
	alpha = Rand->a;
	beta = Rand->b;
	gamma = Rand->g;
	transx = Rand->transx;
	transy = Rand->transy;
	transz = Rand->transz;
	vector<double> txyz;
	vector<vector<double> >tcoords;
	Cmol->new_mcoords.clear();
	for (unsigned i=0; i<Cmol->mcoords.size(); i++){
        com = compute_com(Cmol->xyz, Cmol);
		for (int j=0; j<Cmol->N; j++){
			x=Cmol->mcoords[i][j][0]-com[0];
			y=Cmol->mcoords[i][j][1]-com[1];
			z=Cmol->mcoords[i][j][2]-com[2];

			txyz.push_back(((((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))))+transx));
			txyz.push_back(((((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))))+transy));
			txyz.push_back(((((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)))+transz));

			tcoords.push_back(txyz);
			txyz.clear();
		}
		Cmol->new_mcoords.push_back(tcoords);
		tcoords.clear();
	}
}

void COORD_MC::rototranslate_all(Mol2 *Cmol, double alpha, double beta, double gamma, double transx, double transy, double transz){
    vector<double> txyz (3);
	vector<vector<double> >tcoords;
	Cmol->new_mcoords.clear();
    double x, y, z;
    vector<double> COM;
	for (unsigned i=0; i<Cmol->mcoords.size(); i++){
        COM = compute_com(Cmol->mcoords[i], Cmol);
		for (int j=0; j<Cmol->N; j++){
            x=Cmol->mcoords[i][j][0]-COM[0];
            y=Cmol->mcoords[i][j][1]-COM[1];
            z=Cmol->mcoords[i][j][2]-COM[2];

            txyz[0] = ((((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))))+transx+COM[0]);
            txyz[1] = ((((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))))+transy + COM[1]);
            txyz[2] = ((((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)))+transz + COM[2]);

			tcoords.push_back(txyz);
		}
		Cmol->new_mcoords.push_back(tcoords);
		tcoords.clear();
	}
}

vector<vector<double> >COORD_MC::rototranslate(vector<vector<double> >coordinates, int N, RAND* Rand){
	alpha = Rand->a;
	beta = Rand->b;
	gamma = Rand->g;
	transx = Rand->transx;
	transy = Rand->transy;
	transz = Rand -> transz;
	vector<double>  txyz;
	vector<vector<double> >new_coordinates;
	for(int i=0; i < N ; i++){
		x=coordinates[i][0];
		y=coordinates[i][1];
		z=coordinates[i][2];
		txyz.push_back((((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))))+transx);
		txyz.push_back((((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))))+transy);
		txyz.push_back((((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)))+transz);
		new_coordinates.push_back(txyz);
		txyz.clear();
	}
	return(new_coordinates);
}



double COORD_MC::compute_rmsd(vector<vector<double> >coordinates, vector<vector<double> >new_coord, int N){
    double lrmsd = 0.00;
	for (int i=0; i < N; i++){
		for (int j=0; j<3; j++){
            lrmsd+= (coordinates[i][j]-new_coord[i][j])*(coordinates[i][j]-new_coord[i][j]);
		}
	}
    lrmsd = (sqrt(lrmsd/N));
    return(lrmsd);
}

double COORD_MC::compute_prob(double old_energy, double new_energy, double temp){
//	printf("%.2f  %.2f  %.2f\n", new_energy, old_energy, temp);
// 	k=0.0019858775203792202 kcal/molK
	delta_energy = new_energy - old_energy;
	return(exp((-delta_energy)/(0.00198587752*temp)));
}






#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;


BOOST_PYTHON_MODULE(pyCOORD_MC)
{

    vector<double> (COORD_MC::*cc1)(Mol2 *Cmol) = &COORD_MC::compute_com;
    vector<double> (COORD_MC::*cc2)(vector<vector<double> >coords, Mol2* Cmol)   = &COORD_MC::compute_com;

    vector<vector<double> > (COORD_MC::*rt1)(vector<vector<double> >coordinates, int N, double alpha, double beta, double gamma, double transx, double transy, double transz) = &COORD_MC::rototranslate;
    vector<vector<double> >(COORD_MC::*rt2)(vector<vector<double> >coordinates, Mol2* Lig, double alpha, double beta, double gamma, double transx, double transy, double transz) = &COORD_MC::rototranslate;
    vector<vector<double> >(COORD_MC::*rt3)(vector<vector<double> >coordinates, Mol2* LIG, RAND* Rand) = &COORD_MC::rototranslate;
    vector<vector<double> >(COORD_MC::*rt4)(vector<vector<double> >coordinates, int N, RAND* rand) = &COORD_MC::rototranslate;



    void (COORD_MC::*rta1)(Mol2 *Cmol, RAND *Rand) = &COORD_MC::rototranslate_all;
    void (COORD_MC::*rta2)(Mol2 *Cmol, double alpha, double beta, double gamma, double transx, double transy, double transz) = &COORD_MC::rototranslate_all;



    class_<COORD_MC>("COORD_MC", init< >())
        .def_readwrite("x", & COORD_MC::x)
        .def_readwrite("y", & COORD_MC::y)
        .def_readwrite("z", & COORD_MC::z)
        .def_readwrite("centerx", & COORD_MC::centerx)
        .def_readwrite("centery", & COORD_MC::centery)
        .def_readwrite("centerz", & COORD_MC::centerz)
        .def_readwrite("totalmass", & COORD_MC::totalmass)
        .def_readwrite("coordinates", & COORD_MC::coordinates)
        .def_readwrite("com", & COORD_MC::com)
        .def_readwrite("rmsd", & COORD_MC::rmsd)
        .def_readwrite("delta_energy", & COORD_MC::delta_energy)
        .def_readwrite("alpha", & COORD_MC::alpha)
        .def_readwrite("beta", & COORD_MC::beta)
        .def_readwrite("gamma", & COORD_MC::gamma)
        .def_readwrite("transx", & COORD_MC::transx)
        .def_readwrite("transy", & COORD_MC::transy)
        .def_readwrite("transz", & COORD_MC::transz)
        .def_readwrite("Cmol", & COORD_MC::Cmol)


        .def("compute_com", cc1)
        .def("compute_com", cc2)

        .def("translate", & COORD_MC::translate)
        .def("rotate", & COORD_MC::rotate)

        .def("rototranslate", rt1)
        .def("rototranslate", rt2)
        .def("rototranslate", rt3)
        .def("rototranslate", rt4)


        .def("rototranslate_all", rta1)
        .def("rototranslate_all", rta2)

        .def("compute_rmsd", & COORD_MC::compute_rmsd)
        .def("compute_prob", & COORD_MC::compute_prob)

    ;


}


