#ifndef FULLSEARCH_H
#define FULLSEARCH_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "Mol2.h"
#include "PARSER.h"
#include "Conformer.h"
#include "WRITER.h"
#include "COORD_MC.h"
#include "Grid.h"
#include "Energy2.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/rotor.h>
#include <openbabel/conformersearch.h>


using namespace std;


class FullSearch
{
public:
    PARSER* Input;
    Mol2* RefLig;
    Mol2* Lig;
    Conformer* Conf;
    WRITER* Writer;
    OBMol mol;
    OBForceField* OBff;
    OBRotorList RotorList;
    OBRotorIterator RotorIterator;
    OBRotor *Rotor;
    char info[98];
    Grid* Grids;
    vector<vector<int> > atoms_in_dihedrals;
    double* myxyz;

    FullSearch(PARSER* _Input, Mol2* _Lig, WRITER* Writer);
    FullSearch(PARSER* _Input, Mol2* _Lig, WRITER* Writer, Grid* _Grids);
    OBMol GetMol(const std::string &molfile);
    void copy_to_obmol(vector<vector<double> > vec_xyz);
    vector<vector<double> > copy_from_obmol(OBMol mymol);
    double do_search(void);
};

#endif // FULLSEARCH_H
