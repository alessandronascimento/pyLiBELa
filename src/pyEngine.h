#include <vector>
#include <iostream>
#include <string>
#include <gsl/gsl_rng.h>

#include "pyGrid.h"
#include "pyPARSER.h"
#include "pyCOORD_MC.h"
#include "pyWRITER.h"
#include "pyDocker.h"
#include "pyMol2.h"
#include "pyMC.h"
#include "pyFullSearch.h"
#include "pyConformer.h"
#include "pyFindHB.h"

using namespace std;

class Engine{
public:
	char info[98];
	WRITER* Writer;
	Engine(WRITER* _Writer);
	void Run_docking(PARSER* Input, Mol2* Rec, Mol2* RefLig, Grid* Grids);
	void print_info(char info[98]);

};
