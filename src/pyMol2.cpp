/*
 * Mol2.cpp
 *
 *  Created on: 10/10/2011
 *      Author: Nascimento
 */

#include "pyMol2.h"

using namespace std;

Mol2::Mol2(){
}

Mol2::Mol2(PARSER *Input, string molfile) {
    bool ok;
    if ((molfile.substr(molfile.size()-3, 3) == ".gz") or (molfile.substr(molfile.size()-2, 2) == ".z")){
        ok = this->parse_gzipped_file(Input, molfile);
        if (! ok){
            printf("Could not correctly parse mol2 file %s. Please check.\n", molfile.c_str());
            exit(1);
        }
    }
    else {
        ok = this->parse_mol2file(Input, molfile);
        if (! ok){
            printf("Could not correctly parse mol2 file %s. Please check.\n", molfile.c_str());
            exit(1);
        }
    }
}

bool Mol2::parse_mol2file(PARSER *Input, string molfile) {
	FILE *mol2file;
	int tint;
	float tx, ty, tz;
	vector<double> txyz;
	int tres;
	float tcharge;
	int count=0;
	char tatomtype[10];
    char resname[20];
	string cpstr;
	bool bret = false;

    this->initialize_gaff();


	mol2file = fopen(molfile.c_str(), "r");

	if (mol2file !=NULL){
		str[0]='#';
		while(str[0] !='@'){
			fgets(str, 80, mol2file);

		}
		fgets(str, 80, mol2file);
		this->molname = str;
        this->molname = this->molname.substr(0,this->molname.size()-1);
		fscanf(mol2file, "%d %d %d %d %d", &this->N, &this->Nbonds, &this->Nres, &tint, &tint);

		cpstr = string(str);
		while (cpstr.substr(0,13) != "@<TRIPOS>ATOM"){
			fgets(str, 80, mol2file);
			cpstr = string(str);
		}

		for (int i=0; i<this->N; i++){
            fscanf(mol2file, "%d %s %f %f %f %s %d %s %f\n", &tint, str, &tx, &ty, &tz, tatomtype, &tres, resname, &tcharge);
			txyz.push_back(tx);
			txyz.push_back(ty);
			txyz.push_back(tz);
			this->xyz.push_back(txyz);
			txyz.clear();

			this->charges.push_back(tcharge);
			this->atomnames.push_back(str);

			if (Input->mol2_aa){
				this->amberatoms.push_back(tatomtype);
                atom_param* at = new atom_param;
                this->get_gaff_atomic_parameters(string(tatomtype), at);
                this->radii.push_back(at->radius);
                this->epsilons.push_back(at->epsilon);
                this->epsilons_sqrt.push_back(sqrt(at->epsilon));
                this->masses.push_back(at->mass);
                delete at;
			}
            else{
                if (Input->atomic_model_ff == "AMBER" or Input->atomic_model_ff == "amber"){
                    this->amberatoms.push_back(this->sybyl_2_amber(string(tatomtype)));
                }
                else {
                    this->amberatoms.push_back(this->sybyl_2_gaff(string(tatomtype)));
                }
                atom_param* at = new atom_param;
                this->get_gaff_atomic_parameters(this->amberatoms[i], at);
                this->radii.push_back(at->radius);
                this->epsilons.push_back(at->epsilon);
                this->epsilons_sqrt.push_back(sqrt(at->epsilon));
                this->masses.push_back(at->mass);
                delete at;
                this->sybyl_atoms.push_back(string(tatomtype));
            }

			if (tres > count){
				this->residue_pointer.push_back(i+1);
				count = tres;
				this->resnames.push_back(string(resname));
			}
		}

//		fscanf(mol2file, "%s\n", str);
        cpstr = string(str);
        while (cpstr.substr(0,13) != "@<TRIPOS>BOND"){
            fgets(str, 80, mol2file);
            cpstr = string(str);
        }

		vector<string> bond;
		char s1[6], s2[6], s3[5];
		for (int i=0; i<this->Nbonds; i++){
			fscanf(mol2file, "%d%s%s%s\n", &tint, s1, s2, s3);
			bond.push_back(string(s1));
			bond.push_back(string(s2));
			bond.push_back(string(s3));
			this->bonds.push_back(bond);
			bond.clear();
		}
		bret = true;
	}
	else {
		bret = false;
		printf("Skipping file %s\n", molfile.c_str());
	}
	fclose(mol2file);
	return (bret);
}

bool Mol2::parse_gzipped_file(PARSER* Input, string molfile){
    bool bret = false;
    int tint;
    float tx, ty, tz;
    vector<double> txyz;
    int tres;
    float tcharge;
    int count=0;
    char tatomtype[10];
    char resname[20];
    string cpstr;

/*
 * Here we read the GAFF/AMBER parameters fro LJ potentials
 */

    this->initialize_gaff();

    gzFile mol2file = gzopen(molfile.c_str(), "r");
    if (mol2file != NULL){
        str[0]='#';
        while(str[0] !='@'){
            gzgets(mol2file, str, 80);

        }
        gzgets(mol2file, str, 100);
        this->molname = str;
        this->molname = this->molname.substr(0,this->molname.size()-1);
        gzgets(mol2file, str, 100);
        sscanf(str, "%d %d %d %d %d", &this->N, &this->Nbonds, &this->Nres, &tint, &tint);

        cpstr = string(str);
        while (cpstr.substr(0,13) != "@<TRIPOS>ATOM"){
            gzgets(mol2file, str, 100);
            cpstr = string(str);
        }

        for (int i=0; i<this->N; i++){
            gzgets(mol2file, str, 100);
            sscanf(str, "%d %s %f %f %f %s %d %s %f\n", &tint, str, &tx, &ty, &tz, tatomtype, &tres, resname, &tcharge);
            txyz.push_back(tx);
            txyz.push_back(ty);
            txyz.push_back(tz);
            this->xyz.push_back(txyz);
            txyz.clear();

            this->charges.push_back(tcharge);
            this->atomnames.push_back(str);

            if (Input->mol2_aa){
                this->amberatoms.push_back(tatomtype);
                atom_param* at = new atom_param;
                this->get_gaff_atomic_parameters(string(tatomtype), at);
                this->radii.push_back(at->radius);
                this->epsilons.push_back(at->epsilon);
                this->epsilons_sqrt.push_back(sqrt(at->epsilon));
                this->masses.push_back(at->mass);
                delete at;
            }
            else{
                if (Input->atomic_model_ff == "AMBER" or Input->atomic_model_ff == "amber"){
                    this->amberatoms.push_back(this->sybyl_2_amber(string(tatomtype)));
                }
                else {
                    this->amberatoms.push_back(this->sybyl_2_gaff(string(tatomtype)));
                }
                atom_param* at = new atom_param;
                this->get_gaff_atomic_parameters(this->amberatoms[i], at);
                this->radii.push_back(at->radius);
                this->epsilons.push_back(at->epsilon);
                this->epsilons_sqrt.push_back(sqrt(at->epsilon));
                this->masses.push_back(at->mass);
                delete at;
                this->sybyl_atoms.push_back(string(tatomtype));
            }

            if (tres > count){
                this->residue_pointer.push_back(i+1);
                count = tres;
                this->resnames.push_back(string(resname));
            }
        }

//        gzgets(mol2file, str, 100);
        cpstr = string(str);
        while (cpstr.substr(0,13) != "@<TRIPOS>BOND"){
            gzgets(mol2file, str, 100);
            cpstr = string(str);
        }

        vector<string> bond;
        char s1[6], s2[6], s3[5];
        for (int i=0; i<this->Nbonds; i++){
            gzgets(mol2file, str, 80);
            sscanf(str, "%d%s%s%s\n", &tint, s1, s2, s3);
            bond.push_back(string(s1));
            bond.push_back(string(s2));
            bond.push_back(string(s3));
            this->bonds.push_back(bond);
            bond.clear();
        }

// Checking for internal consistency ....

        if ((int(this->radii.size()) == this->N) and (int(this->epsilons.size()) == this->N)){
            bret = true;
        }
    }
    else {
        printf("Skipping file %s...\n", molfile.c_str());
        bret = false;
    }
    gzclose(mol2file);

    return (bret);
}


Mol2::~Mol2(){
	this->xyz.clear();
	this->charges.clear();
	this->radii.clear();
	this->epsilons.clear();
	this->epsilons_sqrt.clear();
	this->resnames.clear();
	this->bonds.clear();
	this->sybyl_atoms.clear();
//	this->mcoords.clear();
//	this->new_mcoords.clear();
	this->new_xyz.clear();
	this->atomnames.clear();
	this->masses.clear();
	this->amberatoms.clear();
}

void Mol2::initialize_gaff(){
    FILE *gaff_file;
    char str[80];
    char at[3];
    float r, e, m;
    char filename[150];

    char* dir_path = getenv("LIBELA");
    if (dir_path== NULL){
/*
        printf("Environment variable LIBELA is not set.\n");
        printf("Trying to use local folder as LIBELA folder...\n");
*/
        strcpy(filename, "");
        strcat(filename, "param/LJ_parm.dat");
    }
    else {
        strcpy(filename, dir_path);
        strcat(filename, "/param/LJ_parm.dat");
    }

    gaff_file = fopen(filename, "r");
    if (gaff_file!= NULL){
        while (!feof(gaff_file)){
            fgets(str, 80, gaff_file);
            if (str[0] != '#'){
                sscanf(str, "%s %f %f %f", at, &r, &e, &m);
                atom_param v;
                v.type = string(at);
                v.radius = double(r);
                v.epsilon = double(e);
                v.mass = double(m);
                this->gaff_force_field.push_back(v);
            }
        }
    }
    else{
        printf("Error reading file %s\n", filename);
        exit(1);
    }

    fclose(gaff_file);
}


void Mol2::get_gaff_atomic_parameters(string gaff_atom, atom_param* ap){
    bool found=false;
    for (unsigned i=0; i<this->gaff_force_field.size(); i++){
        if (this->gaff_force_field[i].type == gaff_atom){
            ap->type = this->gaff_force_field[i].type;
            ap->epsilon= this->gaff_force_field[i].epsilon;
            ap->radius = this->gaff_force_field[i].radius;
            ap->mass= this->gaff_force_field[i].mass;
            found = true;
        }
    }
    if (!found){
        printf("Could not find atomic parameters for atom %s. Please, check.\n", gaff_atom.c_str());
        exit(1);
    }
}

string Mol2::sybyl_2_gaff(string atom){
    string gaff_atom;
    if (atom == "C.3"){
        gaff_atom = "c3";
    }
    else if (atom =="C.2"){
        gaff_atom = "c2";
    }
    else if (atom =="C.1"){
        gaff_atom = "c1";
    }

    else if (atom =="C.ar"){
        gaff_atom = "ca";
    }

    else if (atom =="C.cat"){
        gaff_atom = "c";
    }

    else if (atom =="N.3"){
        gaff_atom = "n3";
    }

    else if (atom =="N.2"){
        gaff_atom = "n2";
    }

    else if (atom =="N.1"){
        gaff_atom = "n1";
    }

    else if (atom =="N.ar"){
        gaff_atom = "nh";
    }

    else if (atom =="N.am"){
        gaff_atom = "n";
    }

    else if (atom =="N.4"){
        gaff_atom = "n4";
    }

    else if (atom =="N.pl3"){
        gaff_atom = "na";
    }

    else if (atom =="N.p"){
        gaff_atom = "na";
    }

    else if (atom =="O.3"){
        gaff_atom = "oh";
    }

    else if (atom =="O.2"){
        gaff_atom = "o";
    }

    else if (atom =="O.co2"){
        gaff_atom = "o";
    }

    else if (atom =="O.spc" or atom == "O.t3p"){
        gaff_atom = "ow";
    }

    else if (atom =="S.3"){
        gaff_atom = "sh"; //not sure... sh or ss
    }

    else if (atom =="S.2"){
        gaff_atom = "s2";
    }

    else if (atom =="S.O" or atom == "S.o"){
        gaff_atom = "s4";
    }

    else if (atom =="S.O2" or atom == "S.o2"){
        gaff_atom = "s6";
    }

    else if (atom =="P.3"){
        gaff_atom = "p3";
    }

    else if (atom =="F"){
        gaff_atom = "f";
    }

    else if (atom =="H"){
        gaff_atom = "hc";
    }

    else if (atom =="H.spc" or atom=="H.t3p"){
        gaff_atom = "hw";
    }

    else if (atom =="Cl"){
        gaff_atom = "cl";
    }

    else if (atom =="Br"){
        gaff_atom = "br";
    }

    else if (atom =="I"){
        gaff_atom = "i";
    }

    else if (atom =="Mg"){
        gaff_atom = "MG";
    }

    else if (atom =="LP" or atom == "Lp"){
        gaff_atom = "EP";
    }

    else if (atom == "Fe"){
        gaff_atom = "Fe";
    }
    else if (atom == "Zn"){
        gaff_atom = "Zn";
    }
    else if (atom == "Cu"){
        gaff_atom = "Cu";
    }
    else if (atom == "Ca"){
        gaff_atom = "Ca";
    }
    else if (atom == "Si"){
        gaff_atom = "Si";
    }
    else{
        bool found = false;
        for (unsigned i=0; i< this->gaff_force_field.size(); i++){
            if (this->gaff_force_field[i].type == atom){
                gaff_atom = atom;
                found = true;
            }
        }
        if (! found){
            printf("Atom type %s not found among GAFF parameters.\nPlease check Mol2.h source file.\n", atom.c_str());
            exit(1);
        }
    }


    return(gaff_atom);
}

string Mol2::sybyl_2_amber(string atom){
    string amber_atom;
    if (atom == "C.3"){
        amber_atom = "CT";
    }
    else if (atom =="C.2"){
        amber_atom = "C*";
    }
    else if (atom =="C.1"){
        amber_atom = "C*";
    }

    else if (atom =="C.ar"){
        amber_atom = "C*";
    }

    else if (atom =="C.cat"){
        amber_atom = "C*";
    }

    else if (atom =="N.3"){
        amber_atom = "N3";
    }

    else if (atom =="N.2"){
        amber_atom = "N";
    }

    else if (atom =="N.1"){
        amber_atom = "N";
    }

    else if (atom =="N.ar"){
        amber_atom = "N";
    }

    else if (atom =="N.am"){
        amber_atom = "N";
    }

    else if (atom =="N.4"){
        amber_atom = "n4";
    }

    else if (atom =="N.pl3"){
        amber_atom = "N";
    }

    else if (atom =="N.p"){
        amber_atom = "N";
    }

    else if (atom =="O.3"){
        amber_atom = "OH";
    }

    else if (atom =="O.2"){
        amber_atom = "O2";
    }

    else if (atom =="O.co2"){
        amber_atom = "O";
    }

    else if (atom =="O.spc" or atom == "O.t3p"){
        amber_atom = "OW";
    }

    else if (atom =="S.3"){
        amber_atom = "S"; //not sure... sh or ss
    }

    else if (atom =="S.2"){
        amber_atom = "S";
    }

    else if (atom =="S.O" or atom == "S.o"){
        amber_atom = "S";
    }

    else if (atom =="S.O2" or atom == "S.o2"){
        amber_atom = "S";
    }

    else if (atom =="P.3"){
        amber_atom = "P";
    }

    else if (atom =="F"){
        amber_atom = "F";
    }

    else if (atom =="H"){
        amber_atom = "H";
    }

    else if (atom =="H.spc" or atom=="H.t3p"){
        amber_atom = "HW";
    }

    else if (atom =="Cl"){
        amber_atom = "Cl";
    }

    else if (atom =="Br"){
        amber_atom = "Br";
    }

    else if (atom =="I"){
        amber_atom = "I";
    }

    else if (atom =="Mg"){
        amber_atom = "MG";
    }

    else if (atom =="LP" or atom == "Lp"){
        amber_atom = "EP";
    }

    else if (atom == "Fe"){
        amber_atom = "Fe";
    }
    else if (atom == "Zn"){
        amber_atom = "Zn";
    }
    else if (atom == "Cu"){
        amber_atom = "Cu";
    }
    else if (atom == "Ca"){
        amber_atom = "Ca";
    }
    else if (atom == "Si"){
        amber_atom = "Si";
    }
    else{
        bool found = false;
        for (unsigned i=0; i< this->gaff_force_field.size(); i++){
            if (this->gaff_force_field[i].type == atom){
                amber_atom = atom;
                found = true;
            }
        }
        if (! found){
            printf("Atom type %s not found among GAFF parameters.\nPlease check Mol2.h source file.\n", atom.c_str());
            exit(1);
        }
    }
    return(amber_atom);
}

bool Mol2::parse_gzipped_ensemble(PARSER* Input, string molfile, int skipper=1){
    char tstr[80];
    bool bret = false;
    int tint;
    float tx, ty, tz;
    vector<double> txyz;
    int tres;
    float tcharge;
    int count=0;
    char tatomtype[10];
    char resname[20];
    string cpstr;
    vector<vector<double> > tcoord;
    int trajsize=0;

    this->initialize_gaff();

    gzFile mol2file = gzopen(molfile.c_str(), "r");
    if (mol2file != NULL){
        str[0]='#';
        while(str[0] !='@'){
            gzgets(mol2file, str, 100);

        }
        gzgets(mol2file, str, 100);
        this->molname = str;
        this->molname = this->molname.substr(0,this->molname.size()-1);
        gzgets(mol2file, str, 100);
        sscanf(str, "%d %d %d %d %d", &this->N, &this->Nbonds, &this->Nres, &tint, &tint);


        cpstr = string(str);
        while (cpstr.substr(0,6) != "Energy"){
            gzgets(mol2file, str, 100);
            cpstr = string(str);
        }
        sscanf(str, "%s %f\n", tstr, &tx);              //parsing ensemble energy
        this->ensemble_energies.push_back(double(tx));
        trajsize++;

        while (cpstr.substr(0,13) != "RMSD/OVERLAY:"){
            gzgets(mol2file, str, 100);
            cpstr = string(str);
        }
        sscanf(str, "%s %f\n", tstr, &tx);              //parsing ensemble rmsd
        this->ensemble_rmsd.push_back(double(tx));

        cpstr = string(str);
        while (cpstr.substr(0,13) != "@<TRIPOS>ATOM"){
            gzgets(mol2file, str, 100);
            cpstr = string(str);
        }

        for (int i=0; i<this->N; i++){
            gzgets(mol2file, str, 100);
            sscanf(str, "%d %s %f %f %f %s %d %s %f\n", &tint, str, &tx, &ty, &tz, tatomtype, &tres, resname, &tcharge);
            txyz.push_back(tx);
            txyz.push_back(ty);
            txyz.push_back(tz);
            this->xyz.push_back(txyz);
            txyz.clear();

            this->charges.push_back(tcharge);
            this->atomnames.push_back(str);

            if (Input->mol2_aa){
                this->amberatoms.push_back(tatomtype);
                atom_param* at = new atom_param;
                this->get_gaff_atomic_parameters(string(tatomtype), at);
                this->radii.push_back(at->radius);
                this->epsilons.push_back(at->epsilon);
                this->epsilons_sqrt.push_back(sqrt(at->epsilon));
                this->masses.push_back(at->mass);
                delete at;
            }
            else{
                if (Input->atomic_model_ff == "AMBER" or Input->atomic_model_ff == "amber"){
                    this->amberatoms.push_back(this->sybyl_2_amber(string(tatomtype)));
                }
                else {
                    this->amberatoms.push_back(this->sybyl_2_gaff(string(tatomtype)));
                }
                atom_param* at = new atom_param;
                this->get_gaff_atomic_parameters(this->amberatoms[i], at);
                this->radii.push_back(at->radius);
                this->epsilons.push_back(at->epsilon);
                this->epsilons_sqrt.push_back(sqrt(at->epsilon));
                this->masses.push_back(at->mass);
                delete at;
                this->sybyl_atoms.push_back(string(tatomtype));
            }

            if (tres > count){
                this->residue_pointer.push_back(i+1);
                count = tres;
                this->resnames.push_back(string(resname));
            }
        }

        gzgets(mol2file, str, 100);
        if (str[0] != '@'){
            while (str[0] != '@'){
                gzgets(mol2file, str, 100);
            }
        }

        vector<string> bond;
        char s1[6], s2[6], s3[5];
        for (int i=0; i<this->Nbonds; i++){
            gzgets(mol2file, str, 100);
            sscanf(str, "%d%s%s%s\n", &tint, s1, s2, s3);
            bond.push_back(string(s1));
            bond.push_back(string(s2));
            bond.push_back(string(s3));
            this->bonds.push_back(bond);
            bond.clear();
        }

        int n=0;

        this->mcoords.push_back(this->xyz);

        while (! gzeof(mol2file)){
            gzgets(mol2file, str, 100);
            while ((str[0] != 'E' or str[6] != ':') and (!gzeof(mol2file))){
                gzgets(mol2file, str, 100);
            }

            if (!gzeof(mol2file)){
                sscanf(str, "%s %f\n", tstr, &tx);
                trajsize++;
                if (trajsize % skipper == 0){
                    this->ensemble_energies.push_back(double(tx));
                }
            }

            while ((str[0] != 'R' or str[12] != ':') and (!gzeof(mol2file))){
                gzgets(mol2file, str, 100);
            }

            if (!gzeof(mol2file)){
                sscanf(str, "%s %f\n", tstr, &tx);
                if (trajsize % skipper == 0){
                    this->ensemble_rmsd.push_back(double(tx));
                }
            }

            while ((str[0] != '@' or str[9] != 'A') and (!gzeof(mol2file))){
                gzgets(mol2file, str, 100);
            }
            if (!gzeof(mol2file) and (trajsize % skipper == 0)){
                txyz.clear();
                for (int i=0; i<this->N; i++){
                    gzgets(mol2file, str, 100);
                    sscanf(str, "%d %s %f %f %f %s %d %s %f\n", &tint, tstr, &tx, &ty, &tz, tatomtype, &tres, resname, &tcharge);
                    txyz.push_back(tx);
                    txyz.push_back(ty);
                    txyz.push_back(tz);
                    tcoord.push_back(txyz);
                    txyz.clear();
                }
                this->mcoords.push_back(tcoord);
                n++;
                tcoord.clear();
            }
        }
#ifdef DEBUG
        printf("Found %d conformations in file %s\n", n, molfile.c_str());
#endif
        bret = true;
    }

    else {
        printf("Skipping file %s...\n", molfile.c_str());
        bret = false;
    }
    gzclose(mol2file);
    return (bret);
}

vector<vector<double> > Mol2::get_next_xyz(PARSER* Input, gzFile mol2file) {
    char tstr[80];
    int tint;
    float tx, ty, tz;
    vector<double> txyz(3);
    int tres;
    float tcharge;
    char tatomtype[10];
    char resname[20];
    string cpstr;
    vector<vector<double> > tcoord(unsigned(this->N));
    char str[100];                          // making it local

    for (unsigned i=0; i<10; i++){
        str[i] = '#';
    }

    if (mol2file != NULL){
        while ((str[0] != 'E' or str[6] != ':') and (!gzeof(mol2file))){
            gzgets(mol2file, str, 100);
        }

        if (!gzeof(mol2file)){
            sscanf(str, "%s %f\n", tstr, &tx);
            this->energy = double(tx);
        }

        while ((str[0] != 'R' or str[12] != ':') and (!gzeof(mol2file))){
            gzgets(mol2file, str, 100);
        }

        if (!gzeof(mol2file)){
            sscanf(str, "%s %f\n", tstr, &tx);
            this->rmsd = double(tx);
        }

        while ((str[0] != '@' or str[9] != 'A') and (!gzeof(mol2file))){
            gzgets(mol2file, str, 100);
        }
        if (!gzeof(mol2file)){
            for (int i=0; i<this->N; i++){
                gzgets(mol2file, str, 100);
                sscanf(str, "%d %s %f %f %f %s %d %s %f\n", &tint, tstr, &tx, &ty, &tz, tatomtype, &tres, resname, &tcharge);
                txyz[0] = double(tx);
                txyz[1] = double(ty);
                txyz[2] = double(tz);
                tcoord[i] = txyz;
            }
        }
    }
    else {
        printf("Skipping file ...\n");
    }
    return (tcoord);
}

void Mol2::find_longest_axis(){
    vector<int> axis(2);
    double d, dist=0.0;
    for (int i=0; i<this->N-1; i++){
        for (int j=i+1; j<this->N; j++){
            d = this->distance(this->xyz[i], this->xyz[j]);
            if (d > dist){
                dist = d;
                axis[0] = i;
                axis[1] = j;
            }
        }
    }
    this->longest_axis = axis;
    this->radius = dist/2.0;
}

double Mol2::distance(vector<double> atom1, vector<double> atom2) {
    return ( sqrt(((atom2[0]-atom1[0])*(atom2[0]-atom1[0]))+((atom2[1]-atom1[1])*(atom2[1]-atom1[1]))+((atom2[2]-atom1[2])*(atom2[2]-atom1[2]))) );
}


#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(pyMol2)
{
    class_<Mol2>("Mol2", init< >())
        .def_readwrite("N", & Mol2::N)
        .def_readwrite("Nres", & Mol2::Nres)
        .def_readwrite("Natomtypes", & Mol2::Natomtypes)
        .def_readwrite("Nbonds", & Mol2::Nbonds)
        .def_readwrite("molname", & Mol2::molname)
        .def_readwrite("charges", & Mol2::charges)
        .def_readwrite("masses", & Mol2::masses)
        .def_readwrite("amberatoms", & Mol2::amberatoms)
        .def_readwrite("atomtypes_prm", & Mol2::atomtypes_prm)
        .def_readwrite("atomnames", & Mol2::atomnames)
        .def_readwrite("xyz", & Mol2::xyz)
        .def_readwrite("new_xyz", & Mol2::new_xyz)
        .def_readwrite("opt_overlay_xyz", & Mol2::opt_overlay_xyz)
        .def_readwrite("opt_energy_xyz", & Mol2::opt_energy_xyz)
        .def_readwrite("epsilons", & Mol2::epsilons)
        .def_readwrite("epsilons", & Mol2::epsilons_sqrt)
        .def_readwrite("radii", & Mol2::radii)
        .def_readwrite("resnames", & Mol2::resnames)
        .def_readwrite("residue_pointer", & Mol2::residue_pointer)
        .def_readwrite("HBacceptors", & Mol2::HBacceptors)
        .def_readwrite("HBdonors", & Mol2::HBdonors)
        .def_readwrite("line", & Mol2::line)
        .def_readwrite("str", & Mol2::str)
        .def_readwrite("self_obj_function", & Mol2::self_obj_function)
        .def_readwrite("bonds", & Mol2::bonds)
        .def_readwrite("sybyl_atoms", & Mol2::sybyl_atoms)
        .def_readwrite("atm", & Mol2::atm)
        .def_readwrite("conformer_energies", & Mol2::conformer_energies)
        .def_readwrite("energy", & Mol2::energy)
        .def_readwrite("rmsd", & Mol2::rmsd)
        .def_readwrite("longest_axis", & Mol2::longest_axis)
        .def_readwrite("radius", & Mol2::radius)

        .def(init<PARSER*, string>())
        .def("parse_gzipped_file", & Mol2::parse_gzipped_file)
        .def("parse_mol2file", & Mol2::parse_mol2file)
        .def("get_next_xyz", & Mol2::get_next_xyz)
        .def("initialize_gaff", & Mol2::initialize_gaff)
        .def("get_gaff_atomic_parameters", & Mol2::get_gaff_atomic_parameters)
        .def("sybyl_2_gaff", & Mol2::sybyl_2_gaff)
        .def("sybyl_2_amber", & Mol2::sybyl_2_amber)
        .def("find_longest_axis", & Mol2::find_longest_axis)
        .def("distance", & Mol2::distance)
    ;
    class_<atom_param>("atom_param")
        .def_readwrite("type", & Mol2::atom_param::type)
        .def_readwrite("radius", & Mol2::atom_param::radius)
        .def_readwrite("epsilon", & Mol2::atom_param::epsilon)
        .def_readwrite("mass", & Mol2::atom_param::mass)
    ;
}
