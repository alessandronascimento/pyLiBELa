/*
 * Mol2.cpp
 *
 *  Created on: 10/10/2011
 *      Author: Nascimento
 */

#include "pyMol2.h"

using namespace std;

Mol2::Mol2(){
    str = new char[100];
}

Mol2::Mol2(PARSER *Input, string molfile) {
    bool ok;
    str = new char[100];
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

bool Mol2::parse_smiles(PARSER *Input, string smiles_input, string molname){
    bool bret = true;

/*
 * Reading SMILES
 */

    OBMol mol;
    OBConversion conv;
    if (!(conv.SetInFormat("smi") && conv.ReadString(&mol, smiles_input))){
        printf("Skipping smiles %s...\n", smiles_input.c_str());
        bret=false;
    }

/*
 * Create 3D mol from SMILES
 */
    if (mol.NumAtoms() < Input->atom_limit){
        OBBuilder builder;
        builder.Build(mol);
        mol.AddHydrogens(false, true); // adding H atoms: false= not for polar atoms only. true=correct for pH 7.4.

/*
 * Minimize energy
 */

        OBForceField* OBff = OBForceField::FindForceField("MMFF94");
        OBff->Setup(mol);
        OBff->SteepestDescent(100);
        OBff->UpdateCoordinates(mol);

 /*
 * Get data for Mol2 class
 *
 */

        char aname[7];
        OBAtom *atom;
        vector<int> vtemp(2);
        string sybyl_atom;
        string sybyl_gaff;

        this->N = mol.NumAtoms();
        this->molname = molname;
        this->Nres = 0;
        this->residue_pointer.push_back(1);

        this->initialize_gaff2();

        FOR_ATOMS_OF_MOL(atom, mol){
            sprintf(aname, "%s%d", OBElements::GetSymbol(atom->GetAtomicNum()), atom->GetIdx());
            this->atomnames.push_back(string(aname));
            this->charges.push_back(double(atom->GetPartialCharge()));
            if (atom->IsHbondAcceptor()){
                this->HBacceptors.push_back(int(atom->GetIdx()-1));
            }
            else if (atom->IsHbondDonorH()){
                FOR_NBORS_OF_ATOM(nbr, &*atom){
                    vtemp[0] = int(nbr->GetIdx()-1);
                    vtemp[1] = int(atom->GetIdx()-1);
                    this->HBdonors.push_back(vtemp);
                }
            }

            sybyl_atom = this->sybyl_2_gaff(string(atom->GetType()));
            this->amberatoms.push_back(sybyl_atom);
            if (sybyl_atom == ""){
                bret=false;
                printf("bret set to false because of atom %s\n", atom->GetType());
            }
            else{
                this->sybyl_atoms.push_back(this->gaff_2_sybyl(sybyl_atom));
                sybyl_gaff = this->gaff_2_sybyl(sybyl_atom);
                if (sybyl_gaff == ""){
                    bret=false;
                    printf("bret set to false because of atom %s\n", sybyl_atom);
                }
                atom_param* at = new atom_param;
                this->get_gaff_atomic_parameters(sybyl_atom, at);
                this->radii.push_back(at->radius);
                this->epsilons.push_back(at->epsilon);
                this->epsilons_sqrt.push_back(sqrt(at->epsilon));
                this->masses.push_back(at->mass);
                delete at;
            }
        }
        int b1, b2, b3;
        char tmp[10];
        vector<string> bond;
        FOR_BONDS_OF_MOL(b, mol){
            b1 = b->GetBeginAtomIdx();
            sprintf(tmp, "%d", b1);
            bond.push_back(string(tmp));
            b2 = b->GetEndAtomIdx();
            sprintf(tmp, "%d", b2);
            bond.push_back(string(tmp));
            b3 = b->GetBondOrder();
            sprintf(tmp, "%d", b3);
            bond.push_back(string(tmp));
            this->bonds.push_back(bond);
            bond.clear();
        }

        this->xyz = this->copy_from_obmol(mol);
        this->resnames.push_back("LIG");

        if (bret and Input->generate_conformers){
            OBMol RefMol = mol;
            OBff->DiverseConfGen(0.5, Input->conf_search_trials, 50.0, false);
            OBff->GetConformers(mol);

            int generated_conformers=0;
            if (mol.NumConformers() > Input->lig_conformers){
                generated_conformers = Input->lig_conformers;
            }
            else {
                generated_conformers = mol.NumConformers();
            }

            if (mol.NumConformers() > 0){
                for (int i=0; i<generated_conformers; i++){
                    double x[mol.NumAtoms()*3];
                    double* xyz;
                    xyz = x;
                    vector<double> v3;
                    vector<vector<double> > xyz_tmp;
                    mol.SetConformer(i);
                    OBff->Setup(mol);
                    OBff->GetCoordinates(mol);
                    energy = OBff->Energy();
                    if (OBff->GetUnit() == "kJ/mol"){       // Converting to kcal/mol, if needed.
                        energy = energy/4.18;
                    }
                    this->conformer_energies.push_back(energy);

                    OBAlign* align = new OBAlign;
                    align->SetRefMol(RefMol);
                    align->SetTargetMol(mol);
                    align->Align();
                    align->UpdateCoords(&mol);
                    delete align;

                    xyz = mol.GetCoordinates();
                    for (unsigned j=0; j<mol.NumAtoms(); j++){
                        v3.push_back(xyz[3*j]);
                        v3.push_back(xyz[(3*j)+1]);
                        v3.push_back(xyz[(3*j)+2]);
                        xyz_tmp.push_back(v3);
                        v3.clear();
                    }
                    this->mcoords.push_back(xyz_tmp);
                    xyz_tmp.clear();
                }
            }
        }
    }
    else{
        bret = false;
    }
    return bret;
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
    string atom_type;
    bool missing_atom=false;

    this->initialize_gaff2();


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
                    atom_type=this->sybyl_2_amber(string(tatomtype));
                    //                    this->amberatoms.push_back();
                }
                else {
                    atom_type = this->sybyl_2_gaff(string(tatomtype));
                }
                if (atom_type == ""){
                    missing_atom = true;
                    printf("Skipping molecule due to missing atom type: %s.\n", tatomtype);
                }
                else{
                    this->amberatoms.push_back(atom_type);
                    atom_param* at = new atom_param;
                    this->get_gaff_atomic_parameters(this->amberatoms[i], at);
                    this->radii.push_back(at->radius);
                    this->epsilons.push_back(at->epsilon);
                    this->epsilons_sqrt.push_back(sqrt(at->epsilon));
                    this->masses.push_back(at->mass);
                    delete at;
                    this->sybyl_atoms.push_back(string(tatomtype));
                }
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
    return ((bret and (!missing_atom)));
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
    bool missing_atom=false;
    string atom_type;

    /*
 * Here we read the GAFF/AMBER parameters fro LJ potentials
 */

    this->initialize_gaff2();

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
                    atom_type = this->sybyl_2_amber(string(tatomtype));
                }
                else {
                    atom_type = this->sybyl_2_gaff(string(tatomtype));
                }
                if (atom_type == ""){
                    missing_atom = true;
                    printf("Skipping molecule due to missing atom type: %s.\n", tatomtype);
                }
                else {
                    this->amberatoms.push_back(atom_type);
                    atom_param* at = new atom_param;
                    this->get_gaff_atomic_parameters(this->amberatoms[i], at);
                    this->radii.push_back(at->radius);
                    this->epsilons.push_back(at->epsilon);
                    this->epsilons_sqrt.push_back(sqrt(at->epsilon));
                    this->masses.push_back(at->mass);
                    delete at;
                    this->sybyl_atoms.push_back(string(tatomtype));
                }
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

    return (bret and !missing_atom);
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

void Mol2::initialize_gaff2(){
    atom_param ap;
    vector<atom_param> gaff_parameters;

    ap.type = "hc";
    ap.radius = 0.6000;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "ha";
    ap.radius = 1.4735;
    ap.epsilon = 0.0161;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "hn";
    ap.radius = 0.6210;
    ap.epsilon = 0.0100;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "ho";
    ap.radius = 0.3019;
    ap.epsilon = 0.0047;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "hs";
    ap.radius = 0.6112;
    ap.epsilon = 0.0124;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "hp";
    ap.radius = 0.6031;
    ap.epsilon = 0.0144;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "o";
    ap.radius = 1.7107;
    ap.epsilon = 0.1463;
    ap.mass = 16.00;
    gaff_parameters.push_back(ap);

    ap.type = "os";
    ap.radius = 1.7713;
    ap.epsilon = 0.0726;
    ap.mass = 16.00;
    gaff_parameters.push_back(ap);

    ap.type = "oh";
    ap.radius = 1.8200;
    ap.epsilon = 0.0930;
    ap.mass = 16.00;
    gaff_parameters.push_back(ap);

    ap.type = "ow";
    ap.radius = 1.7683;
    ap.epsilon = 0.1520;
    ap.mass = 16.00;
    gaff_parameters.push_back(ap);

    ap.type = "c3";
    ap.radius = 1.9069;
    ap.epsilon = 0.1078;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "c2";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "c1";
    ap.radius = 1.9525;
    ap.epsilon = 0.1596;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "n";
    ap.radius = 1.7852;
    ap.epsilon = 0.1636;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "s";
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.push_back(ap);

    ap.type = "p2";
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.push_back(ap);

    ap.type = "f";
    ap.radius = 1.7029;
    ap.epsilon = 0.0832;
    ap.mass = 19.00;
    gaff_parameters.push_back(ap);

    ap.type = "cl";
    ap.radius = 1.9452;
    ap.epsilon = 0.2638;
    ap.mass = 35.45;
    gaff_parameters.push_back(ap);

    ap.type = "br";
    ap.radius = 2.0275;
    ap.epsilon = 0.3932;
    ap.mass = 79.90;
    gaff_parameters.push_back(ap);

    ap.type = "i";
    ap.radius = 2.1558;
    ap.epsilon = 0.4955;
    ap.mass = 126.9;
    gaff_parameters.push_back(ap);

    ap.type = "n1";
    ap.radius = 1.8372;
    ap.epsilon = 0.1098;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "n2";
    ap.radius = 1.8993;
    ap.epsilon = 0.0941;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "n3";
    ap.radius = 1.8886;
    ap.epsilon = 0.0858;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "na";
    ap.radius = 1.7992;
    ap.epsilon = 0.2042;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "nh";
    ap.radius = 1.7903;
    ap.epsilon = 0.2150;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "n+";
    ap.radius = 1.6028;
    ap.epsilon = 0.7828;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "n9";
    ap.radius = 2.2700;
    ap.epsilon = 0.0095;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "h1";
    ap.radius = 1.3593;
    ap.epsilon = 0.0208;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "h2";
    ap.radius = 1.2593;
    ap.epsilon = 0.0208;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "h3";
    ap.radius = 1.1593;
    ap.epsilon = 0.0208;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "hx";
    ap.radius = 1.0593;
    ap.epsilon = 0.0208;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "h4";
    ap.radius = 1.4235;
    ap.epsilon = 0.0161;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "h5";
    ap.radius = 1.3735;
    ap.epsilon = 0.0161;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "cx";
    ap.radius = 1.9069;
    ap.epsilon = 0.1078;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "cy";
    ap.radius = 1.9069;
    ap.epsilon = 0.1078;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "c";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "cs";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "ca";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "cc";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "cd";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "ce";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "cf";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "cp";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "cq";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "cz";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "cg";
    ap.radius = 1.9525;
    ap.epsilon = 0.1596;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "ch";
    ap.radius = 1.9525;
    ap.epsilon = 0.1596;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "cu";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "cv";
    ap.radius = 1.8606;
    ap.epsilon = 0.0988;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "nb";
    ap.radius = 1.8993;
    ap.epsilon = 0.0941;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "nc";
    ap.radius = 1.8993;
    ap.epsilon = 0.0941;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "nd";
    ap.radius = 1.8993;
    ap.epsilon = 0.0941;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "ne";
    ap.radius = 1.8993;
    ap.epsilon = 0.0941;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "nf";
    ap.radius = 1.8993;
    ap.epsilon = 0.0941;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "no";
    ap.radius = 1.8886;
    ap.epsilon = 0.0858;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "n7";
    ap.radius = 1.9686;
    ap.epsilon = 0.0522;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "n8";
    ap.radius = 2.0486;
    ap.epsilon = 0.0323;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "n4";
    ap.radius = 1.4028;
    ap.epsilon = 3.8748;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "nx";
    ap.radius = 1.4528;
    ap.epsilon = 2.5453;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "ny";
    ap.radius = 1.5028;
    ap.epsilon = 1.6959;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "nz";
    ap.radius = 1.5528;
    ap.epsilon = 1.1450;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "ns";
    ap.radius = 1.8352;
    ap.epsilon = 0.1174;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "nt";
    ap.radius = 1.8852;
    ap.epsilon = 0.0851;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "nu";
    ap.radius = 1.8403;
    ap.epsilon = 0.1545;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "nv";
    ap.radius = 1.8903;
    ap.epsilon = 0.1120;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "s2";
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.push_back(ap);

    ap.type = "s4";
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.push_back(ap);

    ap.type = "s6";
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.push_back(ap);

    ap.type = "sx";
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.push_back(ap);

    ap.type = "sy";
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.push_back(ap);

    ap.type = "sh";
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.push_back(ap);

    ap.type = "ss";
    ap.radius = 1.9825;
    ap.epsilon = 0.2824;
    ap.mass = 32.06;
    gaff_parameters.push_back(ap);

    ap.type = "p3";
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.push_back(ap);

    ap.type = "p4";
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.push_back(ap);

    ap.type = "p5";
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.push_back(ap);

    ap.type = "pb";
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.push_back(ap);

    ap.type = "px";
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.push_back(ap);

    ap.type = "py";
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.push_back(ap);

    ap.type = "pc";
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.push_back(ap);

    ap.type = "pd";
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.push_back(ap);

    ap.type = "pe";
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.push_back(ap);

    ap.type = "pf";
    ap.radius = 2.0732;
    ap.epsilon = 0.2295;
    ap.mass = 30.97;
    gaff_parameters.push_back(ap);

    ap.type = "H";
    ap.radius = 0.6000;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "HO";
    ap.radius = 0.0000;
    ap.epsilon = 0.0000;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "HS";
    ap.radius = 0.6000;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "HC";
    ap.radius = 1.4870;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "H1";
    ap.radius = 1.3870;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "H2";
    ap.radius = 1.2870;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "H3";
    ap.radius = 1.1870;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "HP";
    ap.radius = 1.1000;
    ap.epsilon = 0.0157;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "HA";
    ap.radius = 1.4590;
    ap.epsilon = 0.0150;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "H4";
    ap.radius = 1.4090;
    ap.epsilon = 0.0150;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "H5";
    ap.radius = 1.3590;
    ap.epsilon = 0.0150;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "HW";
    ap.radius = 0.0000;
    ap.epsilon = 0.0000;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "HZ";
    ap.radius = 1.4590;
    ap.epsilon = 0.0150;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "O";
    ap.radius = 1.6612;
    ap.epsilon = 0.2100;
    ap.mass = 16.00;
    gaff_parameters.push_back(ap);

    ap.type = "O2";
    ap.radius = 1.6612;
    ap.epsilon = 0.2100;
    ap.mass = 16.00;
    gaff_parameters.push_back(ap);

    ap.type = "OW";
    ap.radius = 1.7683;
    ap.epsilon = 0.1520;
    ap.mass = 16.00;
    gaff_parameters.push_back(ap);

    ap.type = "OH";
    ap.radius = 1.7210;
    ap.epsilon = 0.2104;
    ap.mass = 16.00;
    gaff_parameters.push_back(ap);

    ap.type = "OS";
    ap.radius = 1.6837;
    ap.epsilon = 0.1700;
    ap.mass = 16.00;
    gaff_parameters.push_back(ap);

    ap.type = "OP";
    ap.radius = 1.8500;
    ap.epsilon = 0.1700;
    ap.mass = 16.00;
    gaff_parameters.push_back(ap);

    ap.type = "C*";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CI";
    ap.radius = 1.9080;
    ap.epsilon = 0.1094;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "C5";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "C4";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CT";
    ap.radius = 1.9080;
    ap.epsilon = 0.1094;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CX";
    ap.radius = 1.9080;
    ap.epsilon = 0.1094;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "C";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "N";
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "N3";
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "S";
    ap.radius = 2.0000;
    ap.epsilon = 0.2500;
    ap.mass = 32.06;
    gaff_parameters.push_back(ap);

    ap.type = "SH";
    ap.radius = 2.0000;
    ap.epsilon = 0.2500;
    ap.mass = 32.06;
    gaff_parameters.push_back(ap);

    ap.type = "P";
    ap.radius = 2.1000;
    ap.epsilon = 0.2000;
    ap.mass = 30.97;
    gaff_parameters.push_back(ap);

    ap.type = "MG";
    ap.radius = 0.7926;
    ap.epsilon = 0.8947;
    ap.mass = 24.305;
    gaff_parameters.push_back(ap);

    ap.type = "C0";
    ap.radius = 1.7131;
    ap.epsilon = 0.45979;
    ap.mass = 40.08;
    gaff_parameters.push_back(ap);

    ap.type = "Zn";
    ap.radius = 1.10;
    ap.epsilon = 0.0125;
    ap.mass = 65.4;
    gaff_parameters.push_back(ap);

    ap.type = "F";
    ap.radius = 1.75;
    ap.epsilon = 0.061;
    ap.mass = 19.00;
    gaff_parameters.push_back(ap);

    ap.type = "Cl";
    ap.radius = 1.948;
    ap.epsilon = 0.265;
    ap.mass = 35.45;
    gaff_parameters.push_back(ap);

    ap.type = "Br";
    ap.radius = 2.22;
    ap.epsilon = 0.320;
    ap.mass = 79.90;
    gaff_parameters.push_back(ap);

    ap.type = "I";
    ap.radius = 2.35;
    ap.epsilon = 0.40;
    ap.mass = 126.9;
    gaff_parameters.push_back(ap);

    ap.type = "EP";
    ap.radius = 0.00;
    ap.epsilon = 0.0000;
    ap.mass = 0.00;
    gaff_parameters.push_back(ap);

    ap.type = "2C";
    ap.radius = 1.9080;
    ap.epsilon = 0.1094;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "3C";
    ap.radius = 1.9080;
    ap.epsilon = 0.1094;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "C8";
    ap.radius = 1.9080;
    ap.epsilon = 0.1094;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CO";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CA";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CB";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CC";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CD";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CK";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CM";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CQ";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CV";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CW";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CR";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CN";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CY";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CZ";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CP";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "CS";
    ap.radius = 1.9080;
    ap.epsilon = 0.0860;
    ap.mass = 12.01;
    gaff_parameters.push_back(ap);

    ap.type = "N2";
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "NA";
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "NB";
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "NC";
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "NT";
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "NY";
    ap.radius = 1.8240;
    ap.epsilon = 0.1700;
    ap.mass = 14.01;
    gaff_parameters.push_back(ap);

    ap.type = "Fe";
    ap.radius = 1.200;
    ap.epsilon = 0.0500;
    ap.mass = 55.00;
    gaff_parameters.push_back(ap);

    ap.type = "Cu";
    ap.radius = 2.200;
    ap.epsilon = 0.2000;
    ap.mass = 63.55;
    gaff_parameters.push_back(ap);

    ap.type = "Ca";
    ap.radius = 1.790;
    ap.epsilon = 0.0140;
    ap.mass = 40.00;
    gaff_parameters.push_back(ap);

    ap.type = "Si";
    ap.radius = 2.220;
    ap.epsilon = 0.3200;
    ap.mass = 28.00;
    gaff_parameters.push_back(ap);

    ap.type = "hw";
    ap.radius = 0.000;
    ap.epsilon = 0.0000;
    ap.mass = 1.008;
    gaff_parameters.push_back(ap);

    ap.type = "K";
    ap.radius = 2.658;
    ap.epsilon = 0.00033;
    ap.mass = 39.098;
    gaff_parameters.push_back(ap);

    this->gaff_force_field = gaff_parameters;
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
    string gaff_atom="";
    if (atom == "C.3"){
        gaff_atom = "c3";
    }
    else if (atom == "C3"){
        gaff_atom = "c3";
    }
    else if (atom =="C.2"){
        gaff_atom = "c2";
    }
    else if (atom =="C2"){
        gaff_atom = "c2";
    }
    else if (atom =="C1"){
        gaff_atom = "c1";
    }
    else if (atom =="C.1"){
        gaff_atom = "c1";
    }

    else if (atom =="C.ar"){
        gaff_atom = "ca";
    }
    else if (atom =="Car"){
        gaff_atom = "ca";
    }

    else if (atom =="C.cat"){
        gaff_atom = "c";
    }
    else if (atom =="Ccat"){
        gaff_atom = "c";
    }
    else if (atom =="C+"){
        gaff_atom = "c";
    }

    else if (atom =="N.3"){
        gaff_atom = "n3";
    }
    else if (atom =="N3"){
        gaff_atom = "n3";
    }


    else if (atom =="N.2"){
        gaff_atom = "n2";
    }
    else if (atom =="N2"){
        gaff_atom = "n2";
    }

    else if (atom =="N.1"){
        gaff_atom = "n1";
    }
    else if (atom =="N1"){
        gaff_atom = "n1";
    }

    else if (atom =="N.ar"){
        gaff_atom = "nh";
    }
    else if (atom =="Nar"){
        gaff_atom = "nh";
    }

    else if (atom =="N.am"){
        gaff_atom = "n";
    }
    else if (atom =="Nam"){
        gaff_atom = "n";
    }

    else if (atom =="N.4"){
        gaff_atom = "n4";
    }
    else if (atom =="N4"){
        gaff_atom = "n4";
    }
    else if (atom =="Ng+"){
        gaff_atom = "n4";
    }

    else if (atom =="N.pl3"){
        gaff_atom = "na";
    }
    else if (atom =="Npl3"){
        gaff_atom = "na";
    }
    else if (atom =="Npl"){
        gaff_atom = "na";
    }

    else if (atom =="N.p"){
        gaff_atom = "na";
    }
    else if (atom =="Np"){
        gaff_atom = "na";
    }

    else if (atom =="O.3"){
        gaff_atom = "oh";
    }
    else if (atom =="O3"){
        gaff_atom = "oh";
    }

    else if (atom =="O.2"){
        gaff_atom = "o";
    }
    else if (atom =="O2"){
        gaff_atom = "o";
    }

    else if (atom =="O.co2"){
        gaff_atom = "o";
    }
    else if (atom =="Oco2"){
        gaff_atom = "o";
    }

    else if (atom =="O-"){
        gaff_atom = "o";
    }

    else if (atom =="O.spc" or atom == "O.t3p"){
        gaff_atom = "ow";
    }

    else if (atom =="S.3"){
        gaff_atom = "sh"; //not sure... sh or ss
    }
    else if (atom =="S3"){
        gaff_atom = "sh"; //not sure... sh or ss
    }

    else if (atom =="S.2"){
        gaff_atom = "s2";
    }
    else if (atom =="S2"){
        gaff_atom = "s2";
    }

    else if (atom =="S.O" or atom == "S.o"){
        gaff_atom = "s4";
    }
    else if (atom =="SO" or atom == "So"){
        gaff_atom = "s4";
    }

    else if (atom =="S.O2" or atom == "S.o2"){
        gaff_atom = "s6";
    }
    else if (atom =="SO2" or atom == "So2"){
        gaff_atom = "s6";
    }

    else if (atom =="P.3"){
        gaff_atom = "p3";
    }
    else if (atom =="P3"){
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
        /*
        if (! found){
            printf("Atom type %s not found among GAFF parameters.\nPlease check Mol2.h source file.\n", atom.c_str());
            exit(1);
        }
*/
    }
    return(gaff_atom);
}

string Mol2::gaff_2_sybyl(string atom){
    string sybyl_atom="";
    if (atom == "c3"){
        sybyl_atom = "C.3";
    }
    else if (atom == "c2"){
        sybyl_atom = "C.2";
    }
    else if (atom == "c1"){
        sybyl_atom = "C.1";
    }
    else if (atom == "ca"){
        sybyl_atom = "C.ar";
    }
    else if (atom == "c"){
        sybyl_atom = "C.cat";
    }
    else if (atom == "n3"){
        sybyl_atom = "N.3";
    }
    else if (atom =="n2"){
        sybyl_atom = "N.2";
    }
    else if (atom =="n1"){
        sybyl_atom = "N.1";
    }
    else if (atom =="nh"){
        sybyl_atom = "N.ar";
    }
    else if (atom =="n"){
        sybyl_atom = "N.am";
    }
    else if (atom =="n4"){
        sybyl_atom = "N.4";
    }
    else if (atom =="na"){
        sybyl_atom = "N.pl";
    }
    else if (atom =="oh"){
        sybyl_atom = "O.3";
    }
    else if (atom =="o"){
        sybyl_atom = "O.2";
    }
    else if (atom =="ow"){
        sybyl_atom = "O.t3p";
    }
    else if (atom =="sh"){
        sybyl_atom = "S.3"; //not sure... sh or ss
    }
    else if (atom =="s2"){
        sybyl_atom = "S.2";
    }
    else if (atom =="s4"){
        sybyl_atom = "S.O";
    }
    else if (atom =="s6"){
        sybyl_atom = "S.O2";
    }
    else if (atom =="p3"){
        sybyl_atom = "P.3";
    }
    else if (atom =="f"){
        sybyl_atom = "F";
    }
    else if (atom =="hc"){
        sybyl_atom = "H";
    }
    else if (atom =="hw"){
        sybyl_atom = "H.t3p";
    }
    else if (atom =="cl"){
        sybyl_atom = "Cl";
    }
    else if (atom =="br"){
        sybyl_atom = "Br";
    }
    else if (atom =="i"){
        sybyl_atom = "I";
    }
    else if (atom =="MG"){
        sybyl_atom = "Mg";
    }
    else if (atom =="EP"){
        sybyl_atom = "LP";
    }
    else if (atom == "Fe"){
        sybyl_atom = "Fe";
    }
    else if (atom == "Zn"){
        sybyl_atom = "Zn";
    }
    else if (atom == "Cu"){
        sybyl_atom = "Cu";
    }
    else if (atom == "Ca"){
        sybyl_atom = "Ca";
    }
    else if (atom == "Si"){
        sybyl_atom = "Si";
    }
    else if (atom =="HO"){
        sybyl_atom = "H";
    }
    return(sybyl_atom);
}


string Mol2::sybyl_2_amber(string atom){
    string amber_atom="";
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
        /*
        if (! found){
            printf("Atom type %s not found among GAFF parameters.\nPlease check Mol2.h source file.\n", atom.c_str());
            exit(1);
        }
*/
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

vector<vector<double> > Mol2::copy_from_obmol(OBMol mymol){
    vector<vector<double > > vec_xyz;
    vector<double> tmp(3);
    double* myxyz = new double[mymol.NumAtoms()*3];
    myxyz = mymol.GetCoordinates();
    for (unsigned i=0; i < mymol.NumAtoms(); i++){
        tmp[0] = (myxyz[3*i]);
        tmp[1] = (myxyz[(3*i)+1]);
        tmp[2] = (myxyz[(3*i)+2]);
        vec_xyz.push_back(tmp);
    }

    tmp.clear();
    return vec_xyz;
}





#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;


BOOST_PYTHON_MODULE(pyMol2)
{

    class_< vector<double> >("vectorDouble")
        .def(vector_indexing_suite<vector<double> >())
    ;
    
    class_< vector<vector<double> > >("vectorvectorDouble")
        .def(vector_indexing_suite<vector<vector<double> > >())
    ;
        class_< vector<int> >("vectorInt")
        .def(vector_indexing_suite<vector<int> >())
    ;
        class_< vector<string> >("vectorString")
        .def(vector_indexing_suite<vector<string> >())
    ;
    
        class_< vector<vector<int> > >("vectorvectorInt")
        .def(vector_indexing_suite<vector<vector<int> > >())
    ;
    
    class_<Mol2>("Mol2", init< >())
        .def(init<PARSER*, string>())
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
        .def_readwrite("epsilons_sqrt", & Mol2::epsilons_sqrt)
        .def_readwrite("radii", & Mol2::radii)
        .def_readwrite("str", & Mol2::str)
        .def_readwrite("resnames", & Mol2::resnames)
        .def_readwrite("residue_pointer", & Mol2::residue_pointer)
        .def_readwrite("HBacceptors", & Mol2::HBacceptors)
        .def_readwrite("HBdonors", & Mol2::HBdonors)
        .def_readwrite("line", & Mol2::line)
        .def_readwrite("self_obj_function", & Mol2::self_obj_function)
        .def_readwrite("bonds", & Mol2::bonds)
        .def_readwrite("sybyl_atoms", & Mol2::sybyl_atoms)
        .def_readwrite("atm", & Mol2::atm)
        .def_readwrite("conformer_energies", & Mol2::conformer_energies)
        .def_readwrite("energy", & Mol2::energy)
        .def_readwrite("rmsd", & Mol2::rmsd)
        .def_readwrite("longest_axis", & Mol2::longest_axis)
        .def_readwrite("radius", & Mol2::radius)
        .def_readwrite("gaff_force_field", & Mol2::gaff_force_field)
        .def("parse_smiles", & Mol2::parse_smiles)
        .def("parse_gzipped_file", & Mol2::parse_gzipped_file)
        .def("parse_mol2file", & Mol2::parse_mol2file)
        .def("get_next_xyz", & Mol2::get_next_xyz)
        .def("initialize_gaff", & Mol2::initialize_gaff)
        .def("initialize_gaff2", & Mol2::initialize_gaff2)
        .def("get_gaff_atomic_parameters", & Mol2::get_gaff_atomic_parameters)
        .def("sybyl_2_gaff", & Mol2::sybyl_2_gaff)
        .def("sybyl_2_amber", & Mol2::sybyl_2_amber)
        .def("gaff_2_sybyl", & Mol2::gaff_2_sybyl)
        .def("find_longest_axis", & Mol2::find_longest_axis)
        .def("distance", & Mol2::distance)
        .def("copy_from_obmol", & Mol2::copy_from_obmol)
    ;

    class_<Mol2::atom_param>("atom_param")
        .def_readwrite("type", & Mol2::atom_param::type)
        .def_readwrite("radius", & Mol2::atom_param::radius)
        .def_readwrite("epsilon", & Mol2::atom_param::epsilon)
        .def_readwrite("mass", & Mol2::atom_param::mass)
    ;
}






