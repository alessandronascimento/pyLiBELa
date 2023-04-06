#include "FullSearch.h"

FullSearch::FullSearch(PARSER* _Input, Mol2* _Lig, WRITER* _Writer)
{
    this->Input = _Input;
    this->Lig = _Lig;
    this->Writer = _Writer;

    if (Input->dock_mode){
        Writer->writeMol2(Lig, Lig->xyz, 0.0, 0.0, "Lig_docked");
        mol = this->GetMol("Lig_docked.mol2.gz");
    }
    else {
        mol = this->GetMol(Input->lig_mol2);
    }

    if (Input->ligand_energy_model == "GAFF" or Input->ligand_energy_model == "gaff"){
        OBff = OBForceField::FindForceField("GAFF");
    }
    else {
        OBff = OBForceField::FindForceField("MMFF94");
    }

    if (Input->verbose){
        OBff->SetLogFile(&cout);
        OBff->SetLogLevel(OBFF_LOGLVL_LOW);
    }

    if (!OBff){
        cout << "Could not find OpenBabel FF parameters!" << endl;
        exit(1);
    }

    OBff->Setup(mol);
    OBff->GetCoordinates(mol);
    RotorList.Setup(mol);
    Rotor = RotorList.BeginRotor(RotorIterator);
    mol.ToInertialFrame();


    vector<int> tmp(4);
    sprintf(info, "Found %d rotatable bonds in ligand %s.", RotorList.Size(), Lig->molname.c_str());
    Writer->print_info(info);
    Writer->print_line();
    for (unsigned i = 0; i < RotorList.Size(); ++i, Rotor = RotorList.NextRotor(RotorIterator)) {
        tmp = Rotor->GetDihedralAtoms();
        atoms_in_dihedrals.push_back(tmp);
        tmp.clear();
    }
}

FullSearch::FullSearch(PARSER* _Input, Mol2* _Lig, WRITER* _Writer, Grid* _Grids)
{
    this->Input = _Input;
    this->Lig = _Lig;
    this->Writer = _Writer;
    this->Grids = _Grids;

    if (Input->dock_mode){
        Writer->writeMol2(Lig, Lig->xyz, 0.0, 0.0, "Lig_docked");
        mol = this->GetMol("Lig_docked.mol2.gz");
    }
    else {
        mol = this->GetMol(Input->lig_mol2);
    }

    if (Input->ligand_energy_model == "GAFF" or Input->ligand_energy_model == "gaff"){
        OBff = OBForceField::FindForceField("GAFF");
    }
    else {
        OBff = OBForceField::FindForceField("MMFF94");
    }

    if (Input->verbose){
        OBff->SetLogFile(&cout);
        OBff->SetLogLevel(OBFF_LOGLVL_LOW);
    }

    if (!OBff){
        cout << "Could not find OpenBabel FF parameters!" << endl;
        exit(1);
    }

    OBff->Setup(mol);
    OBff->GetCoordinates(mol);
    RotorList.Setup(mol);
    Rotor = RotorList.BeginRotor(RotorIterator);
    mol.ToInertialFrame();


    vector<int> tmp(4);
    sprintf(info, "Found %d rotatable bonds in ligand %s.", RotorList.Size(), Lig->molname.c_str());
    Writer->print_info(info);
    Writer->print_line();
    for (unsigned i = 0; i < RotorList.Size(); ++i, Rotor = RotorList.NextRotor(RotorIterator)) {
        tmp = Rotor->GetDihedralAtoms();
        atoms_in_dihedrals.push_back(tmp);
        tmp.clear();
    }
}

OBMol FullSearch::GetMol(const std::string &molfile){
    OBMol mol;

    OBConversion conv;
    OBFormat *format = conv.FormatFromExt(molfile.c_str());
    if (!format || !conv.SetInFormat(format)) {
    printf("Could not find input format for file\n");
    return mol;
  }

    ifstream ifs(molfile.c_str());
    if (!ifs) {
        printf("Could not open %s for reading.\n", molfile.c_str());
        return mol;
    }

    if (!conv.Read(&mol, &ifs)) {
        printf("Could not read molecule from file\n");
        return mol;
    }
    return mol;
}

vector<vector<double> > FullSearch::copy_from_obmol(OBMol mymol){
    vector<vector<double > > vec_xyz;
    vector<double> tmp(3);
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

void FullSearch::copy_to_obmol(vector<vector<double> > vec_xyz){
    double* dxyz = new double[vec_xyz.size()*3];
    for (unsigned i=0; i<vec_xyz.size(); i++){
        dxyz[3*i] = vec_xyz[i][0];
        dxyz[(3*i)+1] = vec_xyz[i][1];
        dxyz[(3*i)+2] = vec_xyz[i][2];
    }
    mol.SetCoordinates(dxyz);
    delete [] dxyz;
}


double FullSearch::do_search(void){
    COORD_MC* Coords = new COORD_MC;
    vector<double> com = Coords->compute_com(Lig->xyz, Lig);
    vector<vector<double> > new_xyz;
    Energy2* Energy = new Energy2(Input);
    double int_energy, energy;

    double x_lower_lim=com[0]-(Input->search_box_x/2);
    double x_upper_lim=com[0]+(Input->search_box_x/2);
    double y_lower_lim=com[1]-(Input->search_box_y/2);
    double y_upper_lim=com[1]+(Input->search_box_y/2);
    double z_lower_lim=com[2]-(Input->search_box_z/2);
    double z_upper_lim=com[2]+(Input->search_box_z/2);

    printf("#%10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s ", "Step", "Energy", "RMSD", "DX", "DY", "DZ", "ALPHA", "BETA", "GAMMA", "NCONF", "ConfEnergy", "Elec", "VDW", "SOLV");
    for (unsigned i=1; i <= RotorList.Size(); i++){
        sprintf(info, "ROT[%3d]", i);
        printf("%10.10s", info);
    }
    printf("\n");

    int count=0;
    for (double x=x_lower_lim; x<=x_upper_lim; x+=Input->translation_step){
        for (double y=y_lower_lim; y<=y_upper_lim; y+=Input->translation_step){
            for (double z=z_lower_lim; z<=z_upper_lim; z+=Input->translation_step){
                for (double alpha=0.; alpha<360.; alpha+=Input->rotation_step){
                    for (double beta=0.0; beta<180.; beta+=Input->rotation_step){
                        for (double gamma=0.0; gamma<360.; gamma+=Input->rotation_step){
                            new_xyz = Coords->rototranslate(Lig->xyz, Lig, alpha, beta, gamma, x, y, z);
                            for (unsigned tangle=0; tangle< RotorList.Size(); tangle++){
                                for (double angle=0.0; angle<360.0; angle+=Input->rotation_step){
                                    this->copy_to_obmol(new_xyz);
                                    mol.SetTorsion(mol.GetAtom(atoms_in_dihedrals[tangle][0]), mol.GetAtom(atoms_in_dihedrals[tangle][1]), mol.GetAtom(atoms_in_dihedrals[tangle][2]), mol.GetAtom(atoms_in_dihedrals[tangle][3]), angle*PI/180.);
                                    OBff->Setup(mol);
                                    int_energy = OBff->Energy();
                                    string unit = OBff->GetUnit();
                                    if (unit == "kJ/mol"){
                                        int_energy = int_energy/4.18;
                                    }
                                    new_xyz = this->copy_from_obmol(mol);
                                    energy = Energy->compute_ene(Grids, Lig, new_xyz);
                                    count++;
                                    if (energy<0){
                                        printf("%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n", count, (energy+int_energy), x, y, z,
                                               alpha, beta, gamma, int_energy, energy);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    delete Energy;
    delete Coords;
}

