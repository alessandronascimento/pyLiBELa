/*
 * Conformer.cpp
 *
 *  Created on: 02/08/2012
 *      Author: Nascimento
 */

#include "pyConformer.h"

using namespace std;
using namespace OpenBabel;

Conformer::Conformer() {
    OBMessageHandler messageHandler;
    messageHandler.SetOutputLevel(obError);
    OpenBabel::obErrorLog = messageHandler;
}

Conformer::~Conformer() {
}

OBMol Conformer::GetMol(const std::string &molfile){
    OBMol mol;

    OBConversion conv;
    OBFormat *format = conv.FormatFromExt(molfile.c_str());
    if (!format || !conv.SetInFormat(format)) {
    printf("Could not find input format for file %s\n", molfile.c_str());
    return mol;
  }

    ifstream ifs(molfile.c_str());
    if (!ifs) {
        printf("Could not open %s for reading. Does this file exist?\n", molfile.c_str());
        return mol;
    }

    if (!conv.Read(&mol, &ifs)) {
        printf("Could not read molecule from file %s. Maybe an OpenBabel issue?\n", molfile.c_str());
        return mol;
    }
    return mol;
}

bool Conformer::generate_conformers_confab(PARSER* Input, Mol2* Lig, string molfile){
    bool file_read;
    OBMol mol;

    mol = this->GetMol(molfile);
    OBMol ref_mol;
    ref_mol = this->GetMol(molfile);


    OBForceField* OBff;

    if (Input->verbose){
        OBff->SetLogFile(&cout);
        OBff->SetLogLevel(OBFF_LOGLVL_LOW);
    }

    if (Input->ligand_energy_model == "GAFF" or Input->ligand_energy_model == "gaff"){
        OBff = OBForceField::FindForceField("GAFF");
    }
    else {
        OBff = OBForceField::FindForceField("MMFF94");
    }

    if (!OBff){
        printf("Could not find OpenBabel FF parameters!\n");
        exit(1);
    }


    // Original conformation energy
    OBff->Setup(mol);
    mol.SetTotalCharge(mol.GetTotalCharge());
    double energy = OBff->Energy();
    if (OBff->GetUnit() == "kJ/mol"){       // Converting to kcal/mol, if needed.
        energy = energy/4.18;
    }

    // Do initial energy minimization prior to conformer generation
    OBff->GetCoordinates(mol);
    OBff->SteepestDescent(Input->conformer_min_steps);
    OBff->UpdateCoordinates(mol);
    energy = OBff->Energy();
    if (OBff->GetUnit() == "kJ/mol"){       // Converting to kcal/mol, if needed.
        energy = energy/4.18;
    }

    // Conformer Search
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
            Lig->conformer_energies.push_back(energy);

            OBAlign* align = new OBAlign;
            align->SetRefMol(ref_mol);
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
            Lig->mcoords.push_back(xyz_tmp);
            xyz_tmp.clear();
        }
        file_read = true;
    }
    else {
        file_read = false;
    }
    return file_read;
}





#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;


BOOST_PYTHON_MODULE(pyConformer)
{

    class_<Conformer>("Conformer", init< >())
	.def("GetMol", & Conformer::GetMol)
	.def("generate_conformers_confab", & Conformer::generate_conformers_confab)

    ;


}

