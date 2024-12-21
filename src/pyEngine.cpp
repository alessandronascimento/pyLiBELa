#include "pyEngine.h"

using namespace std;

Engine::Engine(WRITER* _Writer){
	this->Writer = _Writer;
}

void Engine::Run_docking(PARSER* Input, Mol2* Rec, Mol2* RefLig, Grid* Grids){
    sprintf(info, "Running Dock mode in parallel with up to %d threads...", int(Input->parallel_jobs));
    this->print_info(info);
    vector<string> ligand_list;
    vector<string> ligand_codes;
    
    vector<double> center;
    COORD_MC* Coord = new COORD_MC();
    center = Coord->compute_com(RefLig);
    delete Coord;


/* Here, we read the multifile file and parse the files of the molecules
 * to be docked into a std::vector named ligand_list.
 */
 

    if (Input->use_smiles){
        ifstream multifile(Input->smiles_multifile.c_str());
        if (!multifile.is_open()){
            sprintf(info, "Could not open file %s. Exiting...\n", Input->smiles_multifile.c_str());
            this->print_info(info);
            exit(1);
        }
        string ligand, code;
        multifile >> ligand >> code;
        while((!multifile.eof()) and (ligand != "EOF")){
            ligand_list.push_back(ligand);
            ligand_codes.push_back(code);
            multifile >> ligand >> code;
        }
        multifile.close();
    }
    else {
        if (Input->multifile != ""){
            ifstream multifile(Input->multifile.c_str());
            if (! multifile.is_open()){
                sprintf(info, "Could not open file %s. Exiting...\n", Input->multifile.c_str());
                this->print_info(info);
                exit(1);
            }
            string ligand;
            multifile >> ligand;

            while ((!multifile.eof()) and (ligand != "EOF")){
                ligand_list.push_back(ligand);
                multifile >> ligand;
            }

            multifile.close();
        }
    }
    
	sprintf(info, "Found %d ligands to dock", int(ligand_list.size()));
	this->print_info(info);

/*
 * Here the parallel processing begins. Firstly, a parallel section is created with a pragma
 * indicating the number of parallel threads as defined in the input file.
 * After, the molecule objects are created and conformers are calculated using OpenBabel.
 * Due to OpenBabel thread unsafety issues, the conformer generation is done in serial, using
 * the sync directive #pragma critical. Finally, the docking objects are created and the molecules
 * are actually docked.
 */


#pragma omp parallel num_threads(Input->parallel_jobs)
        {
#pragma omp for schedule(static,1)
            for (int i=0; i< int(ligand_list.size()); i++){
                Mol2* Lig2 = new Mol2;
                bool lig_is_opened = false;
#pragma omp critical
                {
                    if (Input->use_smiles){
                        lig_is_opened = Lig2->parse_smiles(Input, ligand_list[i], ligand_codes[i]);
                    }
                    else{
                        if (ligand_list[i].substr(ligand_list[i].size()-3, 3) == ".gz"){
                            lig_is_opened = Lig2->parse_gzipped_file(Input, ligand_list[i]);
                        }
                        else {
                            lig_is_opened = Lig2->parse_mol2file(Input, ligand_list[i]);
                        }
                        if (lig_is_opened){
                            FindHB* HB = new FindHB;
                            HB->find_ligandHB(ligand_list[i], Lig2);
                            delete HB;
                        }
                        if (Input->generate_conformers){
                            Conformer* Conf = new Conformer;
                            Conf->generate_conformers_confab(Input, Lig2, ligand_list[i]);
                            delete Conf;
                        }
                    }
                }
                Docker* Dock = new Docker(Writer);
                if (lig_is_opened){
                    if (Input->use_grids){
                        Dock->run(Rec, Lig2, RefLig, center, Input, Grids, i+1);
/*                        printf("Elec: %7.3f, VdW: %7.3f, Solv: %7.3f, HB: %7.3f, Total: %7.3f\n", Dock->best_energy_t->elec, Dock->best_energy_t->vdw, 
                            (Dock->best_energy_t->rec_solv+Dock->best_energy_t->lig_solv), 
                            (Dock->best_energy_t->hb_donor+Dock->best_energy_t->hb_acceptor), Dock->best_energy_t->total);
*/                    }
                    else {
                        Dock->run(Rec, Lig2, RefLig, center, Input, i+1);
/*                        printf("Elec: %7.3f, VdW: %7.3f, Solv: %7.3f, HB: %7.3f, Total: %7.3f\n", Dock->best_energy_t->elec, Dock->best_energy_t->vdw, 
                            (Dock->best_energy_t->rec_solv+Dock->best_energy_t->lig_solv), 
                            (Dock->best_energy_t->hb_donor+Dock->best_energy_t->hb_acceptor), Dock->best_energy_t->total);
*/                    }
                    delete Dock;
                }
                delete Lig2;
        }
    }
    ligand_list.clear();
}

void Engine::print_info(char info[98]){
    Writer->print_info(info);
}

#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(pyEngine)
{
    class_<Engine>("Engine", init<WRITER*> ())
//    	.def_readwrite("info", & Engine::info)
    	.def_readwrite("Writer", & Engine::Writer)
    	.def("Run_docking", &Engine::Run_docking)
    	.def("print_info", &Engine::print_info)
    	;
}
