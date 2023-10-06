#include "../src/pyCOORD_MC.cpp"
#include "../src/pyFindHB.cpp"
#include "../src/pyGrid.cpp"
#include "../src/pyMol2.cpp"
#include "../src/pyPARSER.cpp"
#include "../src/pyWRITER.cpp"
#include "../src/pyEnergy2.cpp"
#include "../src/cudaGrid.cuh"
#include "../src/cudaEnergy2.cuh"
#include <cstdlib>
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {


  PARSER* input = new PARSER();

  // ------------- INPUT PARAMS ------------------------------------------
  
  input->generate_conformers = true;
  input->dock_parallel = true;
  input->parallel_jobs = 2;
  input->write_grids = true;
  input->use_grids = true;
  input->write_mol2 = true;
  input->atom_limit = 60;

  input->dielectric_model = "r"; 
  input->scoring_function = 0; 

  input->grid_spacing = 0.5; 
  input->solvation_alpha = 0.1;
  input->solvation_beta = -0.005;

  input->min_tol = 1E-10;
  input->min_delta = 1E-5;
  input->dock_min_tol = 1E-10;
  input->search_box_x = 20.0;
  input->search_box_y = 20.0;
  input->search_box_z = 20.0;
  input->x_dim = 30.0;
  input->y_dim = 30.0;
  input->z_dim = 30.0;
  input->timeout = 30; 
  input->min_timeout = 30;
  input->overlay_optimizer = "mma"; 
  input->energy_optimizer = "mma";

  double delta = 2.5; 
  input->deltaij6 = (delta*delta*delta*delta*delta*delta);
  double delta_es = 2.5; 
  input->deltaij_es6 = pow(delta_es, 6);
  input->deltaij_es3 = (delta_es*delta_es*delta_es);
  
  input->grid_prefix = argv[3]; 

  // ------------------------------------------------------------------------

  Mol2* rec = new Mol2(input, argv[1]);
  Mol2* lig = new Mol2(input, argv[2]);

  FindHB* hb = new FindHB();
  COORD_MC* coord = new COORD_MC();
  WRITER* writer = new WRITER(input);

  hb->find_ligandHB(input->reflig_mol2, lig);
  for (int i = 0; i < rec->residue_pointer.size() - 1; i++) {
    hb->parse_residue(rec->residue_pointer[i] - 1,
                      rec->residue_pointer[i + 1] - 2, rec->resnames[i], rec,
                      lig, 9.0);
  }

  vector<double> center_of_mass = coord->compute_com(lig);
  
  Grid* grid = new Grid(input, writer, rec, center_of_mass);
  grid->write_grids_to_file();

  delete input;
  delete rec;
  delete lig;
  delete hb;
  delete coord;
  delete writer;
  delete grid;

  for (auto i: grid->elec_grid){
    for (auto j: i){
      for (auto k: j){
        printf("%.2f ", k);
      }
      printf("\n");
    } 
    printf("\n");
  }
  return 0;
}
