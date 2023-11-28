#include "../src/cudaEnergy2.cuh"
#include "../src/cudaGrid.cuh"
#include "../src/pyCOORD_MC.cpp"
#include "../src/pyEnergy2.cpp"
#include "../src/pyFindHB.cpp"
#include "../src/pyGrid.cpp"
#include "../src/pyMol2.cpp"
#include "../src/pyPARSER.cpp"
#include "../src/pyWRITER.cpp"
#include <array>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

void write_grid(FILE *file,
                const std::vector<std::vector<std::vector<double>>> &vec,
                bool print = false) {
  for (auto i : vec) {
    for (auto j : i) {
      for (auto k : j) {
        fwrite(&k, sizeof(double), 1, file);
        if (print)
          std::cout << k << ' ';
      }
    }
  }
}

void write_to_file(const Grid *grid) {

  FILE *outgrid;
  double elec, vdwA, vdwB, solv, rec_solv, pbsa, delphi, hb_donor, hb_acceptor;
  int pbsa_flag = 0;

  if (grid->pbsa_loaded) {
    pbsa_flag = 1;
  } else if (grid->delphi_loaded) {
    pbsa_flag = 2;
  }
  outgrid = fopen((grid->Input->grid_prefix + ".grid").c_str(), "wb");
  if (outgrid == NULL) {
    printf("Could not open McGrid file. Please check");
    exit(1);
  }

  write_grid(outgrid, grid->elec_grid);
  write_grid(outgrid, grid->vdwA_grid);
  write_grid(outgrid, grid->vdwB_grid);
  write_grid(outgrid, grid->rec_solv_gauss);
  write_grid(outgrid, grid->solv_gauss);
  write_grid(outgrid, grid->hb_donor_grid);
  write_grid(outgrid, grid->hb_acceptor_grid);

  if (grid->pbsa_loaded)
    write_grid(outgrid, grid->pbsa_grid);
  if (grid->delphi_loaded)
    write_grid(outgrid, grid->delphi_grid);

  fclose(outgrid);
}

int main(int argc, char *argv[]) {

  PARSER *input = new PARSER();
  std::vector<std::string> args;

  // Handles the case where the input args get concatenated into a single string
  if (argc == 2) {
    std::stringstream input{argv[1]};
    std::string temp;
    while (getline(input, temp, ' ')) {
      args.push_back(temp);
    }
  }

  else {
    args.assign(argv + 1, argv + argc);
  }

  assert((args.size() == 3) &&
         "Make sure you passed 3 arguments to the program: "
         "path/to/receptor.mol2, path/to/ligand.mol2, path/to/output_grid");

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
  input->x_dim = 10.0;
  input->y_dim = 10.0;
  input->z_dim = 10.0;
  input->timeout = 30;
  input->min_timeout = 30;
  input->overlay_optimizer = "mma";
  input->energy_optimizer = "mma";

  double delta = 2.5;
  input->deltaij6 = (delta * delta * delta * delta * delta * delta);
  double delta_es = 2.5;
  input->deltaij_es6 = pow(delta_es, 6);
  input->deltaij_es3 = (delta_es * delta_es * delta_es);

  input->grid_prefix = args[2];

  // ------------------------------------------------------------------------

  Mol2 *rec = new Mol2(input, args[0]);
  Mol2 *lig = new Mol2(input, args[1]);

  FindHB *hb = new FindHB();
  COORD_MC *coord = new COORD_MC();
  WRITER *writer = new WRITER(input);

  hb->find_ligandHB(input->reflig_mol2, lig);
  for (int i = 0; i < rec->residue_pointer.size() - 1; i++) {
    hb->parse_residue(rec->residue_pointer[i] - 1,
                      rec->residue_pointer[i + 1] - 2, rec->resnames[i], rec,
                      lig, 9.0);
  }

  std::vector<double> center_of_mass = coord->compute_com(lig);

  Grid *grid = new Grid(input, writer, rec, center_of_mass);
  write_to_file(grid);

  delete input;
  delete rec;
  delete lig;
  delete hb;
  delete coord;
  delete writer;
  delete grid;

  return 0;
}
