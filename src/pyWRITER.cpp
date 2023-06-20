/*
 * WRITER.cpp
 *
 *  Created on: 23/11/2010
 *      Author: Nascimento
 */

#include "pyWRITER.h"

WRITER::WRITER(string prefix){
	WRITER::output_prefix = prefix;
	outputfile=WRITER::output_prefix+".log";
	output = fopen(outputfile.c_str(), "w");
    if ((Input->write_mol2) and (Input->dock_mode)){
        outmol2 = gzopen((prefix+"_dock.mol2.gz").c_str(), "w");
        gzclose(outmol2);
    }
}


WRITER::WRITER(string prefix, PARSER* _Input){
    this->Input = _Input;
    WRITER::output_prefix = prefix;
    outputfile=WRITER::output_prefix+".log";
    output = fopen(outputfile.c_str(), "w");
    if ((Input->write_mol2) and (Input->dock_mode)){
        outmol2 = gzopen((prefix+"_dock.mol2.gz").c_str(), "w");
        gzclose(outmol2);
    }
}

WRITER::WRITER(PARSER* _Input){
    this->Input = _Input;
	WRITER::output_prefix = Input->output;
	outputfile=WRITER::output_prefix+".log";
	output = fopen(outputfile.c_str(), "w");
    if ((Input->write_mol2) and Input->dock_mode){
        outmol2 = gzopen((Input->output+"_dock.mol2.gz").c_str(), "w");
        gzclose(outmol2);
    }
	this->print_welcome();
    this->print_params();
}

void WRITER::print_welcome(void){
	printf("****************************************************************************************************\n");
	printf("****************************************************************************************************\n");
	printf("*                  MCLiBELa - Monte Carlo-based Ligand Binding Energy Landscape                    *\n");
    printf("*                                      Version 1.0  - Build %5d                                  *\n", BUILD);
	printf("*                                                                                                  *\n");
    printf("* University of Sao Paulo                                                                          *\n");
	printf("* More Info:                                                                                       *\n");
	printf("*      http://www.biotechmol.ifsc.usp.br/                                                          *\n");
	printf("*                                                                                                  *\n");
	printf("****************************************************************************************************\n");
	printf("****************************************************************************************************\n");

	fprintf(output, "****************************************************************************************************\n");
	fprintf(output, "****************************************************************************************************\n");
	fprintf(output, "*                  MCLiBELa - Monte Carlo-based Ligand Binding Energy Landscape                    *\n");
    fprintf(output, "*                                      Version 1.0  - Build %5d                                  *\n", BUILD);
	fprintf(output, "*                                                                                                  *\n");
    fprintf(output, "* University of Sao Paulo                                                                          *\n");
	fprintf(output, "* More Info:                                                                                       *\n");
	fprintf(output, "*      http://www.biotechmol.ifsc.usp.br                                                           *\n");
	fprintf(output, "*                                                                                                  *\n");
	fprintf(output, "****************************************************************************************************\n");
	fprintf(output, "****************************************************************************************************\n");

}

void WRITER::print_params(){
	printf("*                                                                                                  *\n");
    if (Input->dock_mode){
        printf("* %-30s %-66.66s*\n", "mode", "Docking");
    }
    printf("* %-30s %-66d*\n", "dock_parallel", Input->dock_parallel);
    if (Input->dock_parallel){
        printf("* %-30s %-66d*\n", "parallel_jobs", Input->parallel_jobs);
    }
    printf("*                                                                                                  *\n");
    printf("* %-30s %-66.66s*\n", "rec_mol2", Input->rec_mol2.c_str());
    printf("* %-30s %-66.66s*\n", "lig_mol2", Input->lig_mol2.c_str());
    printf("* %-30s %-66.66s*\n", "reflig_mol2", Input->reflig_mol2.c_str());
    printf("* %-30s %-66.66s*\n", "multifile", Input->multifile.c_str());
    printf("* %-30s %-66d*\n", "mol2_aa", Input->mol2_aa);
    printf("*                                                                                                  *\n");
    switch (Input->scoring_function){
    case 0:
        printf("* %-30s %-66.66s*\n", "scoring function", "Amber Softcore + Desolvation");
        printf("* %-30s %-66.2f*\n", "deltaij6", Input->deltaij6);
        printf("* %-30s %-66.2f*\n", "deltaij_es6", Input->deltaij_es6);
        break;
    case 1:
        printf("* %-30s %-66.66s*\n", "scoring function", "Amber Softcore");
        printf("* %-30s %-66.2f*\n", "deltaij6", Input->deltaij6);
        printf("* %-30s %-66.2f*\n", "deltaij_es6", Input->deltaij_es6);
        break;
    case 2:
        printf("* %-30s %-66.66s*\n", "scoring function", "Amber FF + Desolvation");
        break;
    case 3:
        printf("* %-30s %-66.66s*\n", "scoring function", "Amber FF");
        break;
    case 4:
        printf("* %-30s %-66.66s*\n", "scoring function", "Amber FF with Gaussian Weighted LJ Potential + Desolvation");
        printf("* %-30s %-66.2f*\n", "LJ_sigma", Input->LJ_sigma);
        break;
    case 5:
        printf("* %-30s %-66.66s*\n", "scoring function", "Amber FF with Gaussian Weighted LJ Potential");
        printf("* %-30s %-66.2f*\n", "LJ_sigma", Input->LJ_sigma);
        break;
    }
    printf("* %-30s %-66.66s*\n", "dielectric_model", Input->dielectric_model.c_str());
    printf("* %-30s %-66.3f*\n", "diel", Input->diel);
    printf("* %-30s %-66.3f*\n", "sigma", Input->sigma);
    printf("* %-30s %-66.3f*\n", "solvation_alpha", Input->solvation_alpha);
    printf("* %-30s %-66.3f*\n", "solvation_beta", Input->solvation_beta);

    printf("* %-30s %-66d*\n", "use grids?", Input->use_grids);
    if (Input->use_grids){
        printf("* %-30s %-66.2f*\n", "grid_spacing", Input->grid_spacing);
        printf("* %-30s %-22.2f%-22.2f%-22.2f*\n", "grid box", Input->x_dim, Input->y_dim, Input->z_dim);
        if (Input->write_grids){
            printf("* %-30s %-66.66s*\n", "write grids", Input->grid_prefix.c_str());
        }
        else {
            printf("* %-30s %-66.66s*\n", "load grids", Input->grid_prefix.c_str());
        }
    }
    printf("* %-30s %-66d*\n", "use PBSA grids?", Input->use_pbsa);
    if (Input->use_pbsa){
        printf("* %-30s %-66.66s*\n", "Electrostatic model", "PBSA");
        printf("* %-30s %-66.66s*\n", "PBSA grid", Input->pbsa_grid.c_str());
    }
    else if (Input->use_delphi){
        printf("* %-30s %-66.66s*\n", "Electrostatic model", "DelPhi");
        printf("* %-30s %-66.66s*\n", "DelPhi grid", Input->delphi_grid.c_str());
        printf("* %-30s %-66d*\n", "Delphi grid gsize", Input->delphi_gsize);

    }
    else{
        printf("* %-30s %-66.66s*\n", "Electrostatic model", "Coulomb");
    }
    printf("*                                                                                                  *\n");
    printf("* %-30s %-66d*\n", "only_score", Input->only_score);
    printf("* %-30s %-66d*\n", "use_score_optimization", Input->use_score_optimization);
    printf("* %-30s %-22.2f%-22.2f%-22.2f*\n", "search_box", Input->search_box_x, Input->search_box_y, Input->search_box_z);
    printf("* %-30s %-66.10f*\n", "minimization_tolerance", Input->min_tol);
    printf("* %-30s %-66.10f*\n", "minimization_delta", Input->min_delta);
    printf("* %-30s %-66.10f*\n", "dock_min_tol", Input->dock_min_tol);
    printf("* %-30s %-66d*\n", "minimization_timeout", Input->min_timeout);
    printf("* %-30s %-66d*\n", "sort_by_energy", Input->sort_by_energy);
    printf("* %-30s %-66.2f*\n", "elec_scale", Input->elec_scale);
    printf("* %-30s %-66.2f*\n", "vdw_scale", Input->vdw_scale);
    printf("* %-30s %-66.66s*\n", "overlay_optimizer", Input->overlay_optimizer.c_str());
    printf("* %-30s %-66.66s*\n", "energy_optimizer", Input->energy_optimizer.c_str());
    printf("* %-30s %-66d*\n", "ignore_h", Input->dock_no_h);
    printf("* %-30s %-66d*\n", "deal", Input->deal);
    printf("* %-30s %-66.2f*\n", "scale_elec_energy", Input->scale_elec_energy);
    printf("* %-30s %-66.2f*\n", "scale_vdw_energy", Input->scale_vdw_energy);
    printf("* %-30s %-66.2f*\n", "Overlay Cutoff", Input->overlay_cutoff);

    printf("*                                                                                                  *\n");
    printf("* %-30s %-66d*\n", "generate_conformers", Input->generate_conformers);
    if (Input->generate_conformers){
        printf("* %-30s %-66d*\n", "number_of_conformers", Input->lig_conformers);
        printf("* %-30s %-66d*\n", "conformers_to_rank", Input->conformers_to_evaluate);
    }
    printf("*                                                                                                  *\n");
    printf("* %-30s %-66.66s*\n", "output_prefix", Input->output.c_str());
    printf("* %-30s %-66d*\n", "write_mol2", Input->write_mol2);
    if (Input->use_writeMol2_score_cutoff){
        printf("* %-30s %-66.2f*\n", "Cutoff in Overlay for Mol2", Input->writeMol2_score_cutoff);
    }
    if (Input->use_writeMol2_energy_cutoff){
        printf("* %-30s %-66.2f*\n", "Cutoff in Energy for Mol2", Input->writeMol2_energy_cutoff);
    }
	printf("*                                                                                                  *\n");
    printf("* %-30s %-66.66s*\n", "Ligand energy model", Input->ligand_energy_model.c_str());
    printf("* %-30s %-66.66s*\n", "Atomic FF model", Input->atomic_model_ff.c_str());
    printf("*                                                                                                  *\n");
    if (Input->eq_mode){
        printf("* %-30s %-66.10f*\n", "MC Temperature", Input->temp);
        printf("* %-30s %-66d*\n", "MC Sampling of Torsions", Input->sample_torsions);
        printf("* %-30s %-66.10f*\n", "MC Maximal COM Translation", Input->cushion);
        printf("* %-30s %-66.10f*\n", "MC Maximal COM Rotation", Input->rotation_step);
        if (Input->mc_full_flex){
            printf("* %-30s %-66d*\n", "MC Fully Flexible Sampling", Input->mc_full_flex);
            printf("* %-30s %-66.10f*\n", "MC atomic displacement", Input->max_atom_displacement);
        }
        printf("* %-30s %-66d*\n", "Entropy bins for rotation", Input->entropy_rotation_bins);
        printf("* %-30s %-66d*\n", "Entropy bins for translation", Input->entropy_translation_bins);
    }
	printf("****************************************************************************************************\n");
	printf("****************************************************************************************************\n");

    //

    fprintf(output, "*                                                                                                  *\n");
    if (Input->dock_mode){
        fprintf(output, "* %-30s %-66.66s*\n", "mode", "Docking");
    }
    fprintf(output, "* %-30s %-66d*\n", "dock_parallel", Input->dock_parallel);
    if (Input->dock_parallel){
        fprintf(output, "* %-30s %-66d*\n", "parallel_jobs", Input->parallel_jobs);
    }
    fprintf(output, "*                                                                                                  *\n");
    fprintf(output, "* %-30s %-66.66s*\n", "rec_mol2", Input->rec_mol2.c_str());
    fprintf(output, "* %-30s %-66.66s*\n", "lig_mol2", Input->lig_mol2.c_str());
    fprintf(output, "* %-30s %-66.66s*\n", "reflig_mol2", Input->reflig_mol2.c_str());
    fprintf(output, "* %-30s %-66.66s*\n", "multifile", Input->multifile.c_str());
    fprintf(output, "* %-30s %-66d*\n", "mol2_aa", Input->mol2_aa);
    fprintf(output, "*                                                                                                  *\n");
    switch (Input->scoring_function){
    case 0:
        fprintf(output, "* %-30s %-66.66s*\n", "scoring function", "Amber Softcore + Dessolvation");
        fprintf(output, "* %-30s %-66.2f*\n", "deltaij6", Input->deltaij6);
        fprintf(output, "* %-30s %-66.2f*\n", "deltaij_es6", Input->deltaij_es6);
        break;
    case 1:
        fprintf(output, "* %-30s %-66.66s*\n", "scoring function", "Amber Softcore");
        fprintf(output, "* %-30s %-66.2f*\n", "deltaij6", Input->deltaij6);
        fprintf(output, "* %-30s %-66.2f*\n", "deltaij_es6", Input->deltaij_es6);
        break;
    case 2:
        fprintf(output, "* %-30s %-66.66s*\n", "scoring function", "Amber FF + Dessolvation");
        break;
    case 3:
        fprintf(output, "* %-30s %-66.66s*\n", "scoring function", "Amber FF");
        break;
    case 4:
        fprintf(output, "* %-30s %-66.66s*\n", "scoring function", "Amber FF with Gaussian Weighted LJ Potential + Desolvation");
        fprintf(output, "* %-30s %-66.2f*\n", "LJ_sigma", Input->LJ_sigma);
        break;
    case 5:
        fprintf(output, "* %-30s %-66.66s*\n", "scoring function", "Amber FF with Gaussian Weighted LJ Potential");
        fprintf(output, "* %-30s %-66.2f*\n", "LJ_sigma", Input->LJ_sigma);
        break;
    }
    if (Input->use_pbsa){
        fprintf(output, "* %-30s %-66.66s*\n", "Electrostatic model", "PBSA");
        fprintf(output, "* %-30s %-66.66s*\n", "PBSA grid", Input->pbsa_grid.c_str());
    }
    else if (Input->use_delphi){
        fprintf(output, "* %-30s %-66.66s*\n", "Electrostatic model", "DelPhi");
        fprintf(output, "* %-30s %-66.66s*\n", "DelPhi grid", Input->delphi_grid.c_str());
        fprintf(output, "* %-30s %-66d*\n", "Delphi grid gsize", Input->delphi_gsize);

    }
    else{
        fprintf(output, "* %-30s %-66.66s*\n", "Electrostatic model", "Coulomb");
    }
    fprintf(output, "* %-30s %-66.66s*\n", "dielectric_model", Input->dielectric_model.c_str());
    fprintf(output, "* %-30s %-66.3f*\n", "diel", Input->diel);
    fprintf(output, "* %-30s %-66.3f*\n", "sigma", Input->sigma);
    fprintf(output, "* %-30s %-66.3f*\n", "solvation_alpha", Input->solvation_alpha);
    fprintf(output, "* %-30s %-66.3f*\n", "solvation_beta", Input->solvation_beta);
    fprintf(output, "* %-30s %-66d*\n", "use grids?", Input->use_grids);
    if (Input->use_grids){
        fprintf(output, "* %-30s %-66.2f*\n", "grid spacing", Input->grid_spacing);
        fprintf(output, "* %-30s %-22.2f%-22.2f%-22.2f*\n", "grid box", Input->x_dim, Input->y_dim, Input->z_dim);
        if (Input->write_grids){
            fprintf(output, "* %-30s %-66.66s*\n", "write grids", Input->grid_prefix.c_str());
        }
        else {
            fprintf(output, "* %-30s %-66.66s*\n", "load grids", Input->grid_prefix.c_str());
        }
    }
    fprintf(output, "*                                                                                                  *\n");
    fprintf(output, "* %-30s %-66d*\n", "only_score", Input->only_score);
    fprintf(output, "* %-30s %-66d*\n", "use_score_optimization", Input->use_score_optimization);
    fprintf(output, "* %-30s %-22.2f%-22.2f%-22.2f*\n", "search_box", Input->search_box_x, Input->search_box_y, Input->search_box_z);
    fprintf(output, "* %-30s %-66.10f*\n", "minimization_tolerance", Input->min_tol);
    fprintf(output, "* %-30s %-66.10f*\n", "minimization_delta", Input->min_delta);
    fprintf(output, "* %-30s %-66.10f*\n", "dock_min_tol", Input->dock_min_tol);
    fprintf(output, "* %-30s %-66d*\n", "minimization_timeout", Input->min_timeout);
    fprintf(output, "* %-30s %-66d*\n", "sort_by_energy", Input->sort_by_energy);
    fprintf(output, "* %-30s %-66.2f*\n", "elec_scale", Input->elec_scale);
    fprintf(output, "* %-30s %-66.2f*\n", "vdw_scale", Input->vdw_scale);
    fprintf(output, "* %-30s %-66.66s*\n", "overlay_optimizer", Input->overlay_optimizer.c_str());
    fprintf(output, "* %-30s %-66.66s*\n", "energy_optimizer", Input->energy_optimizer.c_str());
    fprintf(output, "* %-30s %-66.10f*\n", "Overlay Cuttoff", Input->overlay_cutoff);
    fprintf(output, "* %-30s %-66d*\n", "ignore_h", Input->dock_no_h);
    fprintf(output, "* %-30s %-66d*\n", "deal", Input->deal);
    fprintf(output, "*                                                                                                  *\n");
    fprintf(output, "* %-30s %-66d*\n", "generate_conformers", Input->generate_conformers);
    if (Input->generate_conformers){
        fprintf(output, "* %-30s %-66d*\n", "number_of_conformers", Input->lig_conformers);
        fprintf(output, "* %-30s %-66d*\n", "conformers_to_rank", Input->conformers_to_evaluate);
    }
    fprintf(output, "*                                                                                                  *\n");
    fprintf(output, "* %-30s %-66.66s*\n", "output_prefix", Input->output.c_str());
    fprintf(output, "* %-30s %-66d*\n", "write_mol2", Input->write_mol2);
    if (Input->use_writeMol2_score_cutoff){
        printf("* %-30s %-66.2f*\n", "Cutoff in Overlay for Writting Mol2", Input->writeMol2_score_cutoff);
    }
    if (Input->use_writeMol2_energy_cutoff){
        printf("* %-30s %-66.2f*\n", "Cutoff in Energy for Writting Mol2", Input->writeMol2_energy_cutoff);
    }
	fprintf(output, "*                                                                                                  *\n");
    fprintf(output, "* %-30s %-66.66s*\n", "Ligand energy model", Input->ligand_energy_model.c_str());
    fprintf(output, "* %-30s %-66.66s*\n", "Atomic FF model", Input->atomic_model_ff.c_str());
    fprintf(output, "*                                                                                                  *\n");
    if (Input->eq_mode){
        fprintf(output, "* %-30s %-66d*\n", "MC Sampling of Torsions", Input->sample_torsions);
        fprintf(output, "* %-30s %-66.10f*\n", "MC Maximal COM Translation", Input->cushion);
        fprintf(output, "* %-30s %-66.10f*\n", "MC Maximal COM Rotation", Input->rotation_step);
        if (Input->mc_full_flex){
            fprintf(output, "* %-30s %-66d*\n", "MC Fully Flexible Sampling", Input->mc_full_flex);
            fprintf(output, "* %-30s %-66.10f*\n", "MC atomic displacement", Input->max_atom_displacement);
        }
        fprintf(output, "* %-30s %-66d*\n", "Entropy bins for rotation", Input->entropy_rotation_bins);
        fprintf(output, "* %-30s %-66d*\n", "Entropy bins for translation", Input->entropy_translation_bins);
    }

	fprintf(output, "****************************************************************************************************\n");
	fprintf(output, "****************************************************************************************************\n");
}

void WRITER::write_box(vector<double>center, double min_x, double min_y, double min_z, double max_x, double max_y, double max_z){

	FILE *box;
	box = fopen("box.pdb", "w");

	fprintf (box, "REMARK    CENTER OF THE BOX  %2.3f  %2.3f  %2.3f\n", center[0], center[1], center[2]);
	fprintf (box, "ATOM      1  DUA BOX     1     % 7.3f % 7.3f % 7.3f\n", min_x, min_y, min_z);
	fprintf (box, "ATOM      2  DUB BOX     1     % 7.3f % 7.3f % 7.3f\n", max_x, min_y, min_z);
	fprintf (box, "ATOM      3  DUC BOX     1     % 7.3f % 7.3f % 7.3f\n", max_x, min_y, max_z);
	fprintf (box, "ATOM      4  DUD BOX     1     % 7.3f % 7.3f % 7.3f\n", min_x, min_y, max_z);
	fprintf (box, "ATOM      5  DUE BOX     1     % 7.3f % 7.3f % 7.3f\n", min_x, max_y, min_z);
	fprintf (box, "ATOM      6  DUF BOX     1     % 7.3f % 7.3f % 7.3f\n", max_x, max_y, min_z);
	fprintf (box, "ATOM      7  DUG BOX     1     % 7.3f % 7.3f % 7.3f\n", max_x, max_y, max_z);
	fprintf (box, "ATOM      8  DUH BOX     1     % 7.3f % 7.3f % 7.3f\n", min_x, max_y, max_z);
	fprintf (box, "CONECT    1    2    4    5\n");
	fprintf (box, "CONECT    2    1    3    6\n");
	fprintf (box, "CONECT    3    2    4    7\n");
	fprintf (box, "CONECT    4    1    3    8\n");
	fprintf (box, "CONECT    5    1    6    8\n");
	fprintf (box, "CONECT    6    2    5    7\n");
	fprintf (box, "CONECT    7    3    6    8\n");
	fprintf (box, "CONECT    8    4    5    7\n");
	fprintf (box, "END\n");

	fclose(box);
}

void WRITER::write_pdb(vector<vector<double> >xyz, int N, vector<string> atomnames, double energy, double rmsd, string outname){
	gzFile outpdb;
	outpdb = gzopen((WRITER::output_prefix+".pdb.gz").c_str(), "a");
	gzprintf(outpdb, "REMARK\n");
	gzprintf(outpdb, "REMARK %-9s energy = %9.2f rmsd   = %12.3f\n", outname.c_str(), energy, rmsd);
	for (int i=0; i<N; i++){
		gzprintf(outpdb, "ATOM    %3d%4s  LIG     1    % 8.3f % 7.3f % 7.3f   0.000    0.00  1\n", i+1, atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2]);
	}
	gzprintf(outpdb, "TER\n");
	gzclose(outpdb);
}

void WRITER::write_pdb(Mol2 *Cmol, vector<vector<double> >xyz, double energy, double rmsd, string outname){
	gzFile outpdb;
	outpdb = gzopen((outname+".pdb.gz").c_str(), "a");
	gzprintf(outpdb, "MDL\n");
	gzprintf(outpdb, "REMARK\n");
	gzprintf(outpdb, "REMARK %-9s energy = %9.2f rmsd   = %12.3f\n", outname.c_str(), energy, rmsd);
	int i=0;
	int resn=0;

	while (resn < int(Cmol->residue_pointer.size()-1)){
		while(i < Cmol->residue_pointer[resn+1]){
			gzprintf(outpdb, "ATOM   %4d%4s  %3.3s  %4d    % 8.3f % 7.3f % 7.3f   0.000    0.00  1\n", i+1, Cmol->sybyl_atoms[i].c_str(), Cmol->resnames[resn].c_str(), resn+1, xyz[i][0], xyz[i][1], xyz[i][2]);
			i++;
		}
		resn++;
	}

	while (i < Cmol->N){
		gzprintf(outpdb, "ATOM   %4d%4s  %3.3s  %4d    % 8.3f % 7.3f % 7.3f   0.000    0.00  1\n", i+1, Cmol->sybyl_atoms[i].c_str(), Cmol->resnames[resn].c_str(), resn+1, xyz[i][0], xyz[i][1], xyz[i][2]);
		i++;
	}
	gzprintf(outpdb, "TER\n");
	gzprintf(outpdb, "ENDMDL\n");
	gzclose(outpdb);
}

void WRITER::writeMol2(Mol2* Cmol, vector<vector<double> >xyz, double energy, double rmsd, string outname){
	gzFile outmol2;

	outmol2 = gzopen((outname+".mol2.gz").c_str(), "a");
	gzprintf(outmol2, "\n");
	gzprintf(outmol2, "########## %15.15s: %19.19s\n", "Name", Cmol->molname.c_str());
	gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "Energy Score", energy);
	gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "RMSD", rmsd);
	gzprintf(outmol2, "\n");
	gzprintf(outmol2, "@<TRIPOS>MOLECULE\n");
	gzprintf(outmol2, "%s\n", Cmol->molname.c_str());
	gzprintf(outmol2, "%d %d %d\n", Cmol->N, Cmol->Nbonds, Cmol->Nres);
	gzprintf(outmol2, "SMALL\n");
	gzprintf(outmol2, "USER_CHARGES\n");
	gzprintf(outmol2, "Energy: %7.3f\n", energy);
	gzprintf(outmol2, "RMSD/OVERLAY: %7.3f\n", rmsd);
	gzprintf(outmol2, "@<TRIPOS>ATOM\n");
	int i=0;
	unsigned resn=0;

	while(resn < Cmol->residue_pointer.size()-1){
		while(i < Cmol->residue_pointer[resn+1]-1){
			if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
                gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
			}
			else {
                gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
			}
			i++;
		}
		resn++;
	}
	while(i < Cmol->N){
		if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
            gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
		}
		else{
            gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
		}
		i++;
	}

	gzprintf(outmol2, "@<TRIPOS>BOND\n");

	for (unsigned j=0; j<Cmol->bonds.size(); j++){
		gzprintf(outmol2, "%6d%6.6s%6.6s%5.5s\n", j+1, Cmol->bonds[j][0].c_str(), Cmol->bonds[j][1].c_str(),Cmol->bonds[j][2].c_str());
	}
	gzclose(outmol2);
}

void WRITER::writeMol2(Mol2* Cmol, vector<vector<double> >xyz, energy_result_t* result, double rmsd, string outname){
    gzFile outmol2;

    outmol2 = gzopen((outname+".mol2.gz").c_str(), "a");
    gzprintf(outmol2, "\n");
    gzprintf(outmol2, "########## %15.15s: %19.19s\n", "Name", Cmol->molname.c_str());
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "Energy Score", result->total);
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "Elec Score", result->elec);
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "VDW Score", result->vdw);
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "HB Score", result->hb_donor+result->hb_acceptor);
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "Rest Score", result->restraints);
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "Solv Score", (result->rec_solv + result->lig_solv));
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "RMSD", rmsd);
    gzprintf(outmol2, "\n");
    gzprintf(outmol2, "@<TRIPOS>MOLECULE\n");
    gzprintf(outmol2, "%s\n", Cmol->molname.c_str());
    gzprintf(outmol2, "%d %d %d\n", Cmol->N, Cmol->Nbonds, Cmol->Nres);
    gzprintf(outmol2, "SMALL\n");
    gzprintf(outmol2, "USER_CHARGES\n");
    gzprintf(outmol2, "@<TRIPOS>ATOM\n");
    int i=0;
    unsigned resn=0;

    while(resn < Cmol->residue_pointer.size()-1){
        while(i < Cmol->residue_pointer[resn+1]-1){
            if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
                gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
            }
            else {
                gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
            }
            i++;
        }
        resn++;
    }
    while(i < Cmol->N){
        if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
            gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
        }
        else{
            gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
        }
        i++;
    }

    gzprintf(outmol2, "@<TRIPOS>BOND\n");

    for (unsigned j=0; j<Cmol->bonds.size(); j++){
        gzprintf(outmol2, "%6d%6.6s%6.6s%5.5s\n", j+1, Cmol->bonds[j][0].c_str(), Cmol->bonds[j][1].c_str(),Cmol->bonds[j][2].c_str());
    }
    gzclose(outmol2);
}

void WRITER::writeMol2(Mol2* Cmol, vector<vector<double> >xyz, double energy, double rmsd){
    outmol2 = gzopen((this->output_prefix+"_dock.mol2.gz").c_str(), "a");
	gzprintf(outmol2, "\n");
	gzprintf(outmol2, "########## %15.15s: %19.19s\n", "Name", Cmol->molname.c_str());
	gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "Energy Score", energy);
	gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "RMSD", rmsd);
	gzprintf(outmol2, "\n");
	gzprintf(outmol2, "@<TRIPOS>MOLECULE\n");
	gzprintf(outmol2, "%s\n", Cmol->molname.c_str());
	gzprintf(outmol2, "%d %d %d\n", Cmol->N, Cmol->Nbonds, Cmol->Nres);
	gzprintf(outmol2, "SMALL\n");
	gzprintf(outmol2, "USER_CHARGES\n");
	gzprintf(outmol2, "Energy: %7.3f\n", energy);
	gzprintf(outmol2, "RMSD/OVERLAY: %7.3f\n", rmsd);
	gzprintf(outmol2, "@<TRIPOS>ATOM\n");
	int i=0;
	unsigned resn=0;

	if (int(xyz.size()) != Cmol->N){
		printf("Mismatch in atom number while writting mol2 file. Please check!\n");
		printf("Cmol->N = %d    xyz.size() = %d\n", Cmol->N, int(xyz.size()));
		exit(1);
	}

	while(resn < Cmol->residue_pointer.size()-1){
		while(i < Cmol->residue_pointer[resn+1]-1){
			if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
                gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
			}
			else {
                gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
			}
			i++;
		}
		resn++;
	}
	while(i < Cmol->N){
		if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
            gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
		}
		else{
            gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
		}
		i++;
	}

	gzprintf(outmol2, "@<TRIPOS>BOND\n");

	for (unsigned j=0; j<Cmol->bonds.size(); j++){
		gzprintf(outmol2, "%6d%6.6s%6.6s%5.5s\n", j+1, Cmol->bonds[j][0].c_str(), Cmol->bonds[j][1].c_str(),Cmol->bonds[j][2].c_str());
	}
    gzprintf(outmol2, "\n");
    gzclose(outmol2);
}

void WRITER::writeMol2(Mol2* Cmol, vector<vector<double> >xyz, energy_result_t* result, double rmsd){
    outmol2 = gzopen((this->output_prefix+"_dock.mol2.gz").c_str(), "a");
    gzprintf(outmol2, "\n");
    gzprintf(outmol2, "########## %15.15s: %19.19s\n", "Name", Cmol->molname.c_str());
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "Energy Score", result->total);
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "Elec Score", result->elec);
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "VDW Score", result->vdw);
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "HB Score", result->hb_donor+result->hb_acceptor);
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "Rest Score", result->restraints);
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "Solv Score", (result->rec_solv + result->lig_solv));
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "RMSD", rmsd);
    gzprintf(outmol2, "\n");
    gzprintf(outmol2, "@<TRIPOS>MOLECULE\n");
    gzprintf(outmol2, "%s\n", Cmol->molname.c_str());
    gzprintf(outmol2, "%d %d %d\n", Cmol->N, Cmol->Nbonds, Cmol->Nres);
    gzprintf(outmol2, "SMALL\n");
    gzprintf(outmol2, "USER_CHARGES\n");
    gzprintf(outmol2, "@<TRIPOS>ATOM\n");
    int i=0;
    unsigned resn=0;

    if (int(xyz.size()) != Cmol->N){
        printf("Mismatch in atom number while writting mol2 file. Please check!\n");
        printf("Cmol->N = %d    xyz.size() = %d\n", Cmol->N, int(xyz.size()));
        exit(1);
    }

    while(resn < Cmol->residue_pointer.size()-1){
        while(i < Cmol->residue_pointer[resn+1]-1){
            if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
                gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
            }
            else {
                gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
            }
            i++;
        }
        resn++;
    }
    while(i < Cmol->N){
        if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
            gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
        }
        else{
            gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
        }
        i++;
    }

    gzprintf(outmol2, "@<TRIPOS>BOND\n");

    for (unsigned j=0; j<Cmol->bonds.size(); j++){
        gzprintf(outmol2, "%6d%6.6s%6.6s%5.5s\n", j+1, Cmol->bonds[j][0].c_str(), Cmol->bonds[j][1].c_str(),Cmol->bonds[j][2].c_str());
    }
    gzprintf(outmol2, "\n");
    gzclose(outmol2);
}

void WRITER::writeMol2_Mol_new_xyz(Mol2* Cmol, double energy, double rmsd){
    outmol2 = gzopen((Input->output+"_dock.mol2.gz").c_str(), "a");
	gzprintf(outmol2, "\n");
	gzprintf(outmol2, "########## %15.15s: %19.19s\n", "Name", Cmol->molname.c_str());
	gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "Energy Score", energy);
	gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "RMSD", rmsd);
	gzprintf(outmol2, "\n");
	gzprintf(outmol2, "@<TRIPOS>MOLECULE\n");
	gzprintf(outmol2, "%s\n", Cmol->molname.c_str());
	gzprintf(outmol2, "%d %d %d\n", Cmol->N, Cmol->Nbonds, Cmol->Nres);
	gzprintf(outmol2, "SMALL\n");
	gzprintf(outmol2, "USER_CHARGES\n");
	gzprintf(outmol2, "Energy: %7.3f\n", energy);
	gzprintf(outmol2, "RMSD/OVERLAY: %7.3f\n", rmsd);
	gzprintf(outmol2, "@<TRIPOS>ATOM\n");
	int i=0;
	unsigned resn=0;

	while(resn < Cmol->residue_pointer.size()-1){
		while(i < Cmol->residue_pointer[resn+1]-1){
			if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
                gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), Cmol->new_xyz[i][0], Cmol->new_xyz[i][1], Cmol->new_xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
			}
			else {
                gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), Cmol->new_xyz[i][0], Cmol->new_xyz[i][1], Cmol->new_xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
			}
			i++;
		}
		resn++;
	}
	while(i < Cmol->N){
		if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
            gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f\n", i+1, Cmol->atomnames[i].c_str(), Cmol->new_xyz[i][0], Cmol->new_xyz[i][1], Cmol->new_xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
		}
		else{
            gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f\n", i+1, Cmol->atomnames[i].c_str(), Cmol->new_xyz[i][0], Cmol->new_xyz[i][1], Cmol->new_xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
		}
		i++;
	}

	gzprintf(outmol2, "@<TRIPOS>BOND\n");

	for (unsigned j=0; j<Cmol->bonds.size(); j++){
		gzprintf(outmol2, "%6d%6.6s%6.6s%5.5s\n", j+1, Cmol->bonds[j][0].c_str(), Cmol->bonds[j][1].c_str(),Cmol->bonds[j][2].c_str());
	}
    gzprintf(outmol2, "\n");
    gzclose(outmol2);
}

void WRITER::print_info(char info[98]){
    printf("*%-98.98s*\n", info);
    fprintf(output, "%s\n", info);
    fflush(output);
}

void WRITER::print_info_no_newline(char info[98]){
    printf("*%-98.98s\n", info);
    fprintf(output, "%s\n", info);
    fflush(output);
}

void WRITER::print_line(void){
    printf("****************************************************************************************************\n");
    fprintf(output, "****************************************************************************************************\n");
}

WRITER::~WRITER(){
	fclose(output);
}

void WRITER::print_dock_params(){
    printf("*                                                                                                  *\n");
    printf("* %-30s %-66s*\n", "mode", Input->mode.c_str());
    printf("* %-30.30s %-66.66s*\n", "Receptor File:", Input->rec_mol2.c_str());
    printf("* %-30.30s %-66.66s*\n", "Reference Ligand:", Input->reflig_mol2.c_str());
    printf("* %-30.30s %-66.66s*\n", "Ligand File", Input->lig_mol2.c_str());
    printf("* %-30.30s %-66.66s*\n", "Multifile", Input->multifile.c_str());
    printf("*                                                                                                  *\n");
    printf("* %-30.30s %-66.6f*\n", "VDW softcore term", Input->deltaij6);
    printf("* %-30.30s %-66.6f*\n", "Electrostatic softcore term", Input->deltaij_es6);
    printf("* %-30.30s %-66.6f*\n", "Delta_ij^6", Input->deltaij6);

}


void WRITER::write_pqr(Mol2 *Cmol, string outname){
    gzFile outpdb;
    outpdb = gzopen((outname+".pqr.gz").c_str(), "w");
    int i=0;
    int resn=0;

    while (resn < int(Cmol->residue_pointer.size()-1)){
        while(i < Cmol->residue_pointer[resn+1]-1){
            gzprintf(outpdb, "ATOM %6d%3.3s   %3.3s%2.2s%4d    % 8.3f% 8.3f% 8.3f% 8.4f% 7.4f\n", i+1, Cmol->amberatoms[i].c_str(), Cmol->resnames[resn].c_str(), "A" ,resn+1, Cmol->xyz[i][0], Cmol->xyz[i][1], Cmol->xyz[i][2], Cmol->charges[i], Cmol->radii[i]);
            i++;
        }
        resn++;
    }

    while (i < Cmol->N){
        gzprintf(outpdb, "ATOM %6d%3.3s   %3.3s%2.2s%4d    % 8.3f% 8.3f% 8.3f% 8.4f% 7.4f\n", i+1, Cmol->amberatoms[i].c_str(), Cmol->resnames[resn].c_str(), "A" ,resn+1, Cmol->xyz[i][0], Cmol->xyz[i][1], Cmol->xyz[i][2], Cmol->charges[i], Cmol->radii[i]);
        i++;
    }
    gzprintf(outpdb, "TER\n");
    gzprintf(outpdb, "END\n");
    gzclose(outpdb);
}



#include <boost/python.hpp>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;


//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(wpdboverloads,write_pdb, 5, 6);

BOOST_PYTHON_MODULE(pyWRITER)
{



    void    (WRITER::*wpdb_1)(vector<vector<double> >, int, vector<string>, double, double, string)  = &WRITER::write_pdb;
    void    (WRITER::*wpdb_2)(Mol2*, vector<vector<double> >, double, double, string)   = &WRITER::write_pdb;


    void    (WRITER::*wmol2_1)(Mol2*, vector<vector<double> >, double, double, string)  = &WRITER::writeMol2;
    void    (WRITER::*wmol2_2)(Mol2*, vector<vector<double> >, double, double) = &WRITER::writeMol2;
    void    (WRITER::*wmol2_3)(Mol2*, vector<vector<double> >, energy_result_t*, double) = &WRITER::writeMol2;
    void    (WRITER::*wmol2_4)(Mol2*, vector<vector<double> >, energy_result_t*, double, string)  = &WRITER::writeMol2;



    class_<WRITER>("WRITER", init<string>())
        .def(init<PARSER*>())
        .def(init<string,PARSER*>())
        .def_readwrite("output", &WRITER::output)
        .def_readwrite("outputfile", &WRITER::outputfile)
        .def_readwrite("Cmol", &WRITER::Cmol)
        .def_readwrite("outmol2", &WRITER::outmol2)
        .def_readwrite("Input", &WRITER::Input)



        .def("print_welcome", &WRITER::print_welcome)
        .def("print_params", &WRITER::print_params)
        .def("write_box", &WRITER::write_box)
        .def("write_pdb", wpdb_1)
        .def("write_pdb", wpdb_2)

 //       .def("write_pdb", &WRITER::write_pdb, wpdboverloads())


            
        .def("writeMol2", wmol2_1)
        .def("writeMol2", wmol2_2)
        .def("writeMol2", wmol2_3)
        .def("writeMol2", wmol2_4)


        .def("writeMol2_Mol_new_xyz", &WRITER::writeMol2_Mol_new_xyz)
        .def("print_info", &WRITER::print_info)
        .def("print_info_no_newline", &WRITER::print_info_no_newline)
        .def("print_line", &WRITER::print_line)
        .def("print_dock_params", &WRITER::print_dock_params)
        .def("write_pqr", &WRITER::write_pqr)
    ;  
}
