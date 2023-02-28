
/*
 * PARSER.cpp
 *
 *  Created on: 28/05/2010
 *      Author: Nascimento
 */

#include "PARSER.h"

using namespace std;

PARSER::PARSER(){
	this->flex_lig = false;
	this->flex_rec = false;
	this->dynamic_rate = false;
    this->timeout = 15;
    this->scoring_function = 1;
    this->min_delta = 1.0e-3;
    this->min_tol = 1.0e-4;
    this->min_timeout = 15;
	this->elec_scale = 1.0;
	this->vdw_scale = 1.0;
	this->multifile = "";
    this->overlay_optimizer = "mma";
    this->energy_optimizer = "direct";
    this->deltaij6 = pow(1.75, 6);
    this->deltaij_es6 = pow(1.5, 6);
	this->cushion = 1.5;
	this->sigma = 3.5;
	this->sa_scheme = false;
	this->sa_start_temp = 5000.0;
	this->sa_steps = 500;
	this->sa_mu_t = 1.01;
    this->diel = 2.0;
	this->output = "McLiBELa";
	this->rotation_step = 15.0;
	this->x_dim = 15.0;
	this->y_dim = 15.0;
	this->z_dim = 15.0;
	this->temp = 300.0;
	this->number_steps = 1000;
    this->eq_steps = 0;
	this->dock_no_h = false;
	this->deal_pop_size = 100;
	this->deal_generations = 30;
	this->deal_top_ranked = 20;
	this->deal = false;
    this->generate_conformers = true;
	this->lig_conformers = 10;
	this->conformers_to_evaluate = 1;
    this->conf_search_trials = 10000;
	this->mol2_aa = false;
    this->conformer_min_steps = 1000;
    this->dock_min_tol = 1.0e-4;
	this->sa_mode = false;
	this->dock_mode = false;
	this->eq_mode = false;
    this->mcr_mode = false;
    this->full_search_mode = false;
	this->dock_parallel = false;
    this->parallel_jobs = 1 ;
	this->write_mol2 = true;
    this->grid_spacing = 0.4;
	this->use_grids = false;
    this->search_box_x = 6.0;
    this->search_box_y = 6.0;
    this->search_box_z = 6.0;
	this->load_grid_from_file = false;
	this->write_grids = false;
	this->show_rmsd = false;
	this->sort_by_energy = false;
    this->dielectric_model = "r";
    this->only_score = false;
    this->solvation_alpha = 0.10;
    this->solvation_beta = -0.005;
    this->seed = (rand()/(RAND_MAX + 1.0));
    this->ligsim = false;
    this->mc_stride = 1;
    this->mcr_size = 1;
    this->bi = 2.0;
    this->torsion_step = 10.0;
    this->sample_torsions = false;
    this->verbose = false;
    this->entropy_rotation_bins = 360;
    this->entropy_translation_bins = this->x_dim*2;
    this->ligand_energy_model = "MMFF94";
    this->atomic_model_ff = "GAFF";
    this->pbsa_grid = "";
    this->use_pbsa = false;
    this->use_delphi = false;
    this->delphi_gsize = 75;
    this->mc_full_flex = false;
    this->compute_rotation_entropy = false;
    this->max_atom_displacement = 0.00005;
    this->use_writeMol2_score_cutoff = false;
    this->use_writeMol2_energy_cutoff = false;
    this->writeMol2_score_cutoff = 0.75;
    this->writeMol2_energy_cutoff = 10.0;
    this->use_GW_LJ6 = false;
    this->use_GW_LJ12 = false;
    this->use_GW_Coulomb = false;
    this->use_Erestraints = false;
    this->restraints_weight = 0.0;
    this->translation_step = 0.2;
    this->scale_vdw_energy = 1.0;
    this->scale_elec_energy = 1.0;
    this->overlay_cutoff = 0.0;
    this->use_overlay_cutoff = false;
    this->use_score_optimization = false;
    this->use_only_binding_energy = false;
}

void PARSER::comparing (string param, ifstream &input) {
	if (param == "nsteps") {
		input >> PARSER::number_steps;
	}
	else if (param == "temperature") {
		input >> PARSER::temp;
	}
	else if (param=="sa_start_temp"){
			input >> PARSER::sa_start_temp;
	}
	else if (param == "mode"){
		input >> tmp;
		string mode;
		size_t plus = 0;
		while (plus != string::npos){
			plus = tmp.find("+");
			mode = tmp.substr(0,plus);
			if (mode == "dock" or mode == "DOCK" or mode == "Dock"){
				this->dock_mode = true;
			}
            else if(mode == "full_search"){
                this->full_search_mode = true;
            }
			else if (mode == "sa" or mode == "SA" or mode == "Sa"){
				this->sa_mode = true;
			}
            else if (mode == "equilibrium" or mode == "eq" or mode == "EQ" or mode == "mc" or mode == "MC"){
				this->eq_mode = true;
			}
            else if (mode == "mcr" or mode == "MCR"){
                this->mcr_mode = true;
            }
			tmp=tmp.substr(plus+1, tmp.size());
		}
	}
	else if (param=="sa_steps"){
		input >> PARSER::sa_steps;
	}
	else if (param == "sa_mu_t"){
		input >> this->sa_mu_t;
	}
	else if (param == "cushion") {
		input >> PARSER::cushion;
	}
	else if (param == "deltaij6") {
		input >> PARSER::deltaij6;
    }
    else if (param == "deltaij"){
        double dij;
        input >> dij;
        this->deltaij6 = pow(dij, 6);
    }
	else if (param == "deltaij_es6"){
		input >> PARSER::deltaij_es6;
        this->deltaij_es3 = sqrt(this->deltaij_es6);
	}
    else if (param == "deltaij_es"){
        double dijes;
        input >> dijes;
        this->deltaij_es6 = pow(dijes, 6);
        this->deltaij_es3 = pow(dijes, 3);
    }
    else if (param == "diel"){
		input >> PARSER::diel;
	}
	else if (param == "sigma"){
		input >> PARSER::sigma;
	}
	else if (param == "rec_mol2"){
		input >> this->rec_mol2;
	}
	else if (param == "lig_mol2"){
		input >> this->lig_mol2;
	}
	else if(param == "mol2_aa"){
		input >> this->sa;
		if (sa == "yes" or sa == "YES" or sa == "Yes"){
			this->mol2_aa = true;
		}
		else {
			this->mol2_aa = false;
		}
	}
	else if (param == "output_prefix"){
			input >> PARSER::output;
	}
	else if (param == "rotation_step"){
				input >> PARSER::rotation_step;
	}
    else if (param == "sample_torsions"){
        input >> tmp;
        if (tmp == "yes" or tmp == "YES" or tmp == "Yes"){
            this->sample_torsions = true;
        }
    }
    else if (param == "torsion_step"){
        input >> this->torsion_step;
    }
	else if (param == "timeout"){
					input >> PARSER::timeout;
	}
	else if (param == "grid_box"){
		input >> PARSER::x_dim >> PARSER::y_dim >> PARSER::z_dim;
	}
    else if (param == "lig_traj"){
        input >> PARSER::lig_traj;
        PARSER::flex_lig = true;
    }
    else if (param == "rec_traj"){
        input >> PARSER::rec_traj;
        PARSER::flex_rec = true;
	}
	else if( param == "scoring_function"){
		input >> PARSER::scoring_function;
	}
	else if ( param == "dyn_rate_inf_limit"){
		this->dynamic_rate = true;
		input >> this->dyn_rate_inf_limit;
	}
	else if (param == "dyn_rate_steps"){
		input >> this->dyn_rate_steps;
	}
	else if (param == "dyn_rate_sup_limit"){
		this->dynamic_rate = true;
		input >> this->dyn_rate_sup_limit;
	}
	else if (param == "elec_scale"){
		input >> this->elec_scale;
	}
	else if (param == "vdw_scale"){
		input >> this->vdw_scale;
	}
	else if (param == "reflig_mol2"){
		input >> this->reflig_mol2;
	}
	else if (param == "minimization_tolerance"){
		input >> this->min_tol;
	}
	else if (param == "minimization_delta"){
		input >> this->min_delta;
	}
	else if (param == "minimization_timeout"){
		input >> this->min_timeout;
	}
	else if (param == "multifile"){
		input >> this->multifile;
	}
	else if (param == "overlay_optimizer"){
		input >> this->overlay_optimizer;
	}
	else if (param == "energy_optimizer"){
		input >> this->energy_optimizer;
	}
	else if (param == "ignore_h"){
		input >> this->tmp;
		if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
			this->dock_no_h = true;
		}
		else {
			this->dock_no_h = false;
		}
	}
	else if (param == "deal"){
		input >> this->tmp;
		if (tmp == "yes" or tmp == "YES" or tmp == "Yes"){
			this->deal = true;
		}
		else {
			this->deal = false;
		}
	}
	else if (param == "deal_population_size"){
		input >> this->deal_pop_size;
	}
	else if (param == "deal_generations"){
		input >> this->deal_generations;
	}
	else if (param == "deal_top_ranked"){
		input >> this->deal_top_ranked;
	}

	else if (param == "generate_conformers"){
		input >> tmp;
		if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
			this->generate_conformers = true;
		}
		else {
			this->generate_conformers = false;
		}
        this->flex_lig = true;
	}

	else if (param == "number_of_conformers"){
		input >> this->lig_conformers;
	}

    else if (param == "conformers_to_rank"){
        input >> this->conformers_to_evaluate;
    }

    else if (param == "conformer_min_steps"){
        input >> this->conformer_min_steps;
	}

    else if (param == "conf_search_trials"){
        input >> this->conf_search_trials;
    }

    else if (param == "dock_min_tol"){
		input >> this->dock_min_tol;
	}
	else if (param == "dock_parallel"){
		input >> tmp;
		if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
			this->dock_parallel = true;
		}
	}
	else if (param == "parallel_jobs"){
		input >> (this->parallel_jobs);
	}
	else if (param == "write_mol2"){
		input >> tmp;
		if ( tmp == "no" or tmp == "No" or tmp == "NO"){
			this->write_mol2 = false;
		}
	}
	else if (param == "use_grids"){
		input >> tmp;
		if (tmp == "Yes" or tmp == "Yes" or tmp == "yes"){
			this->use_grids = true;
		}
	}
	else if (param == "grid_spacing"){
		input >> this->grid_spacing;
	}
	else if (param == "search_box"){
		input >> this->search_box_x >> this->search_box_y >> this->search_box_z;
	}
	else if (param == "load_grids"){
		this->load_grid_from_file = true;
		input >> this->grid_prefix;
	}
	else if (param == "write_grids"){
		this->write_grids = true;
		input >> this->grid_prefix;
	}
	else if (param == "compute_rmsd"){
		input >> tmp;
		if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
			this->show_rmsd = true;
		}
	}
	else if (param == "sort_by_energy"){
		input >> tmp;
		if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
			this->sort_by_energy = true;
		}
	}
    else if (param == "dielectric_model"){
        input >> this->dielectric_model ;
    }
    else if (param == "only_score"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "Yes"){
            this->only_score = true;
        }
    }
    else if (param == "solvation_alpha"){
        input >> this->solvation_alpha;
    }
    else if (param == "solvation_beta"){
        input >> this->solvation_beta;
    }
    else if (param == "seed"){
        input >> this->seed;
    }
    else if (param == "ligand_simulation"){
        input >> tmp;
        if(tmp == "yes" or tmp == "1" or tmp == "y")
            this->ligsim = true;
    }
    else if (param == "equilibration_steps"){
        input >> this->eq_steps;
    }
    else if (param == "mc_stride"){
        input >> this->mc_stride;
    }
    else if (param == "mcr_size"){
        input >> this->mcr_size;
    }
    else if (param == "verbose"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            this->verbose = true;
        }
    }
    else if (param == "mcr_coefficients"){
        double bi;
        for (int i=0; i<mcr_size; i++){
            input >> bi;
            mcr_coefficients.push_back(bi);
        }
    }
    else if (param == "entropy_rotation_bins"){
        input >> this->entropy_rotation_bins;
    }

    else if (param == "entropy_translation_bins"){
        input >> this->entropy_translation_bins;
    }
    else if (param == "ligand_energy_model"){
        input >> this->ligand_energy_model;
    }
    else if (param == "atomic_model_ff"){
        input >> this->atomic_model_ff;
    }
    else if (param == "pbsa_grid"){
        input >> this->pbsa_grid;
    }
    else if (param == "use_pbsa"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            this->use_pbsa = true;
        }
    }
    else if (param == "use_delphi"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            this->use_delphi = true;
        }
    }
    else if (param == "delphi_grid"){
        input >> this->delphi_grid;
    }
    else if (param == "delphi_gsize"){
        input >> this->delphi_gsize;
    }
    else if (param == "delphi_cube_grid"){
        input >> this->delphi_cube_grid;
    }
    else if (param == "mc_full_flex"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            mc_full_flex = true;
        }
    }
    else if (param == "compute_rotation_entropy"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            this->compute_rotation_entropy = true;
        }
    }
    else if (param == "max_atom_displacement"){
        input >> this->max_atom_displacement;
    }
    else if (param == "use_writeMol2_score_cutoff"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            use_writeMol2_score_cutoff = true;
        }
    }
    else if (param == "use_writeMol2_energy_cutoff"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            use_writeMol2_energy_cutoff = true;
        }
    }
    else if (param == "writeMol2_score_cutoff"){
        input >> this->writeMol2_score_cutoff;
    }
    else if (param == "writeMol2_energy_cutoff"){
        input >> this->writeMol2_energy_cutoff;
    }
    else if (param == "LJ_sigma"){
        input >> this->LJ_sigma;
    }
    else if (param == "Coulomb_sigma"){
        input >> this->coulomb_sigma;
    }
    else if (param == "use_GW_LJ6"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            this->use_GW_LJ6 = true;
        }
    }
    else if (param == "use_GW_LJ12"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            this->use_GW_LJ12 = true;
        }
    }
    else if (param == "use_GW_Coulomb"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            this->use_GW_Coulomb = true;
        }
    }
    else if (param == "use_docking_restraints"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            this->use_Erestraints = true;
        }
    }
    else if (param == "restraints_weight"){
        input >> this->restraints_weight;
    }
    else if (param == "translation_step"){
        input >> this->translation_step;
    }
    else if (param == "scale_vdw_energy"){
        input >> this->scale_vdw_energy;
    }
    else if (param == "scale_elec_energy"){
        input >> this->scale_elec_energy;
    }
    else if (param == "use_overlay_cutoff"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            this->use_overlay_cutoff = true;
        }
    }
    else if (param == "overlay_cutoff"){
        input >> this->overlay_cutoff;
    }
    else if (param == "use_score_optimization"){
        input >> tmp;
        if (tmp == "yes" or tmp == "YES" or tmp == "Yes"){
            this->use_score_optimization = true;
        }
    }
    else if (param == "use_only_binding_energy"){
        input >> tmp;
        if (tmp == "yes" or tmp == "YES" or tmp == "Yes"){
            this->use_only_binding_energy = true;
        }
    }
	else {
		cout << "Unknown parameter: " << param << endl;
		exit(1);
	}
}


void PARSER::set_parameters(char* arg){

	ifstream input(arg);
	char line[256];
	while (!input.eof()){
		input >> param;
		if ((param.substr(0,1) == "#") or (param.substr(0,2) == "//")){
			input.getline(line, 256);
		}
		else {
			PARSER::comparing (param, input);
		}
	}
	this->check_parameters_sanity();
}

void PARSER::check_parameters_sanity(void){
	if (this->conformers_to_evaluate > this->lig_conformers){
		this->conformers_to_evaluate = this->lig_conformers;
	}
}


#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(pyPARSER)
{
    class_<PARSER>("PARSER", init<void>())
        .def_readwrite("param", & PARSER::param)
        .def_readwrite("number_steps", & PARSER::number_steps)
        .def_readwrite("temp", & PARSER::temp)
        .def_readwrite("sweep_steps", & PARSER::sweep_steps)
        .def_readwrite("eq_steps", & PARSER::eq_steps)
        .def_readwrite("cushion", & PARSER::cushion)
        .def_readwrite("deltaij6", & PARSER::deltaij6)
        .def_readwrite("deltaij_es6", & PARSER::deltaij_es6)
        .def_readwrite("deltaij_es3", & PARSER::deltaij_es3)
        .def_readwrite("diel", & PARSER::diel)
        .def_readwrite("sigma", & PARSER::sigma)
        .def_readwrite("output", & PARSER::output)
        .def_readwrite("rotation_step", & PARSER::rotation_step)
        .def_readwrite("sample_torsions", & PARSER::sample_torsions)
        .def_readwrite("torsion_step", & PARSER::torsion_step)
        .def_readwrite("x_dim", & PARSER::x_dim)
        .def_readwrite("y_dim", & PARSER::y_dim)
        .def_readwrite("z_dim", & PARSER::z_dim)
        .def_readwrite("flex_lig", & PARSER::flex_lig)
        .def_readwrite("timeout", & PARSER::timeout)
        .def_readwrite("tmp", & PARSER::tmp)
        .def_readwrite("scoring_function", & PARSER::scoring_function)
        .def_readwrite("rec_mol2", & PARSER::rec_mol2)
        .def_readwrite("lig_mol2", & PARSER::lig_mol2)
        .def_readwrite("mol2_aa", & PARSER::mol2_aa)
        .def_readwrite("min_delta", & PARSER::min_delta)
        .def_readwrite("min_tol", & PARSER::min_tol)
        .def_readwrite("min_timeout", & PARSER::min_timeout)
        .def_readwrite("mode", & PARSER::mode)
        .def_readwrite("elec_scale", & PARSER::elec_scale)
        .def_readwrite("vdw_scale", & PARSER::vdw_scale)
        .def_readwrite("reflig_mol2", & PARSER::reflig_mol2)
        .def_readwrite("multifile", & PARSER::multifile)
        .def_readwrite("overlay_optimizer", & PARSER::overlay_optimizer)
        .def_readwrite("energy_optimizer", & PARSER::energy_optimizer)
        .def_readwrite("dock_no_h", & PARSER::dock_no_h)
        .def_readwrite("generate_conformers", & PARSER::generate_conformers)
        .def_readwrite("lig_conformers", & PARSER::lig_conformers)
        .def_readwrite("conformers_to_evaluate", & PARSER::conformers_to_evaluate)
        .def_readwrite("conformer_min_steps", & PARSER::conformer_min_steps)
        .def_readwrite("conf_search_trials", & PARSER::conf_search_trials)
        .def_readwrite("dock_min_tol", & PARSER::dock_min_tol)
        .def_readwrite("eq_mode", & PARSER::eq_mode)
        .def_readwrite("sa_mode", & PARSER::sa_mode)
        .def_readwrite("dock_mode", & PARSER::dock_mode)
        .def_readwrite("dock_parallel", & PARSER::dock_parallel)
        .def_readwrite("parallel_jobs", & PARSER::parallel_jobs)
        .def_readwrite("write_mol2", & PARSER::write_mol2)
        .def_readwrite("sa_mu_t", & PARSER::sa_mu_t)
        .def_readwrite("grid_spacing", & PARSER::grid_spacing)
        .def_readwrite("use_grids", & PARSER::use_grids)
        .def_readwrite("search_box_x", & PARSER::search_box_x)
        .def_readwrite("search_box_y", & PARSER::search_box_y)
        .def_readwrite("search_box_z", & PARSER::search_box_z)
        .def_readwrite("load_grid_from_file", & PARSER::load_grid_from_file)
        .def_readwrite("grid_prefix", & PARSER::grid_prefix)
        .def_readwrite("write_grids", & PARSER::write_grids)
        .def_readwrite("show_rmsd", & PARSER::show_rmsd)
        .def_readwrite("sort_by_energy", & PARSER::sort_by_energy)
        .def_readwrite("only_score", & PARSER::only_score)
        .def_readwrite("solvation_alpha", & PARSER::solvation_alpha)
        .def_readwrite("solvation_beta", & PARSER::solvation_beta)
        .def_readwrite("seed", & PARSER::seed)
        .def_readwrite("ligsim", & PARSER::ligsim)
        .def_readwrite("mc_stride", & PARSER::mc_stride)
        .def_readwrite("translation_step", & PARSER::translation_step)
        .def_readwrite("full_search_mode", & PARSER::full_search_mode)
        .def_readwrite("scale_vdw_energy", & PARSER::scale_vdw_energy)
        .def_readwrite("scale_elec_energy", & PARSER::scale_elec_energy)
        .def_readwrite("use_score_optimization", & PARSER::use_score_optimization)
        .def_readwrite("use_only_binding_energy", & PARSER::use_only_binding_energy)
        .def_readwrite("dielectric_model", & PARSER::dielectric_model)
        .def_readwrite("mcr_size", & PARSER::mcr_size)
        .def_readwrite("mcr_coefficients", & PARSER::mcr_coefficients)
        .def_readwrite("mcr_mode", & PARSER::mcr_mode)
        .def_readwrite("bi", & PARSER::bi)
        .def_readwrite("verbose", & PARSER::verbose)
        .def_readwrite("entropy_rotation_bins", & PARSER::entropy_rotation_bins)
        .def_readwrite("entropy_translation_bins", & PARSER::entropy_translation_bins)
        .def_readwrite("ligand_energy_model", & PARSER::ligand_energy_model)
        .def_readwrite("atomic_model_ff", & PARSER::atomic_model_ff)
        .def_readwrite("delphi_grid", & PARSER::delphi_grid)
        .def_readwrite("use_delphi", & PARSER::use_delphi)
        .def_readwrite("delphi_gsize", & PARSER::delphi_gsize)
        .def_readwrite("mc_full_flex", & PARSER::mc_full_flex)
        .def_readwrite("compute_rotation_entropy", & PARSER::compute_rotation_entropy)
        .def_readwrite("max_atom_displacement", & PARSER::max_atom_displacement)
        .def_readwrite("use_writeMol2_score_cutoff", & PARSER::use_writeMol2_score_cutoff)
        .def_readwrite("use_writeMol2_energy_cutoff", & PARSER::use_writeMol2_energy_cutoff)
        .def_readwrite("writeMol2_score_cutoff", & PARSER::writeMol2_score_cutoff)
        .def_readwrite("writeMol2_energy_cutoff", & PARSER::writeMol2_energy_cutoff)
        .def_readwrite("use_Erestraints", & PARSER::use_Erestraints)
        .def_readwrite("restraints_weight", & PARSER::restraints_weight)
        .def_readwrite("use_overlay_cutoff", & PARSER::use_overlay_cutoff)
        .def_readwrite("overlay_cutoff", & PARSER::overlay_cutoff)

        .def("set_parameters", &PARSER::set_parameters)
        .def("comparing", &PARSER::comparing)
        .def("check_parameters_sanity", &PARSER::check_parameters_sanity)
    ;
}


