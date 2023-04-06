/*
 * Docker.cpp
 *
 *  Created on: 24/03/2012
 *      Author: Nascimento
 */

#include "pyDocker.h"

Docker::Docker(WRITER* _Writer){
    this->Writer = _Writer;
}

#ifdef HAS_GUI
Docker::Docker(QtWriter* _QWriter){
    this->QWriter = _QWriter;
}
#endif

void Docker::run(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, unsigned counter) {

    if (!Input->generate_conformers){

        COORD_MC* Coord = new COORD_MC;
        int overlay_status;
        double overlay_fmax;
        vector<double> com_lig = Coord->compute_com(Lig);
        bool match=false;

/*
* Shifting the ligand to match its center of mass with the reference
* ligand center of mass.
*/

        Lig->xyz = Coord->translate(Lig->xyz, Lig->N, com[0]-com_lig[0], com[1]-com_lig[1], com[2]-com_lig[2]);

/*
 *
 * Doing a pre-alignment of the longest axis of the ligands
 *
 */

        Lig->find_longest_axis();
        RefLig->find_longest_axis();

        Optimizer* Opt0 = new Optimizer(Rec, RefLig, Input);
        Optimizer::opt_result_t* opt_result0 = new Optimizer::opt_result_t;
        opt_result0->energy_result = new energy_result_t;

        Opt0->pre_align(Lig, opt_result0);
        Lig->xyz = opt_result0->optimized_xyz;

        delete opt_result0;
        delete Opt0;

/*
 * End of pre-alignment
 */


        Optimizer* Opt = new Optimizer(Rec, RefLig, Input);
        Optimizer::opt_result_t* opt_result = new Optimizer::opt_result_t;
        opt_result->energy_result = new energy_result_t;

        // Optimizing overlay....

        if (! this->minimize_overlay(Input, Opt, Lig, opt_result)){
            sprintf(info, "Overlay optimizer %s is not defined. Exiting...\n", Input->overlay_optimizer.c_str());
        }

        Gaussian* Gauss = new Gaussian;
        double t1, t2, t3, si;
        t1 = Gauss->compute_shape_and_charge_density(Input, RefLig, Lig, RefLig->xyz, opt_result->optimized_xyz);
        t2 = Gauss->compute_shape_and_charge_density(Input, RefLig, RefLig, RefLig->xyz, RefLig->xyz);
        t3 = Gauss->compute_shape_and_charge_density(Input, Lig, Lig, opt_result->optimized_xyz,opt_result->optimized_xyz);
        si = (2*t1) / (t2+t3);
        delete Gauss;

        //Copying new coordinates.
        Lig->xyz = opt_result->optimized_xyz;
        overlay_status = opt_result->optimization_status;
        overlay_fmax = si;

        //        delete opt_result->energy_result;
        delete opt_result;


        //Optimizing Energy...
        Optimizer::opt_result_t* opt_result2 = new Optimizer::opt_result_t;
        opt_result2->energy_result = new energy_result_t;

        if (Input->deal){
            sprintf(info, "Entering Deal....");
            this->print_info(info);
#ifdef HAS_GUI
            Deal* DEAl = new Deal(Rec, Lig, Input);
            sprintf(info, "%-20.20s %-20.20s %-7.3e % -7.3f kcal/mol", Lig->molname.c_str(), Lig->resnames[0].c_str(), fo, DEAl->energies[0]);
            this->print_info(info);
            QWriter->writeMol2(Lig, DEAl->final_coords, DEAl->energies[0], fo);
            delete DEAl;
#endif
        }

        else {
            if (si >= Input->overlay_cutoff){
                match=true;
                if (! this->minimize_energy(Input, Opt, Rec, Lig, RefLig, opt_result2)){
                    sprintf(info, "Energy optimizer %s is not defined. Exiting...\n", Input->energy_optimizer.c_str());
                }

                Gaussian* Gauss2 = new Gaussian;
                t1 = Gauss2->compute_shape_and_charge_density(Input, RefLig, Lig, RefLig->xyz, opt_result2->optimized_xyz);
                t2 = Gauss2->compute_shape_and_charge_density(Input, RefLig, RefLig, RefLig->xyz, RefLig->xyz);
                t3 = Gauss2->compute_shape_and_charge_density(Input, Lig, Lig, opt_result2->optimized_xyz,opt_result2->optimized_xyz);
                si = (2*t1) / (t2+t3);
                delete Gauss2;

                Lig->xyz = opt_result2->optimized_xyz;
            }

            if (match){
                sprintf(info, "%5d %-12.12s %-4.4s %-10.3e %-8.3g %-8.3g %-8.3g %-8.2f %-8.3g %-3d %-2d %-2d %.2f", counter, Lig->molname.c_str(), Lig->resnames[0].c_str(), overlay_fmax,
                        opt_result2->energy_result->elec, opt_result2->energy_result->vdw, opt_result2->energy_result->rec_solv+opt_result2->energy_result->lig_solv,
                        opt_result2->energy_result->hb_donor+opt_result2->energy_result->hb_acceptor, opt_result2->energy_result->total, 0, overlay_status, opt_result2->optimization_status, si);
            }
            else{
                sprintf(info, "%5d %-12.12s %-4.4s %-10.3e  %-8.3g %-8.3g %-8.3g %-8.2f %-8.3g %3d %2d %2d %.2f", counter, Lig->molname.c_str(), Lig->resnames[0].c_str(), overlay_fmax,
                        0.0, 0.0, 0.0, 0.0, 0.0, -1, overlay_status, 0.0, si);
            }

            this->print_info(info);

            if (Input->write_mol2 and match){

#pragma omp critical
                {
                    this->write_mol2(Lig, opt_result2->optimized_xyz, opt_result2->energy_result, overlay_fmax);
                }
            }
        }
        delete Coord;
        delete opt_result2;
        delete Opt;
        //        Lig->xyz = opt_result2->optimized_xyz;
    }
    else {
        this->Dock_conformers(Rec, Lig, RefLig, com, Input, counter);
    }
}

void Docker::run(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, Grid* Grids, unsigned counter) {

    if (!Input->generate_conformers){

        COORD_MC* Coord = new COORD_MC;
        int overlay_status;
        double overlay_fmax;
        vector<double> com_lig = Coord->compute_com(Lig);
        bool match=false;

/*
* Shifting the ligand to match its center of mass with the reference
* ligand center of mass.
*/

        Lig->xyz = Coord->translate(Lig->xyz, Lig->N, com[0]-com_lig[0], com[1]-com_lig[1], com[2]-com_lig[2]);

/*
*
* Doing a pre-alignment of the longest axis of the ligands
*
*/

        Lig->find_longest_axis();
        RefLig->find_longest_axis();

        Optimizer* Opt0 = new Optimizer(Rec, RefLig, Input);
        Optimizer::opt_result_t* opt_result0 = new Optimizer::opt_result_t;
        opt_result0->energy_result = new energy_result_t;

        Opt0->pre_align(Lig, opt_result0);
        Lig->xyz = opt_result0->optimized_xyz;

        delete opt_result0;
        delete Opt0;

/*
* End of pre-alignment
*/

        Optimizer* Opt = new Optimizer(Rec, RefLig, Input, Grids);
        Optimizer::opt_result_t* opt_result = new Optimizer::opt_result_t;
        opt_result->energy_result = new energy_result_t;

        // Optimizing overlay....

        if (! this->minimize_overlay(Input, Opt, Lig, opt_result)){
            sprintf(info, "Overlay optimizer %s is not defined. Exiting...\n", Input->overlay_optimizer.c_str());
        }

        Gaussian* Gauss = new Gaussian;
        double t1, t2, t3, si;
        t1 = Gauss->compute_shape_and_charge_density(Input, RefLig, Lig, RefLig->xyz, opt_result->optimized_xyz);
        t2 = Gauss->compute_shape_and_charge_density(Input, RefLig, RefLig, RefLig->xyz, RefLig->xyz);
        t3 = Gauss->compute_shape_and_charge_density(Input, Lig, Lig, opt_result->optimized_xyz,opt_result->optimized_xyz);
        si = (2*t1) / (t2+t3);
        delete Gauss;


        //Copying new coordinates.
        Lig->xyz = opt_result->optimized_xyz;
        overlay_status = opt_result->optimization_status;
        overlay_fmax = si;

        delete opt_result;

        //Optimizing Energy...

        Optimizer::opt_result_t* opt_result2 = new Optimizer::opt_result_t;
        opt_result2->energy_result = new energy_result_t;

#ifdef DEBUG
        sprintf(info, "Starting Binding Energy Minimization...");
        this->print_info(info);
#endif

        if (Input->deal){
            sprintf(info, "Entering Deal....");
            this->print_info(info);
#ifdef HAS_GUI
            Deal* DEAl = new Deal(Rec, Lig, Input);
            sprintf(info, "%-20.20s %-20.20s %-7.3e % -7.3f kcal/mol", Lig->molname.c_str(), Lig->resnames[0].c_str(), fo, DEAl->energies[0]);
            this->print_info(info);
            Writer->writeMol2(Lig, DEAl->final_coords, DEAl->energies[0], fo);
            delete DEAl;
#endif
        }

        else {
            if (si >= Input->overlay_cutoff){
                match=true;
                if (! this->minimize_energy(Input, Opt, Rec, Lig, RefLig, opt_result2)){
                    sprintf(info, "Energy optimizer %s is not defined. Exiting...\n", Input->energy_optimizer.c_str());
                }

                Gaussian* Gauss2 = new Gaussian;

                t1 = Gauss2->compute_shape_and_charge_density(Input, RefLig, Lig, RefLig->xyz, opt_result2->optimized_xyz);
                t2 = Gauss2->compute_shape_and_charge_density(Input, RefLig, RefLig, RefLig->xyz, RefLig->xyz);
                t3 = Gauss2->compute_shape_and_charge_density(Input, Lig, Lig, opt_result2->optimized_xyz,opt_result2->optimized_xyz);

                si = (2*t1) / (t2+t3);
                delete Gauss2;

                Lig->xyz = opt_result2->optimized_xyz;
            }

            if (match){
                sprintf(info, "%5d %-12.12s %-4.4s %-10.3e  %-8.3g %-8.3g %-8.3g %-8.2f %-8.3g %3d %2d %2d %.2f", counter, Lig->molname.c_str(), Lig->resnames[0].c_str(), overlay_fmax,
                        opt_result2->energy_result->elec, opt_result2->energy_result->vdw, opt_result2->energy_result->rec_solv+opt_result2->energy_result->lig_solv,
                        opt_result2->energy_result->hb_donor+opt_result2->energy_result->hb_acceptor, opt_result2->energy_result->total, 0, overlay_status, opt_result2->optimization_status, si);
            }
            else{
                sprintf(info, "%5d %-12.12s %-4.4s %-10.3e  %-8.3g %-8.3g %-8.3g %-8.2f %-8.3g %3d %2d %2d %.2f", counter, Lig->molname.c_str(), Lig->resnames[0].c_str(), overlay_fmax,
                        0.0, 0.0, 0.0, 0.0, 0.0, -1, overlay_status, 0.0, si);
            }
            this->print_info(info);

            if (Input->write_mol2){
                bool will_write = true;
                if (Input->use_writeMol2_score_cutoff){
                    if (si < Input->writeMol2_score_cutoff){
                        will_write = false;
                    }
                }
                if (Input->use_writeMol2_energy_cutoff){
                    if (opt_result2->energy_result->total > Input->writeMol2_energy_cutoff){
                        will_write = false;
                    }
                }
                if (will_write and match){
                    this->write_mol2(Lig, opt_result2->optimized_xyz, opt_result2->energy_result, overlay_fmax);
                }
            }
        }
        delete Coord;
        delete opt_result2;
        delete Opt;
    }
    else {
        this->Dock_conformers(Rec, Lig, RefLig, com, Input, Grids, counter);
    }
}

Docker::~Docker() {
}

void  Docker::Dock_conformers(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, unsigned counter){
    COORD_MC* Coord = new COORD_MC;
    double best_ene=0.00;
    int best_conf=0, energy_status=0;
    vector<int> overlay_status;
    vector<vector<vector<double> > > new_mcoords;
    vector<double> energies;
    vector<double> overlays;
    double si;
    bool match=false;
    Gaussian* Gauss = new Gaussian;
    double t1, t2, t3;
    t2 = Gauss->compute_shape_and_charge_density(Input, RefLig, RefLig, RefLig->xyz, RefLig->xyz);

    for (unsigned i=0; i<Lig->mcoords.size(); i++){
        vector<double> com_lig = Coord->compute_com(Lig->mcoords[i], Lig);
        if (int(Lig->mcoords[i].size()) == Lig->N){

            Lig->xyz = Coord->translate(Lig->mcoords[i], Lig->N, com[0]-com_lig[0], com[1]-com_lig[1], com[2]-com_lig[2]);

/*
*
* Doing a pre-alignment of the longest axis of the ligands
*
*/

            Lig->find_longest_axis();
            RefLig->find_longest_axis();

            Optimizer* Opt0 = new Optimizer(Rec, RefLig, Input);
            Optimizer::opt_result_t* opt_result0 = new Optimizer::opt_result_t;
            opt_result0->energy_result = new energy_result_t;

            Opt0->pre_align(Lig, opt_result0);
            Lig->xyz = opt_result0->optimized_xyz;

            delete opt_result0;
            delete Opt0;

/*
* End of pre-alignment
*/


            Optimizer* Opt = new Optimizer(Rec, RefLig, Input);
            Optimizer::opt_result_t* opt_result = new Optimizer::opt_result_t;
            opt_result->energy_result = new energy_result_t;

            if (! this->minimize_overlay(Input, Opt, Lig, opt_result)){
                sprintf(info, "Overlay optimizer %s is not defined. Exiting...\n", Input->overlay_optimizer.c_str());
                this->print_info(info);
            }

            t1 = Gauss->compute_shape_and_charge_density(Input, RefLig, Lig, RefLig->xyz, opt_result->optimized_xyz);
            t3 = Gauss->compute_shape_and_charge_density(Input, Lig, Lig, opt_result->optimized_xyz, opt_result->optimized_xyz);
            si = (2*t1) / (t2+t3);

            overlays.push_back(si);
            energies.push_back(opt_result->energy_result->total);
            new_mcoords.push_back(opt_result->optimized_xyz);
            overlay_status.push_back(opt_result->optimization_status);
            delete Opt;
            delete opt_result;
        }
    }
    delete Gauss;

    vector<unsigned> index;

    /*
     * Here the user can choose if the best overlay is the one that maximizes the overlap of the docked ligand
     * with the reference ligand (default) or the on which results in the best binding energy, evaluated after
     * the overlap of the ligands. To choose the binding energy as the parameter, the input "sort_by_energy" must
     * bet set to "yes" in the input file.
     * ASN, Feb/2013.
     */

    if (Input->sort_by_energy){
        index = this->sort_vector(energies);
    }
    else {
        index = this->sort_vector_inv(overlays);
    }

    best_ene= 999999999.;
    vector<vector<double> > new_xyz;
    energy_status=0;

    new_xyz = new_mcoords[index[0]];
    energy_result_t* best_energy_t = new energy_result_t;

    for (int i=0; i<Input->conformers_to_evaluate; i++){
        if (int(new_mcoords[index[i]].size()) == Lig->N){

            Lig->xyz = new_mcoords[index[i]];

            if (overlays[index[i]] >= Input->overlay_cutoff){
                match=true;

                Optimizer* Opt = new Optimizer(Rec, RefLig, Input);

                Optimizer::opt_result_t* opt_result = new Optimizer::opt_result_t;
                opt_result->energy_result = new energy_result_t;

                if (! this->minimize_energy(Input, Opt, Rec, Lig, RefLig, opt_result)){
                    sprintf(info, "Energy optimizer %s is not defined. Exiting...\n", Input->energy_optimizer.c_str());
                    this->print_info(info);
                }

                if (opt_result->energy_result->total < best_ene){
                    best_ene = opt_result->energy_result->total;
                    best_energy_t = opt_result->energy_result;
                    new_xyz = opt_result->optimized_xyz;
                    best_conf=index[i];
                    energy_status = opt_result->optimization_status;
                }
                delete Opt;
                delete opt_result;
            }
        }
    }
    Gaussian* Gauss2 = new Gaussian;
    t1 = Gauss2->compute_shape_and_charge_density(Input, RefLig, Lig, RefLig->xyz, new_xyz);
    //    t2 = Gauss2->compute_shape_and_charge_density(Input, RefLig, RefLig, RefLig->xyz, RefLig->xyz);
    t3 = Gauss2->compute_shape_and_charge_density(Input, Lig, Lig, new_xyz, new_xyz);

    si = (2*t1) / (t2+t3);
    delete Gauss2;

    Lig->xyz = new_xyz;

    if (match){
        sprintf(info, "%5d %-12.12s %-4.4s %-10.3g  %-8.3g %-8.3g %-8.3g %-8.2f %-8.3g %3d %2d %2d %.2f", counter, Lig->molname.c_str(), Lig->resnames[0].c_str(), overlays[best_conf],
                best_energy_t->elec, best_energy_t->vdw, best_energy_t->rec_solv+best_energy_t->lig_solv, best_energy_t->hb_donor+best_energy_t->hb_acceptor, best_energy_t->total, best_conf, overlay_status[best_conf], energy_status, si);
    }
    else {
        sprintf(info, "%5d %-12.12s %-4.4s %-10.3g  %-8.3g %-8.3g %-8.3g %-8.2f %-8.3g %3d %2d %2d %.2f", counter, Lig->molname.c_str(), Lig->resnames[0].c_str(), overlays[best_conf],
                0.0, 0.0, 0.0, 0.0, 0.0, -1, 0, 0, 0.00);
    }
    this->print_info(info);

    if (Input->write_mol2){
        bool will_write = true;
        if (Input->use_writeMol2_score_cutoff){
            if (si < Input->writeMol2_score_cutoff){
                will_write = false;
            }
        }
        if (Input->use_writeMol2_energy_cutoff){
            if (best_energy_t->total > Input->writeMol2_energy_cutoff){
                will_write = false;
            }
        }
        if (will_write and match){
#pragma omp critical
            {
                this->write_mol2(Lig, new_xyz, best_energy_t, si);
            }
        }
    }

    if (Input->show_rmsd){
        double rmsd = Coord->compute_rmsd(RefLig->xyz, new_xyz, RefLig->N);
        sprintf(info, "RMSD: %7.3f", rmsd);
        this->print_info(info);
    }
    new_mcoords.clear();
    new_xyz.clear();
}

void Docker::Dock_conformers(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, Grid* Grids, unsigned counter){
    COORD_MC* Coord = new COORD_MC;
    double best_ene=0.00;
    int best_conf=0;
    double si;
    vector<vector<vector<double> > > new_mcoords;
    vector<double> energies;
    vector<double> overlays;
    vector<int> overlay_status;
    bool match=false;

    Gaussian* Gauss = new Gaussian;
    double t1, t2, t3;
    t2 = Gauss->compute_shape_and_charge_density(Input, RefLig, RefLig, RefLig->xyz, RefLig->xyz);

    for (unsigned i=0; i<Lig->mcoords.size(); i++){
        vector<double> com_lig = Coord->compute_com(Lig->mcoords[i], Lig);
        if (int(Lig->mcoords[i].size()) == Lig->N){

            Lig->xyz = Coord->translate(Lig->mcoords[i], Lig->N, com[0]-com_lig[0], com[1]-com_lig[1], com[2]-com_lig[2]);

/*
*
* Doing a pre-alignment of the longest axis of the ligands
*
*/


            Lig->find_longest_axis();
            RefLig->find_longest_axis();

            Optimizer* Opt0 = new Optimizer(Rec, RefLig, Input);
            Optimizer::opt_result_t* opt_result0 = new Optimizer::opt_result_t;
            opt_result0->energy_result = new energy_result_t;

            Opt0->pre_align(Lig, opt_result0);
            Lig->xyz = opt_result0->optimized_xyz;

            delete opt_result0;
            delete Opt0;


/*
* End of pre-alignment
*/

            Optimizer* Opt = new Optimizer(Rec, RefLig, Input, Grids);
            Optimizer::opt_result_t* opt_result = new Optimizer::opt_result_t;
            opt_result->energy_result = new energy_result_t;

            if (! this->minimize_overlay(Input, Opt, Lig, opt_result)){
                sprintf(info, "Overlay optimizer %s is not defined. Exiting...\n", Input->overlay_optimizer.c_str());
                this->print_info(info);
            }


            t1 = Gauss->compute_shape_and_charge_density(Input, RefLig, Lig, RefLig->xyz, opt_result->optimized_xyz);
            t3 = Gauss->compute_shape_and_charge_density(Input, Lig, Lig, opt_result->optimized_xyz,opt_result->optimized_xyz);
            si = (2*t1) / (t2+t3);

            overlays.push_back(si);
            energies.push_back(opt_result->energy_result->total);
            new_mcoords.push_back(opt_result->optimized_xyz);
            overlay_status.push_back(opt_result->optimization_status);

            delete Opt;
            delete opt_result;
        }
    }

    delete Gauss;
    delete Coord;

    vector<unsigned> index;

/*
* Here the user can choose if the best overlay is the one that maximizes the overlap of the docked ligand
* with the reference ligand (default) or the on which results in the best binding energy, evaluated after
* the overlap of the ligands. To choose the binding energy as the parameter, the input "sort_by_energy" must
* bet set to "yes" in the input file.
* ASN, Feb/2013.
*/

    if (Input->sort_by_energy){
        index = this->sort_vector(energies);
    }
    else {
        index = this->sort_vector_inv(overlays);
    }

    best_ene= 999999999.;
    vector<vector<double> >new_xyz;
    int energy_status=0;

    new_xyz = new_mcoords[index[0]];

    energy_result_t* best_energy_t = new energy_result_t;

    Gaussian* Gauss2 = new Gaussian;
    t2 = Gauss2->compute_shape_and_charge_density(Input, RefLig, RefLig, RefLig->xyz, RefLig->xyz);

    for (int i=0; i<Input->conformers_to_evaluate; i++){
        if (i < int(Lig->mcoords.size())){

            if (overlays[index[i]]>= Input->overlay_cutoff){
                match=true;

                Lig->xyz = new_mcoords[index[i]];

                Optimizer* Opt2 = new Optimizer(Rec, RefLig, Input, Grids);
                Optimizer::opt_result_t* opt_result2 = new Optimizer::opt_result_t;
                opt_result2->energy_result = new energy_result_t;

                if (! this->minimize_energy(Input, Opt2, Rec, Lig, RefLig, opt_result2)){
                    sprintf(info, "Energy optimizer %s is not defined. Exiting...\n", Input->energy_optimizer.c_str());
                    this->print_info(info);
                }

                if (opt_result2->energy_result->total < best_ene){
                    best_ene = opt_result2->energy_result->total;
                    best_energy_t = opt_result2->energy_result;
                    new_xyz = opt_result2->optimized_xyz;
                    best_conf=index[i];
                    energy_status = opt_result2->optimization_status;
                }
                delete opt_result2;
                delete Opt2;
            }
        }
    }
    t1 = Gauss2->compute_shape_and_charge_density(Input, RefLig, Lig, RefLig->xyz, new_xyz);
    t3 = Gauss2->compute_shape_and_charge_density(Input, Lig, Lig, new_xyz,new_xyz);

    si = (2*t1) / (t2+t3);
    delete Gauss2;

    Lig->xyz = new_xyz;

    if (match){
        sprintf(info, "%5d %-12.12s %-4.4s %-10.3g  %-8.3g %-8.3g %-8.3g %-8.2f %-8.3g %3d %2d %2d %.2f", counter, Lig->molname.c_str(), Lig->resnames[0].c_str(), overlays[best_conf],
                best_energy_t->elec, best_energy_t->vdw, (best_energy_t->rec_solv + best_energy_t->lig_solv), (best_energy_t->hb_donor + best_energy_t->hb_acceptor), best_energy_t->total, best_conf, overlay_status[best_conf], energy_status, si);
    }
    else {
        sprintf(info, "%5d %-12.12s %-4.4s %-10.3g  %-8.3g %-8.3g %-8.3g %-8.2f %-8.3g %3d %2d %2d %.2f", counter, Lig->molname.c_str(), Lig->resnames[0].c_str(), overlays[best_conf],
                0.00, 0.00, 0.00, 0.00, 0.00, -1, 0, 0, 0.00);
    }
    this->print_info(info);

    if (Input->write_mol2){
        bool will_write = true;
        if (Input->use_writeMol2_score_cutoff){
            if (si < Input->writeMol2_score_cutoff){
                will_write = false;
            }
        }
        if (Input->use_writeMol2_energy_cutoff){
            if (best_energy_t->total > Input->writeMol2_energy_cutoff){
                will_write = false;
            }
        }
        if (will_write and match){
#pragma omp critical
            {
                this->write_mol2(Lig, new_xyz, best_energy_t, si);
            }
        }
    }

    if (Input->show_rmsd){
        double rmsd = Coord->compute_rmsd(RefLig->xyz, new_xyz, RefLig->N);
        sprintf(info, "RMSD: %7.3f", rmsd);
        this->print_info(info);
    }
}

vector<unsigned> Docker::sort_vector(vector<double> vec){
    vector<unsigned> indexes;
    for (unsigned i=0; i<vec.size(); i++){

#ifdef DEBUG
        printf("vector[%d]: %7.3f\n", i, vec[i]);
#endif

        indexes.push_back(i);
    }

    double tmp;
    unsigned itmp;
    for (unsigned i=0; i< vec.size()-1; i++){
        for(unsigned j=i+1; j<vec.size(); j++){
            if (vec[j] < vec[i]){
                tmp = vec[i];
                itmp = indexes[i];

                vec[i]=vec[j];
                indexes[i] = indexes[j];

                vec[j] = tmp;
                indexes[j] = itmp;
            }
        }
    }
    if (indexes.size() != vec.size()){
        printf("Mismatch in vector sorting. Please check!\n");
        exit(1);
    }
    return indexes;
}

vector<unsigned> Docker::sort_vector_inv(vector<double> vec){
    vector<unsigned> indexes;
    for (unsigned i=0; i<vec.size(); i++){

#ifdef DEBUG
        printf("vector[%d]: %7.3f\n", i, vec[i]);
#endif

        indexes.push_back(i);
    }

    double tmp;
    unsigned itmp;
    for (unsigned i=0; i< vec.size()-1; i++){
        for(unsigned j=i+1; j<vec.size(); j++){
            if (vec[j] > vec[i]){
                tmp = vec[i];
                itmp = indexes[i];

                vec[i]=vec[j];
                indexes[i] = indexes[j];

                vec[j] = tmp;
                indexes[j] = itmp;
            }
        }
    }
    if (indexes.size() != vec.size()){
        printf("Mismatch in vector sorting. Please check!\n");
        exit(1);
    }
    return indexes;
}

bool Docker::minimize_overlay(PARSER* Input, Optimizer* Opt, Mol2* Lig, Optimizer::opt_result_t* opt_result){
    bool ret = false;
    if (Input->overlay_optimizer == "lbfgs" or Input->overlay_optimizer == "lbfgs2"){
        Opt->minimize_overlay_nlopt_lbfgs(Lig, opt_result);
        ret = true;
    }
    else if (Input->overlay_optimizer == "ln_auglag"){
        Opt->minimize_overlay_nlopt_ln_auglag(Lig, opt_result);
        ret = true;
    }
    else if (Input->overlay_optimizer == "ld_auglag"){
        Opt->minimize_overlay_nlopt_ld_auglag(Lig, opt_result);
        ret = true;
    }
    else if (Input->overlay_optimizer == "mma"){
        Opt->minimize_overlay_nlopt_mma(Lig, opt_result);
        ret = true;
    }
    else if (Input->overlay_optimizer == "subplex"){
        Opt->minimize_overlay_nlopt_subplex(Lig, opt_result);
        ret = true;
    }
    else if (Input->overlay_optimizer == "cobyla"){
        Opt->minimize_overlay_nlopt_cobyla(Lig, opt_result);
        ret = true;
    }
    else if (Input->overlay_optimizer == "crs"){
        Opt->minimize_overlay_nlopt_crs(Lig, opt_result);
        ret = true;
    }
    else if (Input->overlay_optimizer == "direct"){
        Opt->minimize_overlay_nlopt_direct(Lig, opt_result);
        ret = true;
    }
    else if (Input->overlay_optimizer == "none"){
        opt_result->optimized_xyz = Lig->xyz;
        opt_result->optimization_status = 0;
        opt_result->f_min = 0.0;
    }
    else {
        ret = false;
    }
    return(ret);
}

bool Docker::minimize_energy(PARSER* Input, Optimizer* Opt, Mol2* Rec, Mol2* Lig, Optimizer::opt_result_t* opt_result){
    bool ret = false;
    if (Input->energy_optimizer == "lbfgs" or Input->energy_optimizer == "lbfgs2"){
        Opt->minimize_energy_nlopt_lbfgs(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "adaptative"){
        Opt->minimize_energy_adaptative(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "cobyla"){
        Opt->minimize_energy_nlopt_cobyla(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "crs"){
        Opt->minimize_energy_nlopt_crs(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "direct"){
        Opt->minimize_energy_nlopt_direct(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "direct_only"){
        Opt->minimize_energy_nlopt_direct_only(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "ld_auglag"){
        Opt->minimize_energy_nlopt_ld_auglag(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "ln_auglag"){
        Opt->minimize_energy_nlopt_ln_auglag(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "mma"){
        Opt->minimize_energy_nlopt_mma(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "simplex"){
        Opt->minimize_energy_nlopt_simplex(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "stogo"){
        Opt->minimize_energy_nlopt_stogo(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "isres"){
        Opt->minimize_energy_nlopt_isres(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "esch"){
        Opt->minimize_energy_nlopt_esch(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "subplex"){
        Opt->minimize_energy_nlopt_subplex(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "SA"){
        Opt->Simulated_Annealing(Lig, opt_result);
        ret = true;
    }
    else if (Input->energy_optimizer == "none"){
        if (! Input->only_score){
            opt_result->energy_result->total = 0.00;
            opt_result->energy_result->vdw = 0.00;
            opt_result->energy_result->elec = 0.00;
            opt_result->energy_result->rec_solv = 0.00;
            opt_result->energy_result->lig_solv = 0.00;
            opt_result->energy_result->hb_donor = 0.00;
            opt_result->energy_result->hb_acceptor = 0.00;
            opt_result->optimization_status = 0;
            opt_result->optimized_xyz = Lig->xyz;
        }
        else {
            Energy2* Ene = new Energy2(Input);
            if (Input->use_grids){
                Ene->compute_ene(Opt->Grids, Lig, Lig->xyz, opt_result->energy_result);
            }
            else {
                Ene->compute_ene(Rec, Lig, Lig->xyz, opt_result->energy_result);
            }
            opt_result->optimization_status = 0;
            opt_result->optimized_xyz = Lig->xyz;
            delete Ene;
        }
        ret = true;
    }
    else{
        ret = false;
    }
    return (ret);
}

bool Docker::minimize_energy(PARSER* Input, Optimizer* Opt, Mol2* Rec, Mol2* Lig, Mol2* Reflig, Optimizer::opt_result_t* opt_result){
    bool ret = false;

    // Using the similarity index to define whether energy optimization will or will no be
    // done. In this case, si has to be greater than or equal to overlay_cutoff
    //
    if (Input->use_overlay_cutoff){
        Gaussian* Gauss = new Gaussian;
        double t1, t2, t3, si;
        t1 = Gauss->compute_shape_and_charge_density(Input, Reflig, Lig, Reflig->xyz, Lig->xyz);
        t2 = Gauss->compute_shape_and_charge_density(Input, Reflig, Reflig, Reflig->xyz, Reflig->xyz);
        t3 = Gauss->compute_shape_and_charge_density(Input, Lig, Lig, Lig->xyz,Lig->xyz);
        si = (2*t1) / (t2+t3);
        delete Gauss;

        if (si < Input->overlay_cutoff){
            opt_result->energy_result->total = 0.00;
            opt_result->energy_result->vdw = 0.00;
            opt_result->energy_result->elec = 0.00;
            opt_result->energy_result->rec_solv = 0.00;
            opt_result->energy_result->lig_solv = 0.00;
            opt_result->energy_result->hb_donor = 0.00;
            opt_result->energy_result->hb_acceptor = 0.00;
            opt_result->optimization_status = 0;
            opt_result->optimized_xyz = Lig->xyz;
            ret = true;
        }
        else {
            if (Input->energy_optimizer == "lbfgs" or Input->energy_optimizer == "lbfgs2"){
                Opt->minimize_energy_nlopt_lbfgs(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "adaptative"){
                Opt->minimize_energy_adaptative(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "cobyla"){
                Opt->minimize_energy_nlopt_cobyla(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "crs"){
                Opt->minimize_energy_nlopt_crs(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "direct"){
                Opt->minimize_energy_nlopt_direct(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "direct_only"){
                Opt->minimize_energy_nlopt_direct_only(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "ld_auglag"){
                Opt->minimize_energy_nlopt_ld_auglag(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "ln_auglag"){
                Opt->minimize_energy_nlopt_ln_auglag(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "mma"){
                Opt->minimize_energy_nlopt_mma(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "simplex"){
                Opt->minimize_energy_nlopt_simplex(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "stogo"){
                Opt->minimize_energy_nlopt_stogo(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "isres"){
                Opt->minimize_energy_nlopt_isres(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "esch"){
                Opt->minimize_energy_nlopt_esch(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "subplex"){
                Opt->minimize_energy_nlopt_subplex(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "SA"){
                Opt->Simulated_Annealing(Lig, opt_result);
                ret = true;
            }
            else if (Input->energy_optimizer == "none"){
                if (! Input->only_score){
                    opt_result->energy_result->total = 0.00;
                    opt_result->energy_result->vdw = 0.00;
                    opt_result->energy_result->elec = 0.00;
                    opt_result->energy_result->rec_solv = 0.00;
                    opt_result->energy_result->lig_solv = 0.00;
                    opt_result->energy_result->hb_donor = 0.00;
                    opt_result->energy_result->hb_acceptor = 0.00;
                    opt_result->optimization_status = 0;
                    opt_result->optimized_xyz = Lig->xyz;
                }
                else {
                    Energy2* Ene = new Energy2(Input);
                    if (Input->use_grids){
                        Ene->compute_ene(Opt->Grids, Lig, Lig->xyz, opt_result->energy_result);
                    }
                    else {
                        Ene->compute_ene(Rec, Lig, Lig->xyz, opt_result->energy_result);
                    }
                    opt_result->optimization_status = 0;
                    opt_result->optimized_xyz = Lig->xyz;
                    delete Ene;
                }
                ret = true;
            }
            else{
                ret = false;
            }
        }
    }
    else {
        //
        // Not using the overlay cutoff. Works similarly to the next function.
        //
        if (Input->energy_optimizer == "lbfgs" or Input->energy_optimizer == "lbfgs2"){
            Opt->minimize_energy_nlopt_lbfgs(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "adaptative"){
            Opt->minimize_energy_adaptative(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "cobyla"){
            Opt->minimize_energy_nlopt_cobyla(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "crs"){
            Opt->minimize_energy_nlopt_crs(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "direct"){
            Opt->minimize_energy_nlopt_direct(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "direct_only"){
            Opt->minimize_energy_nlopt_direct_only(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "ld_auglag"){
            Opt->minimize_energy_nlopt_ld_auglag(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "ln_auglag"){
            Opt->minimize_energy_nlopt_ln_auglag(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "mma"){
            Opt->minimize_energy_nlopt_mma(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "simplex"){
            Opt->minimize_energy_nlopt_simplex(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "stogo"){
            Opt->minimize_energy_nlopt_stogo(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "isres"){
            Opt->minimize_energy_nlopt_isres(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "esch"){
            Opt->minimize_energy_nlopt_esch(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "subplex"){
            Opt->minimize_energy_nlopt_subplex(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "SA"){
            Opt->Simulated_Annealing(Lig, opt_result);
            ret = true;
        }
        else if (Input->energy_optimizer == "none"){
            if (! Input->only_score){
                opt_result->energy_result->total = 0.00;
                opt_result->energy_result->vdw = 0.00;
                opt_result->energy_result->elec = 0.00;
                opt_result->energy_result->rec_solv = 0.00;
                opt_result->energy_result->lig_solv = 0.00;
                opt_result->energy_result->hb_donor = 0.00;
                opt_result->energy_result->hb_acceptor = 0.00;
                opt_result->optimization_status = 0;
                opt_result->optimized_xyz = Lig->xyz;
            }
            else {
                Energy2* Ene = new Energy2(Input);
                if (Input->use_grids){
                    Ene->compute_ene(Opt->Grids, Lig, Lig->xyz, opt_result->energy_result);
                }
                else {
                    Ene->compute_ene(Rec, Lig, Lig->xyz, opt_result->energy_result);
                }
                opt_result->optimization_status = 0;
                opt_result->optimized_xyz = Lig->xyz;
                delete Ene;
            }
            ret = true;
        }
        else{
            ret = false;
        }
    }
    return (ret);
}


#ifdef HAS_GUI

Docker::Docker(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, QtWriter* Writer, unsigned counter) {

    if (!Input->generate_conformers){

        COORD_MC* Coord = new COORD_MC;
        int overlay_status;
        double overlay_fmax;
        vector<double> com_lig = Coord->compute_com(Lig);

        /*
         * Shifting the ligand to match its center of mass with the reference
         * ligand center of mass.
         */

        Lig->xyz = Coord->translate(Lig->xyz, Lig->N, com[0]-com_lig[0], com[1]-com_lig[1], com[2]-com_lig[2]);

        /*
*
* Doing a pre-alignment of the longest axis of the ligands
*
*/

        Lig->find_longest_axis();
        RefLig->find_longest_axis();

        Optimizer* Opt0 = new Optimizer(Rec, RefLig, Input);
        Optimizer::opt_result_t* opt_result0 = new Optimizer::opt_result_t;
        opt_result0->energy_result = new energy_result_t;

        Opt0->pre_align(Lig, opt_result0);
        Lig->xyz = opt_result0->optimized_xyz;

        delete opt_result0;
        delete Opt0;

        /*
* End of pre-alignment
*/

        Optimizer* Opt = new Optimizer(Rec, RefLig, Input);
        Optimizer::opt_result_t* opt_result = new Optimizer::opt_result_t;
        opt_result->energy_result = new energy_result_t;

        // Optimizing overlay....

        if (! this->minimize_overlay(Input, Opt, Lig, opt_result)){
            sprintf(info, "Overlay optimizer %s is not defined. Exiting...\n", Input->overlay_optimizer.c_str());
        }

        //Copying new coordinates.
        Lig->xyz = opt_result->optimized_xyz;
        overlay_status = opt_result->optimization_status;
        overlay_fmax = opt_result->f_min;

        //        delete opt_result->energy_result;
        delete opt_result;

        //Optimizing Energy...
        Optimizer::opt_result_t* opt_result2 = new Optimizer::opt_result_t;
        opt_result2->energy_result = new energy_result_t;

        if (Input->deal){
            sprintf(info, "Entering Deal....");
            this->print_info(info);
#ifdef HAS_GUI
            Deal* DEAl = new Deal(Rec, Lig, Input);
            sprintf(info, "%-20.20s %-20.20s %-7.3e % -7.3f kcal/mol", Lig->molname.c_str(), Lig->resnames[0].c_str(), fo, DEAl->energies[0]);
            this->print_info(info);
            Writer->writeMol2(Lig, DEAl->final_coords, DEAl->energies[0], fo);
            delete DEAl;
#endif
        }

        else {
            if (! this->minimize_energy(Input, Opt, Rec, Lig, opt_result2)){
                sprintf(info, "Energy optimizer %s is not defined. Exiting...\n", Input->energy_optimizer.c_str());
            }

            Lig->xyz = opt_result2->optimized_xyz;

            sprintf(info, "%5d %-10.10s %-10.10s %-11.3e %-8.3g %-8.3g %-8.3g %8.3g %8.3g %-3d %-2d %-2d", counter, Lig->molname.c_str(), Lig->resnames[0].c_str(), overlay_fmax,
                    opt_result2->energy_result->elec, opt_result2->energy_result->vdw, opt_result2->energy_result->rec_solv+opt_result2->energy_result->lig_solv, opt_result2->energy_result->hb_donor+opt_result2->energy_result->hb_acceptor,
                    opt_result2->energy_result->total, 0, overlay_status, opt_result2->optimization_status);
            this->print_info(info);

            if (Input->write_mol2){
                this->write_mol2(Lig, opt_result2->optimized_xyz, opt_result2->energy_result, overlay_fmax);
            }

        }
        delete Coord;
        delete opt_result2;
        delete Opt;
    }
    else {
        this->Dock_conformers(Rec, Lig, RefLig, com, Input, Writer, counter);
    }
}

Docker::Docker(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, QtWriter* Writer, Grid* Grids, unsigned counter) {

    if (!Input->generate_conformers){

        COORD_MC* Coord = new COORD_MC;
        int overlay_status;
        double overlay_fmax;
        vector<double> com_lig = Coord->compute_com(Lig);

        /*
         * Shifting the ligand to match its center of mass with the reference
         * ligand center of mass.
         */

        Lig->xyz = Coord->translate(Lig->xyz, Lig->N, com[0]-com_lig[0], com[1]-com_lig[1], com[2]-com_lig[2]);

        /*
*
* Doing a pre-alignment of the longest axis of the ligands
*
*/

        Lig->find_longest_axis();
        RefLig->find_longest_axis();

        Optimizer* Opt0 = new Optimizer(Rec, RefLig, Input);
        Optimizer::opt_result_t* opt_result0 = new Optimizer::opt_result_t;
        opt_result0->energy_result = new energy_result_t;

        Opt0->pre_align(Lig, opt_result0);
        Lig->xyz = opt_result0->optimized_xyz;

        delete opt_result0;
        delete Opt0;

        /*
* End of pre-alignment
*/

        Optimizer* Opt = new Optimizer(Rec, RefLig, Input, Grids);
        Optimizer::opt_result_t* opt_result = new Optimizer::opt_result_t;
        opt_result->energy_result = new energy_result_t;

        // Optimizing overlay....

        if (! this->minimize_overlay(Input, Opt, Lig, opt_result)){
            sprintf(info, "Overlay optimizer %s is not defined. Exiting...\n", Input->overlay_optimizer.c_str());
        }

        //Copying new coordinates.
        Lig->xyz = opt_result->optimized_xyz;
        overlay_status = opt_result->optimization_status;
        overlay_fmax = opt_result->f_min;

        delete opt_result;

        //Optimizing Energy...

        Optimizer::opt_result_t* opt_result2 = new Optimizer::opt_result_t;
        opt_result2->energy_result = new energy_result_t;

#ifdef DEBUG
        sprintf(info, "Starting Binding Energy Minimization...");
        this->print_info(info);
#endif

        if (Input->deal){
            sprintf(info, "Entering Deal....");
            this->print_info(info);
#ifdef HAS_GUI
            Deal* DEAl = new Deal(Rec, Lig, Input);
            sprintf(info, "%-20.20s %-20.20s %-7.3e % -7.3f kcal/mol", Lig->molname.c_str(), Lig->resnames[0].c_str(), fo, DEAl->energies[0]);
            this->print_info(info);
            Writer->writeMol2(Lig, DEAl->final_coords, DEAl->energies[0], fo);
            delete DEAl;
#endif
        }

        else {
            if (! this->minimize_energy(Input, Opt, Rec, Lig, opt_result2)){
                sprintf(info, "Energy optimizer %s is not defined. Exiting...\n", Input->energy_optimizer.c_str());
            }

            Gaussian* Gauss = new Gaussian;
            double t1, t2, t3, si;
            t1 = Gauss->compute_shape_and_charge_density(Input, RefLig, Lig, RefLig->xyz, opt_result2->optimized_xyz);
            t2 = Gauss->compute_shape_and_charge_density(Input, RefLig, RefLig, RefLig->xyz, RefLig->xyz);
            t3 = Gauss->compute_shape_and_charge_density(Input, Lig, Lig, opt_result2->optimized_xyz,opt_result2->optimized_xyz);

            si = (2*t1) / (t2+t3);
            delete Gauss;

            Lig->xyz = opt_result2->optimized_xyz;

            sprintf(info, "%5d %-12.12s %-4.4s %-10.3e  %-8.3g %-8.3g %-8.3g %-8.3g %-8.3g %3d %2d %2d %.2f", counter, Lig->molname.c_str(), Lig->resnames[0].c_str(), overlay_fmax,
                    opt_result2->energy_result->elec, opt_result2->energy_result->vdw, opt_result2->energy_result->rec_solv+opt_result2->energy_result->lig_solv,
                    opt_result2->energy_result->hb_donor+opt_result2->energy_result->hb_acceptor, opt_result2->energy_result->total, 0, overlay_status, opt_result2->optimization_status, si);
            this->print_info(info);

            if (Input->write_mol2){
                this->write_mol2(Lig, opt_result2->optimized_xyz, opt_result2->energy_result, si);
            }

        }
        delete Coord;
        delete opt_result2;
        delete Opt;
    }
    else {
        this->Dock_conformers(Rec, Lig, RefLig, com, Input, Writer, Grids, counter);
    }
}


void  Docker::Dock_conformers(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, QtWriter* Writer, unsigned counter){
    COORD_MC* Coord = new COORD_MC;
    double best_ene=0.00;
    int best_conf=0, energy_status=0;
    vector<int> overlay_status;
    vector<vector<vector<double> > > new_mcoords;
    vector<double> energies;
    vector<double> overlays;
    double si;

    for (unsigned i=0; i<Lig->mcoords.size(); i++){
        vector<double> com_lig = Coord->compute_com(Lig->mcoords[i], Lig);
        if (int(Lig->mcoords[i].size()) == Lig->N){

            Lig->xyz = Coord->translate(Lig->mcoords[i], Lig->N, com[0]-com_lig[0], com[1]-com_lig[1], com[2]-com_lig[2]);

            /*
*
* Doing a pre-alignment of the longest axis of the ligands
*
*/

            Lig->find_longest_axis();
            RefLig->find_longest_axis();

            Optimizer* Opt0 = new Optimizer(Rec, RefLig, Input);
            Optimizer::opt_result_t* opt_result0 = new Optimizer::opt_result_t;
            opt_result0->energy_result = new energy_result_t;

            Opt0->pre_align(Lig, opt_result0);
            Lig->xyz = opt_result0->optimized_xyz;

            delete opt_result0;
            delete Opt0;

            /*
* End of pre-alignment
*/

            Optimizer* Opt = new Optimizer(Rec, RefLig, Input);
            Optimizer::opt_result_t* opt_result = new Optimizer::opt_result_t;
            opt_result->energy_result = new energy_result_t;

            if (! this->minimize_overlay(Input, Opt, Lig, opt_result)){
                sprintf(info, "Overlay optimizer %s is not defined. Exiting...\n", Input->overlay_optimizer.c_str());
            }

            overlays.push_back(opt_result->f_min);
            energies.push_back(opt_result->energy_result->total);
            new_mcoords.push_back(opt_result->optimized_xyz);
            overlay_status.push_back(opt_result->optimization_status);
            delete Opt;
            delete opt_result;
        }
    }

    vector<unsigned> index;


    if (Input->sort_by_energy){
        index = this->sort_vector(energies);
    }
    else {
        index = this->sort_vector_inv(overlays);
    }

    best_ene= 999999999.;
    vector<vector<double> > new_xyz;
    energy_status=0;

    new_xyz = new_mcoords[index[0]];
    energy_result_t* best_energy_t = new energy_result_t;

    for (int i=0; i<Input->conformers_to_evaluate; i++){
        if (int(new_mcoords[index[i]].size()) == Lig->N){

            Lig->xyz = new_mcoords[index[i]];

            Optimizer* Opt = new Optimizer(Rec, RefLig, Input);

            Optimizer::opt_result_t* opt_result = new Optimizer::opt_result_t;
            opt_result->energy_result = new energy_result_t;

            if (! this->minimize_energy(Input, Opt, Rec, Lig, opt_result)){
                sprintf(info, "Energy optimizer %s is not defined. Exiting...\n", Input->energy_optimizer.c_str());
            }

            if (opt_result->energy_result->total < best_ene){
                best_ene = opt_result->energy_result->total;
                best_energy_t = opt_result->energy_result;
                new_xyz = opt_result->optimized_xyz;
                best_conf=index[i];
                energy_status = opt_result->optimization_status;
            }
            delete Opt;
            delete opt_result;
        }
    }
    Gaussian* Gauss = new Gaussian;
    double t1, t2, t3;
    t1 = Gauss->compute_shape_and_charge_density(Input, RefLig, Lig, RefLig->xyz, new_xyz);
    t2 = Gauss->compute_shape_and_charge_density(Input, RefLig, RefLig, RefLig->xyz, RefLig->xyz);
    t3 = Gauss->compute_shape_and_charge_density(Input, Lig, Lig, new_xyz, new_xyz);

    si = (2*t1) / (t2+t3);
    delete Gauss;

    Lig->xyz = new_xyz;

    sprintf(info, "%5d %-12.12s %-4.4s %-10.3e  %-8.3g %-8.3g %-8.3g %-8.3g %-8.3g %3d %2d %2d %.2f", counter, Lig->molname.c_str(), Lig->resnames[0].c_str(), overlays[best_conf],
            best_energy_t->elec, best_energy_t->vdw, best_energy_t->rec_solv+best_energy_t->lig_solv, best_energy_t->hb_donor+best_energy_t->hb_acceptor, best_energy_t->total, best_conf, overlay_status[best_conf], energy_status, si);
    this->print_info(info);

    if (Input->write_mol2){
        bool will_write = true;
        if (Input->use_writeMol2_score_cutoff){
            if (si < Input->writeMol2_score_cutoff){
                will_write = false;
            }
        }
        if (Input->use_writeMol2_energy_cutoff){
            if (best_energy_t->total > Input->writeMol2_energy_cutoff){
                will_write = false;
            }
        }
        if (will_write){
#pragma omp critical
            {
                this->write_mol2(Lig, new_xyz, best_energy_t, si);
            }
        }
    }

    if (Input->show_rmsd){
        double rmsd = Coord->compute_rmsd(RefLig->xyz, new_xyz, RefLig->N);
        sprintf(info, "RMSD: %7.3f", rmsd);
        this->print_info(info);
    }
}

void  Docker::Dock_conformers(Mol2* Rec, Mol2* Lig, Mol2* RefLig, vector<double> com, PARSER* Input, QtWriter* Writer, Grid* Grids, unsigned counter){
    COORD_MC* Coord = new COORD_MC;
    double best_ene=0.00;
    int best_conf=0;
    double si;
    vector<vector<vector<double> > > new_mcoords;
    vector<double> energies;
    vector<double> overlays;
    vector<int> overlay_status;

    for (unsigned i=0; i<Lig->mcoords.size(); i++){
        vector<double> com_lig = Coord->compute_com(Lig->mcoords[i], Lig);
        if (int(Lig->mcoords[i].size()) == Lig->N){

            Lig->xyz = Coord->translate(Lig->mcoords[i], Lig->N, com[0]-com_lig[0], com[1]-com_lig[1], com[2]-com_lig[2]);

            /*
*
* Doing a pre-alignment of the longest axis of the ligands
*
*/

            Lig->find_longest_axis();
            RefLig->find_longest_axis();

            Optimizer* Opt0 = new Optimizer(Rec, RefLig, Input);
            Optimizer::opt_result_t* opt_result0 = new Optimizer::opt_result_t;
            opt_result0->energy_result = new energy_result_t;

            Opt0->pre_align(Lig, opt_result0);
            Lig->xyz = opt_result0->optimized_xyz;

            delete opt_result0;
            delete Opt0;

            /*
* End of pre-alignment
*/

            Optimizer* Opt = new Optimizer(Rec, RefLig, Input, Grids);
            Optimizer::opt_result_t* opt_result = new Optimizer::opt_result_t;
            opt_result->energy_result = new energy_result_t;

            if (! this->minimize_overlay(Input, Opt, Lig, opt_result)){
                sprintf(info, "Overlay optimizer %s is not defined. Exiting...\n", Input->overlay_optimizer.c_str());
            }

            overlays.push_back(opt_result->f_min);
            energies.push_back(opt_result->energy_result->total);
            new_mcoords.push_back(opt_result->optimized_xyz);
            overlay_status.push_back(opt_result->optimization_status);

            delete Opt;
            delete opt_result;
        }
    }

    delete Coord;

    vector<unsigned> index;

    /*
     * Here the user can choose if the best overlay is the one that maximizes the overlap of the docked ligand
     * with the reference ligand (default) or the on which results in the best binding energy, evaluated after
     * the overlap of the ligands. To choose the binding energy as the parameter, the input "sort_by_energy" must
     * bet set to "yes" in the input file.
     * ASN, Feb/2013.
     */

    if (Input->sort_by_energy){
        index = this->sort_vector(energies);
    }
    else {
        index = this->sort_vector_inv(overlays);
    }

    best_ene= 999999999.;
    vector<vector<double> >new_xyz;
    int energy_status=0;

    new_xyz = new_mcoords[index[0]];

    energy_result_t* best_energy_t = new energy_result_t;

    for (int i=0; i<Input->conformers_to_evaluate; i++){
        if (i < int(Lig->mcoords.size())){

            Lig->xyz = new_mcoords[index[i]];

            Optimizer* Opt = new Optimizer(Rec, RefLig, Input, Grids);
            Optimizer::opt_result_t* opt_result = new Optimizer::opt_result_t;
            opt_result->energy_result = new energy_result_t;

            if (! this->minimize_energy(Input, Opt, Rec, Lig, opt_result)){
                sprintf(info, "Energy optimizer %s is not defined. Exiting...\n", Input->energy_optimizer.c_str());
            }

            if (opt_result->energy_result->total < best_ene){
                best_ene = opt_result->energy_result->total;
                best_energy_t = opt_result->energy_result;
                new_xyz = opt_result->optimized_xyz;
                best_conf=index[i];
                energy_status = opt_result->optimization_status;
            }
            delete Opt;
            delete opt_result;
        }
    }
    Gaussian* Gauss = new Gaussian;
    double t1, t2, t3;
    t1 = Gauss->compute_shape_and_charge_density(Input, RefLig, Lig, RefLig->xyz, new_xyz);
    t2 = Gauss->compute_shape_and_charge_density(Input, RefLig, RefLig, RefLig->xyz, RefLig->xyz);
    t3 = Gauss->compute_shape_and_charge_density(Input, Lig, Lig, new_xyz, new_xyz);

    si = (2*t1) / (t2+t3);
    delete Gauss;

    Lig->xyz = new_xyz;

    sprintf(info, "%5d %-12.12s %-4.4s %-10.3g  %-8.3g %-8.3g %-8.3g %-8.3g %-8.3g %3d %2d %2d %.2f", counter, Lig->molname.c_str(), Lig->resnames[0].c_str(), overlays[best_conf],
            best_energy_t->elec, best_energy_t->vdw, best_energy_t->rec_solv+best_energy_t->lig_solv, best_energy_t->hb_donor+best_energy_t->hb_acceptor, best_energy_t->total, best_conf, overlay_status[best_conf], energy_status, si);
    this->print_info(info);

    if (Input->write_mol2){
        bool will_write = true;
        if (Input->use_writeMol2_score_cutoff){
            if (si < Input->writeMol2_score_cutoff){
                will_write = false;
            }
        }
        if (Input->use_writeMol2_energy_cutoff){
            if (best_energy_t->total > Input->writeMol2_energy_cutoff){
                will_write = false;
            }
        }
        if (will_write){
#pragma omp critical
            {
                this->write_mol2(Lig, new_xyz, best_energy_t, si);
            }
        }
    }
}

#endif

void Docker::print_info(char info[98]){
#ifdef HAS_GUI
    QWriter->print_info(info);
#else
    Writer->print_info(info);
#endif
}

void Docker::print_line(){
#ifdef HAS_GUI
    QWriter->print_line();
#else
    Writer->print_line();
#endif
}


void Docker::write_mol2(Mol2* Lig, vector<vector<double> > new_xyz, double ene, double si){
#ifdef HAS_GUI
    QWriter->writeMol2(Lig, new_xyz, ene, si);
#else
    Writer->writeMol2(Lig, new_xyz, ene, si);
#endif
}

void Docker::write_mol2(Mol2* Lig, vector<vector<double> > new_xyz, energy_result_t* result, double si){
#ifdef HAS_GUI
    QWriter->writeMol2(Lig, new_xyz, result, si);
#else
    Writer->writeMol2(Lig, new_xyz, result, si);
#endif
}
