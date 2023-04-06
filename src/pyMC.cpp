#include "pyMC.h"
#include "pyEnergy2.cpp"
//#include "pyCOORD_MC.cpp"

using namespace OpenBabel;

MC::MC(WRITER* _Writer)
{
    info= new char[98];

    srand(rand());
    r = gsl_rng_alloc (gsl_rng_ranlxs2);
    Writer = _Writer;
}

MC::MC(Mol2* Lig, PARSER* Input, WRITER* _Writer){

    info= new char[98];
    srand(rand());
    r = gsl_rng_alloc (gsl_rng_ranlxs2);
    Writer = _Writer;

    if (Input->write_mol2){
        gzFile outmol2;
        string outname = Input->output + "_" + Lig->molname + "_MC";
        outmol2 = gzopen((outname+".mol2.gz").c_str(), "w");
        gzclose(outmol2);
        if (Input->ligsim){
            gzFile outmol2_lig;
            outname = Input->output + "_" + Lig->molname + "_MC.ligsim";
            outmol2_lig = gzopen((outname+".mol2.gz").c_str(), "w");
            gzclose(outmol2_lig);
        }
    }


/*! Here we will use OpenBabel API to generate a OBMol and use it to get and set torsion angles for all
 * rotatable bonds.
*/


    mol = this->GetMol(Input->lig_mol2);


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

    this->copy_to_obmol(Lig->xyz);

    OBff->Setup(mol);
    OBff->GetCoordinates(mol);
    OBff->SteepestDescent(Input->conformer_min_steps, 1.0e-10);
    OBff->GetCoordinates(mol);
    RotorList.Setup(mol);
    Rotor = RotorList.BeginRotor(RotorIterator);
    mol.ToInertialFrame();


    vector<int> tmp(4);
    Writer->print_line();
    sprintf(info, "Found %lu rotatable bonds in ligand %s.", RotorList.Size(), Lig->molname.c_str());
    Writer->print_info(info);
    Writer->print_line();
    for (unsigned i = 0; i < RotorList.Size(); ++i, Rotor = RotorList.NextRotor(RotorIterator)) {
        tmp = Rotor->GetDihedralAtoms();
        atoms_in_dihedrals.push_back(tmp);
        tmp.clear();
    }

// Preparing Writer pointer and GLS for random number sorting...

    myxyz = new double[mol.NumAtoms()*3];
}


MC::~MC(){
    gsl_rng_free (r);
    delete info;
}

void MC::write_conformers(Mol2* Lig){
    for (unsigned i=0; i< Lig->mcoords.size(); i++){
        Writer->writeMol2(Lig, Lig->mcoords[i], 0.0, 0.0, "teste");
    }
}

void MC::run(Grid* Grids, Mol2* RefLig, Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T){
    this->xyz = xyz;
    int nReject = 0;

    double sum_x = 0.0;
    double sum_xsquared = 0.0;
    long double sum_Boltzmann_ene = 0.0L;
    long double sum_Boltzmann2_ene = 0.0L;
    long double sum_Boltzmann2_ene_squared = 0.0L;
    bool ligand_is_in = false;
    long double independent_average = 0.0L;

    this->average_energy = 0.0;
    this->energy_standard_deviation = 0.0;
    this->Boltzmann_weighted_average_energy = 0.0L;
    this->MCR_Boltzmann_weighted_average = 0.0L;
    this->MCR_Boltzmann_weighted_stdev = 0.0L;

    double k=0.0019858775203792202;

    this->MaxMin.push_back(99999.0);         // min x
    this->MaxMin.push_back(-99999.0);        // max x
    this->MaxMin.push_back(99999.0);         // min y
    this->MaxMin.push_back(-99999.0);        // max y
    this->MaxMin.push_back(99999.0);         // min z
    this->MaxMin.push_back(-99999.0);        // max z

    gsl_rng_set(r, Input->seed);

    sprintf(info, "Starting Monte Carlo equilibration simulation with %5d steps...", (Input->eq_steps));
    Writer->print_info(info);
    gzFile mc_output;
    mc_output = gzopen((Input->output + "_" + Lig->molname + "_mc.dat.gz").c_str(), "w");

    Energy2* Energy = new Energy2(Input);
    COORD_MC* Coord = new COORD_MC;
    vector<double> original_com = Coord->compute_com(Lig);

    vector<double> com = Coord->compute_com(Lig);

    vector<double> rot_angles(3);
    for (unsigned i=0; i<3; i++){
        rot_angles[i]= 0.0;
    }

    int n_rot = int(RotorList.Size());

    McEntropy* Entropy = new McEntropy(Input, Coord, original_com, n_rot);

    energy_result_t* energy_t = new energy_result_t;

    double energy=0.0, new_energy=0.0, p=0.0, rnumber=0.0, rmsd=0.0;
    step_t* step = new step_t;

    while (!ligand_is_in){
        if (Input->sample_torsions){
            this->take_step_torsion(Input, Lig, step);
        }
        else if (Input->mc_full_flex){
            this->take_step_full_flex(Input, Lig, step);
        }
        else {
            this->take_step(Input, Lig, step);
        }
        com = Coord->compute_com(step->xyz, Lig);
        ligand_is_in = this->ligand_is_inside_box(Input, step, original_com, com);
    }
    ligand_is_in = false;

    Energy->compute_ene(Grids, Lig, step->xyz, energy_t);
    if (! Input->use_only_binding_energy){
        energy = energy_t->total+step->internal_energy;
    }

    sprintf(info, "#%10s %10s %10s", "Step", "Energy", "RMSD");
    gzprintf(mc_output,"###################################################################################################################################################################\n");
    gzprintf(mc_output,"#BoxsideX= %10.4f BoxsideY= %10.4f BoxsideY= %10.4f Temperature= %10.4f \n", Input->x_dim, Input->y_dim, Input->z_dim, Input->temp);
    gzprintf(mc_output,"#NROT = %10.10d\n", int(RotorList.Size()));
    gzprintf(mc_output,"###################################################################################################################################################################\n");
    if (n_rot == 0){
        gzprintf(mc_output, "#%10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s\n", "Step", "Energy", "RMSD", "DX", "DY", "DZ", "ALPHA", "BETA", "GAMMA", "NCONF", "ConfEnergy", "Elec", "VDW", "SOLV");
    }
    else {
        gzprintf(mc_output, "#%10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s ", "Step", "Energy", "RMSD", "DX", "DY", "DZ", "ALPHA", "BETA", "GAMMA", "NCONF", "ConfEnergy", "Elec", "VDW", "SOLV");
        for (unsigned i=1; i <= RotorList.Size(); i++){
            sprintf(info, "ROT[%3d]", i);
            gzprintf(mc_output, "%10.10s", info);
        }
        gzprintf(mc_output, "\n");
    }

    int count=0;
    int eqcount = 0;


    //Equilibration implementation

    while (eqcount <= Input->eq_steps){

        while (!ligand_is_in){
            if (Input->sample_torsions){
                this->take_step_torsion(Input, Lig, step);
            }
            else if (Input->mc_full_flex){
                this->take_step_full_flex(Input, Lig, step);
            }
            else {
                this->take_step(Input, Lig, step);
            }
            com = Coord->compute_com(step->xyz, Lig);
            ligand_is_in = this->ligand_is_inside_box(Input, step, original_com, com);
        }
        ligand_is_in = false;

        Energy->compute_ene(Grids, Lig, step->xyz, energy_t);
        if (! Input->use_only_binding_energy){
            new_energy = energy_t->total+step->internal_energy;
        }
        else {
            new_energy = energy_t->total;
        }

        if (new_energy <= energy){
            Lig->mcoords = Lig->new_mcoords;
            this->xyz = step->xyz;
            energy = new_energy;
            rmsd = Coord->compute_rmsd(Lig->xyz, step->xyz, Lig->N);

            this->increment_angles(&rot_angles, step);
            this->MaxMinCM(com[0],com[1],com[2],this->MaxMin);
            eqcount++;
        }
        else{
            p = this->Boltzmman(energy, new_energy, T, Input->bi);
            rnumber = gsl_rng_uniform(r);
            if (p > rnumber){
                Lig->mcoords = Lig->new_mcoords;
                this->xyz = step->xyz;
                energy = new_energy;
                rmsd = Coord->compute_rmsd(Lig->xyz, step->xyz, Lig->N);
                this->increment_angles(&rot_angles, step);
                this->MaxMinCM(com[0],com[1],com[2],this->MaxMin);
                eqcount++;
            }
        }
    }

    Writer->print_line();
    sprintf(info, "Equilibration done with %5d steps. Current system energy: %9.3f kcal/mol.", Input->eq_steps, energy);
    Writer->print_info(info);
    Writer->print_line();
    sprintf(info, "%10s %10s %10s", "Step", "Energy", "RMSD");
    Writer->print_info(info);

    while (count <= Input->number_steps){

        while (!ligand_is_in){

            if (Input->sample_torsions){
                this->take_step_torsion(Input, Lig, step);
            }
            else if (Input->mc_full_flex){
                this->take_step_full_flex(Input, Lig, step);
            }
            else {
                this->take_step(Input, Lig, step);
            }
            com = Coord->compute_com(step->xyz, Lig);
            ligand_is_in = this->ligand_is_inside_box(Input, step, original_com, com);
        }
        ligand_is_in = false;

        Energy->compute_ene(Grids, Lig, step->xyz, energy_t);
        if (Input->use_only_binding_energy){
            new_energy = energy_t->total;
        }
        else {
            new_energy = energy_t->total+step->internal_energy;
        }

        if (new_energy <= energy){
            Lig->mcoords = Lig->new_mcoords;
            this->xyz = step->xyz;
            energy = new_energy;
            rmsd = Coord->compute_rmsd(Lig->xyz, step->xyz, Lig->N);
            this->increment_angles(&rot_angles, step);
            this->MaxMinCM(com[0],com[1],com[2],this->MaxMin);
            count++;
            sum_x += energy;
            sum_xsquared += (energy*energy);
            sum_Boltzmann_ene += energy*exp(((-(Input->bi-1.0))*energy)/(k*T));
            sum_Boltzmann2_ene += exp(((-(Input->bi-1.0))*energy)/(k*T));
            sum_Boltzmann2_ene_squared += exp(((-(Input->bi-1.0))*energy)/(k*T)) * exp(((-(Input->bi-1.0))*energy)/(k*T));

            Entropy->update(com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->torsion_angles);

            if (count % Input->mc_stride == 0){

                sprintf(info, "%10d %10.4f %10.4f",count, energy, rmsd);
                Writer->print_info(info);

                if (Input->write_mol2){
                    Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, string(Input->output + "_" + Lig->molname + "_MC"));
                }

                if (n_rot == 0){
                    gzprintf(mc_output, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f %10.6f %10.6f %10.6f\n", count, energy, rmsd, com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->nconf, step->internal_energy, energy_t->elec, energy_t->vdw, (energy_t->lig_solv + energy_t->rec_solv));
                }

                else {
                    gzprintf(mc_output, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f %10.6f %10.6f %10.6f ", count, energy, rmsd, com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->nconf, step->internal_energy, energy_t->elec, energy_t->vdw, (energy_t->lig_solv + energy_t->rec_solv));
                    for (unsigned i=1; i <= RotorList.Size(); i++){
                        gzprintf(mc_output, "%10.3f", step->torsion_angles[i-1]);
                    }
                    gzprintf(mc_output, "\n");
                }
                independent_average += long(energy);
            }
        }
        else{
            p = this->Boltzmman(energy, new_energy, T, Input->bi);
            rnumber = gsl_rng_uniform(r);
            if (p > rnumber){
                Lig->mcoords = Lig->new_mcoords;
                this->xyz = step->xyz;
                energy = new_energy;
                rmsd = Coord->compute_rmsd(Lig->xyz, step->xyz, Lig->N);
                this->increment_angles(&rot_angles, step);
                this->MaxMinCM(com[0],com[1],com[2],this->MaxMin);
                count++;
                sum_x += energy;
                sum_xsquared += (energy*energy);
                sum_Boltzmann_ene += energy*exp(((-(Input->bi-1.0))*energy)/(k*T));
                sum_Boltzmann2_ene += exp(((-(Input->bi-1.0))*energy)/(k*T));
                sum_Boltzmann2_ene_squared += exp(((-(Input->bi-1.0))*energy)/(k*T)) * exp(((-(Input->bi-1.0))*energy)/(k*T));

                Entropy->update(com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->torsion_angles);

                if (count % Input->mc_stride == 0){

                    sprintf(info, "%10d %10.4f %10.4f",count, energy, rmsd);
                    Writer->print_info(info);

                    if (Input->write_mol2){
                        Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, string(Input->output + "_" + Lig->molname + "_MC"));
                    }
                    if (n_rot == 0){
                        gzprintf(mc_output, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f %10.6f %10.6f %10.6f\n", count, energy, rmsd, com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->nconf, step->internal_energy, energy_t->elec, energy_t->vdw, (energy_t->lig_solv+energy_t->rec_solv));
                    }
                    else {
                        gzprintf(mc_output, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f %10.6f %10.6f %10.6f ", count, energy, rmsd, com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->nconf, step->internal_energy, energy_t->elec, energy_t->vdw, (energy_t->lig_solv+energy_t->rec_solv));
                        for (unsigned i=1; i <= RotorList.Size(); i++){
                            gzprintf(mc_output, "%10.3f", step->torsion_angles[i-1]);
                        }
                        gzprintf(mc_output, "\n");
                    }
                    independent_average += long(energy);
                }
            }
            else{
                nReject++;
                sum_x += energy;
                sum_xsquared += (energy*energy);
                sum_Boltzmann_ene += energy*exp(((-(Input->bi-1.0))*energy)/(k*T));
                sum_Boltzmann2_ene += exp(((-(Input->bi-1.0))*energy)/(k*T));
                sum_Boltzmann2_ene_squared += exp(((-(Input->bi-1.0))*energy)/(k*T)) * exp(((-(Input->bi-1.0))*energy)/(k*T));
            }
        }
    }

    gzclose(mc_output);
    delete step;
    delete Energy;
    delete Coord;
    delete energy_t;
    Writer->print_line();
    sprintf(info, "Finished MC simulation for complex. Acceptance rate: %5.3f", double((Input->number_steps*1.0)/(Input->number_steps+nReject)));
    Writer->print_info(info);
    Writer->print_line();
    sprintf(info, "Conformer energies:");
    Writer->print_info(info);

    for (unsigned i=0; i< Lig->conformer_energies.size(); i++){
        sprintf(info, "%5d %10.3f kcal/mol", i+1, Lig->conformer_energies[i]);
        Writer->print_info(info);
    }
    Writer->print_line();

    this->XSize = this->MaxMin[1] - this->MaxMin[0];
    this->YSize = this->MaxMin[3] - this->MaxMin[2];
    this->ZSize = this->MaxMin[5] - this->MaxMin[4];
    sprintf(info, "Max Dimensions:");
    Writer->print_info(info);
    sprintf(info, "%10.3f %10.3f %10.3f", this->XSize, this->YSize, this->ZSize);
    Writer->print_info(info);
    Writer->print_line();

    double avg_xsquared = sum_xsquared / (Input->number_steps+nReject);
    long double avg_Boltzmann2_ene_squared = sum_Boltzmann2_ene_squared / (Input->number_steps+nReject);
    independent_average = independent_average / (1.0L*Input->number_steps/Input->mc_stride);

    this->average_energy = sum_x/(Input->number_steps+nReject);
    this->energy_standard_deviation = sqrt((avg_xsquared - (this->average_energy*this->average_energy))/(Input->number_steps+nReject-1.0));
    this->Boltzmann_weighted_average_energy = sum_Boltzmann_ene/sum_Boltzmann2_ene;
    this->MCR_Boltzmann_weighted_average = sum_Boltzmann2_ene/(Input->number_steps+nReject);
    this->MCR_Boltzmann_weighted_stdev = sqrt((avg_Boltzmann2_ene_squared - (this->MCR_Boltzmann_weighted_average*this->MCR_Boltzmann_weighted_average)) / (Input->number_steps+nReject-1.0));

    this->average_bound_energy = this->average_energy;

    sprintf(info, "Average Monte Carlo energy: %10.3f +/- %10.3f @ %7.2g K", this->average_energy, this->energy_standard_deviation, T);
    Writer->print_info(info);
    sprintf(info, "Boltzmann-weighted average energy: %10.4Lg @ %7.2g K", this->Boltzmann_weighted_average_energy, T);
    Writer->print_info(info);
    sprintf(info, "Average Monte Carlo energy over independent steps: %10.3Lg @ %7.2g K", independent_average, T);
    Writer->print_info(info);

    Writer->print_line();

    McEntropy::entropy_t* McEnt = new McEntropy::entropy_t;
    McEnt->Srot = 0.0; McEnt->Storsion = 0.0; McEnt->Strans = 0.0;

    McEntropy::entropy_t* Max_Ent = new McEntropy::entropy_t;

    Entropy->get_results(McEnt, Max_Ent, Input->number_steps);

    sprintf(info, "First-Order Approximation Translation Entropy (TS): %10.4g kcal/mol @ %7.2g K", McEnt->Strans*T, T);
    Writer->print_info(info);
    sprintf(info, "First-Order Approximation Rotation Entropy (TS):    %10.4g kcal/mol @ %7.2g K", McEnt->Srot*T, T);
    Writer->print_info(info);
    sprintf(info, "First-Order Approximation Torsion Entropy (TS):     %10.4g kcal/mol @ %7.2g K", McEnt->Storsion*T, T);
    Writer->print_info(info);
    sprintf(info, "First-Order Approximation Total Entropy (S):        %10.4g kcal/(mol.K)@ %7.2g K", McEnt->S, T);
    Writer->print_info(info);
    sprintf(info, "First-Order Approximation -TS (-TS):                %10.4g kcal/mol @ %7.2g K", -McEnt->TS, T);
    Writer->print_info(info);
    sprintf(info, "First-Order Approximation -TS @ 300K:               %10.4g kcal/mol @ %7.2g K", -McEnt->S*300., 300.);
    Writer->print_info(info);

    this->boundTS=McEnt->TS;

    Writer->print_line();

    sprintf(info, "Maximal Entropies Computed for this System:");
    Writer->print_info(info);
    sprintf(info, "First-Order Approximation Translation Entropy (TS): %10.4g kcal/mol @ %7.2g K", Max_Ent->Strans*T, T);
    Writer->print_info(info);
    sprintf(info, "First-Order Approximation Rotation Entropy (TS):    %10.4g kcal/mol @ %7.2g K", Max_Ent->Srot*T, T);
    Writer->print_info(info);
    sprintf(info, "First-Order Approximation Torsion Entropy (TS):     %10.4g kcal/mol @ %7.2g K", Max_Ent->Storsion*T, T);
    Writer->print_info(info);
    sprintf(info, "First-Order Approximation Total Entropy (S):        %10.4g kcal/(mol.K)@ %7.2g K", Max_Ent->S, T);
    Writer->print_info(info);
    sprintf(info, "First-Order Approximation -TS (-TS):                %10.4g kcal/mol @ %7.2g K", -Max_Ent->TS, T);
    Writer->print_info(info);
    sprintf(info, "First-Order Approximation -TS @ 300K:               %10.4g kcal/mol @ %7.2g K", -Max_Ent->S*300., 300.);
    Writer->print_info(info);

    Writer->print_line();

    sprintf(info, "Entropy loss (-TdS): %10.4g kcal/mol (%10.4f %s) @ %7.2g K", (-McEnt->TS - (-Max_Ent->TS)), ((-McEnt->TS/-Max_Ent->TS)*100), "%", T);
    Writer->print_info(info);

    Writer->print_line();

    //Lig->xyz = this->xyz;

    delete McEnt;
    delete Max_Ent;
    delete Entropy;
}


void MC::ligand_run(Mol2* RefLig, Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T){

    if (Input->ligsim){

        double sum_x = 0.0;
        double sum_xsquared = 0.0;
        long double sum_Boltzmann_ene = 0.0L;
        long double sum_Boltzmann2_ene = 0.0L;
        long double sum_Boltzmann2_ene_squared = 0.0L;
        bool ligand_is_in = false;
        long double independent_average = 0.0L;

        this->average_energy = 0.0;
        this->energy_standard_deviation = 0.0;
        this->Boltzmann_weighted_average_energy = 0.0L;
        this->MCR_Boltzmann_weighted_average = 0.0L;
        this->MCR_Boltzmann_weighted_stdev = 0.0L;

        double k=0.0019858775203792202;

        this->MaxMin.push_back(99999.0);         // min x
        this->MaxMin.push_back(-99999.0);        // max x
        this->MaxMin.push_back(99999.0);         // min y
        this->MaxMin.push_back(-99999.0);        // max y
        this->MaxMin.push_back(99999.0);         // min z
        this->MaxMin.push_back(-99999.0);        // max z

        this->xyz = xyz;
        int nReject = 0;

        gsl_rng_set(r, Input->seed);

        sprintf(info, "Ligand Monte Carlo simulation with %5d steps", (Input->number_steps));
        Writer->print_info(info);

        gzFile mc_output_lig;
        mc_output_lig = gzopen(string(Input->output + "_" + Lig->molname + "_mc.ligsim.dat.gz").c_str(), "w");

        Energy2* Energy = new Energy2(Input);
        COORD_MC* Coord = new COORD_MC;
        vector<double> original_com = Coord->compute_com(Lig->xyz, Lig);

        int count=0;
        vector<double> com(3);
        vector<double> rot_angles(3);
        for (unsigned i=0; i<3; i++){
            rot_angles[i] = 0.0;
        }
        com = Coord->compute_com(Lig);

        int n_rot = int(RotorList.Size());

        McEntropy* Entropy = new McEntropy(Input, Coord, original_com, n_rot);

        double energy, new_energy, p, rnumber, rmsd;
        step_t* step = new step_t;

        while (!ligand_is_in){
//            if (Input->generate_conformers){
//                this->take_step_flex(Input, Lig, step);
//            }
//            else
            if (Input->sample_torsions){
                this->take_step_torsion(Input, Lig, step);
            }
            else if(Input->mc_full_flex){
                this->take_step_full_flex(Input, Lig, step);
            }
            else {
                this->take_step(Input, Lig, step);
            }
            com = Coord->compute_com(step->xyz, Lig);
            ligand_is_in = this->ligand_is_inside_box(Input, step, original_com, com);
        }
        ligand_is_in = false;
        energy = step->internal_energy;

        gzprintf(mc_output_lig,"###################################################################################################################################################################\n");
        gzprintf(mc_output_lig,"#BoxsideX= %10.4f BoxsideY= %10.4f BoxsideY= %10.4f Temperature= %10.4f \n", Input->x_dim, Input->y_dim, Input->z_dim, Input->temp);
        gzprintf(mc_output_lig,"#NROT = %10.10d\n", int(RotorList.Size()));
        gzprintf(mc_output_lig,"###################################################################################################################################################################\n");
        if (n_rot == 0){
            gzprintf(mc_output_lig, "#%10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s\n", "Step", "Energy", "RMSD", "DX", "DY", "DZ", "ALPHA", "BETA", "GAMMA", "NCONF", "ConfEnergy");
        }
        else {
            gzprintf(mc_output_lig, "#%10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s ", "Step", "Energy", "RMSD", "DX", "DY", "DZ", "ALPHA", "BETA", "GAMMA", "NCONF", "ConfEnergy");
            for (unsigned i=1; i <= RotorList.Size(); i++){
                sprintf(info, "ROT[%3d]", i);
                gzprintf(mc_output_lig, "%10.10s", info);
            }
            gzprintf(mc_output_lig, "\n");
        }

        Writer->print_line();
        sprintf(info, "%10s %10s %10s", "Step", "Energy", "RMSD");
        Writer->print_info(info);

        while (count <= Input->number_steps){

            while (!ligand_is_in){

//                if (Input->generate_conformers){
//                    this->take_step_flex(Input, Lig, step);
//                }
//                else
                if (Input->sample_torsions){
                    this->take_step_torsion(Input, Lig, step);
                }
                else if (Input->mc_full_flex){
                    this->take_step_full_flex(Input, Lig, step);
                }
                else {
                    this->take_step(Input, Lig, step);
                }
                com = Coord->compute_com(step->xyz, Lig);
                ligand_is_in = this->ligand_is_inside_box(Input, step, original_com, com);
            }
            ligand_is_in = false;

            new_energy = step->internal_energy;

            if (new_energy <= energy){
                Lig->mcoords = Lig->new_mcoords;
                this->xyz = step->xyz;
                energy = new_energy;
                rmsd = Coord->compute_rmsd(Lig->xyz, step->xyz, Lig->N);
                this->increment_angles(&rot_angles, step);
                this->MaxMinCM(com[0], com[1], com[2], this->MaxMin);

                count++;
                sum_x += energy;
                sum_xsquared += (energy*energy);
                sum_Boltzmann_ene += energy*exp(((-(Input->bi-1.0))*energy)/(k*T));
                sum_Boltzmann2_ene += exp((-(Input->bi-1.0))*energy/(k*T));
                sum_Boltzmann2_ene_squared += (exp((-(Input->bi-1.0))*energy/(k*T)))*(exp((-(Input->bi-1.0))*energy/(k*T)));

                Entropy->update(com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->torsion_angles);

                if (count % Input->mc_stride == 0){
                    sprintf(info, "%10d %10.4f %10.4f",count, energy, rmsd);
                    Writer->print_info(info);

                    if (Input->write_mol2){
                        Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, string(Input->output + "_" + Lig->molname + "_MC.ligsim"));
                    }
                    if (n_rot == 0 ){
                        gzprintf(mc_output_lig, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f\n", count, energy, rmsd, com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->nconf, step->internal_energy);
                    }
                    else {
                        gzprintf(mc_output_lig, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f ", count, energy, rmsd, com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->nconf, step->internal_energy);
                        for (unsigned i=1; i <= RotorList.Size(); i++){
                            gzprintf(mc_output_lig, "%10.3f", step->torsion_angles[i-1]);
                        }
                        gzprintf(mc_output_lig, "\n");
                    }
                    independent_average += long(energy);
                }
            }
            else{
                p = this->Boltzmman(energy, new_energy, T, Input->bi);
                rnumber = gsl_rng_uniform(r);
                if (p > rnumber){
                    Lig->mcoords = Lig->new_mcoords;
                    this->xyz = step->xyz;
                    energy = new_energy;
                    rmsd = Coord->compute_rmsd(Lig->xyz, step->xyz, Lig->N);
                    this->increment_angles(&rot_angles, step);
                    this->MaxMinCM(com[0],com[1],com[2],this->MaxMin);

                    count++;
                    sum_x += energy;
                    sum_xsquared += (energy*energy);
                    sum_Boltzmann_ene += energy*exp(((-(Input->bi-1.0))*energy)/(k*T));
                    sum_Boltzmann2_ene += exp((-(Input->bi-1.0))*energy/(k*T));
                    sum_Boltzmann2_ene_squared += (exp((-(Input->bi-1.0))*energy/(k*T)))*(exp((-(Input->bi-1.0))*energy/(k*T)));

                    Entropy->update(com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->torsion_angles);

                    if (count % Input->mc_stride == 0){

                        sprintf(info, "%10d %10.4f %10.4f",count, energy, rmsd);
                        Writer->print_info(info);

                        if (Input->write_mol2){
                            Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, string(Input->output + "_" + Lig->molname + "_MC.ligsim"));
                        }

                        if (n_rot == 0){
                            gzprintf(mc_output_lig, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f\n", count, energy, rmsd, com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->nconf, step->internal_energy);
                        }
                        else {
                            gzprintf(mc_output_lig, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f ", count, energy, rmsd, com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->nconf, step->internal_energy);
                            for (unsigned i=1; i <= RotorList.Size(); i++){
                                gzprintf(mc_output_lig, "%10.3f", step->torsion_angles[i-1]);
                            }
                            gzprintf(mc_output_lig, "\n");
                        }
                        independent_average += long(energy);
                    }
                }
                else{
                    nReject++;
                    sum_x += energy;
                    sum_xsquared += (energy*energy);
                    sum_Boltzmann_ene += energy*exp(((-(Input->bi-1.0))*energy)/(k*T));
                    sum_Boltzmann2_ene += exp((-(Input->bi-1.0))*energy/(k*T));
                    sum_Boltzmann2_ene_squared += (exp((-(Input->bi-1.0))*energy/(k*T)))*(exp((-(Input->bi-1.0))*energy/(k*T)));
                }
            }
        }

        gzclose(mc_output_lig);
        delete step;
        delete Energy;
        delete Coord;
        Writer->print_line();
        sprintf(info, "Finished MC simulation for ligand. Acceptance rate: %5.3f", (Input->number_steps*1.0)/(Input->number_steps+nReject));
        Writer->print_info(info);
        Writer->print_line();

        this->XSize = this->MaxMin[1] - this->MaxMin[0];
        this->YSize = this->MaxMin[3] - this->MaxMin[2];
        this->ZSize = this->MaxMin[5] - this->MaxMin[4];
        sprintf(info, "Max Dimensions:");
        Writer->print_info(info);
        sprintf(info, "%10.3f %10.3f %10.3f", this->XSize, this->YSize, this->ZSize);
        Writer->print_info(info);
        Writer->print_line();

        double avg_xsquared = sum_xsquared / (Input->number_steps+nReject);
        long double avg_Boltzmann2_ene_squared = sum_Boltzmann2_ene_squared / (Input->number_steps+nReject);
        independent_average = independent_average / (1.0L*Input->number_steps/Input->mc_stride);

        this->average_energy = double(sum_x/(Input->number_steps+nReject));
        this->energy_standard_deviation = sqrt((avg_xsquared - (this->average_energy*this->average_energy))/(Input->number_steps+nReject-1.0));
        this->Boltzmann_weighted_average_energy = sum_Boltzmann_ene/sum_Boltzmann2_ene;
        this->MCR_Boltzmann_weighted_average = sum_Boltzmann2_ene/(Input->number_steps+nReject);
        this->MCR_Boltzmann_weighted_stdev = sqrt((avg_Boltzmann2_ene_squared - (this->MCR_Boltzmann_weighted_average*this->MCR_Boltzmann_weighted_average))/(Input->number_steps+nReject-1.0));

        this->average_freeligand_energy=this->average_energy;

        sprintf(info, "Average Monte Carlo ligand energy: %10.3f +- %10.3g @ %7.2g K", this->average_energy, this->energy_standard_deviation, T);
        Writer->print_info(info);
        sprintf(info, "Boltzmann-weighted average ligand energy: %10.3Lg @ %7.2g K", this->Boltzmann_weighted_average_energy, T);
        Writer->print_info(info);
        sprintf(info, "Average Monte Carlo ligand energy over independent steps: %10.3Lg @ %7.2g K", independent_average, T);
        Writer->print_info(info);

        Writer->print_line();

        McEntropy::entropy_t* McEnt = new McEntropy::entropy_t;

        McEntropy::entropy_t* Max_Ent = new McEntropy::entropy_t;

        McEnt->Srot = 0.0; McEnt->Storsion = 0.0; McEnt->Strans = 0.0;
        Entropy->get_results(McEnt, Max_Ent, Input->number_steps);

        sprintf(info, "First-Order Approximation Ligand Translation Entropy (TS): %10.4g kcal/mol @ %7.2g K", McEnt->Strans*T, T);
        Writer->print_info(info);
        sprintf(info, "First-Order Approximation Ligand Rotation Entropy (TS):    %10.4g kcal/mol @ %7.2g K", McEnt->Srot*T, T);
        Writer->print_info(info);
        sprintf(info, "First-Order Approximation Ligand Torsion Entropy (TS):     %10.4g kcal/mol @ %7.2g K", McEnt->Storsion*T, T);
        Writer->print_info(info);
        sprintf(info, "First-Order Approximation Ligand Total Entropy (S):        %10.4g kcal/(mol.K)@ %7.2g K", McEnt->S, T);
        Writer->print_info(info);
        sprintf(info, "First-Order Approximation Ligand -TS (-TS):                %10.4g kcal/mol @ %7.2g K", -McEnt->TS, T);
        Writer->print_info(info);
        sprintf(info, "First-Order Approximation Ligand -TS @ 300K:               %10.4g kcal/mol @ %7.2g K", -McEnt->S*300., 300.);
        Writer->print_info(info);

        this->freeTS=McEnt->TS;

        Writer->print_line();

        sprintf(info, "Maximal Entropies Computed for this System:");
        Writer->print_info(info);
        sprintf(info, "First-Order Approximation Ligand Translation Entropy (TS): %10.4g kcal/mol @ %7.2g K", Max_Ent->Strans*T, T);
        Writer->print_info(info);
        sprintf(info, "First-Order Approximation Ligand Rotation Entropy (TS):    %10.4g kcal/mol @ %7.2g K", Max_Ent->Srot*T, T);
        Writer->print_info(info);
        sprintf(info, "First-Order Approximation Ligand Torsion Entropy (TS):     %10.4g kcal/mol @ %7.2g K", Max_Ent->Storsion*T, T);
        Writer->print_info(info);
        sprintf(info, "First-Order Approximation Ligand Total Entropy (S):        %10.4g kcal/(mol.K)@ %7.2g K", Max_Ent->S, T);
        Writer->print_info(info);
        sprintf(info, "First-Order Approximation Ligand -TS (-TS):                %10.4g kcal/mol @ %7.2g K", -Max_Ent->TS, T);
        Writer->print_info(info);
        sprintf(info, "First-Order Approximation Ligand -TS @ 300K:               %10.4g kcal/mol @ %7.2g K", -Max_Ent->S*300., 300.);
        Writer->print_info(info);

        Writer->print_line();

        sprintf(info, "Entropy loss (-TdS): %10.4g kcal/mol (%10.4f %s) @ %7.2g K", (-McEnt->TS - (-Max_Ent->TS)), ((-McEnt->TS/-Max_Ent->TS)*100.), "%", T);
        Writer->print_info(info);

        Writer->print_line();

        delete McEnt;
        delete Max_Ent;
        delete Entropy;
    }
}


void MC::run(Mol2* Rec, Mol2* RefLig, Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T){
    this->xyz = xyz;
    int nReject = 0;

    double sum_x = 0.0;
    double sum_xsquared = 0.0;
    long double sum_Boltzmann_ene = 0.0;
    long double sum_Boltzmann_ene_squared = 0.0;

    double k=0.0019858775203792202;

    this->MaxMin.push_back(99999.0);         // min x
    this->MaxMin.push_back(-99999.0);        // max x
    this->MaxMin.push_back(99999.0);         // min y
    this->MaxMin.push_back(-99999.0);        // max y
    this->MaxMin.push_back(99999.0);         // min z
    this->MaxMin.push_back(-99999.0);        // max z

    gsl_rng_set(r, Input->seed);

    if (Input->eq_mode){
        sprintf(info, "Starting Monte Carlo equilibration simulation with %5d steps...", (Input->eq_steps));
        Writer->print_info(info);
        gzFile mc_output;
        mc_output = gzopen((Input->output + "_mc.dat.gz").c_str(), "w");

        Energy2* Energy = new Energy2(Input);
        COORD_MC* Coord = new COORD_MC;
        double energy, new_energy, p, rnumber, rmsd;
        step_t* step = new step_t;

        if (Input->generate_conformers){
            this->take_step_flex(Input, Lig, step);
        }
        else if (Input->sample_torsions){
            this->take_step_torsion(Input, Lig, step);
        }
        else {
            this->take_step(Input, Lig, step);
        }

        energy = (Energy->compute_ene(Rec, Lig, step->xyz)+Lig->conformer_energies[step->nconf]);



        sprintf(info, "#%10s %10s %10s", "Step", "Energy", "RMSD");
        gzprintf(mc_output,"###################################################################################################################################################################\n");
        gzprintf(mc_output,"#BoxsideX= %10.4f BoxsideY= %10.4f BoxsideY= %10.4f Temperature= %10.4f \n", Input->x_dim, Input->y_dim, Input->z_dim, Input->temp);
        gzprintf(mc_output,"#NROT = %10.10d\n", int(RotorList.Size()));
        gzprintf(mc_output,"###################################################################################################################################################################\n");
        if (! Input->sample_torsions){
            gzprintf(mc_output, "#%10.10s 10.10s 10.10s 10.10s 10.10s 10.10s 10.10s 10.10s 10.10s 10.10s 10.10s\n", "Step", "Energy", "RMSD", "DX", "DY", "DZ", "ALPHA", "BETA", "GAMMA", "NCONF", "ConfEnergy");
        }
        else {
            gzprintf(mc_output, "#%10.10s 10.10s 10.10s 10.10s 10.10s 10.10s 10.10s 10.10s 10.10s 10.10s 10.10s ", "Step", "Energy", "RMSD", "DX", "DY", "DZ", "ALPHA", "BETA", "GAMMA", "NCONF", "ConfEnergy");
            for (unsigned i=1; i <= RotorList.Size(); i++){
                sprintf(info, "Torsion[%3d]", i);
                gzprintf(mc_output, "%10.10s", info);
            }
            gzprintf(mc_output, "\n");
        }

        int count=0;
        int eqcount = 0;
        vector<double> com = Coord->compute_com(Lig);
        vector<double> rot_angles(3);
        for (unsigned i=0; i< rot_angles.size(); i++){
            rot_angles[i] = 0.0;
        }

        //Equilibration implementation
        while (eqcount <= Input->eq_steps){

//            if (Input->generate_conformers){
//                this->take_step_flex(Input, Lig, step);
//            }
//            else
            if (Input->sample_torsions){
                this->take_step_torsion(Input, Lig, step);
            }
            else if (Input->mc_full_flex){
                this->take_step_full_flex(Input, Lig, step);
            }
            else {
                this->take_step(Input, Lig, step);
            }
            new_energy = (Energy->compute_ene(Rec, Lig, step->xyz)+Lig->conformer_energies[step->nconf]);

            if (new_energy <= energy){
                Lig->mcoords = Lig->new_mcoords;
                this->xyz = step->xyz;
                energy = new_energy;
                rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                com[0] += step->dx;
                com[1] += step->dy;
                com[2] += step->dz;
                this->increment_angles(&rot_angles, step);
                this->MaxMinCM(com[0],com[1],com[2],this->MaxMin);
                eqcount++;
            }
            else{
                p = this->Boltzmman(energy, new_energy, T, Input->bi);
                rnumber = gsl_rng_uniform(r);
                if (p > rnumber){
                    Lig->mcoords = Lig->new_mcoords;
                    this->xyz = step->xyz;
                    energy = new_energy;
                    rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                    com[0] += step->dx;
                    com[1] += step->dy;
                    com[2] += step->dz;
                    this->increment_angles(&rot_angles, step);
                    this->MaxMinCM(com[0],com[1],com[2],this->MaxMin);
                    eqcount++;
                }
            }
        }

        Writer->print_line();
        sprintf(info, "Equilibration done with %5d steps. Current system energy: %9.3f kcal/mol.", Input->eq_steps, energy);
        Writer->print_info(info);
        Writer->print_line();
        sprintf(info, "%10s %10s %10s", "Step", "Energy", "RMSD");
        Writer->print_info(info);

        while (count <= Input->number_steps){

//            if (Input->generate_conformers){
//                this->take_step_flex(Input, Lig, step);
//            }
//            else {
            if (Input->sample_torsions){
                this->take_step_torsion(Input, Lig, step);
            }
            else if (Input->mc_full_flex){
                this->take_step_full_flex(Input, Lig, step);
            }
            else {
                this->take_step(Input, Lig, step);
            }

            new_energy = (Energy->compute_ene(Rec, Lig, step->xyz)+Lig->conformer_energies[step->nconf]);

            if (new_energy <= energy){
                Lig->mcoords = Lig->new_mcoords;
                this->xyz = step->xyz;
                energy = new_energy;
                rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                com[0] += step->dx;
                com[1] += step->dy;
                com[2] += step->dz;
                this->increment_angles(&rot_angles, step);
                this->MaxMinCM(com[0],com[1],com[2],this->MaxMin);
                count++;
                sum_x += energy;
                sum_xsquared += (energy*energy);
                sum_Boltzmann_ene += exp(((-Input->bi-1.0)*energy)/(k*T));
                sum_Boltzmann_ene_squared = (exp(((-Input->bi-1.0)*energy)/(k*T)))*(exp(((-Input->bi-1.0)*energy)/(k*T)));
            }
            else{
                p = this->Boltzmman(energy, new_energy, T, Input->bi);
                rnumber = gsl_rng_uniform(r);
                if (p > rnumber){
                    Lig->mcoords = Lig->new_mcoords;
                    this->xyz = step->xyz;
                    energy = new_energy;
                    rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                    com[0] += step->dx;
                    com[1] += step->dy;
                    com[2] += step->dz;
                    this->increment_angles(&rot_angles, step);
                    this->MaxMinCM(com[0],com[1],com[2],this->MaxMin);
                    count++;
                    sum_x += energy;
                    sum_xsquared += (energy*energy);
                    sum_Boltzmann_ene += exp(((-Input->bi-1.0)*energy)/(k*T));
                    sum_Boltzmann_ene_squared = (exp(((-Input->bi-1.0)*energy)/(k*T)))*(exp(((-Input->bi-1.0)*energy)/(k*T)));
                }
                else{
                    nReject++;
                    sum_x += energy;
                    sum_xsquared += (energy*energy);
                }
            }

            if ((Input->write_mol2) and (count % Input->mc_stride == 0)){
                Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, string(Input->output + "_" + Lig->molname + "_MC"));
            }
            if (! Input->sample_torsions){
                gzprintf(mc_output, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f\n", count, energy, rmsd, com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->nconf, step->internal_energy);
            }
            else {
                gzprintf(mc_output, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f ", count, energy, rmsd, com[0], com[1], com[2], rot_angles[0], rot_angles[1], rot_angles[2], step->nconf, step->internal_energy);
                for (unsigned i=1; i <= RotorList.Size(); i++){
                    gzprintf(mc_output, "%10.3f", step->torsion_angles[i-1]);
                }
                gzprintf(mc_output, "\n");
            }
            sprintf(info, "%10d %10.3f %10.3f", count, energy, rmsd);
            count++;
            if (count % Input->mc_stride == 0){
                sprintf(info, "%10d %10.4f %10.4f",count, energy, rmsd);
                Writer->print_info(info);
            }
        }

        gzclose(mc_output);
        delete step;
        delete Energy;
        delete Coord;
        Writer->print_line();
        sprintf(info, "Finished MC simulation after %d accepted steps, %d total steps, acceptance rate %5.3f", Input->number_steps, (Input->number_steps+nReject),double((Input->number_steps*1.0)/(Input->number_steps+nReject)));
        Writer->print_info(info);
        Writer->print_line();
        sprintf(info, "Conformer energies:");
        Writer->print_info(info);

        for (unsigned i=0; i< Lig->conformer_energies.size(); i++){
            sprintf(info, "%5d %10.3f kcal/mol", i, Lig->conformer_energies[i]);
            Writer->print_info(info);
        }
        Writer->print_line();

        this->XSize = this->MaxMin[1] - this->MaxMin[0];
        this->YSize = this->MaxMin[3] - this->MaxMin[2];
        this->ZSize = this->MaxMin[5] - this->MaxMin[4];
        sprintf(info, "Max Dimensions:");
        Writer->print_info(info);
        sprintf(info, "%10.3f %10.3f %10.3f", this->XSize, this->YSize, this->ZSize);
        Writer->print_info(info);
        Writer->print_line();

        this->average_energy = double(sum_x/(Input->number_steps+nReject));
        this->energy_standard_deviation = (sum_xsquared - (sum_x*sum_x))/(Input->number_steps+nReject-1);
        this->energy_standard_deviation = sqrt(this->energy_standard_deviation);
        this->Boltzmann_weighted_average_energy = sum_Boltzmann_ene/(Input->number_steps+nReject);
        this->MCR_Boltzmann_weighted_stdev = (sum_Boltzmann_ene_squared - (sum_Boltzmann_ene*sum_Boltzmann_ene))/(Input->number_steps+nReject-1);
        this->MCR_Boltzmann_weighted_stdev  = sqrt(this->MCR_Boltzmann_weighted_stdev);

        sprintf(info, "Average Monte Carlo energy: %10.3f kcal/mol +- (%10.3f kcal/mol) @ %7.2f K", this->average_energy, this->energy_standard_deviation, T);
        Writer->print_info(info);
        sprintf(info, "Boltzmann-weighted average energy: %10.3Lf kcal/mol @ %7.2f K", this->Boltzmann_weighted_average_energy, T);
        Writer->print_info(info);
        Writer->print_line();
    }
}


void MC::MaxMinCM(double XCOM, double YCOM, double ZCOM, vector <double> Max){

    if (XCOM < Max[0]){
        this->MaxMin[0]=XCOM;
    }
    if (XCOM > Max[1]){
        this->MaxMin[1]=XCOM;
    }
    if (YCOM < Max[2]){
        this->MaxMin[2]=YCOM;
    }
    if (YCOM > Max[3]){
        this->MaxMin[3]=YCOM;
    }
    if (ZCOM < Max[4]){
        this->MaxMin[4]=ZCOM;
    }
    if (ZCOM > Max[5]){
        this->MaxMin[5]=ZCOM;
    }
}


void MC::take_step(PARSER* Input, Mol2* Lig, step_t* step){
    COORD_MC* Coord = new COORD_MC;
    double rnumber; //transx, transy, transz, a, b, g, rnumber;

    rnumber = gsl_rng_uniform(r);
    step->dx = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dy = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dz = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));

    rnumber = gsl_rng_uniform(r);
    step->dalpha = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dbeta = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dgamma = -Input->rotation_step + (rnumber*(2*Input->rotation_step));

    step->nconf = 0;

    step->xyz = Coord->rototranslate(this->xyz, Lig, step->dalpha, step->dbeta,step->dgamma, step->dx, step->dy, step->dz);


    delete Coord;
}

void MC::take_step_flex(PARSER* Input, Mol2* Lig, step_t* step){
    COORD_MC* Coord = new COORD_MC;
    double rnumber; //transx, transy, transz, a, b, g, rnumber;
    int ln=0;

    rnumber = gsl_rng_uniform(r);
    step->dx = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dy = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dz = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));

    rnumber = gsl_rng_uniform(r);
    step->dalpha = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dbeta = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dgamma = -Input->rotation_step + (rnumber*(2*Input->rotation_step));

    Coord->rototranslate_all(Lig, step->dalpha, step->dbeta, step->dgamma, step->dx, step->dy, step->dz);

    ln = gsl_rng_uniform_int(r, Lig->mcoords.size());
    step->nconf = ln;

    delete Coord;
    step->xyz = Lig->new_mcoords[ln];
    step->internal_energy = Lig->conformer_energies[ln];
}

void MC::take_step_torsion(PARSER* Input, Mol2* Lig, step_t* step){

    COORD_MC* Coord = new COORD_MC;
    double rnumber;
    step->torsion_angles.clear();

// Do rotation and translation

    rnumber = gsl_rng_uniform(r);
    step->dx = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dy = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dz = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));

    rnumber = gsl_rng_uniform(r);
    step->dalpha = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dbeta = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dgamma = -Input->rotation_step + (rnumber*(2*Input->rotation_step));

    step->xyz = Coord->rototranslate(this->xyz, Lig, step->dalpha, step->dbeta,step->dgamma, step->dx, step->dy, step->dz);


// Copy coordinates to OBMol

    this->copy_to_obmol(step->xyz);

    delete Coord;

// Do torsion search

    double current_angle, new_angle;
    for (unsigned i=0; i< RotorList.Size(); i++){
        rnumber = gsl_rng_uniform(r);
        current_angle = mol.GetTorsion(mol.GetAtom(atoms_in_dihedrals[i][0]), mol.GetAtom(atoms_in_dihedrals[i][1]), mol.GetAtom(atoms_in_dihedrals[i][2]), mol.GetAtom(atoms_in_dihedrals[i][3]));
        new_angle = (current_angle + (-Input->torsion_step + (rnumber*(2*Input->torsion_step))));
        new_angle = this->check_angle(new_angle);
        mol.SetTorsion(mol.GetAtom(atoms_in_dihedrals[i][0]), mol.GetAtom(atoms_in_dihedrals[i][1]), mol.GetAtom(atoms_in_dihedrals[i][2]), mol.GetAtom(atoms_in_dihedrals[i][3]), new_angle*PI/180.);
        step->torsion_angles.push_back(new_angle);
    }


// copy coordinates and internal energy to type step_t

    step->xyz = this->copy_from_obmol(mol);
    OBff->Setup(mol);
    step->internal_energy = OBff->Energy();
    string unit = OBff->GetUnit();
    if (unit == "kJ/mol"){
        step->internal_energy = step->internal_energy/4.18;
    }
    step->nconf = 0;
}



double MC::Boltzmman(double ene, double new_ene, double t, double b){
    double de = (new_ene - ene);
    double x = (-(b-1.0)*de) / (0.0019858775203792202*t); // k=0.0019858775203792202 kcal/(mol.K)
    return(exp(x));
}

vector<vector<double> > MC::copy_from_obmol(OBMol mymol){
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

void MC::copy_to_obmol(vector<vector<double> > vec_xyz){
    double* dxyz = new double[vec_xyz.size()*3];
    for (unsigned i=0; i<vec_xyz.size(); i++){
        dxyz[3*i] = vec_xyz[i][0];
        dxyz[(3*i)+1] = vec_xyz[i][1];
        dxyz[(3*i)+2] = vec_xyz[i][2];
    }
    mol.SetCoordinates(dxyz);
    delete [] dxyz;
}

OBMol MC::GetMol(const std::string &molfile){
    OBMol mol;
    OBConversion conv;
    OBFormat *format = conv.FormatFromExt(molfile.c_str());
    if (!format || !conv.SetInFormat(format)) {
    printf("Could not find input format for file %s\n", molfile.c_str());
    return mol;
  }

    ifstream ifs(molfile.c_str());
    if (!ifs) {
        printf("Could not open %s for reading.\n", molfile.c_str());
        return mol;
    }

    if (!conv.Read(&mol, &ifs)) {
        printf("Could not read molecule from file %s\n", molfile.c_str());
        return mol;
    }
    return mol;
}

double MC::check_angle(double angle){
    if ( angle > 360.0){
        while (angle > 360.0){
            angle -= 360.0;
        }
    }
    else if (angle < 0.0){
        while (angle < 0.0){
            angle += 360.0;
        }
    }
    return angle;
}

void MC::increment_angles(vector<double> *angles, step_t* step){
        angles->at(0) += step->dalpha;
        angles->at(1) += step->dbeta;
        angles->at(2) += step->dgamma;

        for (unsigned i=0; i < 3; i++){
            if (angles->at(i) > 360.){
                angles->at(i) -= 360.;
            }
            else if (angles->at(i) < -360.){
                angles->at(i) += 360.;
            }

            if (angles->at(i) < 0.0){       //make all angles positive, in the range 0-360 degres.
                angles->at(i) += 360.;
            }
        }

        if (angles->at(1) > 180.){
            angles->at(1) -= 180.;
        }
        else if (angles->at(1) < -180.){
            angles->at(1) += 180.;
        }

        if (angles->at(1) < 0.0){       //beta is within the range of 0-180 degres.
            angles->at(1) += 180.;
        }
}

bool MC::ligand_is_inside_box(PARSER* Input, step_t* step, vector<double> original_com, vector<double> current_com){
    bool ret = false;
    if (abs(current_com[0]-original_com[0]) <= (Input->search_box_x/2.0)){
        if (abs(current_com[1]-original_com[1]) <= (Input->search_box_y/2.0)){
            if (abs(current_com[2]-original_com[2]) <= (Input->search_box_z/2.0)){
                ret=true;
            }
        }
    }
    return ret;
}

void MC::take_step_full_flex(PARSER* Input, Mol2* Lig, step_t* step){
    COORD_MC* Coord = new COORD_MC;
    double rnumber;
    double current_angle;
    step->torsion_angles.clear();
    vector<vector<double> > ref_xyz = this->xyz;

    // Do rigid body rotation and translation

    rnumber = gsl_rng_uniform(r);
    step->dx = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dy = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dz = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));

    rnumber = gsl_rng_uniform(r);
    step->dalpha = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dbeta = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dgamma = -Input->rotation_step + (rnumber*(2*Input->rotation_step));

    step->xyz = Coord->rototranslate(this->xyz, Lig, step->dalpha, step->dbeta,step->dgamma, step->dx, step->dy, step->dz);

    // Now, let's do a random shift in the internal atomic coordinates

    double dx, dy, dz;
    for (int i=0; i< Lig->N; i++){
        rnumber = gsl_rng_uniform(r);
        dx = -(Input->max_atom_displacement) + (1.0 * (rnumber*(2*Input->max_atom_displacement)));
        rnumber = gsl_rng_uniform(r);
        dy = -(Input->max_atom_displacement) + (1.0 * (rnumber*(2*Input->max_atom_displacement)));
        rnumber = gsl_rng_uniform(r);
        dz = -(Input->max_atom_displacement) + (1.0 * (rnumber*(2*Input->max_atom_displacement)));

        step->xyz[i][0] = step->xyz[i][0] + dx;
        step->xyz[i][1] = step->xyz[i][1] + dy;
        step->xyz[i][2] = step->xyz[i][2] + dz;
    }

// copy coordinates and internal energy to type step_t

    this->copy_to_obmol(step->xyz);

// compute torsion angles with Babel API

    for (unsigned i=0; i< RotorList.Size(); i++){
            current_angle = mol.GetTorsion(mol.GetAtom(atoms_in_dihedrals[i][0]), mol.GetAtom(atoms_in_dihedrals[i][1]), mol.GetAtom(atoms_in_dihedrals[i][2]), mol.GetAtom(atoms_in_dihedrals[i][3]));
            current_angle = this->check_angle(current_angle);
            step->torsion_angles.push_back(current_angle);
    }

    delete Coord;

    OBff->Setup(mol);
    step->internal_energy = OBff->Energy();
    string unit = OBff->GetUnit();
    if (unit == "kJ/mol"){
        step->internal_energy = step->internal_energy/4.18;
    }
    step->nconf = 0;
}


/*
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;


BOOST_PYTHON_MODULE(pyMC)
{

    void (MC::*r1)(Mol2 *Rec, Mol2* Reflig , Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T) =&MC::run;
    void (MC::*r2)(Grid* Grids, Mol2* RefLig, Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T) =&MC::run;

    class_<MC>("MC", init<WRITER*>())
        .def(init<Mol2*, PARSER*, WRITER*>())
        .def_readwrite("XSize", & MC::XSize)
        .def_readwrite("YSize", & MC::YSize)
        .def_readwrite("ZSize", & MC::ZSize)
        .def_readwrite("MaxMin", & MC::ZSize)
        .def_readwrite("xyz", & MC::xyz)
        .def_readwrite("average_energy", & MC::average_energy)
        .def_readwrite("energy_standard_deviation", & MC::energy_standard_deviation)
        .def_readwrite("Boltzmann_weighted_average_energy", & MC::Boltzmann_weighted_average_energy)
        .def_readwrite("MCR_Boltzmann_weighted_average", & MC::MCR_Boltzmann_weighted_average)
        .def_readwrite("MCR_Boltzmann_weighted_stdev", & MC::MCR_Boltzmann_weighted_stdev)
        .def_readwrite("average_bound_energy", & MC::average_bound_energy)
        .def_readwrite("average_freeligand_energy", & MC::average_freeligand_energy)
        .def_readwrite("average_deltaE", & MC::average_deltaE)
        .def_readwrite("boundTS", & MC::boundTS)
        .def_readwrite("freeTS", & MC::freeTS)

        .def_readwrite("r", & MC::r)
        .def_readwrite("Writer", & MC::Writer)
        .def_readwrite("info", & MC::info)
        .def_readwrite("myxyz", & MC::myxyz)

        .def_readwrite("mol", & MC::mol)
        .def_readwrite("OBff", & MC::OBff)
        .def_readwrite("RotorList", & MC::RotorList)
        .def_readwrite("RotorIterator", & MC::RotorIterator)
        .def_readwrite("Rotor", & MC::Rotor)

        .def("run", r1)
        .def("run", r2)

        .def("Boltzmman", & MC::Boltzmman)
        .def("take_step", & MC::take_step)
        .def("take_step_flex", & MC::take_step_flex)
        .def("take_step_torsion", & MC::take_step_torsion)
        .def("take_step_full_flex", & MC::take_step_full_flex)
        .def("write_conformers", & MC::write_conformers)
        .def("MaxMinCM", & MC::MaxMinCM)
        .def("ligand_run", & MC::ligand_run)

        .def("copy_to_obmol", & MC::copy_to_obmol)
        .def("copy_from_obmol", & MC::copy_from_obmol)

        .def("GetMol", & MC::GetMol)
        .def("check_angle", & MC::check_angle)
        .def("increment_angles", & MC::increment_angles)
        .def("ligand_is_inside_box", & MC::ligand_is_inside_box)


    ;

    class_<MC::step_t>("step_t")
        .def_readwrite("xyz", & MC::step_t::xyz)
        .def_readwrite("dx", & MC::step_t::dx)
        .def_readwrite("dy ", & MC::step_t::dy)
        .def_readwrite("dz", & MC::step_t::dz)
        .def_readwrite("dalpha", & MC::step_t::dalpha)
        .def_readwrite("dbeta", & MC::step_t::dbeta)
        .def_readwrite("dgamma", & MC::step_t::dgamma)
        .def_readwrite("nconf", & MC::step_t::nconf)
        .def_readwrite("torsion_angles", & MC::step_t::torsion_angles)
        .def_readwrite("internal_energy", & MC::step_t::internal_energy)
    ;
}
*/
