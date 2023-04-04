#include "pyEnergy2.h"
#include "pyGrid.cpp"

Energy2::Energy2(PARSER* _Input)
{
    this->Input = _Input;
}

double Energy2::distance(double x1, double x2, double y1, double y2, double z1, double z2)
{
    return ( sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) );
}

double Energy2::distance_squared(double x1, double x2, double y1, double y2, double z1, double z2)
{
    return (((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1)));
}

double Energy2::angle(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3){
    //
    // Law of cosines: c** = a** + b** - 2ab cos(C)
    //
    double ab = distance(x1, x2, y1, y2, z1, z2);
    if (ab == 0.0)
        ab=0.0001 ;
    double ac = distance(x1, x3, y1, y3, z1, z3);
    if (ac == 0.0)
        ac=0.0001 ;
    double bc = distance(x2, x3, y2, y3, z2, z3);
    if (bc == 0.0)
        bc=0.0001 ;
    double angle = acos(((ab*ab)+(bc*bc)-(ac*ac))/(2*ab*bc));
    angle = angle * 180.0 / PI;
    return (angle);
}


bool Energy2::atom_is_acceptor(int a, Mol2* Lig){
    bool accept = false;
    for (unsigned i=0; i<Lig->HBacceptors.size(); i++){
        if (Lig->HBacceptors[i] == a){
            accept = true;
        }
    }
    return accept;
}

int Energy2::H_is_donor(int a, Mol2* Lig){
    int donor = -1;
    if (Lig->HBdonors.size() > 0){
        for (unsigned i=0; i<Lig->HBdonors.size(); i++){
            if (Lig->HBdonors[i][1] == a){
                donor = Lig->HBdonors[i][0];
            }
        }
    }
    return donor;
}

double Energy2::compute_ene(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz){
    double ene = 0.0;
    switch(Input->scoring_function){
    case 0:
        ene = this->compute_energy_softcore_solvation(Rec, Lig, lig_xyz);
        break;
    case 1:
        ene = this->compute_energy_softcore(Rec, Lig, lig_xyz);
        break;
    case 2:
        ene = this->compute_energy_hardcore_solvation(Rec, Lig, lig_xyz);
        break;
    case 3:
        ene = this->compute_energy_hardcore(Rec, Lig, lig_xyz);
        break;
    default:
        ene = this->compute_energy_hardcore_solvation(Rec, Lig, lig_xyz);
        break;
    }
    return(ene);
}

double Energy2::compute_ene(Grid* Grids, Mol2* Lig, vector<vector<double> > lig_xyz){
    double ene=0.0;
    switch(Input->scoring_function){
    case 0:
        ene = this->compute_ene_from_grids_softcore_solvation(Grids, Lig, lig_xyz);
        break;
    case 1:
        ene = this->compute_ene_from_grids_softcore(Grids, Lig, lig_xyz);
        break;
    case 2:
        ene = this->compute_ene_from_grids_hardcore_solvation(Grids, Lig, lig_xyz);
        break;
    case 3:
        ene = this->compute_ene_from_grids_hardcore(Grids, Lig, lig_xyz);
        break;
    case 4:
        ene = this->compute_ene_from_grids_hardcore_solvation(Grids, Lig, lig_xyz);
        break;
    case 5:
        ene = this->compute_ene_from_grids_hardcore(Grids, Lig, lig_xyz);
        break;
    }
    return(ene);
}

double Energy2::compute_ene(Grid* Grids, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result){
    double e=0.0;
    switch(Input->scoring_function){
    case 0:
        e = this->compute_ene_from_grids_softcore_solvation(Grids, Lig, lig_xyz, energy_result);
        break;
    case 1:
        e = this->compute_ene_from_grids_softcore(Grids, Lig, lig_xyz, energy_result);
        break;
    case 2:
        e = this->compute_ene_from_grids_hardcore_solvation(Grids, Lig, lig_xyz, energy_result);
        break;
    case 3:
        e = this->compute_ene_from_grids_hardcore(Grids, Lig, lig_xyz, energy_result);
        break;
    case 4:
        e = this->compute_ene_from_grids_hardcore_solvation(Grids, Lig, lig_xyz, energy_result);
        break;
    case 5:
        e = this->compute_ene_from_grids_hardcore(Grids, Lig, lig_xyz, energy_result);
        break;
    }
    return (e);
}

double Energy2::compute_ene(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result){
    double e;
    switch(Input->scoring_function){
    case 0:
        e = this->compute_energy_softcore_solvation(Rec, Lig, lig_xyz, energy_result);
        break;
    case 1:
        e = this->compute_energy_softcore(Rec, Lig, lig_xyz, energy_result);
        break;
    case 2:
        e = this->compute_energy_hardcore_solvation(Rec, Lig, lig_xyz, energy_result);
        break;
    case 3:
        e = this->compute_energy_hardcore(Rec, Lig, lig_xyz, energy_result);
        break;
    default:
        e = this->compute_energy_softcore_solvation(Rec, Lig, lig_xyz, energy_result);
        break;
    }
    return (e);
}

double Energy2::compute_energy_softcore_solvation(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz){
    double rec_solv = 0.00;
    double lig_solv = 0.00;
    double vdw=0.0;
    double elec=0.0;
    double acoef, bcoef, lig_solv_affinity, lig_solv_distf, rec_solv_distf, rec_solv_affinity, dij, dij2, deff, dij6;//, rij, eij;
    double sqrt2 = sqrt(2);

    for (int i=0; i< Rec->N; i++){
        rec_solv_affinity = (Input->solvation_alpha*((Rec->charges[i])*(Rec->charges[i])))+Input->solvation_beta;
        for (int j=0; j< Lig->N; j++){

            //! distances

            dij2 = this->distance_squared(Rec->xyz[i][0], lig_xyz[j][0],Rec->xyz[i][1], lig_xyz[j][1], Rec->xyz[i][2], lig_xyz[j][2]);
            dij6 = dij2*dij2*dij2;
            dij = sqrt(dij2);
            deff = pow(((dij*dij*dij)+Input->deltaij_es3), (1.0/3.0));

            //! electrostactic energy (softcore)
            if (Input->dielectric_model == "constant"){
                elec+= (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*deff);

                lig_solv_affinity = (Input->solvation_alpha*((Lig->charges[j])*(Lig->charges[j])))+Input->solvation_beta;
                lig_solv_distf = ((4.0/3.0)*PI*(Rec->radii[i]*Rec->radii[i]*Rec->radii[i]));
                lig_solv_distf = (lig_solv_distf * exp(-(deff*deff)/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
                rec_solv_distf = ((4.0/3.0)*PI*(Lig->radii[j]*Lig->radii[j]*Lig->radii[j]));
                rec_solv_distf = (rec_solv_distf * exp(-(deff*deff)/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
            }
            else {
                deff = pow(((dij6)+Input->deltaij_es6), (1.0/3.0));
                elec+= (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*(pow((dij6+Input->deltaij_es6),(1.0/3.0))));

                lig_solv_affinity = (Input->solvation_alpha*((Lig->charges[j])*(Lig->charges[j])))+Input->solvation_beta;
                lig_solv_distf = ((4.0/3.0)*PI*(Rec->radii[i]*Rec->radii[i]*Rec->radii[i]));
                lig_solv_distf = (lig_solv_distf * exp(-deff/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
                rec_solv_distf = ((4.0/3.0)*PI*(Lig->radii[j]*Lig->radii[j]*Lig->radii[j]));
                rec_solv_distf = (rec_solv_distf * exp(-deff/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
            }

            //! VDW energy (softcore)

            acoef = (Rec->epsilons_sqrt[i] * pow(2*Rec->radii[i], 6)) * (Lig->epsilons_sqrt[j]*pow(2*Lig->radii[j], 6));
            bcoef = (sqrt2*Rec->epsilons_sqrt[i]*pow(2*Rec->radii[i], 3)) * (sqrt2*Lig->epsilons_sqrt[j]*pow(2*Lig->radii[j], 3));
            vdw+= ((acoef/pow((dij6+Input->deltaij6),2)) - (bcoef/(dij6 + Input->deltaij6)));

            //! Solvation energy

            rec_solv+= rec_solv_affinity*rec_solv_distf;
            lig_solv+= lig_solv_affinity*lig_solv_distf;
        }
    }
    //    printf("Elec: %7.3f    VDW: %7.3f    RecSolv: %7.3f    LigSolv: %7.3f    Total: %7.3f\n", elec, vdw, rec_solv, lig_solv, elec+vdw+rec_solv+lig_solv);
    return((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*vdw)+rec_solv+lig_solv);
}

double Energy2::compute_energy_softcore(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz){
    double vdw=0.0;
    double elec=0.0;
    double dij, dij2, deff, acoef, bcoef, dij6;
    double sqrt2 = sqrt(2);
    for (int i=0; i < Rec->N; i++){
        for (int j=0; j< Lig->N; j++){

            //!distances
            dij2 = distance_squared(Rec->xyz[i][0], lig_xyz[j][0], Rec->xyz[i][1], lig_xyz[j][1], Rec->xyz[i][2], lig_xyz[j][2]);
            dij = sqrt(dij2);
            dij6 = dij2*dij2*dij2;
            deff = pow(((dij*dij*dij)+Input->deltaij_es3), (1.0/3.0));

            //! electrostactic energy (softcore)
            if (Input->dielectric_model == "constant"){
                elec+= (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*deff);
            }
            else {
                elec+= (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*(pow((dij6+Input->deltaij_es6),(1.0/3.0))));
            }

            //! VDW energy (softcore)

            acoef = (Rec->epsilons_sqrt[i] * pow(2*Rec->radii[i], 6)) * (Lig->epsilons_sqrt[j]*pow(2*Lig->radii[j], 6));
            bcoef = (sqrt2*Rec->epsilons_sqrt[i]*pow(2*Rec->radii[i], 3)) * (sqrt2*Lig->epsilons_sqrt[j]*pow(2*Lig->radii[j], 3));
            vdw+= ((acoef/pow((dij6+Input->deltaij6),2)) - (bcoef/(dij6 + Input->deltaij6)));
        }
    }
    //    printf("Elec: %7.3f    VDW: %7.3f    Total: %7.3f\n", elec, vdw, elec+vdw);
    return((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*vdw));
}

double Energy2::compute_energy_hardcore_solvation(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz){
    double bcoef, acoef, lig_solv_affinity, rec_solv_affinity, lig_solv_distf, rec_solv_distf, rij, eij, dij2, dij, deff;
    double vdw=0.0;
    double elec=0.00;
    double rec_solv=0.0;
    double lig_solv=0.0;
    for (int i=0; i<Rec->N; i++){
        rec_solv_affinity = (Input->solvation_alpha*((Rec->charges[i])*(Rec->charges[i])))+Input->solvation_beta;
        for (int j=0; j < Lig->N; j++){
            rij = Rec->radii[i] + Lig->radii[j];
            eij = sqrt(Rec->epsilons[i]*Lig->epsilons[j]);
            dij2 = this->distance_squared(Rec->xyz[i][0], lig_xyz[j][0], Rec->xyz[i][1], lig_xyz[j][1], Rec->xyz[i][2], lig_xyz[j][2]);
            dij = sqrt(dij2);

            if (Input->dielectric_model == "constant"){
                elec += (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*dij);
            }
            else if (Input->dielectric_model == "4r") {
                elec += (332.0*Rec->charges[i]*Lig->charges[j])/(4*dij2);
            }
            else {                                          // Input->dielectric_model = "r"
                elec += (332.0*Rec->charges[i]*Lig->charges[j])/(dij2);
            }

            bcoef = 2*eij*(rij*rij*rij*rij*rij*rij);
            acoef = eij*(rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij);
            vdw+= ((acoef/(dij2*dij2*dij2*dij2*dij2*dij2)) - (bcoef/(dij2*dij2*dij2)));

            deff = dij2;

            lig_solv_affinity = (Input->solvation_alpha*((Lig->charges[j])*(Lig->charges[j])))+Input->solvation_beta;
            lig_solv_distf = ((4.0/3.0)*PI*(Rec->radii[i]*Rec->radii[i]*Rec->radii[i]));
            lig_solv_distf = (lig_solv_distf * exp(-(deff)/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
            rec_solv_distf = ((4.0/3.0)*PI*(Lig->radii[j]*Lig->radii[j]*Lig->radii[j]));
            rec_solv_distf = (rec_solv_distf * exp(-(deff)/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
            rec_solv+= rec_solv_affinity*rec_solv_distf;
            lig_solv+= lig_solv_affinity*lig_solv_distf;
        }
    }
    //    printf("Elec: %7.3f    VDW: %7.3f    RecSolv: %7.3f    LigSolv: %7.3f    Total: %7.3f\n", elec, vdw, rec_solv, lig_solv, elec+vdw+rec_solv+lig_solv);
    return ((Input->scale_vdw_energy*vdw)+(Input->scale_elec_energy*elec)+rec_solv+lig_solv);
}

double Energy2::compute_energy_hardcore(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz){
    double vdw=0.00;
    double elec=0.00;
    double rij, eij, dij, dij2, acoef, bcoef;
    for (int i=0; i< Rec->N; i++){
        for (int j=0; j< Lig->N; j++){
            rij = Rec->radii[i]+Lig->radii[j];
            eij = sqrt(Rec->epsilons[i]*Lig->epsilons[j]);
            dij2 = this->distance_squared(Rec->xyz[i][0], lig_xyz[j][0], Rec->xyz[i][1], lig_xyz[j][1], Rec->xyz[i][2], lig_xyz[j][2]); //rec_crd[i][0], lig_crd[j][0],rec_crd[i][1], lig_crd[j][1], rec_crd[i][2], lig_crd[j][2]);
            dij = sqrt(dij2);
            if (Input->dielectric_model == "constant"){
                elec += (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*dij);
            }
            else if (Input->dielectric_model == "4r") {
                elec += (332.0*Rec->charges[i]*Lig->charges[j])/(4*dij2);
            }
            else {                                          // Input->dielectric_model = "r"
                elec += (332.0*Rec->charges[i]*Lig->charges[j])/(dij2);
            }
            bcoef = 2*eij*(rij*rij*rij*rij*rij*rij);
            acoef = eij*(rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij);
            vdw+= ((acoef/(dij2*dij2*dij2*dij2*dij2*dij2)) - (bcoef/(dij2*dij2*dij2)));
        }
    }
    //    printf("Elec: %7.3f    VDW: %7.3f    Total: %7.3f\n", elec, vdw,elec+vdw);
    return ((Input->scale_vdw_energy*vdw)+(Input->scale_elec_energy*elec));
}

double Energy2::compute_ene_from_grids_softcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz){
    double elec =0.0, vdwA=0.0, vdwB = 0.00, rec_solv=0.00, lig_solv=0.00, lig_affinity, hb_donor=0.0, hb_acceptor=0.0;
    int a1=0, b1=0, c1=0, a2=0, b2=0, c2=0, donor_index;
    double af, bf, cf, angle_term;
    for (int i=0; i< Lig->N; i++){
        af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
        bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
        cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
        a1 = floor(af);
        b1 = floor(bf);
        c1 = floor(cf);
        a2 = ceil(af);
        if (a2 <= a1){
            a2++;
        }
        b2 = ceil(bf);
        if (b2 <= b1){
            b2++;
        }
        c2 = ceil(cf);
        if (c2 <= c1){
            c2++;
        }

        if (a1 > 0 and b1 > 0 and c1 > 0 and a2 < Grids->npointsx and b2 < Grids->npointsy and c2 < Grids->npointsz){

            GridInterpol* GI = new GridInterpol;
            this->trilinear_interpolation(Grids, af, bf, cf, a1, b1, c1, a2, b2, c2, GI);

            if (Input->use_pbsa and Grids->pbsa_loaded){
                elec += Lig->charges[i]* GI->pbsa;
            }
            else if (Input->use_delphi and Grids->delphi_loaded){
                elec += Lig->charges[i] * GI->delphi;
            }
            else {
                elec += Lig->charges[i]* GI->elec;
            }
            vdwA += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 6) * GI->vdwA;
            vdwB += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 3) * GI->vdwB;

            lig_affinity = (Input->solvation_alpha*Lig->charges[i]*Lig->charges[i]) + Input->solvation_beta;
            rec_solv += GI->solv_gauss * (4.0/3.0) * PI * pow(Lig->radii[i], 3);
            lig_solv += lig_affinity * GI->rec_solv_gauss;

            if (this->atom_is_acceptor(i, Lig)){
                hb_donor += GI->hb_donor;
            }
            else {
                if (Lig->HBacceptors.size() > 0){
                    donor_index = this->H_is_donor(i, Lig);
                    if (donor_index >= 0){       // atom is a H donor
                        double x = (af*Grids->grid_spacing)+Grids->xbegin;
                        double y = (bf*Grids->grid_spacing)+Grids->ybegin;
                        double z = (cf*Grids->grid_spacing)+Grids->zbegin;
                        angle_term = this->angle(Lig->xyz[donor_index][0], Lig->xyz[donor_index][1], Lig->xyz[donor_index][2], Lig->xyz[i][0], Lig->xyz[i][1], Lig->xyz[i][2], x, y, z);
                        angle_term = cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0);
                        hb_acceptor += GI->hb_acceptor*angle_term;
                    }
                }
            }
            delete GI;
        }
        else {
            elec += 999999.9;
            vdwA += 999999.9;
            vdwB += 999999.9;
#ifdef DEBUG
            //            printf("Ligand %s is slipping from computation box. Inaccurate energies will be estimated.\n", Lig->molname.c_str());
#endif
        }
    }
    return((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*(vdwA-vdwB))+rec_solv+lig_solv+hb_donor+hb_acceptor);
}

double Energy2::compute_ene_from_grids_softcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz){
    double elec =0.0, vdwA=0.0, vdwB = 0.00, hb_donor=0.0, hb_acceptor=0.0;
    int a1=0, b1=0, c1=0, a2=0, b2=0, c2=0, donor_index;
    double af, bf, cf, angle_term;
    for (int i=0; i< Lig->N; i++){
        af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
        bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
        cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
        a1 = floor(af);
        b1 = floor(bf);
        c1 = floor(cf);
        a2 = ceil(af);
        if (a2 <= a1){
            a2++;
        }
        b2 = ceil(bf);
        if (b2 <= b1){
            b2++;
        }
        c2 = ceil(cf);
        if (c2 <= c1){
            c2++;
        }

        if (a1 > 0 and b1 > 0 and c1 > 0 and a2 < Grids->npointsx and b2 < Grids->npointsy and c2 < Grids->npointsz){

            GridInterpol* GI = new GridInterpol;
            this->trilinear_interpolation(Grids, af, bf, cf, a1, b1, c1, a2, b2, c2, GI);

            if (Input->use_pbsa and Grids->pbsa_loaded){
                elec += Lig->charges[i]* GI->pbsa;
            }
            else if (Input->use_delphi and Grids->delphi_loaded){
                elec += Lig->charges[i] * GI->delphi;
            }
            else {
                elec += Lig->charges[i]* GI->elec;
            }
            vdwA += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 6) * GI->vdwA;
            vdwB += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 3) * GI->vdwB;

            if (this->atom_is_acceptor(i, Lig)){
                hb_donor += GI->hb_donor;
            }
            else {
                if (Lig->HBacceptors.size() > 0){
                    donor_index = this->H_is_donor(i, Lig);
                    if (donor_index >= 0){       // atom is a H donor
                        double x = (af*Grids->grid_spacing)+Grids->xbegin;
                        double y = (bf*Grids->grid_spacing)+Grids->ybegin;
                        double z = (cf*Grids->grid_spacing)+Grids->zbegin;
                        angle_term = this->angle(Lig->xyz[donor_index][0], Lig->xyz[donor_index][1], Lig->xyz[donor_index][2], Lig->xyz[i][0], Lig->xyz[i][1], Lig->xyz[i][2], x, y, z);
                        angle_term = cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0);
                        hb_acceptor += GI->hb_acceptor*angle_term;
                    }
                }
            }

            delete GI;
        }
        else {
            elec += 999999.9;
            vdwA += 999999.9;
            vdwB += 999999.9;
#ifdef DEBUG
            //            printf("Ligand %s is slipping from computation box. Inaccurate energies will be estimated.\n", Lig->molname.c_str());
#endif
        }
    }
    return((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*(vdwA-vdwB))+ hb_donor + hb_acceptor);
}

double Energy2::compute_ene_from_grids_hardcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz){
    double elec =0.0, vdwA=0.0, vdwB = 0.00, rec_solv=0.00, lig_solv=0.00, hb_donor=0.0, hb_acceptor=0.0, lig_affinity, angle_term;
    int a1=0, b1=0, c1=0, a2=0, b2=0, c2=0, donor_index;
    double af, bf, cf;
    double sqrt2=sqrt(2.0);
    for (int i=0; i< Lig->N; i++){
        af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
        bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
        cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
        a1 = floor(af);
        b1 = floor(bf);
        c1 = floor(cf);
        a2 = ceil(af);
        if (a2 <= a1){
            a2++;
        }
        b2 = ceil(bf);
        if (b2 <= b1){
            b2++;
        }
        c2 = ceil(cf);
        if (c2 <= c1){
            c2++;
        }


        if (a1 > 0 and b1 > 0 and c1 > 0 and a2 < Grids->npointsx and b2 < Grids->npointsy and c2 < Grids->npointsz){

            GridInterpol* GI = new GridInterpol;
            this->trilinear_interpolation(Grids, af, bf, cf, a1, b1, c1, a2, b2, c2, GI);

            if (Input->use_pbsa and Grids->pbsa_loaded){
                elec += Lig->charges[i]* GI->pbsa;
            }
            else if (Input->use_delphi and Grids->delphi_loaded){
                elec += Lig->charges[i] * GI->delphi;
            }
            else {
                elec += Lig->charges[i]* GI->elec;
            }

            vdwA += Lig->epsilons_sqrt[i]*64.0*pow(Lig->radii[i],6) * GI->vdwA;
            vdwB += sqrt2*Lig->epsilons_sqrt[i]*8.0*pow(Lig->radii[i], 3) * GI->vdwB;

            lig_affinity = (Input->solvation_alpha*Lig->charges[i]*Lig->charges[i]) + Input->solvation_beta;
            rec_solv +=  GI->solv_gauss * (4.0/3.0) * PI * pow(Lig->radii[i], 3);
            lig_solv += lig_affinity*GI->rec_solv_gauss;

            if (this->atom_is_acceptor(i, Lig)){
                hb_donor += GI->hb_donor;
            }
            else {
                if (Lig->HBacceptors.size() > 0){
                    donor_index = this->H_is_donor(i, Lig);
                    if (donor_index >= 0){       // atom is a H donor
                        double x = (af*Grids->grid_spacing)+Grids->xbegin;
                        double y = (bf*Grids->grid_spacing)+Grids->ybegin;
                        double z = (cf*Grids->grid_spacing)+Grids->zbegin;
                        angle_term = this->angle(Lig->xyz[donor_index][0], Lig->xyz[donor_index][1], Lig->xyz[donor_index][2], Lig->xyz[i][0], Lig->xyz[i][1], Lig->xyz[i][2], x, y, z);
                        angle_term = cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0);
                        hb_acceptor += GI->hb_acceptor*angle_term;
                    }
                }
            }
            delete GI;
        }
        else {
            elec += 999999.9;
            vdwA += 999999.9;
            vdwB += 999999.9;
#ifdef DEBUG
            printf("Ligand %s is slipping from computation box. Inaccurate energies will be estimated.\n", Lig->molname.c_str());
#endif
        }
    }
    return((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*(vdwA-vdwB))+rec_solv+lig_solv + (hb_donor+hb_acceptor));
}

double Energy2::compute_ene_from_grids_hardcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz){
    double elec =0.0, vdwA=0.0, vdwB = 0.00, hb_acceptor=0.0, hb_donor=0.0;
    int a1=0, b1=0, c1=0, a2=0, b2=0, c2=0, donor_index;
    double af, bf, cf, angle_term;
    double sqrt2=sqrt(2.0);
    for (int i=0; i< Lig->N; i++){
        af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
        bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
        cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
        a1 = floor(af);
        b1 = floor(bf);
        c1 = floor(cf);
        a2 = ceil(af);
        if (a2 <= a1){
            a2++;
        }
        b2 = ceil(bf);
        if (b2 <= b1){
            b2++;
        }
        c2 = ceil(cf);
        if (c2 <= c1){
            c2++;
        }

        if (a1 > 0 and b1 > 0 and c1 > 0 and a2 < Grids->npointsx and b2 < Grids->npointsy and c2 < Grids->npointsz){

            GridInterpol* GI = new GridInterpol;
            this->trilinear_interpolation(Grids, af, bf, cf, a1, b1, c1, a2, b2, c2, GI);

            if (Input->use_pbsa and Grids->pbsa_loaded){
                elec += Lig->charges[i]* GI->pbsa;
            }
            else if (Input->use_delphi and Grids->delphi_loaded){
                elec += Lig->charges[i] * GI->delphi;
            }
            else {
                elec += Lig->charges[i]* GI->elec;
            }


            vdwA += Lig->epsilons_sqrt[i]*64.0*pow(Lig->radii[i],6) * GI->vdwA;
            vdwB += sqrt2*Lig->epsilons_sqrt[i]*8.0*pow(Lig->radii[i], 3) * GI->vdwB;

            if (this->atom_is_acceptor(i, Lig)){
                hb_donor += GI->hb_donor;
            }
            else {
                if (Lig->HBacceptors.size() > 0){
                    donor_index = this->H_is_donor(i, Lig);
                    if (donor_index >= 0){       // atom is a H donor
                        double x = (af*Grids->grid_spacing)+Grids->xbegin;
                        double y = (bf*Grids->grid_spacing)+Grids->ybegin;
                        double z = (cf*Grids->grid_spacing)+Grids->zbegin;
                        angle_term = this->angle(Lig->xyz[donor_index][0], Lig->xyz[donor_index][1], Lig->xyz[donor_index][2], Lig->xyz[i][0], Lig->xyz[i][1], Lig->xyz[i][2], x, y, z);
                        angle_term = cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0);
                        hb_acceptor += GI->hb_acceptor*angle_term;
                    }
                }
            }

            delete GI;
        }
        else {
            elec += 999999.9;
            vdwA += 999999.9;
            vdwB += 999999.9;
#ifdef DEBUG
            //            printf("Ligand %s is slipping from computation box. Inaccurate energies will be estimated.\n", Lig->molname.c_str());
#endif
        }
    }
    return((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*(vdwA-vdwB)) + (hb_acceptor+hb_donor));
}

double Energy2::compute_ene_from_grids_hardcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result){
    double elec =0.0, vdwA=0.0, vdwB = 0.00, rec_solv=0.00, lig_solv=0.00, lig_affinity, hb_donor=0.0, hb_acceptor=0.0;
    int a1=0, b1=0, c1=0, a2=0, b2=0, c2=0, donor_index;
    double af, bf, cf, angle_term;
    double sqrt2=sqrt(2);
    for (int i=0; i< Lig->N; i++){
        af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
        bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
        cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
        a1 = floor(af);
        b1 = floor(bf);
        c1 = floor(cf);
        a2 = ceil(af);
        if (a2 <= a1){
            a2++;
        }
        b2 = ceil(bf);
        if (b2 <= b1){
            b2++;
        }
        c2 = ceil(cf);
        if (c2 <= c1){
            c2++;
        }

        if (a1 > 0 and b1 > 0 and c1 > 0 and a2 < Grids->npointsx and b2 < Grids->npointsy and c2 < Grids->npointsz){

            GridInterpol* GI = new GridInterpol;
            this->trilinear_interpolation(Grids, af, bf, cf, a1, b1, c1, a2, b2, c2, GI);

            if (Input->use_pbsa and Grids->pbsa_loaded){
                elec += Lig->charges[i]* GI->pbsa;
            }
            else if (Input->use_delphi and Grids->delphi_loaded){
                elec += Lig->charges[i] * GI->delphi;
            }
            else {
                elec += Lig->charges[i]* GI->elec;
            }

            vdwA += Lig->epsilons_sqrt[i]*64.0*pow(Lig->radii[i],6) * GI->vdwA;
            vdwB += sqrt2*Lig->epsilons_sqrt[i]*8.0*pow(Lig->radii[i], 3) * GI->vdwB;

            lig_affinity = (Input->solvation_alpha*Lig->charges[i]*Lig->charges[i]) + Input->solvation_beta;
            rec_solv += GI->solv_gauss * (4.0/3.0) * PI * pow(Lig->radii[i], 3);
            lig_solv += lig_affinity * GI->rec_solv_gauss;

            if (this->atom_is_acceptor(i, Lig)){
                hb_donor += GI->hb_donor;
            }
            else {
                if (Lig->HBacceptors.size() > 0){
                    donor_index = this->H_is_donor(i, Lig);
                    if (donor_index >= 0){       // atom is a H donor
                        double x = (af*Grids->grid_spacing)+Grids->xbegin;
                        double y = (bf*Grids->grid_spacing)+Grids->ybegin;
                        double z = (cf*Grids->grid_spacing)+Grids->zbegin;
                        angle_term = this->angle(Lig->xyz[donor_index][0], Lig->xyz[donor_index][1], Lig->xyz[donor_index][2], Lig->xyz[i][0], Lig->xyz[i][1], Lig->xyz[i][2], x, y, z);
                        angle_term = cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0);
                        hb_acceptor += GI->hb_acceptor*angle_term;
                    }
                }
            }
            delete GI;
        }
        else {
            elec += 999999.9;
            vdwA += 999999.9;
            vdwB += 999999.9;
        }
    }
    energy_result->vdw = Input->scale_vdw_energy*(vdwA-vdwB);
    energy_result->elec = Input->scale_elec_energy*elec;
    energy_result->rec_solv = rec_solv;
    energy_result->lig_solv = lig_solv;
    energy_result->hb_donor = hb_donor;
    energy_result->hb_acceptor = hb_acceptor;
    energy_result->total = (Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*(vdwA-vdwB))+rec_solv+lig_solv+hb_donor+hb_acceptor;

    return((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*(vdwA-vdwB))+rec_solv+lig_solv+hb_donor+hb_acceptor);
}

double Energy2::compute_ene_from_grids_hardcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result){
    double elec =0.0, vdwA=0.0, vdwB = 0.00, hb_donor=0.0, hb_acceptor=0.0;
    int a1=0, b1=0, c1=0, a2=0, b2=0, c2=0, donor_index;
    double af, bf, cf, angle_term;
    double sqrt2= sqrt(2);
    for (int i=0; i< Lig->N; i++){
        af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
        bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
        cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
        a1 = floor(af);
        b1 = floor(bf);
        c1 = floor(cf);
        a2 = ceil(af);
        if (a2 <= a1){
            a2++;
        }
        b2 = ceil(bf);
        if (b2 <= b1){
            b2++;
        }
        c2 = ceil(cf);
        if (c2 <= c1){
            c2++;
        }

        if (a1 > 0 and b1 > 0 and c1 > 0 and a2 < Grids->npointsx and b2 < Grids->npointsy and c2 < Grids->npointsz){
            GridInterpol* GI = new GridInterpol;
            this->trilinear_interpolation(Grids, af, bf, cf, a1, b1, c1, a2, b2, c2, GI);

            if (Input->use_pbsa and Grids->pbsa_loaded){
                elec += Lig->charges[i]* GI->pbsa;
            }
            else if (Input->use_delphi and Grids->delphi_loaded){
                elec += Lig->charges[i] * GI->delphi;
            }
            else {
                elec += Lig->charges[i]* GI->elec;
            }

            vdwA += Lig->epsilons_sqrt[i]*64.0*pow(Lig->radii[i],6);
            vdwB += sqrt2*Lig->epsilons_sqrt[i]*8.0*pow(Lig->radii[i], 3);
//            vdwA += sqrt(Lig->epsilons[i]*pow((2*Lig->radii[i]), 12))* GI->vdwA;
//            vdwB += sqrt(2.0 * Lig->epsilons[i]*pow((2*Lig->radii[i]), 6)) * GI->vdwB;

            if (this->atom_is_acceptor(i, Lig)){
                hb_donor += GI->hb_donor;
            }
            else {
                if (Lig->HBacceptors.size() > 0){
                    donor_index = this->H_is_donor(i, Lig);
                    if (donor_index >= 0){       // atom is a H donor
                        double x = (af*Grids->grid_spacing)+Grids->xbegin;
                        double y = (bf*Grids->grid_spacing)+Grids->ybegin;
                        double z = (cf*Grids->grid_spacing)+Grids->zbegin;
                        angle_term = this->angle(Lig->xyz[donor_index][0], Lig->xyz[donor_index][1], Lig->xyz[donor_index][2], Lig->xyz[i][0], Lig->xyz[i][1], Lig->xyz[i][2], x, y, z);
                        angle_term = cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0);
                        hb_acceptor += GI->hb_acceptor*angle_term;
                    }
                }
            }
            delete GI;
        }
        else {
            elec += 999999.9;
            vdwA += 999999.9;
            vdwB += 999999.9;
        }
    }
    energy_result->vdw = Input->scale_vdw_energy*(vdwA-vdwB);
    energy_result->elec = Input->scale_elec_energy*elec;
    energy_result->rec_solv = 0.00;
    energy_result->lig_solv = 0.00;
    energy_result->hb_donor = hb_donor;
    energy_result->hb_acceptor = hb_acceptor;
    energy_result->total = (Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*(vdwA-vdwB)) + hb_donor + hb_acceptor;

    return((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*(vdwA-vdwB))+(hb_donor+hb_acceptor));
}

double Energy2::compute_ene_from_grids_softcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result){
    double elec =0.0, vdwA=0.0, vdwB = 0.00, rec_solv=0.00, lig_solv=0.00, lig_affinity, hb_donor=0.0, hb_acceptor=0.0;
    int a1=0, b1=0, c1=0, a2=0, b2=0, c2=0, donor_index;
    double af, bf, cf, angle_term;
    for (int i=0; i< Lig->N; i++){
        af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
        bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
        cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
        a1 = floor(af);
        b1 = floor(bf);
        c1 = floor(cf);
        a2 = ceil(af);
        if (a2 <= a1){
            a2++;
        }
        b2 = ceil(bf);
        if (b2 <= b1){
            b2++;
        }
        c2 = ceil(cf);
        if (c2 <= c1){
            c2++;
        }

        if (a1 > 0 and b1 > 0 and c1 > 0 and a2 < Grids->npointsx and b2 < Grids->npointsy and c2 < Grids->npointsz){
            GridInterpol* GI = new GridInterpol;
            this->trilinear_interpolation(Grids, af, bf, cf, a1, b1, c1, a2, b2, c2, GI);

            if (Input->use_pbsa and Grids->pbsa_loaded){
                elec += Lig->charges[i]* GI->pbsa;
            }
            else if (Input->use_delphi and Grids->delphi_loaded){
                elec += Lig->charges[i] * GI->delphi;
            }
            else {
                elec += Lig->charges[i]* GI->elec;
            }
            vdwA += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 6) * GI->vdwA;
            vdwB += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 3) * GI->vdwB;

            lig_affinity = (Input->solvation_alpha*Lig->charges[i]*Lig->charges[i]) + Input->solvation_beta;
            rec_solv += GI->solv_gauss * (4.0/3.0) * PI * pow(Lig->radii[i], 3);
            lig_solv += lig_affinity * GI->rec_solv_gauss;

            if (this->atom_is_acceptor(i, Lig)){
                hb_donor += GI->hb_donor;
            }
            else {
                if (Lig->HBdonors.size() > 0){
                    donor_index = this->H_is_donor(i, Lig);
                    if (donor_index >= 0){       // atom is a H donor
                        double x = (af*Grids->grid_spacing)+Grids->xbegin;
                        double y = (bf*Grids->grid_spacing)+Grids->ybegin;
                        double z = (cf*Grids->grid_spacing)+Grids->zbegin;
                        angle_term = this->angle(Lig->xyz[donor_index][0], Lig->xyz[donor_index][1], Lig->xyz[donor_index][2], Lig->xyz[i][0], Lig->xyz[i][1], Lig->xyz[i][2], x, y, z);
                        angle_term = cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0);
                        hb_acceptor += GI->hb_acceptor*angle_term;
                    }
                }
            }
            delete GI;
        }
        else {
            elec += 999999.9;
            vdwA += 999999.9;
            vdwB += 999999.9;
        }
    }
    energy_result->elec = Input->scale_elec_energy*elec;
    energy_result->vdw = Input->scale_vdw_energy*(vdwA-vdwB);
    energy_result->rec_solv = rec_solv;
    energy_result->lig_solv = lig_solv;
    energy_result->hb_donor = hb_donor;
    energy_result->hb_acceptor = hb_acceptor;
    energy_result->total = (Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*(vdwA-vdwB))+rec_solv+lig_solv+hb_donor+hb_acceptor;

    return ((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*(vdwA-vdwB))+rec_solv+lig_solv+hb_donor+hb_acceptor);
}

double Energy2::compute_ene_from_grids_softcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result){
    double elec =0.0, vdwA=0.0, vdwB = 0.00, hb_donor=0.0, hb_acceptor=0.0;
    int a1=0, b1=0, c1=0, a2=0, b2=0, c2=0, donor_index;
    double af, bf, cf, angle_term;
    for (int i=0; i< Lig->N; i++){
        af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
        bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
        cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
        a1 = floor(af);
        b1 = floor(bf);
        c1 = floor(cf);
        a2 = ceil(af);
        if (a2 <= a1){
            a2++;
        }
        b2 = ceil(bf);
        if (b2 <= b1){
            b2++;
        }
        c2 = ceil(cf);
        if (c2 <= c1){
            c2++;
        }

        if (a1 > 0 and b1 > 0 and c1 > 0 and a2 < Grids->npointsx and b2 < Grids->npointsy and c2 < Grids->npointsz){
            GridInterpol* GI = new GridInterpol;
            this->trilinear_interpolation(Grids, af, bf, cf, a1, b1, c1, a2, b2, c2, GI);

            if (Input->use_pbsa and Grids->pbsa_loaded){
                elec += Lig->charges[i]* GI->pbsa;
            }
            else if (Input->use_delphi and Grids->delphi_loaded){
                elec += Lig->charges[i] * GI->delphi;
            }
            else {
                elec += Lig->charges[i]* GI->elec;
            }
            vdwA += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 6) * GI->vdwA;
            vdwB += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 3) * GI->vdwB;

            if (this->atom_is_acceptor(i, Lig)){
                hb_donor += GI->hb_donor;
            }
            else {
                if (Lig->HBacceptors.size() > 0){
                    donor_index = this->H_is_donor(i, Lig);
                    if (donor_index >= 0){       // atom is a H donor
                        double x = (af*Grids->grid_spacing)+Grids->xbegin;
                        double y = (bf*Grids->grid_spacing)+Grids->ybegin;
                        double z = (cf*Grids->grid_spacing)+Grids->zbegin;
                        angle_term = this->angle(Lig->xyz[donor_index][0], Lig->xyz[donor_index][1], Lig->xyz[donor_index][2], Lig->xyz[i][0], Lig->xyz[i][1], Lig->xyz[i][2], x, y, z);
                        angle_term = cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0) * cos(angle_term * PI / 180.0);
                        hb_acceptor += GI->hb_acceptor*angle_term;
                    }
                }
            }

            delete GI;
        }
        else {
            elec += 999999.9;
            vdwA += 999999.9;
            vdwB += 999999.9;
        }
    }
    energy_result->elec = Input->scale_elec_energy*elec;
    energy_result->vdw = Input->scale_vdw_energy*(vdwA-vdwB);
    energy_result->rec_solv = 0.00;
    energy_result->lig_solv = 0.00;
    energy_result->hb_donor = hb_donor;
    energy_result->hb_acceptor = hb_acceptor;
    energy_result->total = (Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*(vdwA-vdwB)) + hb_donor + hb_acceptor;

    return ((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*(vdwA-vdwB)) + hb_donor + hb_acceptor);
}

double Energy2::compute_energy_softcore_solvation(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result){
    double rec_solv = 0.00;
    double lig_solv = 0.00;
    double vdw=0.0;
    double elec=0.0;
    double acoef, bcoef, lig_solv_affinity, lig_solv_distf, rec_solv_distf, rec_solv_affinity, dij, dij2, deff, dij6;
    double sqrt2 = sqrt(2);

    for (int i=0; i< Rec->N; i++){
        rec_solv_affinity = (Input->solvation_alpha*((Rec->charges[i])*(Rec->charges[i])))+Input->solvation_beta;
        for (int j=0; j< Lig->N; j++){

            //! distances

            dij2 = this->distance_squared(Rec->xyz[i][0], lig_xyz[j][0],Rec->xyz[i][1], lig_xyz[j][1], Rec->xyz[i][2], lig_xyz[j][2]);
            dij6 = dij2*dij2*dij2;
            dij = sqrt(dij2);
            deff = pow(((dij*dij*dij)+Input->deltaij_es3), (1.0/3.0));

            //! electrostactic energy (softcore)
            if (Input->dielectric_model == "constant"){
                elec+= (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*deff);

                lig_solv_affinity = (Input->solvation_alpha*((Lig->charges[j])*(Lig->charges[j])))+Input->solvation_beta;
                lig_solv_distf = ((4.0/3.0)*PI*(Rec->radii[i]*Rec->radii[i]*Rec->radii[i]));
                lig_solv_distf = (lig_solv_distf * exp(-(deff*deff)/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
                rec_solv_distf = ((4.0/3.0)*PI*(Lig->radii[j]*Lig->radii[j]*Lig->radii[j]));
                rec_solv_distf = (rec_solv_distf * exp(-(deff*deff)/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
            }
            else {
                deff = pow((dij6+Input->deltaij_es6), (1.0/3.0));
                elec+= (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*(pow((dij6+Input->deltaij_es6),(1.0/3.0))));

                lig_solv_affinity = (Input->solvation_alpha*((Lig->charges[j])*(Lig->charges[j])))+Input->solvation_beta;
                lig_solv_distf = ((4.0/3.0)*PI*(Rec->radii[i]*Rec->radii[i]*Rec->radii[i]));
                lig_solv_distf = (lig_solv_distf * exp(-deff/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
                rec_solv_distf = ((4.0/3.0)*PI*(Lig->radii[j]*Lig->radii[j]*Lig->radii[j]));
                rec_solv_distf = (rec_solv_distf * exp(-deff/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);

            }

            //! VDW energy (softcore)

            acoef = (Rec->epsilons_sqrt[i] * pow(2*Rec->radii[i], 6)) * (Lig->epsilons_sqrt[j]*pow(2*Lig->radii[j], 6));
            bcoef = (sqrt2*Rec->epsilons_sqrt[i]*pow(2*Rec->radii[i], 3)) * (sqrt2*Lig->epsilons_sqrt[j]*pow(2*Lig->radii[j], 3));
            vdw+= ((acoef/pow((dij6+Input->deltaij6),2)) - (bcoef/(dij6 + Input->deltaij6)));

            //! Solvation energy

            rec_solv+= rec_solv_affinity*rec_solv_distf;
            lig_solv+= lig_solv_affinity*lig_solv_distf;
        }
    }
    energy_result->elec = Input->scale_elec_energy*elec;
    energy_result->vdw = Input->scale_vdw_energy*vdw;
    energy_result->rec_solv = rec_solv;
    energy_result->lig_solv = lig_solv;
    energy_result->total = ((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*vdw)+rec_solv+lig_solv);

    return((Input->elec_scale*elec)+(Input->scale_vdw_energy*vdw)+rec_solv+lig_solv);
}

double Energy2::compute_energy_softcore(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result){
    double vdw=0.0;
    double elec=0.0;
    double acoef, bcoef,dij, dij2, deff, dij6;
    double sqrt2 = sqrt(2);

    for (int i=0; i< Rec->N; i++){
        for (int j=0; j< Lig->N; j++){

            //! distances

            dij2 = this->distance_squared(Rec->xyz[i][0], lig_xyz[j][0],Rec->xyz[i][1], lig_xyz[j][1], Rec->xyz[i][2], lig_xyz[j][2]);
            dij6 = dij2*dij2*dij2;
            dij = sqrt(dij2);
            deff = pow(((dij*dij*dij)+Input->deltaij_es3), (1.0/3.0));

            //! electrostactic energy (softcore)
            if (Input->dielectric_model == "constant"){
                elec+= (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*deff);
            }
            else{
                deff = pow((dij6+Input->deltaij_es6), (1.0/3.0));
                elec+= (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*(pow((dij6+Input->deltaij_es6),(1.0/3.0))));
            }

            //! VDW energy (softcore)

            acoef = (Rec->epsilons_sqrt[i] * pow(2*Rec->radii[i], 6)) * (Lig->epsilons_sqrt[j]*pow(2*Lig->radii[j], 6));
            bcoef = (sqrt2*Rec->epsilons_sqrt[i]*pow(2*Rec->radii[i], 3)) * (sqrt2*Lig->epsilons_sqrt[j]*pow(2*Lig->radii[j], 3));
            vdw+= ((acoef/pow((dij6+Input->deltaij6),2)) - (bcoef/(dij6 + Input->deltaij6)));
        }
    }
    energy_result->elec = Input->scale_elec_energy*elec;
    energy_result->vdw = Input->scale_vdw_energy*vdw;
    energy_result->rec_solv = 0.00;
    energy_result->lig_solv = 0.00;
    energy_result->total = (Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*vdw);

    return((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*vdw));
}

double Energy2::compute_energy_hardcore_solvation(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result){
    double bcoef, acoef, lig_solv_affinity, rec_solv_affinity, lig_solv_distf, rec_solv_distf, rij, eij, dij2, dij, deff;
    double vdw=0.0;
    double elec=0.00;
    double rec_solv=0.0;
    double lig_solv=0.0;
    for (int i=0; i<Rec->N; i++){
        rec_solv_affinity = (Input->solvation_alpha*((Rec->charges[i])*(Rec->charges[i])))+Input->solvation_beta;
        for (int j=0; j < Lig->N; j++){
            rij = Rec->radii[i] + Lig->radii[j];
            eij = sqrt(Rec->epsilons[i]*Lig->epsilons[j]);
            dij2 = this->distance_squared(Rec->xyz[i][0], lig_xyz[j][0], Rec->xyz[i][1], lig_xyz[j][1], Rec->xyz[i][2], lig_xyz[j][2]);
            dij = sqrt(dij2);

            if (Input->dielectric_model == "constant"){
                elec += (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*dij);
            }
            else if (Input->dielectric_model == "4r") {
                elec += (332.0*Rec->charges[i]*Lig->charges[j])/(4*dij2);
            }
            else {                                          // Input->dielectric_model = "r"
                elec += (332.0*Rec->charges[i]*Lig->charges[j])/(dij2);
            }

            bcoef = 2*eij*(rij*rij*rij*rij*rij*rij);
            acoef = eij*(rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij);
            vdw+= ((acoef/(dij2*dij2*dij2*dij2*dij2*dij2)) - (bcoef/(dij2*dij2*dij2)));

            deff = dij2;

            lig_solv_affinity = (Input->solvation_alpha*((Lig->charges[j])*(Lig->charges[j])))+Input->solvation_beta;
            lig_solv_distf = ((4.0/3.0)*PI*(Rec->radii[i]*Rec->radii[i]*Rec->radii[i]));
            lig_solv_distf = (lig_solv_distf * exp(-(deff)/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
            rec_solv_distf = ((4.0/3.0)*PI*(Lig->radii[j]*Lig->radii[j]*Lig->radii[j]));
            rec_solv_distf = (rec_solv_distf * exp(-(deff)/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
            rec_solv+= rec_solv_affinity*rec_solv_distf;
            lig_solv+= lig_solv_affinity*lig_solv_distf;
        }
    }
    energy_result->elec = Input->scale_elec_energy*elec;
    energy_result->vdw = Input->scale_vdw_energy*vdw;
    energy_result->rec_solv = rec_solv;
    energy_result->lig_solv = lig_solv;
    energy_result->total = ((Input->scale_vdw_energy*vdw)+(Input->scale_elec_energy*elec)+rec_solv+lig_solv);

    return ((Input->scale_vdw_energy*vdw)+(Input->scale_elec_energy*elec)+rec_solv+lig_solv);
}

double Energy2::compute_energy_hardcore(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result){
    double bcoef, acoef, rij, eij, dij2, dij;
    double vdw=0.0;
    double elec=0.00;
    for (int i=0; i<Rec->N; i++){
        for (int j=0; j < Lig->N; j++){
            rij = Rec->radii[i] + Lig->radii[j];
            eij = sqrt(Rec->epsilons[i]*Lig->epsilons[j]);
            dij2 = this->distance_squared(Rec->xyz[i][0], lig_xyz[j][0], Rec->xyz[i][1], lig_xyz[j][1], Rec->xyz[i][2], lig_xyz[j][2]);
            dij = sqrt(dij2);

            if (Input->dielectric_model == "constant"){
                elec += (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*dij);
            }
            else if (Input->dielectric_model == "4r") {
                elec += (332.0*Rec->charges[i]*Lig->charges[j])/(4*dij2);
            }
            else {                                          // Input->dielectric_model = "r"
                elec += (332.0*Rec->charges[i]*Lig->charges[j])/(dij2);
            }

            bcoef = 2*eij*(rij*rij*rij*rij*rij*rij);
            acoef = eij*(rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij);
            vdw+= ((acoef/(dij2*dij2*dij2*dij2*dij2*dij2)) - (bcoef/(dij2*dij2*dij2)));
        }
    }
    energy_result->elec = Input->scale_elec_energy*elec;
    energy_result->vdw = Input->scale_vdw_energy*vdw;
    energy_result->rec_solv = 0.00;
    energy_result->lig_solv = 0.00;
    energy_result->total = (Input->scale_vdw_energy*vdw)+(Input->scale_elec_energy*elec);

    return((Input->scale_elec_energy*elec)+(Input->scale_vdw_energy*vdw));
}

double Energy2::trilinear_interpolation(vector<vector<vector<double> > > grid, double x, double y, double z, int x0, int y0, int z0, int x1, int y1, int z1){
    double xd, yd, zd;
    double c00, c10, c01, c11, c0, c1, c;
    xd=(x-x0)/(x1-x0);
    yd=(y-y0)/(y1-y0);
    zd=(z-z0)/(z1-z0);

    c00=(grid[x0][y0][z0]*(1-xd)) + (grid[x1][y0][z0]*xd);
    c10=(grid[x0][y1][z0]*(1-xd)) + (grid[x1][y1][z0]*xd);
    c01=(grid[x0][y0][z1]*(1-xd)) + (grid[x1][y0][z1]*xd);
    c11=(grid[x0][y1][z1]*(1-xd)) + (grid[x1][y1][z1]*xd);

    c0=(c00*(1-yd))+(c10*yd);
    c1=(c01*(1-yd))+(c11*yd);

    c=(c0*(1-zd))+(c1*zd);

    return (c);
}

void Energy2::trilinear_interpolation(Grid* Grids, double x, double y, double z, int x0, int y0, int z0, int x1, int y1, int z1, GridInterpol* GI){
    double xd, yd, zd;
    double c00, c10, c01, c11, c0, c1;

    xd=double((x-x0)/(x1-x0));
    yd=double((y-y0)/(y1-y0));
    zd=double((z-z0)/(z1-z0));

    // PBSA Grid

    if (Input->use_pbsa and Grids->pbsa_loaded){

        c00=(Grids->pbsa_grid[x0][y0][z0]*(1-xd)) + (Grids->pbsa_grid[x1][y0][z0]*xd);
        c10=(Grids->pbsa_grid[x0][y1][z0]*(1-xd)) + (Grids->pbsa_grid[x1][y1][z0]*xd);
        c01=(Grids->pbsa_grid[x0][y0][z1]*(1-xd)) + (Grids->pbsa_grid[x1][y0][z1]*xd);
        c11=(Grids->pbsa_grid[x0][y1][z1]*(1-xd)) + (Grids->pbsa_grid[x1][y1][z1]*xd);

        c0=(c00*(1-yd))+(c10*yd);
        c1=(c01*(1-yd))+(c11*yd);

        GI->pbsa=(c0*(1-zd))+(c1*zd);
    }

    // DelPhi Grid

    else if (Input->use_delphi and Grids->delphi_loaded){

        c00=(Grids->delphi_grid[x0][y0][z0]*(1-xd)) + (Grids->delphi_grid[x1][y0][z0]*xd);
        c10=(Grids->delphi_grid[x0][y1][z0]*(1-xd)) + (Grids->delphi_grid[x1][y1][z0]*xd);
        c01=(Grids->delphi_grid[x0][y0][z1]*(1-xd)) + (Grids->delphi_grid[x1][y0][z1]*xd);
        c11=(Grids->delphi_grid[x0][y1][z1]*(1-xd)) + (Grids->delphi_grid[x1][y1][z1]*xd);

        c0=(c00*(1-yd))+(c10*yd);
        c1=(c01*(1-yd))+(c11*yd);

        GI->delphi=(c0*(1-zd))+(c1*zd);
    }

    else {
        c00=(Grids->elec_grid[x0][y0][z0]*(1-xd)) + (Grids->elec_grid[x1][y0][z0]*xd);
        c10=(Grids->elec_grid[x0][y1][z0]*(1-xd)) + (Grids->elec_grid[x1][y1][z0]*xd);
        c01=(Grids->elec_grid[x0][y0][z1]*(1-xd)) + (Grids->elec_grid[x1][y0][z1]*xd);
        c11=(Grids->elec_grid[x0][y1][z1]*(1-xd)) + (Grids->elec_grid[x1][y1][z1]*xd);

        c0=(c00*(1-yd))+(c10*yd);
        c1=(c01*(1-yd))+(c11*yd);

        GI->elec=(c0*(1-zd))+(c1*zd);
    }

    // VDW attractive part

    c00=(Grids->vdwA_grid[x0][y0][z0]*(1-xd)) + (Grids->vdwA_grid[x1][y0][z0]*xd);
    c10=(Grids->vdwA_grid[x0][y1][z0]*(1-xd)) + (Grids->vdwA_grid[x1][y1][z0]*xd);
    c01=(Grids->vdwA_grid[x0][y0][z1]*(1-xd)) + (Grids->vdwA_grid[x1][y0][z1]*xd);
    c11=(Grids->vdwA_grid[x0][y1][z1]*(1-xd)) + (Grids->vdwA_grid[x1][y1][z1]*xd);

    c0=(c00*(1-yd))+(c10*yd);
    c1=(c01*(1-yd))+(c11*yd);

    GI->vdwA=(c0*(1-zd))+(c1*zd);

    // VDW repulsive part

    c00=(Grids->vdwB_grid[x0][y0][z0]*(1-xd)) + (Grids->vdwB_grid[x1][y0][z0]*xd);
    c10=(Grids->vdwB_grid[x0][y1][z0]*(1-xd)) + (Grids->vdwB_grid[x1][y1][z0]*xd);
    c01=(Grids->vdwB_grid[x0][y0][z1]*(1-xd)) + (Grids->vdwB_grid[x1][y0][z1]*xd);
    c11=(Grids->vdwB_grid[x0][y1][z1]*(1-xd)) + (Grids->vdwB_grid[x1][y1][z1]*xd);

    c0=(c00*(1-yd))+(c10*yd);
    c1=(c01*(1-yd))+(c11*yd);

    GI->vdwB=(c0*(1-zd))+(c1*zd);

    // HB donor grid

    c00=(Grids->hb_donor_grid[x0][y0][z0]*(1-xd)) + (Grids->hb_donor_grid[x1][y0][z0]*xd);
    c10=(Grids->hb_donor_grid[x0][y1][z0]*(1-xd)) + (Grids->hb_donor_grid[x1][y1][z0]*xd);
    c01=(Grids->hb_donor_grid[x0][y0][z1]*(1-xd)) + (Grids->hb_donor_grid[x1][y0][z1]*xd);
    c11=(Grids->hb_donor_grid[x0][y1][z1]*(1-xd)) + (Grids->hb_donor_grid[x1][y1][z1]*xd);

    c0=(c00*(1-yd))+(c10*yd);
    c1=(c01*(1-yd))+(c11*yd);

    GI->hb_donor=(c0*(1-zd))+(c1*zd);

    // HB acceptor grid

    c00=(Grids->hb_acceptor_grid[x0][y0][z0]*(1-xd)) + (Grids->hb_acceptor_grid[x1][y0][z0]*xd);
    c10=(Grids->hb_acceptor_grid[x0][y1][z0]*(1-xd)) + (Grids->hb_acceptor_grid[x1][y1][z0]*xd);
    c01=(Grids->hb_acceptor_grid[x0][y0][z1]*(1-xd)) + (Grids->hb_acceptor_grid[x1][y0][z1]*xd);
    c11=(Grids->hb_acceptor_grid[x0][y1][z1]*(1-xd)) + (Grids->hb_acceptor_grid[x1][y1][z1]*xd);

    c0=(c00*(1-yd))+(c10*yd);
    c1=(c01*(1-yd))+(c11*yd);

    GI->hb_acceptor=(c0*(1-zd))+(c1*zd);


    if ((Input->scoring_function == 0) or (Input->scoring_function == 2) or (Input->scoring_function == 4)){

        //Receptor Desolvation
        c00=(Grids->rec_solv_gauss[x0][y0][z0]*(1-xd)) + (Grids->rec_solv_gauss[x1][y0][z0]*xd);
        c10=(Grids->rec_solv_gauss[x0][y1][z0]*(1-xd)) + (Grids->rec_solv_gauss[x1][y1][z0]*xd);
        c01=(Grids->rec_solv_gauss[x0][y0][z1]*(1-xd)) + (Grids->rec_solv_gauss[x1][y0][z1]*xd);
        c11=(Grids->rec_solv_gauss[x0][y1][z1]*(1-xd)) + (Grids->rec_solv_gauss[x1][y1][z1]*xd);

        c0=(c00*(1-yd))+(c10*yd);
        c1=(c01*(1-yd))+(c11*yd);

        GI->rec_solv_gauss=(c0*(1-zd))+(c1*zd);

        // Ligand Desolvation

        c00=(Grids->solv_gauss[x0][y0][z0]*(1-xd)) + (Grids->solv_gauss[x1][y0][z0]*xd);
        c10=(Grids->solv_gauss[x0][y1][z0]*(1-xd)) + (Grids->solv_gauss[x1][y1][z0]*xd);
        c01=(Grids->solv_gauss[x0][y0][z1]*(1-xd)) + (Grids->solv_gauss[x1][y0][z1]*xd);
        c11=(Grids->solv_gauss[x0][y1][z1]*(1-xd)) + (Grids->solv_gauss[x1][y1][z1]*xd);

        c0=(c00*(1-yd))+(c10*yd);
        c1=(c01*(1-yd))+(c11*yd);

        GI->solv_gauss=(c0*(1-zd))+(c1*zd);
    }
    else {
        GI->rec_solv_gauss =0.0;
        GI->solv_gauss=0.0;
    }
}

double Energy2::evaluate_forces_hardcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz){
    double elec =0.0, vdwA=0.0, vdwB = 0.00;
    int a1=0, b1=0, c1=0, a2=0, b2=0, c2=0;
    double af, bf, cf, d;
    for (int i=0; i< Lig->N; i++){
        af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
        bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
        cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
        a1 = floor(af);
        b1 = floor(bf);
        c1 = floor(cf);
        a2 = ceil(af);
        b2 = ceil(bf);
        c2 = ceil(cf);

        if (a1 > 0 and b1 > 0 and c1 > 0 and a2 < int(Grids->elec_grid.size()) and b2 < int(Grids->elec_grid[0].size()) and c2 < int(Grids->elec_grid[0][0].size())){

            GridInterpol* GI = new GridInterpol;
            this->trilinear_interpolation(Grids, af, bf, cf, a1, b1, c1, a2, b2, c2, GI);

            d = this->distance(Lig->xyz[i][0], af, Lig->xyz[i][1], bf, Lig->xyz[i][2], cf);
            elec += Lig->charges[i]*GI->elec/d;
            vdwA += (GI->vdwA * sqrt(Lig->epsilons[i]*pow((2*Lig->radii[i]), 12)))/(12*d);
            vdwB += (sqrt(2.0 * Lig->epsilons[i]*pow((2*Lig->radii[i]), 6)) * GI->vdwB) / (6*d);
        }
    }
    return (elec+vdwA-vdwB);
}



#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(pyEnergy2)
{

    double (Energy2::*cess1)(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz) = &Energy2::compute_energy_softcore_solvation;
    double (Energy2::*cess2)(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result) = &Energy2::compute_energy_softcore_solvation;

    double (Energy2::*ces1)(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz) = &Energy2::compute_energy_softcore;
    double (Energy2::*ces2)(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result) = &Energy2::compute_energy_softcore;

    double (Energy2::*cehs1)(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz) = &Energy2::compute_energy_hardcore_solvation;
    double (Energy2::*cehs2)(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result) = &Energy2::compute_energy_hardcore_solvation;

    double (Energy2::*ceh1)(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz) = &Energy2::compute_energy_hardcore;
    double (Energy2::*ceh2)(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result) = &Energy2::compute_energy_hardcore;

    double (Energy2::*cefghs1)(Grid* Grids, Mol2* Lig, vector<vector<double> >xyz) = &Energy2::compute_ene_from_grids_hardcore_solvation;
    double (Energy2::*cefghs2)(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result) = &Energy2::compute_ene_from_grids_hardcore_solvation;

    double (Energy2::*cefgh1)(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz) = &Energy2::compute_ene_from_grids_hardcore;
    double (Energy2::*cefgh2)(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result) = &Energy2::compute_ene_from_grids_hardcore;

    double (Energy2::*cefgss1)(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz) = &Energy2::compute_ene_from_grids_softcore_solvation;
    double (Energy2::*cefgss2)(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result) = &Energy2::compute_ene_from_grids_softcore_solvation;

    double (Energy2::*cefgs1)(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz) = &Energy2::compute_ene_from_grids_softcore;
    double (Energy2::*cefgs2)(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result) = &Energy2::compute_ene_from_grids_softcore;

    double (Energy2::*ce1)(Mol2* Rec,   Mol2* Lig, vector<vector<double> >lig_xyz) = &Energy2::compute_ene;
    double (Energy2::*ce2)(Grid *Grids, Mol2* Lig, vector<vector<double> >lig_xyz) = &Energy2::compute_ene;
    double (Energy2::*ce3)(Grid* Grids, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result) = &Energy2::compute_ene;
    double (Energy2::*ce4)(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result) = &Energy2::compute_ene;


    double (Energy2::*ti1)(vector<vector<vector<double> > > grid, double x, double y, double z, int x0, int y0, int z0, int x1, int y1, int z1) = &Energy2::trilinear_interpolation;
    void (Energy2::*ti2)(Grid* Grids, double x, double y, double z, int x0, int y0, int z0, int x1, int y1, int z1, Energy2::GridInterpol* GI) = &Energy2::trilinear_interpolation;



    class_<Energy2>("Energy2", init<PARSER*>())
        .def("atom_is_acceptor",  & Energy2::atom_is_acceptor)
        .def("H_is_donor",  & Energy2::H_is_donor)
        .def("angle",  & Energy2::angle)
        .def("distance",  & Energy2::distance)
        .def("distance_squared",  & Energy2::distance_squared)

        .def("compute_energy_softcore_solvation", cess1)
        .def("compute_energy_softcore", ces1)
        .def("compute_energy_hardcore_solvation", cehs1)
        .def("compute_energy_hardcore", ceh1)
        .def("compute_ene_from_grids_softcore_solvation", cefgss1)
        .def("compute_ene_from_grids_softcore", cefgs1)
        .def("compute_ene_from_grids_hardcore_solvation", cefghs1)
        .def("compute_ene_from_grids_hardcore", cefgh1)


        .def("compute_energy_softcore_solvation", cess2)
        .def("compute_energy_softcore", ces2)
        .def("compute_energy_hardcore_solvation", cehs2)
        .def("compute_energy_hardcore", ceh2)

        .def("compute_ene", ce1)
        .def("compute_ene", ce2)
        .def("compute_ene", ce3)
        .def("compute_ene", ce4)

        .def("compute_ene_from_grids_hardcore_solvation", cefghs2)
        .def("compute_ene_from_grids_hardcore", cefgh2)
        .def("compute_ene_from_grids_softcore_solvation", cefgss2)
        .def("compute_ene_from_grids_softcore", cefgs2)

        .def("trilinear_interpolation", ti1)
        .def("trilinear_interpolation", ti2)

        .def("evaluate_forces_hardcore_solvation",  & Energy2::evaluate_forces_hardcore_solvation)
    ;

    class_<Energy2::GridInterpol>("GridInterpol")
            .def_readwrite("vdwA", & Energy2::GridInterpol::vdwA)
            .def_readwrite("vdwB", & Energy2::GridInterpol::vdwB)
            .def_readwrite("elec", & Energy2::GridInterpol::elec)
            .def_readwrite("pbsa", & Energy2::GridInterpol::pbsa)
            .def_readwrite("delphi", & Energy2::GridInterpol::delphi)
            .def_readwrite("solv_gauss", & Energy2::GridInterpol::solv_gauss)
            .def_readwrite("rec_solv_gauss", & Energy2::GridInterpol::rec_solv_gauss)
            .def_readwrite("hb_donor", & Energy2::GridInterpol::hb_donor)
            .def_readwrite("hb_acceptor", & Energy2::GridInterpol::hb_acceptor)
    ;

}

