/*
 * Optimizer.cpp
 *
 *  Created on: 21/03/2012
 *      Author: Nascimento
 */

#include "pyOptimizer.h"

Optimizer::Optimizer(Mol2* _Rec, Mol2* _RefLig, PARSER* _Parser) {
    this->RefLig = _RefLig;
    this->Rec = _Rec;
    this->Parser = _Parser;
}
Optimizer::Optimizer(Mol2* _Rec, Mol2* _RefLig, PARSER* _Parser, Grid* _Grids) {
    this->RefLig = _RefLig;
    this->Rec = _Rec;
    this->Parser = _Parser;
    this->Grids = _Grids;
}

Optimizer::~Optimizer() {
}

double Optimizer::distance(vector<double> atom1, vector<double> atom2) {
    return ( sqrt(((atom2[0]-atom1[0])*(atom2[0]-atom1[0]))+((atom2[1]-atom1[1])*(atom2[1]-atom1[1]))+((atom2[2]-atom1[2])*(atom2[2]-atom1[2]))) );
}
double Optimizer::distance_squared(vector<double> atom1, vector<double> atom2) {
    return ( (((atom2[0]-atom1[0])*(atom2[0]-atom1[0]))+((atom2[1]-atom1[1])*(atom2[1]-atom1[1]))+((atom2[2]-atom1[2])*(atom2[2]-atom1[2]))) );
}

vector<vector<double> > Optimizer::update_coords(const std::vector<double> &x, Mol2* Lig2){
    COORD_MC* Coord = new COORD_MC;
    return(Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4], x[5]));
}

double Optimizer::evaluate_rmsd(Mol2* Lig1, Mol2* Lig2){
    double rmsd=0.0;
    for (int i=0; i<Lig1->N; i++){
        for (int j=0; j< 3; j++){
            rmsd += (Lig1->xyz[i][j] - Lig2->xyz[i][j]) * (Lig1->xyz[i][j] - Lig2->xyz[i][j]);
        }
    }
    rmsd = sqrt(rmsd);
    return(rmsd);
}

double Optimizer::evaluate_energy(Mol2* Lig2, vector<vector<double> > new_xyz){
    Energy2* Ene = new Energy2(Parser);
    energy_result_t* energy_results = new energy_result_t;
    double energy=0.0;
    double e_restraints = 0.0;
    double d2;
    if (Parser->use_grids){
        energy = Ene->compute_ene(Grids, Lig2, new_xyz, energy_results);
    }
    else {
        energy = Ene->compute_ene(Rec, Lig2, new_xyz, energy_results);
    }
    if (Parser->use_Erestraints){
        for (int i=0; i<Lig2->N; i++){
            d2 = 0.0;
            for (int j=0; j< 3; j++){
                d2 += (Lig2->xyz[i][j]-new_xyz[i][j])*(Lig2->xyz[i][j]-new_xyz[i][j]);
            }
            e_restraints += 0.5*Parser->restraints_weight*d2;
        }
    }
    energy += e_restraints;
    delete Ene;
    delete energy_results;
    return(energy);
}

void Optimizer::evaluate_energy2(Mol2* Lig2, vector<vector<double> > new_xyz, energy_result_t* energy_result){
    Energy2* Ene = new Energy2(Parser);
    if (Parser->use_grids){
        Ene->compute_ene(Grids, Lig2, new_xyz, energy_result);
    }
    else{
        Ene->compute_ene(Rec, Lig2, new_xyz, energy_result);
    }
    double e_restraints = 0.0;
    double d2;
    if (Parser->use_Erestraints){
        for (int i=0; i<Lig2->N; i++){
            d2 = 0.0;
            for (int j=0; j< 3; j++){
                d2 += (Lig2->xyz[i][j]-new_xyz[i][j])*(Lig2->xyz[i][j]-new_xyz[i][j]);
            }
            e_restraints += 0.5*Parser->restraints_weight*d2;
        }
    }
    energy_result->restraints = e_restraints;
    energy_result->total += energy_result->restraints;
}

double Optimizer::objective_energy_function(const std::vector<double> &x, std::vector<double> &grad, void *data){
    vector<vector<double> > new_xyz;
    COORD_MC* Coord = new COORD_MC;
    Gaussian* Gauss = new Gaussian;
    double f, f2, t1, t2, t3, si;

    Mol2* Lig2 = (Mol2*) data;
    new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4], x[5]);
    f = evaluate_energy(Lig2, new_xyz);
    if (Parser->use_score_optimization){
        t1 = (Gauss->compute_shape_and_charge_density(Parser, RefLig, RefLig, RefLig->xyz));
        t3 = (Gauss->compute_shape_and_charge_density(Parser, Lig2, Lig2, Lig2->xyz));
        t2 = (Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz));
        si = (2*t2)/(t1+t3);
        f=f*si;
    }

    if(!grad.empty()){
        new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0]+Parser->min_delta, x[1], x[2], x[3], x[4], x[5]);
        f2 = evaluate_energy(Lig2, new_xyz);
        if (Parser->use_score_optimization){
            t2 = (Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz));
            si = (2*t2)/(t1+t3);
            f2=f2*si;
        }
        grad[0] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1]+Parser->min_delta, x[2], x[3], x[4], x[5]);
        f2 = evaluate_energy(Lig2, new_xyz);
        if (Parser->use_score_optimization){
            t2 = (Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz));
            si = (2*t2)/(t1+t3);
            f2=f2*si;
        }
        grad[1] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2]+Parser->min_delta, x[3], x[4], x[5]);
        f2 = evaluate_energy(Lig2, new_xyz);
        if (Parser->use_score_optimization){
            t2 = (Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz));
            si = (2*t2)/(t1+t3);
            f2=f2*si;
        }
        grad[2] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3]+Parser->min_delta, x[4], x[5]);
        f2 = evaluate_energy(Lig2, new_xyz);
        if (Parser->use_score_optimization){
            t2 = (Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz));
            si = (2*t2)/(t1+t3);
            f2=f2*si;
        }
        grad[3] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4]+Parser->min_delta, x[5]);
        f2 = evaluate_energy(Lig2, new_xyz);
        if (Parser->use_score_optimization){
            t2 = (Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz));
            si = (2*t2)/(t1+t3);
            f2=f2*si;
        }
        grad[4] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4], x[5]+Parser->min_delta);
        f2 = evaluate_energy(Lig2, new_xyz);
        if (Parser->use_score_optimization){
            t2 = (Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz));
            si = (2*t2)/(t1+t3);
            f2=f2*si;
        }
        grad[5] = (f2-f)/Parser->min_delta;
    }
    delete Coord;
    delete Gauss;
    return (f);
}

double Optimizer::objective_overlay_function(const std::vector<double> &x, std::vector<double> &grad, void *data){
    Gaussian* Gauss = new Gaussian;
    COORD_MC* Coord = new COORD_MC;
    vector<vector<double> > new_xyz;
    double f, f2;
    Mol2* Lig2 = (Mol2*) data;

    new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4], x[5]);

    f = (Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz));

    if (!grad.empty()){
        new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0] + Parser->min_delta, x[1], x[2], x[3], x[4], x[5]);
        f2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz);
        grad[0] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(Lig2->xyz, Lig2,x[0], x[1] + Parser->min_delta, x[2], x[3], x[4], x[5]);
        f2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz);
        grad[1] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2]+ Parser->min_delta, x[3], x[4], x[5]);
        f2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz);
        grad[2] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3]+ Parser->min_delta, x[4], x[5]);
        f2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz);
        grad[3] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4]+Parser->min_delta, x[5]);
        f2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz);
        grad[4] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4], x[5]+Parser->min_delta);
        f2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz);
        grad[5] = (f2-f)/Parser->min_delta;
    }
    delete Gauss;
    delete Coord;
    return (f);
}

double Optimizer::objective_prealign_function(const std::vector<double> &x, std::vector<double> &grad, void *data){
    COORD_MC* Coord = new COORD_MC;
    vector<vector<double> > new_xyz;
    double f=0.0, f2=0.0;
    Mol2* Lig2 = (Mol2*) data;

    new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4], x[5]);

    f += distance_squared(new_xyz[Lig2->longest_axis[0]], RefLig->xyz[RefLig->longest_axis[0]]);
    f += distance_squared(new_xyz[Lig2->longest_axis[1]], RefLig->xyz[RefLig->longest_axis[1]]);

    f2 += distance_squared(new_xyz[Lig2->longest_axis[0]], RefLig->xyz[RefLig->longest_axis[1]]);
    f2 += distance_squared(new_xyz[Lig2->longest_axis[1]], RefLig->xyz[RefLig->longest_axis[0]]);

    if (f2 < f){
        f = f2;
    }

    delete Coord;
    return (f);
}

void  Optimizer::minimize_energy_nlopt_cobyla(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LN_COBYLA,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);
    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void  Optimizer::minimize_energy_nlopt_lbfgs(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LD_LBFGS,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);
    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

}

void  Optimizer::minimize_energy_nlopt_ln_auglag(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LN_AUGLAG,6);
    nlopt::opt local_opt(nlopt::LD_MMA,6);

    local_opt.set_xtol_rel(Parser->dock_min_tol);
    local_opt.set_maxtime(Parser->min_timeout);
    opt->set_local_optimizer(local_opt);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);
    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_energy_nlopt_ld_auglag(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LD_AUGLAG,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] = -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);
    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

}

void  Optimizer::minimize_energy_nlopt_mma(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LD_MMA,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);

    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

}

void Optimizer::minimize_energy_nlopt_subplex(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LN_SBPLX,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);

    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

}

void Optimizer::minimize_energy_nlopt_simplex(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LN_NELDERMEAD,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);

    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

}

void Optimizer::minimize_energy_nlopt_crs(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GN_CRS2_LM,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);

    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

}

void Optimizer::minimize_energy_nlopt_direct(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GN_DIRECT_L_RAND,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres;

    try {
        nres = opt->optimize(x,f_minimum);
    }
    catch(...){
        if (Parser->verbose){
            printf("DIRECT optimization of molecule %s stopped with an exception.\n", Lig2->molname.c_str());
        }
    }

    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);

    delete opt;

    Lig2->xyz = xyz;
    nlopt::opt *opt2 = new nlopt::opt(nlopt::LN_AUGLAG,6);
    opt2->set_lower_bounds(lb);
    opt2->set_upper_bounds(ub);

    opt2->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt2->set_xtol_rel(Parser->dock_min_tol);
    opt2->set_maxtime(Parser->min_timeout);
    vector<double> x2(6);
    x2[0] = 0.0;
    x2[1] = 0.0;
    x2[2] = 0.0;
    x2[3] = 0.0;
    x2[4] = 0.0;
    x2[5] = 0.0;

    f_minimum=0.00;
    nres = opt2->optimize(x2,f_minimum);
    xyz = Optimizer::update_coords(x2, Lig2);
    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);

    if (f_minimum > 9999.9) {
        opt_result->f_min = 9999.9;
    }
    else {
        opt_result->f_min = f_minimum;
    }
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

    delete opt2;

}

void Optimizer::minimize_energy_nlopt_direct_only(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GN_DIRECT_L_RAND,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);

    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);
    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_energy_nlopt_stogo(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GD_STOGO_RAND,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] = -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);

    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

}

void Optimizer::minimize_energy_nlopt_isres(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GN_ISRES,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);

    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

}

void Optimizer::minimize_energy_nlopt_esch(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GN_ESCH,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);

    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

}



void Optimizer::minimize_overlay_nlopt_cobyla(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LN_COBYLA,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x);
    lb[4] = -(Parser->search_box_y);
    lb[5] = -(Parser->search_box_z);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x;
    ub[4] = Parser->search_box_y;
    ub[5] = Parser->search_box_z;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
    opt->set_xtol_rel(Parser->min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double fo;
    nlopt::result nres = opt->optimize(x,fo);
    vector<vector<double> > xyz = update_coords(x, Lig2);
    double f_minimum;
    if (Parser->energy_optimizer != "none"){
        f_minimum = evaluate_energy(Lig2, xyz);
    }
    else {
        f_minimum = 0.00;
    }
    delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_lbfgs(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LD_LBFGS,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x);
    lb[4] = -(Parser->search_box_y);
    lb[5] = -(Parser->search_box_z);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x;
    ub[4] = Parser->search_box_y;
    ub[5] = Parser->search_box_z;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
    opt->set_xtol_rel(Parser->min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double fo;
    nlopt::result nres = opt->optimize(x,fo);
    vector<vector<double> > xyz = update_coords(x, Lig2);
    double f_minimum;
    if (Parser->energy_optimizer != "none"){
        f_minimum = evaluate_energy(Lig2, xyz);
    }
    else {
        f_minimum = 0.00;
    }
    delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_ln_auglag(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LN_AUGLAG,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] = -(Parser->search_box_x);
    lb[4] = -(Parser->search_box_y);
    lb[5] = -(Parser->search_box_z);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x;
    ub[4] = Parser->search_box_y;
    ub[5] = Parser->search_box_z;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
    opt->set_xtol_rel(Parser->min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double fo;
    nlopt::result nres = opt->optimize(x,fo);
    vector<vector<double> > xyz = update_coords(x, Lig2);
    double f_minimum;
    if (Parser->energy_optimizer != "none"){
        f_minimum = evaluate_energy(Lig2, xyz);
    }
    else {
        f_minimum = 0.00;
    }
    delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_ld_auglag(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LD_AUGLAG,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x);
    lb[4] = -(Parser->search_box_y);
    lb[5] = -(Parser->search_box_z);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x;
    ub[4] = Parser->search_box_y;
    ub[5] = Parser->search_box_z;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
    opt->set_xtol_rel(Parser->min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double fo;
    nlopt::result nres = opt->optimize(x,fo);
    vector<vector<double> > xyz = update_coords(x, Lig2);
    double f_minimum;
    if (Parser->energy_optimizer != "none"){
        f_minimum = evaluate_energy(Lig2, xyz);
    }
    else {
        f_minimum = 0.00;
    }
    delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_mma(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LD_MMA,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] = -(Parser->search_box_x);
    lb[4] = -(Parser->search_box_y);
    lb[5] = -(Parser->search_box_z);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x;
    ub[4] = Parser->search_box_y;
    ub[5] = Parser->search_box_z;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
    opt->set_xtol_rel(Parser->min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double fo;
    nlopt::result nres = opt->optimize(x,fo);
    vector<vector<double> > xyz = update_coords(x, Lig2);
    double f_minimum;
    if (Parser->energy_optimizer != "none"){
        f_minimum = evaluate_energy(Lig2, xyz);
    }
    else {
        f_minimum = 0.00;
    }

    delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_subplex(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LN_SBPLX,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] = -(Parser->search_box_x);
    lb[4] = -(Parser->search_box_y);
    lb[5] = -(Parser->search_box_z);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x;
    ub[4] = Parser->search_box_y;
    ub[5] = Parser->search_box_z;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
    opt->set_xtol_rel(Parser->min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double fo;
    nlopt::result nres = opt->optimize(x,fo);
    vector<vector<double> > xyz = update_coords(x, Lig2);
    double f_minimum;
    if (Parser->energy_optimizer != "none"){
        f_minimum = evaluate_energy(Lig2, xyz);
    }
    else {
        f_minimum = 0.00;
    }
    delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_crs(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GN_CRS2_LM,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] = -(Parser->search_box_x);
    lb[4] = -(Parser->search_box_y);
    lb[5] = -(Parser->search_box_z);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x;
    ub[4] = Parser->search_box_y;
    ub[5] = Parser->search_box_z;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
    opt->set_xtol_rel(Parser->min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double fo;
    nlopt::result nres = opt->optimize(x,fo);
    vector<vector<double> > xyz = update_coords(x, Lig2);
    double f_minimum;
    if (Parser->energy_optimizer != "none"){
        f_minimum = evaluate_energy(Lig2, xyz);
    }
    else {
        f_minimum = 0.00;
    }
    delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_direct(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GN_DIRECT_L_RAND,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] = -(Parser->search_box_x);
    lb[4] = -(Parser->search_box_y);
    lb[5] = -(Parser->search_box_z);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x;
    ub[4] = Parser->search_box_y;
    ub[5] = Parser->search_box_z;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
    opt->set_xtol_rel(Parser->min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double fo;
    nlopt::result nres = opt->optimize(x,fo);
    vector<vector<double> > xyz = update_coords(x, Lig2);
    double f_minimum;
    if (Parser->energy_optimizer != "none"){
        f_minimum = evaluate_energy(Lig2, xyz);
    }
    else {
        f_minimum = 0.00;
    }
    delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void Optimizer::Simulated_Annealing(Mol2* Lig, opt_result_t* opt_result){
    gsl_rng* r;
    r = gsl_rng_alloc(gsl_rng_ranlxs2);
    srand(time(NULL));
    gsl_rng_set(r, (rand() % 100));
    vector<vector<double> > new_xyz;

    SA* Anneal = new SA();
    new_xyz = Anneal->optimize(Lig, Parser, Grids, r);
    delete Anneal;
    evaluate_energy2(Lig, new_xyz, opt_result->energy_result);

    opt_result->f_min = opt_result->energy_result->total;
    opt_result->optimized_xyz = new_xyz;
    opt_result->optimization_status = 0;
    gsl_rng_free(r);
}

void Optimizer::minimize_energy_adaptative(Mol2* Lig2, opt_result_t* opt_result){
    double tol = 1.0e-4;
    double deltaij = 2.75;
    double deltaij_es = 1.65;
    Parser->deltaij6 = pow(deltaij,6);
    Parser->deltaij_es6 = pow(deltaij_es,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    int count = 1;
    while (deltaij >= 0.0){
        cout << "Adaptative cycle " << count << " ..." << endl;
        cout << "deltaij = " << deltaij << endl;
        cout << "deltaij_es = " << deltaij_es << endl;

        nlopt::opt *opt = new nlopt::opt(nlopt::LD_AUGLAG,6);
        opt->set_lower_bounds(lb);
        opt->set_upper_bounds(ub);

        opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
        opt->set_xtol_rel(tol);
        opt->set_maxtime(Parser->min_timeout);

        vector<double> x(6);
        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 0.0;
        x[3] = 0.0;
        x[4] = 0.0;
        x[5] = 0.0;

        double f_minimum=0.00;
        nlopt::result nres = opt->optimize(x,f_minimum);
        vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
        Optimizer::evaluate_energy2(Lig2, xyz, opt_result->energy_result);
        delete opt;

        opt_result->f_min = f_minimum;
        opt_result->optimization_status = nres;
        opt_result->optimized_xyz = xyz;


        deltaij -= 0.25;
        deltaij_es -=0.15;

        Parser->deltaij6 = pow(deltaij, 6);
        Parser->deltaij_es6 = pow(deltaij_es, 6);
        count++;
    }
}

double Optimizer::superpose_function(const std::vector<double> &x, std::vector<double> &grad, void *data){
    unique_ptr<COORD_MC> Coord(new COORD_MC);
    vector<vector<double> > new_xyz;
    double f, f2;
    align_t* align_data= (align_t*) data;

    new_xyz = Coord->rototranslate(align_data->current_xyz, RefLig, x[0], x[1], x[2], x[3], x[4], x[5]);

    f = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));

    if (!grad.empty()){
        new_xyz = Coord->rototranslate(align_data->current_xyz, RefLig, x[0] + Parser->min_delta, x[1], x[2], x[3], x[4], x[5]);
        f2 = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));
        grad[0] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(align_data->current_xyz, RefLig, x[0], x[1] + Parser->min_delta, x[2], x[3], x[4], x[5]);
        f2 = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));
        grad[1] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(align_data->current_xyz, RefLig, x[0], x[1], x[2]+ Parser->min_delta, x[3], x[4], x[5]);
        f2 = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));
        grad[2] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(align_data->current_xyz, RefLig, x[0], x[1], x[2], x[3]+ Parser->min_delta, x[4], x[5]);
        f2 = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));
        grad[3] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(align_data->current_xyz, RefLig, x[0], x[1], x[2], x[3], x[4]+Parser->min_delta, x[5]);
        f2 = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));
        grad[4] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(align_data->current_xyz, RefLig, x[0], x[1], x[2], x[3], x[4], x[5]+Parser->min_delta);
        f2 = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));
        grad[5] = (f2-f)/Parser->min_delta;
    }
    return (f);
}

void Optimizer::minimize_alignment_nlopt_simplex(align_t* align_data, align_result_t* opt_result, vector<double> current_com){
    unique_ptr<nlopt::opt> opt(new nlopt::opt(nlopt::LN_NELDERMEAD,6));

    unique_ptr<COORD_MC> Coord(new COORD_MC);

    vector<double> original_com = Coord->compute_com(RefLig);

    double dx = original_com[0] - current_com[0];
    double dy = original_com[1] - current_com[1];
    double dz = original_com[2] - current_com[2];

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = 0.0;
    lb[2] = -180.0;
    lb[3] = -(Parser->search_box_x);
    lb[4] = -(Parser->search_box_y);
    lb[5] = -(Parser->search_box_z);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 180.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x;
    ub[4] = Parser->search_box_y;
    ub[5] = Parser->search_box_z;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::superpose_function, align_data);
    opt->set_xtol_rel(Parser->min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = dx;
    x[4] = dy;
    x[5] = dz;

    double fo;
    nlopt::result nres = opt->optimize(x,fo);

    if (x[0] < 0.0){
        x[0] = x[0] + 360.0;
    }

    if (x[2] < 0.0 ){
        x[2] = x[2] + 360.0;
    }

    opt_result->rmsd = fo;
    opt_result->translation[0] = x[3];
    opt_result->translation[1] = x[4];
    opt_result->translation[2] = x[5];
    opt_result->rotation[0] = x[0];
    opt_result->rotation[1] = x[1];
    opt_result->rotation[2] = x[2];
    opt_result->opt_status = nres;
}

void Optimizer::pre_align(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LN_NELDERMEAD,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] = -3.0;
    lb[4] = -3.0;
    lb[5] = -3.0;
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = 3.0;
    ub[4] = 3.0;
    ub[5] = 3.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_prealign_function, Lig2);
    opt->set_xtol_rel(Parser->min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double fo;
    nlopt::result nres = opt->optimize(x,fo);
    vector<vector<double> > xyz = update_coords(x, Lig2);
    double f_minimum = 0.0;
    delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

Mol2* get_Rec(){
    return Optimizer::Rec;
}

//#define BOOST_PYTHON_STATIC_LIB
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/return_arg.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(pyOptimizer)
{ 

//    class_< vector< vector<double> > >("vectorvectorDouble")
//            .def(vector_indexing_suite< vector< vector<double> >  >())
//        ;
        boost::python::class_<Optimizer>("Optimizer", init<Mol2*, Mol2*, PARSER*>())
            .def(init<Mol2*, Mol2*, PARSER*, Grid*>())
//            .add_static_property("Rec", make_getter(Optimizer::Rec), make_setter(Optimizer::Rec))
              .def("get_Rec", &get_Rec, return_self<const>())
//            .def("set_ref_lig", &Optimizer::set_ref_lig)
//            .def("set_parser", &Optimizer::set_parser)
//            .def("set_grids", &Optimizer::set_grids)
//            .def("set_wirter", &Optimizer::set_writer)
//            .def("run", &Optimizer::run)
       //     .def_readwrite("Rec", &Optimizer::Rec)
       //     .def_readwrite("RefLig", &Optimizer::RefLig)
       //     .def_readwrite("Parser", &Optimizer::Parser)
       //     .def_readwrite("Grids", &Optimizer::Grids)
            .def("evaluate_rmsd", & Optimizer::evaluate_rmsd).staticmethod("evaluate_rmsd")
            .def("Simulated_Annealing", & Optimizer::Simulated_Annealing).staticmethod("Simulated_Annealing")
            .def("update_coords", & Optimizer::update_coords).staticmethod("update_coords")
            .def("distance", & Optimizer::distance).staticmethod("distance")
            .def("distance_squared", & Optimizer::distance_squared).staticmethod("distance_squared")

            .def("evaluate_energy", & Optimizer::evaluate_energy).staticmethod("evaluate_energy")
            .def("evaluate_energy2", & Optimizer::evaluate_energy2).staticmethod("evaluate_energy2")

            .def("objective_energy_function", & Optimizer::objective_energy_function).staticmethod("objective_energy_function")
            .def("objective_overlay_function", & Optimizer::objective_overlay_function).staticmethod("objective_overlay_function")
            .def("objective_prealign_function", & Optimizer::objective_prealign_function).staticmethod("objective_prealign_function")
            .def("superpose_function", & Optimizer::superpose_function).staticmethod("superpose_function")

            .def("minimize_energy_nlopt_cobyla", & Optimizer::minimize_energy_nlopt_cobyla).staticmethod("minimize_energy_nlopt_cobyla")
            .def("minimize_energy_nlopt_lbfgs", & Optimizer::minimize_energy_nlopt_lbfgs).staticmethod("minimize_energy_nlopt_lbfgs")
            .def("minimize_energy_nlopt_ln_auglag", & Optimizer::minimize_energy_nlopt_ln_auglag).staticmethod("minimize_energy_nlopt_ln_auglag")

            .def("minimize_energy_nlopt_ld_auglag", & Optimizer::minimize_energy_nlopt_ld_auglag).staticmethod("minimize_energy_nlopt_ld_auglag")
            .def("minimize_energy_nlopt_mma", & Optimizer::minimize_energy_nlopt_mma).staticmethod("minimize_energy_nlopt_mma")
            .def("minimize_energy_nlopt_subplex", & Optimizer::minimize_energy_nlopt_subplex).staticmethod("minimize_energy_nlopt_subplex")
            .def("minimize_energy_nlopt_simplex", & Optimizer::minimize_energy_nlopt_simplex).staticmethod("minimize_energy_nlopt_simplex")
            .def("minimize_energy_nlopt_crs", & Optimizer::minimize_energy_nlopt_crs).staticmethod("minimize_energy_nlopt_crs")
            .def("minimize_energy_nlopt_direct", & Optimizer::minimize_energy_nlopt_direct).staticmethod("minimize_energy_nlopt_direct")
            .def("minimize_energy_nlopt_direct_only", & Optimizer::minimize_energy_nlopt_direct_only).staticmethod("minimize_energy_nlopt_direct_only")
            .def("minimize_energy_nlopt_stogo", & Optimizer::minimize_energy_nlopt_stogo).staticmethod("minimize_energy_nlopt_stogo")
            .def("minimize_energy_nlopt_isres", & Optimizer::minimize_energy_nlopt_isres).staticmethod("minimize_energy_nlopt_isres")
            .def("minimize_energy_nlopt_esch", & Optimizer::minimize_energy_nlopt_esch).staticmethod("minimize_energy_nlopt_esch")
            .def("minimize_energy_adaptative", & Optimizer::minimize_energy_adaptative).staticmethod("minimize_energy_nlopt_adaptive")

            .def("minimize_overlay_nlopt_cobyla", & Optimizer::minimize_overlay_nlopt_cobyla).staticmethod("minimize_energy_nlopt_cobyla")
            .def("minimize_overlay_nlopt_subplex", & Optimizer::minimize_overlay_nlopt_subplex).staticmethod("minimize_energy_nlopt_subplex")
            .def("minimize_overlay_nlopt_lbfgs", & Optimizer::minimize_overlay_nlopt_lbfgs).staticmethod("minimize_energy_nlopt_lbfgs")
            .def("minimize_overlay_nlopt_ln_auglag", & Optimizer::minimize_overlay_nlopt_ln_auglag).staticmethod("minimize_energy_nlopt_ln_auglag")
            .def("minimize_overlay_nlopt_ld_auglag", & Optimizer::minimize_overlay_nlopt_ld_auglag).staticmethod("minimize_energy_nlopt_ld_auglag")
            .def("minimize_overlay_nlopt_mma", & Optimizer::minimize_overlay_nlopt_mma).staticmethod("minimize_energy_nlopt_mma")
            .def("minimize_overlay_nlopt_crs", & Optimizer::minimize_overlay_nlopt_crs).staticmethod("minimize_energy_nlopt_crs")
            .def("minimize_overlay_nlopt_direct", & Optimizer::minimize_overlay_nlopt_direct).staticmethod("minimize_energy_nlopt_direct")

            .def("minimize_alignment_nlopt_simplex", & Optimizer::minimize_alignment_nlopt_simplex).staticmethod("minimize_alignment_nlopt_simplex")

            .def("pre_align", & Optimizer::pre_align).staticmethod("pre_align")
            ;

            class_<Optimizer::opt_result_t>("opt_result_t")
                .def_readwrite("optimization_status", & Optimizer::opt_result_t::optimization_status)
                .def_readwrite("f_min", & Optimizer::opt_result_t::f_min)
                .def_readwrite("energy_result", & Optimizer::opt_result_t::energy_result)
                .def_readwrite("optimized_xyz", & Optimizer::opt_result_t::optimized_xyz)
            ;

            class_<Optimizer::align_result_t>("align_result_t")
                .def_readwrite("rmsd", & Optimizer::align_result_t::rmsd)
                .def_readwrite("translation", & Optimizer::align_result_t::translation)
                .def_readwrite("rotation", & Optimizer::align_result_t::rotation)
                .def_readwrite("opt_status", & Optimizer::align_result_t::opt_status)
            ;

            class_<Optimizer::align_t>("align_t")
                .def_readwrite("ref_xyz", & Optimizer::align_t::ref_xyz)
                .def_readwrite("current_xyz", & Optimizer::align_t::current_xyz)

           ;
/*    
        try {
                    // Initialize Optimizer module
                    Optimizer::Rec = new Mol2();
                    Optimizer::RefLig = new Mol2();
                    Optimizer::Parser = new PARSER();
                    Optimizer::Writer= new WRITER();
                    Optimizer::Grids = new Grid(Parser,Writer);

                } catch (std::exception const & e) {
                    PyErr_SetString(PyExc_RuntimeError, e.what());
                    return;
                }
*/                
}
