#include "pyMcEntropy.h"

McEntropy::McEntropy(PARSER* _Input, COORD_MC* _Coord, vector<double> _com, int _n_rot)
{
    this->k=0.0019858775203792202;          // Boltzmann constant in kcal/(mol.K)
    this->Input = _Input;
    this->Coord = _Coord;
    this->n_rot = _n_rot;
    this->com = _com;
    this->rot_bins = Input->entropy_rotation_bins;
    this->trans_bins = Input->entropy_translation_bins;
    this->translation_window = int(round(Input->x_dim*2));
    translation_step = translation_window*1.0/trans_bins;
    rotation_step = 360.0/rot_bins;
    vector<double> vtmp;

    for (int i=0; i< this->trans_bins; i++){
        hist_x.push_back(0.0);
        hist_y.push_back(0.0);
        hist_z.push_back(0.0);
    }

    for (int i=0; i< this->rot_bins; i++){
        hist_alpha.push_back(0.0);
        hist_gamma.push_back(0.0);
        vtmp.push_back(0.0);
    }

    for (int i=0; i<(this->rot_bins/2); i++){
        hist_beta.push_back(0.0);
    }

    for (int i=0; i<this->n_rot; i++){
          hist_torsions.push_back(vtmp);
    }
}

void McEntropy::update(double x, double y, double z, double alpha, double beta, double gamma, vector<double> torsion){

    unsigned dx, dy, dz, a, b, g, angle;
    dx = unsigned(round((x-com[0]+(translation_window*1.0/2.))/(translation_step)));
    dy = unsigned(round((y-com[1]+(translation_window*1.0/2.))/(translation_step)));
    dz = unsigned(round((z-com[2]+(translation_window*1.0/2.))/(translation_step)));
    hist_x[dx] += 1.0;
    hist_y[dy] += 1.0;
    hist_z[dz] += 1.0;

    a = unsigned(round(alpha/rotation_step));
    hist_alpha[a] += 1.0;
    b = unsigned(round(beta/rotation_step));
    hist_beta[b] += 1.0;
    g = unsigned(round(gamma/rotation_step));
    hist_gamma[g] += 1.0;

    //    hist_beta[int(round(beta/(rotation_step/2.0)))] += 1.0;

    for (unsigned i=0; i< unsigned(this->n_rot); i++){
        angle = unsigned(round(torsion[i]/rotation_step));
        hist_torsions[i][angle] += 1.0;
    }
}

void McEntropy::update_trajectory(double x, double y, double z, double alpha, double beta, double gamma, vector<double> torsion){

    int dx, dy, dz, a, b, g, angle;
    dx = int(round((x+(translation_window*1.0/2.))/(translation_step)));
    dy = int(round((y+(translation_window*1.0/2.))/(translation_step)));
    dz = int(round((z+(translation_window*1.0/2.))/(translation_step)));
    hist_x[dx] += 1.0;
    hist_y[dy] += 1.0;
    hist_z[dz] += 1.0;

    a = int(round(alpha/rotation_step));
    hist_alpha[a] += 1.0;
    b = int(round(beta/rotation_step));
    hist_beta[b] += 1.0;
    g = int(round(gamma/rotation_step));
    hist_gamma[g] += 1.0;

    //    hist_beta[int(round(beta/(rotation_step/2.0)))] += 1.0;

    for (int i=0; i< this->n_rot; i++){
        angle = int(round(torsion[i]/rotation_step));
        hist_torsions[i][angle] += 1.0;
    }
}

void McEntropy::get_results(entropy_t* entropy, entropy_t* max_entropy, int count){
    for (int i=0; i< trans_bins; i++){
        hist_x[i] = hist_x[i]/count;
        if (hist_x[i] > 0.0){
            entropy->Strans += hist_x[i] * log(hist_x[i]);
        }

        hist_y[i] = hist_y[i]/count;
        if (hist_y[i] > 0.0){
            entropy->Strans += hist_y[i] * log(hist_y[i]);
        }

        hist_z[i] = hist_z[i]/count;
        if (hist_z[i] > 0.0){
            entropy->Strans += hist_z[i] * log(hist_z[i]);
        }
    }

    for (int i=0; i< rot_bins; i++){
        hist_alpha[i] = hist_alpha[i]/count;
        if (hist_alpha[i]> 0.0){
            entropy->Srot += hist_alpha[i] * log(hist_alpha[i]);
        }

        hist_gamma[i] = hist_gamma[i]/count;
        if (hist_gamma[i]> 0.0){
            entropy->Srot += hist_gamma[i] * log(hist_gamma[i]);
        }

        for (int j=0; j< n_rot; j++){
            hist_torsions[j][i] = hist_torsions[j][i]/count;
            if (hist_torsions[j][i] > 0.0){
                entropy->Storsion += hist_torsions[j][i] * log(hist_torsions[j][i]);
            }
        }
    }

    for (int i=0; i< rot_bins/2; i++){
        hist_beta[i] = hist_beta[i]/count;
        if (hist_beta[i]> 0.0){
            entropy->Srot += hist_beta[i] * log(hist_beta[i]);
        }
    }

    entropy->Strans = -k*entropy->Strans;
    entropy->Srot = -k*entropy->Srot;
    entropy->Storsion = -k*entropy->Storsion;
    entropy->S = entropy->Strans + entropy->Srot + entropy->Storsion;
    entropy->TS = Input->temp*entropy->S;

    // Computing results for maximal entropy

    max_entropy->Srot = 0.0;
    double rot_bin_prob = 1.0/rot_bins;
    for (int i=0; i< this->rot_bins; i++){
        max_entropy->Srot += rot_bin_prob * log (rot_bin_prob); // for alpha
        max_entropy->Srot += rot_bin_prob * log (rot_bin_prob); // for gamma;
    }
    for (int i=0; i< this->rot_bins/2; i++){
        max_entropy->Srot += (2*rot_bin_prob) * log (2*rot_bin_prob); // for beta;
    }

    max_entropy->Strans = 0.0;
    double trans_bin_prob = 1.0/trans_bins;
    for (int i=0; i<this->trans_bins; i++){
        max_entropy->Strans +=  3*(trans_bin_prob * log(trans_bin_prob));
    }

    max_entropy->Storsion=0.0;

    for (int j=0; j< this->n_rot; j++){
        for (int i=0; i< this->rot_bins; i++){
            max_entropy->Storsion += rot_bin_prob * log(rot_bin_prob);
        }
    }
    max_entropy->Strans = -k*max_entropy->Strans;
    max_entropy->Srot = -k*max_entropy->Srot;
    max_entropy->Storsion = -k*max_entropy->Storsion;
    max_entropy->S = max_entropy->Strans + max_entropy->Srot + max_entropy->Storsion;
    max_entropy->TS = Input->temp*max_entropy->S;
}

void McEntropy::get_results(entropy_t* entropy, int count){
    for (int i=0; i< trans_bins; i++){
        hist_x[i] = hist_x[i]/count;
        if (hist_x[i] > 0.0){
            entropy->Strans += hist_x[i] * log(hist_x[i]);
        }

        hist_y[i] = hist_y[i]/count;
        if (hist_y[i] > 0.0){
            entropy->Strans += hist_y[i] * log(hist_y[i]);
        }

        hist_z[i] = hist_z[i]/count;
        if (hist_z[i] > 0.0){
            entropy->Strans += hist_z[i] * log(hist_z[i]);
        }
    }

    for (int i=0; i< rot_bins; i++){
        hist_alpha[i] = hist_alpha[i]/count;
        if (hist_alpha[i]> 0.0){
            entropy->Srot += hist_alpha[i] * log(hist_alpha[i]);
        }

        hist_gamma[i] = hist_gamma[i]/count;
        if (hist_gamma[i]> 0.0){
            entropy->Srot += hist_gamma[i] * log(hist_gamma[i]);
        }


        if (Input->sample_torsions){
            for (int j=0; j< n_rot; j++){
                hist_torsions[j][i] = hist_torsions[j][i]/count;
                if (hist_torsions[j][i] > 0.0){
                    entropy->Storsion += hist_torsions[j][i] * log(hist_torsions[j][i]);
                }
            }
        }
    }

    for (int i=0; i< rot_bins/2; i++){
        hist_beta[i] = hist_beta[i]/count;
        if (hist_beta[i]> 0.0){
            entropy->Srot += hist_beta[i] * log(hist_beta[i]);
        }
    }

    entropy->Strans = -k*entropy->Strans;
    entropy->Srot = -k*entropy->Srot;
    entropy->Storsion = -k*entropy->Storsion;
    entropy->S = entropy->Strans + entropy->Srot + entropy->Storsion;
    entropy->TS = Input->temp*entropy->S;
}

McEntropy::~McEntropy(){

    this->hist_x.clear();
    this->hist_y.clear();
    this->hist_z.clear();

    this->hist_alpha.clear();
    this->hist_beta.clear();
    this->hist_gamma.clear();

    this->hist_torsions.clear();

}





#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;


BOOST_PYTHON_MODULE(pyMcEntropy)
{

    void (McEntropy::*gr1)(McEntropy::entropy_t* entropy, McEntropy::entropy_t *max_entropy, int count) =&McEntropy::get_results;
    void (McEntropy::*gr2)(McEntropy::entropy_t* entropy, int count) =&McEntropy::get_results;


    class_<McEntropy>("McEntropy", init<PARSER*, COORD_MC*, vector<double>, int >())
        .def_readwrite("rot_bins", & McEntropy::rot_bins)
        .def_readwrite("trans_bins", & McEntropy::trans_bins)
        .def_readwrite("translation_window", & McEntropy::translation_window)
        .def_readwrite("n_rot", & McEntropy::n_rot)
        .def_readwrite("translation_step", & McEntropy::translation_step)
        .def_readwrite("rotation_step", & McEntropy::rotation_step)
        .def_readwrite("hist_x", & McEntropy::hist_x)
        .def_readwrite("hist_y", & McEntropy::hist_y)
        .def_readwrite("hist_z", & McEntropy::hist_z)
        .def_readwrite("hist_alpha", & McEntropy::hist_alpha)
        .def_readwrite("hist_beta", & McEntropy::hist_beta)
        .def_readwrite("hist_gamma", & McEntropy::hist_gamma)
        .def_readwrite("com", & McEntropy::com)
        .def_readwrite("hist_torsions", & McEntropy::hist_torsions)

        .def_readwrite("Input", & McEntropy::Input)
        .def_readwrite("Coord", & McEntropy::Coord)
        .def_readwrite("k", & McEntropy::k)

        .def("update", & McEntropy::update)
        .def("update_trajectory", & McEntropy::update_trajectory)

        .def("get_results", gr1)
        .def("get_results", gr2)


    ;

    class_<McEntropy::entropy_t>("entropy_t")
        .def_readwrite("Strans", & McEntropy::entropy_t::Strans)
        .def_readwrite("Srot", & McEntropy::entropy_t::Srot)
        .def_readwrite("Storsion", & McEntropy::entropy_t::Storsion)
        .def_readwrite("S", & McEntropy::entropy_t::S)
        .def_readwrite("TS", & McEntropy::entropy_t::TS)
    ;
}
