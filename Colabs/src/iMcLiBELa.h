#ifndef IMCLIBELA_H
#define IMCLIBELA_H

#define ORGANIZATION "Universidade de Sao Paulo"
#define NAME "iMcLiBELa"
#define LONG_NAME "iMcLiBELa - Monte Carlo-Based Ligand Binding Energy Landscape"
#define VERSION "1.0"
#define PI 3.14159265359


struct energy_result_t{
    double vdw;
    double elec;
    double rec_solv;
    double lig_solv;
    double hb_donor;
    double hb_acceptor;
    double restraints;
    double total;
};

#endif
