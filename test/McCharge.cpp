#include <iostream>
#include <fstream>
#include "../LiBELa/Mol2.cpp"
#include "../LiBELa/PARSER.cpp"

using namespace std;

int main(int argc, char *argv[]){

	double alpha = 0.2;
	double beta = -0.005;
	string infile;
	int c;

	if (argc < 2){
      		printf("Usage %s -i <input.mol2> -a <alpha> -b <beta> [-h]\n", argv[0]);
      		exit(1);
	}

    	while ((c = getopt(argc, argv, "i:a:b:h")) != -1)
      		switch (c){
      			case 'a':
          			alpha= double(atof(optarg));
          			break;
      			case 'i':
          			infile = string(optarg);
          			break;
      			case 'h':
          			printf("Usage %s -i <inputfile> -a <alpha> -b <beta> [-h]\n", argv[0]);
          			break;
          			exit(1);
      			case '?':
          			printf("Usage %s -i <inputfile> -a <alpha> -b <beta> [-h]\n", argv[0]);
				break;
          			exit(1);
      			case 'b':
          			beta = double(atof(optarg));
          			break;
      		}

    
    PARSER* Input = new PARSER;
//    printf("Reading file %s...\n", argv[1]);
    Mol2* mol1 = new Mol2(Input, infile.c_str());
	double net_charge=0.0;
	double solvation=0.0;
	for (int i=0; i< mol1->N; i++){
		net_charge += mol1->charges[i];
		solvation += (alpha*mol1->charges[i]*mol1->charges[i])+beta;
	} 
	delete mol1;
	printf("%-40.40s % 7.5f % 7.5f %7d\n", infile.c_str(), solvation, net_charge, int(round(net_charge)));
	return 0;
}

