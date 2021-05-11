//
// Created by Magda on 25.04.2021.
//

#ifndef AT_2D_CLASS_LATTICE_H
#define AT_2D_CLASS_LATTICE_H

#include "Node.h"

#include <iostream>
#include <random>
#include <chrono>

class Lattice {
public:
	Lattice(int, double, double, double, double, double, double, double, double);
	~Lattice();
	void printS();
	void printSigma();
	void printBoth();
	void init(double);
	void monteCarloStep();
	double getE();
	double gSigma();
	double gS();
	double dis();
	double magnetizationS();
	double magnetizationSigma();

	Node** lattice;
	int n;
	double J_1;
	double J_2;
	double K_1;
	double K_2;
	double M_0;
	double R_1;
	double R_2;
	double T;
	// tables of the nearest and next to the nearest neighbours
	int* iP1;
	int* iN1;
	int* iP2;
	int* iN2;
};


#endif //AT_2D_CLASS_LATTICE_H
