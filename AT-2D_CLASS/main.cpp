#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "Node.h"
#include "Lattice.h"
#include "func.h"

int main(int argc, char **argv) {
	// model parameters
	int n = 100; // size
	double d = 0.0; // dissonance between S and sigma in the initial state
	double T = 0.01;
	int MCS = 1000; // number of Monte Carlo Steps
	int repeat = 10;
	// parameters from command line T and A
	if (argc != 3) {
		return -1;
	}
	T = atof(argv[1]);
	double A = atof(argv[2]);
	// double	A = 0;
	////////////////////////////////
	double J_1 = -1.0;
	double J_2 = 1.0;
	double K_1 = 0;
	double K_2 = 0;
	double R_1 = 1.0 + A;
	double R_2 = 0;
	double M_0 = 2.0 * R_1;
	///////////////////////////////
	int skip = 0;

	// files management
	std::stringstream stream;
	stream << std::fixed << std::setprecision(2) << T;
	std::string s = stream.str();

	std::string time = "TIMEmodel2-L=" + str_i(n) + "-J_1=" + DoubleToString(J_1) + "-J_2=" + DoubleToString(J_2)
					   + "-K_1=" + DoubleToString(K_1) + "-K_2=" + DoubleToString(K_2)
					   + "-T=" + s + "-M_0=" + DoubleToString(M_0)
					   + "-R_1=" + DoubleToString(R_1) + "-R_2=" + DoubleToString(R_2)
					   + "-d=" + DoubleToString(d) + ".txt";
	std::ofstream avergeTimeFile(time.c_str());
	std::string diss = "DISmodel2-L=" + str_i(n) + "-J_1=" + DoubleToString(J_1) + "-J_2=" +
					   DoubleToString(J_2) + "-K_1=" + DoubleToString(K_1) + "-K_2=" + DoubleToString(K_2) +
					   "-T=" + s + "-M_0=" + DoubleToString(M_0)
					   + "-R_1=" + DoubleToString(R_1) + "-R_2=" + DoubleToString(R_2)
					   + "-d=" + DoubleToString(d) + ".txt";

	std::ofstream dis_file(diss.c_str());

	for (int r = 0; r < repeat; r++) { // repetitions loop
		double counter = 0;
		double mS = 0;
		double m2S = 0;
		double mSigma = 0;
		double m2Sigma = 0;
		Lattice square(n, J_1, J_2, K_1, K_2, M_0, R_1, R_2, T); // creating a lattice
		square.init(d);
		int stepKeep = 0;
		for (int step = 0; step < MCS; step++) { // Monte Carlo loop
			//dis_file << dis(S, sigma, n) << std::endl;
			if (step%10 == 0) {
				std::cout << "\t" << step;
			}
			stepKeep = step + 1;
			square.monteCarloStep();
			// magnetization
			if (step >= skip) {
				counter++;

				mS += abs(square.magnetizationS());
				m2S += square.magnetizationS() * square.magnetizationS();
				mSigma += abs(square.magnetizationSigma());
				m2Sigma += square.magnetizationSigma() * square.magnetizationSigma();
			}
			// stop condition
			if (square.gS() * square.gS() == 1.0 && square.gSigma() * square.gSigma() == 1.0) {
				break;
			}
		} // end MC loop
		avergeTimeFile << stepKeep << "\t" << square.gS() << "\t" << square.gSigma() << std::endl;
		// susceptibility
		double xS = 0;
		double xSigma = 0;
		xS = ((double) square.n * (double) square.n) * 1 / T * (m2S / counter - mS * mS / (counter * counter));
		xSigma = ((double) square.n * (double) square.n) * 1 / T * (m2Sigma / counter - mSigma * mSigma / (counter * counter));
		//std::cout << m2S;
		dis_file << square.dis() << "\t" << xS << "\t" << xSigma << std::endl;
	} // end reps loop
	dis_file.close();
	avergeTimeFile.close();
	std::cout << std::endl;
	return 0;
}
