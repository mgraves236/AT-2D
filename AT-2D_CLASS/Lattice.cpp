//
// Created by Magda on 25.04.2021.
//

#include "Lattice.h"

// Mersenne Twister
auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937 mt_rand(seed);
std::uniform_real_distribution<double> doubleDist(0.0, 1.0);


Lattice::Lattice(int n, double J_1, double J_2, double K_1,
				 double K_2, double M_0, double R_1, double R_2, double T) {
	this->n = n;
	lattice = new Node*[n];
	for (int i = 0; i < n; i++) {
		lattice[i] = new Node[n];
	}
	this->J_1 = J_1;
	this->J_2 = J_2;
	this->K_1 = K_1;
	this->K_2 = K_2;
	this->R_1 = R_1;
	this->R_2 = R_2;
	this->M_0 = M_0;
	this->T = T;
	// next iN previous iP
	iN1 = new int[n];
	iN2 = new int[n];
	iP1 = new int[n];
	iP2 = new int[n];
	for (int i = 0; i < n; i++) {
		iP1[i] = i - 1;
		iN1[i] = i + 1;
		iP2[i] = i - 2;
		iN2[i] = i + 2;
	}
	// boundary conditions
	iP1[0] = n - 1;
	iP2[0] = n - 2;
	iP2[1] = n - 1;
	iN1[n - 1] = 0;
	iN2[n - 1] = 1;
	iN2[n - 2] = 0;
}

Lattice::~Lattice() {
	for(int i = 0; i < n; i++) {
		delete lattice[i];
	}
	delete lattice;
	delete iN1;
	delete iN2;
	delete iP1;
	delete iP2;
}

void Lattice::printS() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << lattice[i][j].S << " ";
		}
		std::cout << std::endl;
	}
}

void Lattice::printSigma() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << lattice[i][j].sigma << " ";
		}
		std::cout << std::endl;
	}
}

void Lattice::printBoth() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (lattice[i][j].S == 1 && lattice[i][j].sigma == 1) {
				std::cout << "1" << " ";
			} else if (lattice[i][j].S == 1 && lattice[i][j].sigma == -1) {
				std::cout << "2" << " ";
			} else if (lattice[i][j].S == -1 && lattice[i][j].sigma == 1) {
				std::cout << "3" << " ";
			} else if (lattice[i][j].S == -1 && lattice[i][j].sigma == -1) {
				std::cout << "4" << " ";
			}
		}
		std::cout << std::endl;
	}
}

void Lattice::init(double d) {
	// chooosing S spins
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double rand = doubleDist(mt_rand);
			if (rand <= .5) {
				lattice[i][j].S = 1;
			} else {
				lattice[i][j].S = 1;
			}
		}
	}
	// choosing sigma
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double rand = doubleDist(mt_rand);
			if (d < rand) {
				lattice[i][j].sigma = lattice[i][j].S;
 			} else {
				lattice[i][j].sigma = - lattice[i][j].S;
			}
		}
	}
}

void Lattice::monteCarloStep() {
	std::uniform_int_distribution<int> intDist(0, n - 1);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			int randRow = intDist(mt_rand);
			int randCol = intDist(mt_rand);
			// change S
			double E0 = getE();
			lattice[randRow][randCol].S = - lattice[randRow][randCol].S;
			double E = getE();
			lattice[randRow][randCol].S = - lattice[randRow][randCol].S;
			if (E - E0 <= 0) {
				lattice[randRow][randCol].S = - lattice[randRow][randCol].S;
			} else if (T != 0.0) {
				double rand = doubleDist(mt_rand);
				if (rand < exp(- (E - E0) / T)) {
					lattice[randRow][randCol].S = - lattice[randRow][randCol].S;
				}
			}
			// change sigma
			E0 = getE();
			lattice[randRow][randCol].sigma = - lattice[randRow][randCol].sigma;
			E = getE();
			lattice[randRow][randCol].sigma = - lattice[randRow][randCol].sigma;
			if (E - E0 <= 0) {
				lattice[randRow][randCol].sigma = - lattice[randRow][randCol].sigma;
			} else if (T != 0.0) {
				double rand = doubleDist(mt_rand);
				if (rand <  exp(- (E - E0) / T)) {
					lattice[randRow][randCol].sigma = - lattice[randRow][randCol].sigma;
				}
			}
		}
	}
}

double Lattice::getE() {
	double sumJ1 = 0;
	double sumJ2 = 0;
	double sumK1 = 0;
	double sumK2 = 0;
	double sumR1 = 0;
	double sumR2 = 0;
	double sumM0 = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			sumJ1 += lattice[i][j].S *
					(lattice[iN1[i]][j].S + lattice[i][iN1[j]].S +
					 lattice[iP1[i]][j].S + lattice[i][iP1[j]].S);
			sumJ2 += lattice[i][j].S *
					 (lattice[iN2[i]][j].S + lattice[i][iN2[j]].S +
					  lattice[iP2[i]][j].S + lattice[i][iP2[j]].S);
			sumK1 += lattice[i][j].sigma *
					 (lattice[iN1[i]][j].sigma + lattice[i][iN1[j]].sigma +
					  lattice[iP1[i]][j].sigma + lattice[i][iP1[j]].sigma);
			sumK2 += lattice[i][j].sigma *
					 (lattice[iN2[i]][j].sigma + lattice[i][iN2[j]].sigma +
					  lattice[iP2[i]][j].sigma + lattice[i][iP2[j]].sigma);
			sumR1 += lattice[i][j].sigma *
					 (lattice[iN1[i]][j].S + lattice[i][iN1[j]].S +
					  lattice[iP1[i]][j].S + lattice[i][iP1[j]].S);
			sumR2 += lattice[i][j].sigma *
					 (lattice[iN2[i]][j].S + lattice[i][iN2[j]].S +
					  lattice[iP2[i]][j].S + lattice[i][iP2[j]].S);
			sumM0 += lattice[i][j].sigma * lattice[i][j].S;
		}
	}
	double E = (-J_1 * sumJ1 - J_2 * sumJ2
				-K_1 * sumK1 - K_2 * sumK2
				-R_1 * sumR1 - R_2 * sumR2) * 0.5
						 -M_0 * sumM0;
	return E;
}

double Lattice::gS() {
	int sum = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			sum += lattice[i][j].S *
				   (lattice[iN1[i]][j].S + lattice[i][iN1[j]].S +
					lattice[iP1[i]][j].S + lattice[i][iP1[j]].S);
		}
	}
	return sum * 0.25 / ((double) (n * n));
}

double Lattice::gSigma() {
	int sum = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			sum += lattice[i][j].sigma *
				   (lattice[iN1[i]][j].sigma + lattice[i][iN1[j]].sigma +
					lattice[iP1[i]][j].sigma + lattice[i][iP1[j]].sigma);
		}
	}
	return sum * 0.25 / ((double) (n * n));
}

double Lattice::dis() {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			sum += 1 - lattice[i][j].S * lattice[i][j].sigma;
		}
	}
	return 1 /((double) 2) * 1 / ((double) n * n) * sum;
}

double Lattice::magnetizationS() {
	double m = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			m += lattice[i][j].S;
		}
	}
	return m / (double) (n * n);
}

double Lattice::magnetizationSigma() {
	double m = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			m += lattice[i][j].sigma;
		}
	}
	return m / (double) (n * n);
}
