//
// Created by Magda on 28.04.2021.
//

#include "func.h"

std::string str_i(int x) {
	std::stringstream out;
	out << x;
	return out.str();
}

std::string DoubleToString(double a) {
	std::ostringstream temp;
	temp << a;
	return temp.str();
}
