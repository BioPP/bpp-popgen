/*
 * File AbstractODataSet.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 24 2004
 */

#include "AbstractODataSet.h"

// From STL
#include <fstream>
using namespace std;

AbstractODataSet::~AbstractODataSet() {}

void AbstractODataSet::write(const string & path, const DataSet & data_set, bool overwrite) const throw (Exception) {
	ofstream output(path.c_str(), overwrite ? (ios::out) : (ios::out|ios::app));
	write(output, data_set);
	output.close();
}
