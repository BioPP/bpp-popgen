/*
 * File AbstractIDataSet.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 24 2004
 */

#include "AbstractIDataSet.h"

// From STL
#include <fstream>
using namespace std;

AbstractIDataSet::~AbstractIDataSet() {}

void AbstractIDataSet::read(const string & path, DataSet & data_set) throw (Exception) {
	ifstream input(path.c_str(), ios::in);
	read(input, data_set);
	input.close();
}

DataSet * AbstractIDataSet::read(istream & is) throw (Exception) {
	DataSet * data_set = new DataSet();
	read(is, * data_set);
	return data_set;
}

DataSet * AbstractIDataSet::read(const string & path) throw (Exception) {
	DataSet * data_set = new DataSet();
	read(path, * data_set);
	return data_set;
}
