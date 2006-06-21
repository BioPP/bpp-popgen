/*
 * File Genetix.cpp
 * Authors : Sylvain Gaillard <yragael2001@yahoo.fr>
 *           Khalid Belkhir
 * Last modification : Monday August 02 2004
 *
*/
/*
Copyright or © or Copr. CNRS, (November 17, 2004)


This software is a computer program whose purpose is to provide classes
for population genetics analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/
#include "Genetix.h"

Genetix::Genetix() {}

Genetix::~Genetix() {}

const string Genetix::getFormatName() {
	return "Genetix ver 4.05";
}

const string Genetix::getFormatDescription() {
	return "Genetix is a software for populations genetic for Windows(tm)";
}

void Genetix::read(istream & is, DataSet & data_set) throw (Exception) {
	if (!is)
		throw IOException("Genetix::read: fail to open stream.");
	// Loci number
	string temp = FileTools::getNextLine(is);
	unsigned int loc_nbr;
	stringstream(temp) >> loc_nbr;
	data_set.initAnalyzedLoci(loc_nbr);

	// Groups number
	temp = FileTools::getNextLine(is);
	unsigned int grp_nbr;
	stringstream(temp) >> grp_nbr;

	// Loci data
	for (unsigned int i = 0 ; i < loc_nbr ; i++) {
		// Locus name
		string name = FileTools::getNextLine(is);
		name = TextTools::removeSurroundingWhiteSpaces(name);
		LocusInfo tmp_loc(name);
		// Alleles
		stringstream values(FileTools::getNextLine(is));
		unsigned int nbr_al;
		values >> nbr_al;
		for (unsigned int j = 0 ; j < nbr_al ; j++) {
			string al_id;
			values >> al_id;
			BasicAlleleInfo tmp_al(al_id);
			tmp_loc.addAlleleInfo(tmp_al);
		}
		data_set.setLocusInfo(i, tmp_loc);
	}
	
	// Groups
	for (unsigned int i = 0 ; i < grp_nbr ; i++) {
		data_set.addEmptyGroup(i);
		// Group name ... Now used khalid
		temp = FileTools::getNextLine(is);
        data_set.setGroupName(i, temp);
        
		// Number of individuals
		unsigned int ind_nbr;
		temp = FileTools::getNextLine(is);
		stringstream tmp(temp);
		tmp >> ind_nbr;
		for (unsigned int j = 0 ; j < ind_nbr ; j++) {
			temp = FileTools::getNextLine(is);
			string ind_name(temp.begin(), temp.begin()+11);
			temp = string(temp.begin()+11, temp.end());
			data_set.addEmptyIndividualToGroup(i, TextTools::removeSurroundingWhiteSpaces(ind_name) + string("_") + TextTools::toString(i+1) + string("_") + TextTools::toString(j+1));
			data_set.initIndividualGenotypeInGroup(i, j);
			StringTokenizer alleles(temp, string(" "));
			//cout << alleles.numberOfRemainingTokens() << endl;
			for (unsigned int k = 0 ; k < loc_nbr ; k++) {
				string tmp_string = alleles.nextToken();
				vector<string> tmp_alleles;
				tmp_alleles.push_back(string(tmp_string.begin(), tmp_string.begin() + 3));
				tmp_alleles.push_back(string(tmp_string.begin() + 3, tmp_string.begin() + 6));
				if (tmp_alleles[0] != string("000") && tmp_alleles[1] != string("000"))
					data_set.setIndividualMonolocusGenotypeByAlleleIdInGroup(i, j, k, tmp_alleles);
			}
		}
	}
}

void Genetix::read(const string & path, DataSet & data_set) throw (Exception) {
	AbstractIDataSet::read(path, data_set);
}
	
DataSet * Genetix::read(istream & is) throw (Exception) {
	return AbstractIDataSet::read(is);
}

DataSet * Genetix::read(const string & path) throw (Exception) {
	return AbstractIDataSet::read(path);
}
