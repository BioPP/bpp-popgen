/*
 * File Genetix.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday August 02 2004
 *
 * Copyright (C) 2004 Sylvain Gaillard and the
 *                    PopGenLib Development Core Team
 *
 * This file is part of PopGenLib.
 *
 * PopGenLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopGenLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopGenLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
		// Group name ... not used
		temp = FileTools::getNextLine(is);
		data_set.addEmptyGroup(i);
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
