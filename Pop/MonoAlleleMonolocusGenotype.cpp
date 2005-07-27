/*
 * File MonoAlleleMonolocusGenotype.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
 *
*/
/*
Copyright or © or Copr. CNRS, (November 17, 2004)


This software is a computer program whose purpose is to provide classes
for sequences analysis.

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
#include "MonoAlleleMonolocusGenotype.h"

//** Class constructor: *******************************************************/
MonoAlleleMonolocusGenotype::MonoAlleleMonolocusGenotype(unsigned int allele_index) {
	_allele_index = allele_index;
}

MonoAlleleMonolocusGenotype::MonoAlleleMonolocusGenotype(vector<unsigned int> allele_index) throw (BadIntegerException) {
	if (allele_index.size() != 1)
		throw BadIntegerException("MonoAlleleMonolocusGenotype::MonoAlleleMonolocusGenotype: allele_index must conaines one value.", allele_index.size());
	_allele_index = allele_index[0];
}

MonoAlleleMonolocusGenotype::MonoAlleleMonolocusGenotype(const MonoAlleleMonolocusGenotype & mmg) {
	_allele_index = mmg.getAlleleIndex()[0];
}

//** Class destructor: ********************************************************/
MonoAlleleMonolocusGenotype::~MonoAlleleMonolocusGenotype() {}

//** Other methodes: **********************************************************/
MonoAlleleMonolocusGenotype & MonoAlleleMonolocusGenotype::operator= (const MonoAlleleMonolocusGenotype & mmg) {
	_allele_index = mmg.getAlleleIndex()[0];
	return * this;
}

bool MonoAlleleMonolocusGenotype::operator== (const MonoAlleleMonolocusGenotype & mmg) const {
	return (_allele_index == mmg.getAlleleIndex()[0]);
}

vector<unsigned int> MonoAlleleMonolocusGenotype::getAlleleIndex() const {
	vector<unsigned int> index;
	index.push_back(_allele_index);
	return index;
}

Clonable * MonoAlleleMonolocusGenotype::clone() const {
	return new MonoAlleleMonolocusGenotype(* this);
}
