/*
 * File MonoAlleleMonolocusGenotype.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday July 05 2004
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
