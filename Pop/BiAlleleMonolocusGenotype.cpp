/*
 * File BiAlleleMonolocusGenotype.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday July 05 2004
 */

#include "BiAlleleMonolocusGenotype.h"

//** Class constructor: *******************************************************/
BiAlleleMonolocusGenotype::BiAlleleMonolocusGenotype(unsigned int first_allele_index,
		unsigned int second_allele_index) {
	_allele_index.push_back(first_allele_index);
	_allele_index.push_back(second_allele_index);
}

BiAlleleMonolocusGenotype::BiAlleleMonolocusGenotype(vector<unsigned int> allele_index) throw (BadIntegerException) {
	if (allele_index.size() != 2)
		throw BadIntegerException("BiAlleleMonolocusGenotype::BiAlleleMonolocusGenotype: allele_index must containes two values.", allele_index.size());
	_allele_index.push_back(allele_index[0]);
	_allele_index.push_back(allele_index[1]);
}

BiAlleleMonolocusGenotype::BiAlleleMonolocusGenotype(const BiAlleleMonolocusGenotype & bmg) {
	for (unsigned int i=0 ; i<2 ; i++)
		_allele_index.push_back(bmg.getAlleleIndex()[i]);
}

//** Class destructor: ********************************************************/
BiAlleleMonolocusGenotype::~BiAlleleMonolocusGenotype() {
	_allele_index.clear();
}

//** Other methodes: **********************************************************/
BiAlleleMonolocusGenotype & BiAlleleMonolocusGenotype::operator= (const BiAlleleMonolocusGenotype & bmg) {
	for (unsigned int i=0 ; i<2 ; i++)
		_allele_index.push_back(bmg.getAlleleIndex()[i]);
	return * this;
}

bool BiAlleleMonolocusGenotype::operator== (const BiAlleleMonolocusGenotype & bmg) const {
	return ((_allele_index[0] == bmg.getAlleleIndex()[0] && _allele_index[1] == bmg.getAlleleIndex()[1])
			|| (_allele_index[0] == bmg.getAlleleIndex()[1] && _allele_index[1] == bmg.getAlleleIndex()[0]));
}

unsigned int BiAlleleMonolocusGenotype::getFirstAlleleIndex() const {
	return _allele_index[0];
}

unsigned int BiAlleleMonolocusGenotype::getSecondAlleleIndex() const {
	return _allele_index[1];
}

bool BiAlleleMonolocusGenotype::isHomozygous() const {
	return (_allele_index[0] == _allele_index[1]);
}

vector<unsigned int> BiAlleleMonolocusGenotype::getAlleleIndex() const {
	return _allele_index;
}

Clonable * BiAlleleMonolocusGenotype::clone() const {
	return new BiAlleleMonolocusGenotype(* this);
}
