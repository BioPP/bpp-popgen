/*
 * File DiploidMonolocusGenotype.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 22 2004
 */

#include "DiploidMonolocusGenotype.h"

//** Class constructor: *******************************************************/
DiploidMonolocusGenotype::DiploidMonolocusGenotype(unsigned int first_allele_key,
		unsigned int second_allele_key) {
	_allele_keys.push_back(first_allele_key);
	_allele_keys.push_back(second_allele_key);
}

DiploidMonolocusGenotype::DiploidMonolocusGenotype(const DiploidMonolocusGenotype & dmg) {
	for (unsigned int i=0 ; i<2 ; i++)
		_allele_keys.push_back(dmg.getAlleleKey(i));
}

//** Class destructor: ********************************************************/
DiploidMonolocusGenotype::~DiploidMonolocusGenotype() {}

//** Other methodes: **********************************************************/
unsigned int DiploidMonolocusGenotype::getAlleleKey(unsigned int allele_key) const
throw (IndexOutOfBoundsException) {
	if (allele_key == 0 || allele_key == 1)
		return _allele_keys[allele_key];
	else
		throw IndexOutOfBoundsException("DiploidMonolocusGenotype::getAlleleKey: allele_key out of bounds.",
				allele_key, 0, 1);
}

bool DiploidMonolocusGenotype::isHomozygous() const {
	if (_allele_keys[0] == _allele_keys[1])
		return true;
	else
		return false;
}

unsigned int DiploidMonolocusGenotype::getAlleleKey() const {
	return _allele_keys[0];
}

unsigned int DiploidMonolocusGenotype::getPloidy() const {
	return LocusInfo::DIPLOID;
}

unsigned int DiploidMonolocusGenotype::getSize() const {
	return 2;
}

Clonable * DiploidMonolocusGenotype::clone() const {
	return new DiploidMonolocusGenotype(* this);
}
