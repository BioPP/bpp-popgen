/*
 * File DiploidMonolocusGenotype.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 07 2004
 */

#include "DiploidMonolocusGenotype.h"

//** Class constructor: *******************************************************/
DiploidMonolocusGenotype::DiploidMonolocusGenotype(unsigned int first_allele_key,
		unsigned int second_allele_key) {
	_allelekeys.push_back(first_allele_key);
	_allelekeys.push_back(second_allele_key);
}

DiploidMonolocusGenotype::DiploidMonolocusGenotype(const DiploidMonolocusGenotype & dmg) {
	for (unsigned int i=0 ; i<2 ; i++)
		_allelekeys.push_back(dmg.getAlleleKey(i));
}

//** Class destructor: ********************************************************/
DiploidMonolocusGenotype::~DiploidMonolocusGenotype() {}

//** Other methodes: **********************************************************/
unsigned int DiploidMonolocusGenotype::getAlleleKey(unsigned int allele_key) const
throw (IndexOutOfBoundsException) {
	if (allele_key == 0 || allele_key == 1)
		return _allelekeys[allele_key];
	else
		throw IndexOutOfBoundsException("DiploidMonolocusGenotype::getAlleleKey: allele_key out of bounds.",
				allele_key, 0, 1);
}

bool DiploidMonolocusGenotype::isHomozygous() const {
	if (_allelekeys[0] == _allelekeys[1])
		return true;
	else
		return false;
}

unsigned int DiploidMonolocusGenotype::getAlleleKey() const {
	return _allelekeys[0];
}

unsigned int DiploidMonolocusGenotype::getPloidy() const {
	return LocusInfo::DIPLOID;
}

Clonable * DiploidMonolocusGenotype::clone() const {
	return new DiploidMonolocusGenotype(* this);
}
