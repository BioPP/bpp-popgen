/*
 * File HaploidMonolocusGenotype.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 22 2004
 */

#include "HaploidMonolocusGenotype.h"

//** Class constructor: *******************************************************/
HaploidMonolocusGenotype::HaploidMonolocusGenotype(unsigned int allele_key) {
	_allele_key = allele_key;
}

HaploidMonolocusGenotype::HaploidMonolocusGenotype(const HaploidMonolocusGenotype & hmg) {
	_allele_key = hmg.getAlleleKey();
}

//** Class destructor: ********************************************************/
HaploidMonolocusGenotype::~HaploidMonolocusGenotype() {}

//** Other methodes: **********************************************************/
unsigned int HaploidMonolocusGenotype::getAlleleKey() const {
	return _allele_key;
}

unsigned int HaploidMonolocusGenotype::getPloidy() const {
	return LocusInfo::HAPLOID;
}

unsigned int HaploidMonolocusGenotype::getSize() const {
	return 1;
}

Clonable * HaploidMonolocusGenotype::clone() const {
	return new HaploidMonolocusGenotype(* this);
}
