/*
 * File HaploidMonolocusGenotype.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 07 2004
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

Clonable * HaploidMonolocusGenotype::clone() const {
	return new HaploidMonolocusGenotype(* this);
}
