/*
 * File Genotype.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 07 2004
 */

#include "Genotype.h"

//** Class constructor: *******************************************************/
Genotype::Genotype(AnalyzedLoci * analyzed_loci) {
	_analyzedLoci = analyzed_loci;

	// Set the _loci size to the right number of loci
	_loci.resize(analyzed_loci->getNumberOfLoci());

	// Set all the _loci pointers to NULL
	for (unsigned int i=0 ; i<_loci.size() ; i++)
		_loci[i] = NULL;
}

//** Class destructor: *******************************************************/
Genotype::~Genotype() {
	for (unsigned int i=0 ; i<_loci.size() ; i++)
		delete _loci[i];
}

//** Other methodes: *********************************************************/
void Genotype::setMonolocusGenotype(unsigned int locus_index,
		const MonolocusGenotype & monogen) throw (IndexOutOfBoundsException) {
	if (locus_index < _loci.size())
		_loci[locus_index] = dynamic_cast<MonolocusGenotype *>(monogen.clone());
	else
		throw IndexOutOfBoundsException("Genotype::setMonolocusGenotype: locus_index out of bounds.",
				locus_index, 0, _loci.size());
}

void Genotype::setMonolocusGenotypeByAlleleKey(unsigned int locus_index,
		const vector<unsigned int> allele_keys) throw (Exception) {
	if (locus_index < _loci.size()) {
		if (_analyzedLoci->getPloidyByLocusIndex(locus_index) == LocusInfo::HAPLODIPLOID)
			setMonolocusGenotype(locus_index, HaplodiploidMonolocusGenotype(allele_keys[0]));
		if (_analyzedLoci->getPloidyByLocusIndex(locus_index) == LocusInfo::HAPLOID)
			setMonolocusGenotype(locus_index, HaploidMonolocusGenotype(allele_keys[0]));
		if (_analyzedLoci->getPloidyByLocusIndex(locus_index) == LocusInfo::DIPLOID)
			setMonolocusGenotype(locus_index, DiploidMonolocusGenotype(allele_keys[0], allele_keys[1]));
	}
	else
		throw IndexOutOfBoundsException("Genotype::setMonolocusGenotype: locus_index out of bounds.",
				locus_index, 0, _loci.size());
}

void Genotype::setMonolocusGenotypeByAlleleId(unsigned int locus_index,
		const vector<unsigned int> allele_id) throw (Exception) {
	vector<unsigned int> allele_keys;
	for (unsigned int i = 0 ; i < allele_id.size() ; i++)
		allele_keys.push_back(_analyzedLoci->getLocusInfoByIndex(locus_index)->getAlleleInfoKey(allele_id[i]));
	setMonolocusGenotypeByAlleleKey(locus_index, allele_keys);
}
