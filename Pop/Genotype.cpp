/*
 * File Genotype.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 22 2004
 */

#include "Genotype.h"

//** Class constructor: *******************************************************/
Genotype::Genotype(const AnalyzedLoci * analyzed_loci) throw (NullPointerException) {
	if (analyzed_loci == NULL)
		throw NullPointerException("Genotype::Genotype: analyzed_loci is NULL.");
	_analyzedLoci = analyzed_loci;

	// Set the _loci size to the right number of loci
	_loci.resize(analyzed_loci->getNumberOfLoci());

	// Set all the _loci pointers to NULL
	for (unsigned int i=0 ; i<_loci.size() ; i++)
		_loci[i] = NULL;
}

Genotype::Genotype(const Genotype & genotype) {
	_analyzedLoci = genotype._analyzedLoci;
	for (unsigned int i=0 ; i<genotype._loci.size() ; i++)
		_loci[i] = dynamic_cast<MonolocusGenotype *>(genotype._loci[i]->clone());
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
	if (allele_keys.size() == 0)
		throw Exception("Genotype::setMonolocusGenotypeByAlleleKey: no key in allele_keys.");
	if (allele_keys.size() != 1 && (_analyzedLoci->getPloidyByLocusIndex(locus_index) == LocusInfo::HAPLODIPLOID || _analyzedLoci->getPloidyByLocusIndex(locus_index) == LocusInfo::HAPLOID))
		throw Exception("Genotype::setMonolocusGenotypeByAlleleKey: too many value in allele_keys for a HAPLOID or HAPLODIPLOID locus.");
	if (allele_keys.size() != 2 && _analyzedLoci->getPloidyByLocusIndex(locus_index) == LocusInfo::DIPLOID)
		throw Exception("Genotype::setMonolocusGenotypeByAlleleKey: allele_keys doesn't have the right size for a DIPLOID locus.");
	
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
	if (allele_id.size() == 0)
		throw Exception("Genotype::setMonolocusGenotypeByAlleleId: no id in allele_id.");
	if (allele_id.size() != 1 && (_analyzedLoci->getPloidyByLocusIndex(locus_index) == LocusInfo::HAPLODIPLOID || _analyzedLoci->getPloidyByLocusIndex(locus_index) == LocusInfo::HAPLOID))
		throw Exception("Genotype::setMonolocusGenotypeByAlleleId: too many value in allele_id for a HAPLOID or HAPLODIPLOID locus.");
	if (allele_id.size() != 2 && _analyzedLoci->getPloidyByLocusIndex(locus_index) == LocusInfo::DIPLOID)
		throw Exception("Genotype::setMonolocusGenotypeByAlleleId: allele_id doesn't have the right size for a DIPLOID locus.");
	
	vector<unsigned int> allele_keys;
	for (unsigned int i = 0 ; i < allele_id.size() ; i++)
		allele_keys.push_back(_analyzedLoci->getLocusInfoByIndex(locus_index)->getAlleleInfoKey(allele_id[i]));
	try {
		setMonolocusGenotypeByAlleleKey(locus_index, allele_keys);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Genotype::setMonolocusGenotypeByAlleleId: locus_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (Exception & e) {
		throw Exception("Genotype::setMonolocusGenotypeByAlleleId: allele_keys.size() doesn't match ploidy.");
	}
}

unsigned int Genotype::getPloidy(unsigned int locus_index) const throw (IndexOutOfBoundsException) {
	if (locus_index >= _loci.size())
		throw IndexOutOfBoundsException("Genotype::getPloidy: locus_index out of bounds", locus_index, 0, _loci.size());
	return _loci[locus_index]->getPloidy();
}

const MonolocusGenotype * Genotype::getMonolocusGenotype(unsigned int locus_index) const throw (IndexOutOfBoundsException) {
	if (locus_index >= _loci.size())
		throw IndexOutOfBoundsException("Genotype::getMonolocusGenotype: locus_index out of bounds", locus_index, 0, _loci.size());
	return _loci[locus_index];
}

const AnalyzedLoci * Genotype::getAnalyzedLoci() {
	return _analyzedLoci;
}
