/*
 * File MultilocusGenotype.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday July 05 2004
 */

#include "MultilocusGenotype.h"

//** Class constructor: *******************************************************/
MultilocusGenotype::MultilocusGenotype(unsigned int loci_number) throw (BadIntegerException) {
	if (loci_number < 1)
		throw BadIntegerException("MultilocusGenotype::MultilocusGenotype: loci_number must be > 0.", loci_number);

	// Set the _loci size to the right number of loci
	_loci.resize(loci_number);

	// Set all the _loci pointers to NULL
	for (unsigned int i=0 ; i<loci_number ; i++)
		_loci[i] = NULL;
}

MultilocusGenotype::MultilocusGenotype(const MultilocusGenotype & genotype) {
	for (unsigned int i=0 ; i<genotype._loci.size() ; i++)
		_loci[i] = dynamic_cast<MonolocusGenotype *>(genotype._loci[i]->clone());
}

//** Class destructor: *******************************************************/
MultilocusGenotype::~MultilocusGenotype() {
	for (unsigned int i=0 ; i<_loci.size() ; i++)
		delete _loci[i];
}

//** Other methodes: *********************************************************/
void MultilocusGenotype::setMonolocusGenotype(unsigned int locus_index,
		const MonolocusGenotype & monogen) throw (IndexOutOfBoundsException) {
	if (locus_index < _loci.size())
		_loci[locus_index] = dynamic_cast<MonolocusGenotype *>(monogen.clone());
	else
		throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotype: locus_index out of bounds.",
				locus_index, 0, _loci.size());
}

void MultilocusGenotype::setMonolocusGenotypeByAlleleKey(unsigned int locus_index,
		const vector<unsigned int> allele_keys) throw (Exception) {
	if (allele_keys.size() < 1)
		throw Exception("MultilocusGenotype::setMonolocusGenotypeByAlleleKey: no key in allele_keys.");
	
	if (locus_index < _loci.size()) {
		if (allele_keys.size() == 1)
			setMonolocusGenotype(locus_index, MonoAlleleMonolocusGenotype(allele_keys));
		if (allele_keys.size() > 1)
			setMonolocusGenotype(locus_index, BiAlleleMonolocusGenotype(allele_keys));
	}
	else
		throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotype: locus_index out of bounds.",
				locus_index, 0, _loci.size());
}

void MultilocusGenotype::setMonolocusGenotypeByAlleleId(unsigned int locus_index,
		const vector<string> allele_id, const LocusInfo & locus_info) throw (Exception) {
	vector<unsigned int> allele_keys;
	for (unsigned int i = 0 ; i < allele_id.size() ; i++) {
		try {
			allele_keys.push_back(locus_info.getAlleleInfoKey(allele_id[i]));
		}
		catch (AlleleNotFoundException & anfe) {
			throw AlleleNotFoundException("MultilocusGenotype::setMonolocusGenotypeByAlleleId: id not found.", anfe.getIdentifier());
		}
	}
	try {
		setMonolocusGenotypeByAlleleKey(locus_index, allele_keys);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotypeByAlleleId: locus_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

void MultilocusGenotype::setMonolocusGenotypeAsMissing(unsigned int locus_index) throw (IndexOutOfBoundsException) {
	if (locus_index >= _loci.size())
		throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotypeAsMissing: locus_index out of bounds.", locus_index, 0, _loci.size());
	if (_loci[locus_index] != NULL)
		delete _loci[locus_index];
	_loci[locus_index] = NULL;
}

bool MultilocusGenotype::isMonolocusGenotypeMissing(unsigned int locus_index) throw (IndexOutOfBoundsException) {
	if (locus_index >= _loci.size())
		throw IndexOutOfBoundsException("MultilocusGenotype::isMonolocusGenotypeMissing: locus_index out of bounds.", locus_index, 0, _loci.size());
	return _loci[locus_index] == NULL;
}

const MonolocusGenotype * MultilocusGenotype::getMonolocusGenotype(unsigned int locus_index) const throw (IndexOutOfBoundsException) {
	if (locus_index >= _loci.size())
		throw IndexOutOfBoundsException("MultilocusGenotype::getMonolocusGenotype: locus_index out of bounds", locus_index, 0, _loci.size());
	return _loci[locus_index];
}
