/*
 * File MultilocusGenotype.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
 *
 * Copyright (C) 2004 Sylvain Gaillard and the
 *                    PopGenLib Development Core Team
 *
 * This file is part of PopGenLib.
 *
 * PopGenLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopGenLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopGenLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
	for (unsigned int i=0 ; i<genotype.size() ; i++) {
		if (! genotype.isMonolocusGenotypeMissing(i))
			_loci.push_back(dynamic_cast<MonolocusGenotype *>((genotype.getMonolocusGenotype(i))->clone()));
		else
			_loci.push_back(NULL);
	}
}

//** Class destructor: *******************************************************/
MultilocusGenotype::~MultilocusGenotype() {
	for (unsigned int i=0 ; i<_loci.size() ; i++)
		delete _loci[i];
	_loci.clear();
}

//** Other methodes: *********************************************************/
void MultilocusGenotype::setMonolocusGenotype(unsigned int locus_position,
		const MonolocusGenotype & monogen) throw (IndexOutOfBoundsException) {
	if (locus_position < _loci.size())
		_loci[locus_position] = dynamic_cast<MonolocusGenotype *>(monogen.clone());
	else
		throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotype: locus_position out of bounds.",
				locus_position, 0, _loci.size());
}

void MultilocusGenotype::setMonolocusGenotypeByAlleleKey(unsigned int locus_position,
		const vector<unsigned int> allele_keys) throw (Exception) {
	if (allele_keys.size() < 1)
		throw Exception("MultilocusGenotype::setMonolocusGenotypeByAlleleKey: no key in allele_keys.");
	
	if (locus_position < _loci.size()) {
		if (allele_keys.size() == 1)
			setMonolocusGenotype(locus_position, MonoAlleleMonolocusGenotype(allele_keys));
		if (allele_keys.size() > 1)
			setMonolocusGenotype(locus_position, BiAlleleMonolocusGenotype(allele_keys));
	}
	else
		throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotype: locus_position out of bounds.",
				locus_position, 0, _loci.size());
}

void MultilocusGenotype::setMonolocusGenotypeByAlleleId(unsigned int locus_position,
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
		setMonolocusGenotypeByAlleleKey(locus_position, allele_keys);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotypeByAlleleId: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

void MultilocusGenotype::setMonolocusGenotypeAsMissing(unsigned int locus_position) throw (IndexOutOfBoundsException) {
	if (locus_position >= _loci.size())
		throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotypeAsMissing: locus_position out of bounds.", locus_position, 0, _loci.size());
	if (_loci[locus_position] != NULL)
		delete _loci[locus_position];
	_loci[locus_position] = NULL;
}

bool MultilocusGenotype::isMonolocusGenotypeMissing(unsigned int locus_position) const throw (IndexOutOfBoundsException) {
	if (locus_position >= _loci.size())
		throw IndexOutOfBoundsException("MultilocusGenotype::isMonolocusGenotypeMissing: locus_position out of bounds.", locus_position, 0, _loci.size());
	return _loci[locus_position] == NULL;
}

const MonolocusGenotype * MultilocusGenotype::getMonolocusGenotype(unsigned int locus_position) const throw (IndexOutOfBoundsException) {
	if (locus_position >= _loci.size())
		throw IndexOutOfBoundsException("MultilocusGenotype::getMonolocusGenotype: locus_position out of bounds", locus_position, 0, _loci.size());
	return _loci[locus_position];
}

unsigned int MultilocusGenotype::size() const {
	return _loci.size();
}

unsigned int MultilocusGenotype::countNonMissingLoci() const {
	unsigned int count = 0;
	for (unsigned int i = 0 ; i < _loci.size() ; i++)
		if (_loci[i] != NULL) count++;
	return count;
}

unsigned int MultilocusGenotype::countHomozygousLoci() const {
	unsigned int count = 0;
	for (unsigned int i = 0 ; i < _loci.size() ; i++) {
		try {
			if (dynamic_cast<BiAlleleMonolocusGenotype *>(_loci[i])->isHomozygous()) count++;
		}
		catch (...) {}
	}
	return count;
}

unsigned int MultilocusGenotype::countHeterozygousLoci() const {
	unsigned int count = 0;
	for (unsigned int i = 0 ; i < _loci.size() ; i++) {
		try {
			if (!(dynamic_cast<BiAlleleMonolocusGenotype *>(_loci[i])->isHomozygous())) count++;
		}
		catch (...) {}
	}
	return count;
}
