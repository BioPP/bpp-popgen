/*
 * File AnalyzedLoci.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 01 2004
 */

#include "AnalyzedLoci.h"

//** Constructors: ***********************************************************/
AnalyzedLoci::AnalyzedLoci(unsigned int number_of_loci) {
	_loci.resize(number_of_loci);
	for (unsigned int i = 0 ; i < _loci.size() ; i++)
		_loci[i] = NULL;
}

//** Destructor: *************************************************************/
AnalyzedLoci::~AnalyzedLoci() {
	for (unsigned int i = 0 ; i < _loci.size() ; i++)
		delete _loci[i];
}

//** Other methodes: *********************************************************/
// LocusInfo
void AnalyzedLoci::setLocusInfo(unsigned int locus_index, const LocusInfo & locus)
throw (IndexOutOfBoundsException) {
	if (locus_index >= 0 && locus_index < _loci.size())
		_loci[locus_index] = new LocusInfo(locus);
	else
		throw IndexOutOfBoundsException("AnalyzedLoci::setLocusInfo: locus_index out of bounds",
				locus_index, 0, _loci.size());
}

const LocusInfo * AnalyzedLoci::getLocusInfo(const string & locus_name)
throw (LocusInfoNotFoundException) {
	for (unsigned int i = 0 ; i < _loci.size() ; i++)
		if (_loci[i] != NULL && _loci[i]->getName() == locus_name)
			return _loci[i];
	throw LocusInfoNotFoundException("AnalyzedLoci::getLocusInfo: locus not found.",
			locus_name);
}

const LocusInfo * AnalyzedLoci::getLocusInfo(unsigned int locus_index)
throw (NullPointerException) {
	if (_loci[locus_index] != NULL)
		return _loci[locus_index];
	else
		throw NullPointerException("AnalyzedLoci::getLocusInfo: no locus defined here.");
}

// AlleleInfo
void AnalyzedLoci::addAlleleInfo(const string & locus_name, const AlleleInfo& allele)
throw (Exception) {
	bool locus_found = false;
	for (vector<LocusInfo *>::iterator it = _loci.begin() ; it != _loci.end() ; it++) {
		if ((*it)->getName() == locus_name) {
			locus_found = true;
			try {
				(*it)->addAlleleInfo(allele);
			}
			catch (Exception e) {
				throw e;
			}
		}
	}
	if (!locus_found)
		throw LocusInfoNotFoundException("AnalyzedLoci::addAlleleInfo: locus not found.",
				locus_name);
}

void AnalyzedLoci::addAlleleInfo(unsigned int locus_index, const AlleleInfo & allele)
throw (Exception) {
	if (locus_index >= 0 && locus_index < _loci.size()) {
		try {
			_loci[locus_index]->addAlleleInfo(allele);
		}
		catch (Exception e) {
			throw e;
		}
	}
	else
		throw IndexOutOfBoundsException("AnalyzedLoci::addAlleleInfo: locus_index out of bounds.",
				locus_index, 0, _loci.size());
}
