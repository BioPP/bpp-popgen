/*
 * File AnalyzedLoci.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday July 21 2004
 */

#include "AnalyzedLoci.h"

//** Constructors: ***********************************************************/
AnalyzedLoci::AnalyzedLoci(unsigned int number_of_loci) {
	_loci.resize(number_of_loci);
	for (unsigned int i = 0 ; i < _loci.size() ; i++)
		_loci[i] = NULL;
}

AnalyzedLoci::AnalyzedLoci(const AnalyzedLoci & analyzed_loci) {
	for (unsigned int i = 0 ; i < analyzed_loci.getNumberOfLoci() ; i++)
		_loci.push_back(new LocusInfo(* analyzed_loci.getLocusInfoAtPosition(i)));
}

//** Destructor: *************************************************************/
AnalyzedLoci::~AnalyzedLoci() {
	for (unsigned int i = 0 ; i < _loci.size() ; i++)
		delete _loci[i];
}

//** Other methodes: *********************************************************/
// LocusInfo
void AnalyzedLoci::setLocusInfo(unsigned int locus_position, const LocusInfo & locus)
throw (IndexOutOfBoundsException) {
	if (locus_position >= 0 && locus_position < _loci.size())
		_loci[locus_position] = new LocusInfo(locus);
	else
		throw IndexOutOfBoundsException("AnalyzedLoci::setLocusInfo: locus_position out of bounds",
				locus_position, 0, _loci.size());
}

unsigned int AnalyzedLoci::getLocusInfoPosition(const string & locus_name) const
throw (BadIdentifierException) {
	for (unsigned int i = 0 ; i < _loci.size() ; i++)
		if (_loci[i] != NULL && _loci[i]->getName() == locus_name)
			return i;
	throw BadIdentifierException("AnalyzedLoci::getLocusInfoPosition: locus not found.", locus_name);
}

const LocusInfo * AnalyzedLoci::getLocusInfoByName(const string & locus_name) const
throw (BadIdentifierException) {
	for (unsigned int i = 0 ; i < _loci.size() ; i++)
		if (_loci[i] != NULL && _loci[i]->getName() == locus_name)
			return _loci[i];
	throw BadIdentifierException("AnalyzedLoci::getLocusInfo: locus not found.",
			locus_name);
}

const LocusInfo * AnalyzedLoci::getLocusInfoAtPosition(unsigned int locus_position) const
throw (Exception) {
	if (locus_position >= _loci.size())
		throw IndexOutOfBoundsException("AnalyzedLoci::getLocusInfoAtPosition: locus_position out of bounds.", locus_position, 0, _loci.size());
	if (_loci[locus_position] != NULL)
		return _loci[locus_position];
	else
		throw NullPointerException("AnalyzedLoci::getLocusInfo: no locus defined here.");
}

// AlleleInfo
void AnalyzedLoci::addAlleleInfoByLocusName(const string & locus_name,
		const AlleleInfo& allele)
throw (Exception) {
	bool locus_found = false;
	for (vector<LocusInfo *>::iterator it = _loci.begin() ;
			it != _loci.end() ; it++) {
		if ((*it)->getName() == locus_name) {
			locus_found = true;
			try {
				(*it)->addAlleleInfo(allele);
			}
			catch (BadIdentifierException & bie) {
				throw BadIdentifierException("AnalyzedLoci::addAlleleInfoByLocusName: allele id already in use.", bie.getIdentifier());
			}
		}
	}
	if (!locus_found)
		throw LocusNotFoundException("AnalyzedLoci::addAlleleInfoByLocusName: locus_name not found.",
				locus_name);
}

void AnalyzedLoci::addAlleleInfoByLocusPosition(unsigned int locus_position,
		const AlleleInfo & allele)
throw (Exception) {
	if (locus_position >= 0 && locus_position < _loci.size()) {
		try {
			_loci[locus_position]->addAlleleInfo(allele);
		}
		catch (BadIdentifierException & bie) {
			throw BadIdentifierException("AnalyzedLoci::addAlleleInfoByLocusPosition: allele id is already in use.", bie.getIdentifier());
		}
	}
	else
		throw IndexOutOfBoundsException("AnalyzedLoci::addAlleleInfoByLocusPosition: locus_position out of bounds.",
				locus_position, 0, _loci.size());
}

// General
unsigned int AnalyzedLoci::getNumberOfLoci() const {
	return _loci.size();
}

vector<unsigned int> AnalyzedLoci::getNumberOfAlleles() const {
	vector<unsigned int> allele_count;
	for (unsigned int i = 0 ; i < _loci.size() ; i++)
		allele_count.push_back(_loci[i]->getNumberOfAlleles());
	return allele_count;
}

unsigned int AnalyzedLoci::getPloidyByLocusName(const string & locus_name) const
throw (LocusNotFoundException) {
	for (unsigned int i = 0 ; i < _loci.size() ; i++)
		if (_loci[i] != NULL && _loci[i]->getName() == locus_name)
			return _loci[i]->getPloidy();
	throw LocusNotFoundException("AnalyzedLoci::getLocusInfo: locus_name not found.",
			locus_name);
}

unsigned int AnalyzedLoci::getPloidyByLocusPosition(unsigned int locus_position) const
throw (IndexOutOfBoundsException) {
	if (locus_position >= _loci.size())
		throw IndexOutOfBoundsException("AnalyzedLoci::getPloidyByLocusPosition: locus_position out of bounds.", locus_position, 0, _loci.size());
	return _loci[locus_position]->getPloidy();
}
