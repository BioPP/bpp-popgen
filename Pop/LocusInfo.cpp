/*
 * File LocusInfo.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Saturday July 03 2004
 */

// From Utils
#include <Utils/TextTools.h>

#include "LocusInfo.h"
#include "GeneralExceptions.h"

unsigned int LocusInfo::HAPLODIPLOID = 0;
unsigned int LocusInfo::HAPLOID = 1;
unsigned int LocusInfo::DIPLOID = 2;

//** Class constructor: *******************************************************/
LocusInfo::LocusInfo(const string &name, const unsigned int ploidy) {
	_name = name;
	_ploidy = ploidy;
}

//** Class destructor: *******************************************************/
LocusInfo::~LocusInfo() {
	for (unsigned int i = 0 ; i < _alleles.size() ; i++)
		delete _alleles[i];
	_alleles.clear();
}

//** Other methodes: *********************************************************/
// Name
string LocusInfo::getName() const {
	return _name;
}

// Ploidie
unsigned int LocusInfo::getPloidy() const {
	return _ploidy;
}

// AlleleInfos
void LocusInfo::addAlleleInfo(const AlleleInfo &allele) throw (BadIdentifierException) {
	// Check if the allele id is not already in use
	for (unsigned int i = 0 ; i < _alleles.size() ; i++)
		if (_alleles[i]->getId() == allele.getId())
			throw BadIdentifierException("LocusInfo::addAlleleInfo: Id already in use.",allele.getId());
	_alleles.push_back(dynamic_cast<AlleleInfo *>(allele.clone()));
}

AlleleInfo * LocusInfo::getAlleleInfoById(const string & id) const
throw (AlleleNotFoundException) {
	for (unsigned int i = 0 ; i < _alleles.size() ; i++)
		if (_alleles[i]->getId() == id)
			return _alleles[i];
	throw AlleleNotFoundException("LocusInfo::getAlleleInfoById: AlleleInfo id unknown.", id);
}

AlleleInfo * LocusInfo::getAlleleInfoByKey(unsigned int key) const throw (IndexOutOfBoundsException) {
	if (key >= _alleles.size())
		throw IndexOutOfBoundsException("LocusInfo::getAlleleInfoByKey: key out of bounds.", key, 0, _alleles.size());
	return _alleles[key];
}

unsigned int LocusInfo::getAlleleInfoKey(const string & id) const
throw (AlleleNotFoundException) {
	for (unsigned int i = 0 ; i < _alleles.size() ; i++)
		if (_alleles[i]->getId() == id)
			return i;
	throw AlleleNotFoundException("LocusInfo::getAlleleInfoKey: AlleleInfo id not found.", id);
}

unsigned int LocusInfo::getNumberOfAlleles() const {
	return _alleles.size();
}

void LocusInfo::clear() {
	for (unsigned int i = 0 ; i < _alleles.size() ; i++)
		delete _alleles[i];
	_alleles.clear();
}
