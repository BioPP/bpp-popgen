/*
 * File LocusInfo.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 01 2004
 */

// From Utils
#include <Utils/TextTools.h>

#include "LocusInfo.h"
#include "LocusInfoContainerExceptions.h"

unsigned int LocusInfo::HAPLODIPLOID = 0;
unsigned int LocusInfo::HAPLOID = 1;
unsigned int LocusInfo::DIPLOID = 2;

//** Class constructor: *******************************************************/
LocusInfo::LocusInfo(const string &name, const unsigned int ploidie) {
	_name = name;
	_ploidie = ploidie;
}

//** Class destructor: *******************************************************/
LocusInfo::~LocusInfo() {
	for (map<unsigned int, AlleleInfo *>::iterator it = _alleles.begin() ;
			it != _alleles.end() ; it++)
		delete it->second;
}

//** Other methodes: *********************************************************/
// Name
string LocusInfo::getName() const {
	return _name;
}

// Ploidie
unsigned int LocusInfo::getPloidie() const {
	return _ploidie;
}

// AlleleInfos
void LocusInfo::addAlleleInfo(const AlleleInfo &allele) throw (Exception) {
	// Check if the allele id is not already in use
	for (map<unsigned int, AlleleInfo *>::const_iterator it = _alleles.begin() ;
			it != _alleles.end() ; it++)
		if (it->second->getId() == allele.getId()) {
			string mes = "LocusInfo::addAlleleInfo: Id already in use (";
			mes += TextTools::toString(allele.getId());
			mes += ").";
			throw Exception(mes);
		}
	if (_alleles.size() == 0)
		_alleles.insert(make_pair(1, dynamic_cast<AlleleInfo *>(allele.clone())));
	else
		_alleles.insert(make_pair(_alleles.rbegin()->first + 1,
				dynamic_cast<AlleleInfo *>(allele.clone())));
}

AlleleInfo * LocusInfo::getAlleleInfoById(unsigned int id) const
throw (AlleleInfoNotFoundException) {
	for (map<unsigned int, AlleleInfo *>::const_iterator it = _alleles.begin() ;
			it != _alleles.end() ; it++) {
		if (it->second->getId() == id)
			return it->second;
	}
	throw AlleleInfoNotFoundException("LocusInfo::getAlleleInfoById: AlleleInfo id unknown.", id);
}

AlleleInfo * LocusInfo::getAlleleInfoByKey(unsigned int key) const throw (BadIntegerException) {
	map<unsigned int, AlleleInfo *>::const_iterator it = _alleles.find(key);
	if (it == _alleles.end())
		throw BadIntegerException("LocusInfo::getAlleleInfoByKey: Unknown key.", key);
	return it->second;
}

unsigned int LocusInfo::getAlleleInfoKey(unsigned int id) const
throw (AlleleInfoNotFoundException) {
	for (map<unsigned int, AlleleInfo *>::const_iterator it = _alleles.begin() ;
			it != _alleles.end() ; it++) {
		if (it->second->getId() == id)
			return it->first;
	}
	throw AlleleInfoNotFoundException("LocusInfo::getAlleleInfoPosition: AlleleInfo id unknown.",
			id);
}

unsigned int LocusInfo::getNumberOfAlleleInfos() const {
	return _alleles.size();
}

vector<unsigned int> LocusInfo::getAlleleInfosKeys() const throw (Exception) {
	vector<unsigned int> keys;
	for (map<unsigned int, AlleleInfo *>::const_iterator it = _alleles.begin() ;
			it != _alleles.end() ; it++)
		keys.push_back(it->first);
	if (keys.size() == 0) throw Exception("LocusInfo::getAlleleInfosKeys: No AlleleInfo in this LocusInfo.");
	return keys;
}

void LocusInfo::clear() {
	for (map<unsigned int, AlleleInfo *>::const_iterator it = _alleles.begin() ;
			it != _alleles.end() ; it++)
		delete it->second;
	_alleles.clear();
}
