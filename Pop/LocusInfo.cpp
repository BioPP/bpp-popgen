/*
 * File LocusInfo.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
 *
 * Copyright (C) 2004 Sylvain Gaillard and the
 *                    PopLib Development Core Team
 *
 * This file is part of PopLib.
 *
 * PopLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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

LocusInfo::LocusInfo(const LocusInfo & locus_info) {
	_name = locus_info.getName();
	_ploidy = locus_info.getPloidy();
	for (unsigned int i = 0 ; i < locus_info.getNumberOfAlleles() ; i++) {
		Clonable * tmp_allele = locus_info.getAlleleInfoByKey(i)->clone();
		_alleles.push_back(dynamic_cast<AlleleInfo *>(tmp_allele));
	}
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
