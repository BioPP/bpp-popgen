/*
 * File MonoAlleleMonolocusGenotype.cpp
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

#include "MonoAlleleMonolocusGenotype.h"

//** Class constructor: *******************************************************/
MonoAlleleMonolocusGenotype::MonoAlleleMonolocusGenotype(unsigned int allele_index) {
	_allele_index = allele_index;
}

MonoAlleleMonolocusGenotype::MonoAlleleMonolocusGenotype(vector<unsigned int> allele_index) throw (BadIntegerException) {
	if (allele_index.size() != 1)
		throw BadIntegerException("MonoAlleleMonolocusGenotype::MonoAlleleMonolocusGenotype: allele_index must conaines one value.", allele_index.size());
	_allele_index = allele_index[0];
}

MonoAlleleMonolocusGenotype::MonoAlleleMonolocusGenotype(const MonoAlleleMonolocusGenotype & mmg) {
	_allele_index = mmg.getAlleleIndex()[0];
}

//** Class destructor: ********************************************************/
MonoAlleleMonolocusGenotype::~MonoAlleleMonolocusGenotype() {}

//** Other methodes: **********************************************************/
MonoAlleleMonolocusGenotype & MonoAlleleMonolocusGenotype::operator= (const MonoAlleleMonolocusGenotype & mmg) {
	_allele_index = mmg.getAlleleIndex()[0];
	return * this;
}

bool MonoAlleleMonolocusGenotype::operator== (const MonoAlleleMonolocusGenotype & mmg) const {
	return (_allele_index == mmg.getAlleleIndex()[0]);
}

vector<unsigned int> MonoAlleleMonolocusGenotype::getAlleleIndex() const {
	vector<unsigned int> index;
	index.push_back(_allele_index);
	return index;
}

Clonable * MonoAlleleMonolocusGenotype::clone() const {
	return new MonoAlleleMonolocusGenotype(* this);
}
