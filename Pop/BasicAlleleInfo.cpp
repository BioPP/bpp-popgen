/*
 * File BasicAlleleInfo.cpp
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

#include "BasicAlleleInfo.h"

//** Class constructor: *******************************************************/
BasicAlleleInfo::BasicAlleleInfo(const string & id) {
	_id = id;
}

BasicAlleleInfo::BasicAlleleInfo(const BasicAlleleInfo &allele) {
	_id = allele.getId();
}
//** Class destructor: *******************************************************/
BasicAlleleInfo::~BasicAlleleInfo() {}

//** Other methodes: *********************************************************/
BasicAlleleInfo & BasicAlleleInfo::operator= (const BasicAlleleInfo & allele) {
	_id = allele.getId();
	return * this;
}

bool BasicAlleleInfo::operator== (const BasicAlleleInfo & allele) const {
	return (_id == allele.getId());
}

bool BasicAlleleInfo::operator!= (const BasicAlleleInfo & allele) const {
	return !(_id == allele.getId());
}

Clonable * BasicAlleleInfo::clone() const {
	return new BasicAlleleInfo(* this);
}

void BasicAlleleInfo::setId(const string & allele_id) {
	_id = allele_id;
}

string BasicAlleleInfo::getId() const {
	return _id;
}
