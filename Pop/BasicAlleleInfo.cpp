/*
 * File BasicAlleleInfo.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday June 30 2004
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
