/*
 * File AlleleInfo.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 01 2004
 */

#include "AlleleInfo.h"

//** Class constructor: *******************************************************/
AlleleInfo::AlleleInfo(unsigned int id) {
	_id = id;
}

AlleleInfo::AlleleInfo(const AlleleInfo &allele) {
	_id = allele.getId();
}
//** Class destructor: *******************************************************/
AlleleInfo::~AlleleInfo() {}

//** Other methodes: *********************************************************/
Clonable * AlleleInfo::clone() const {
	return new AlleleInfo(* this);
}

bool AlleleInfo::operator== (const AlleleInfo & allele) const {
	if (_id == allele.getId())
		return true;
	else
		return false;
}

unsigned int AlleleInfo::getId() const {
	return _id;
}
