/*
 * File LocusInfoContainerExceptions.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 01 2004
 */

#include "LocusInfoContainerExceptions.h"

LocusInfoNotFoundException::LocusInfoNotFoundException(const char *text,
		const string &locus_name): Exception("LocusInfoNotFoundException: " +
			string(text) + "(" + locus_name + ")"),
			_locus_name(locus_name) {}
			
LocusInfoNotFoundException::LocusInfoNotFoundException(const string &text,
		const string &locus_name): Exception("LocusInfoNotFoundException: " +
			text + "(" + locus_name + ")"),
			_locus_name(locus_name) {}

LocusInfoNotFoundException::~LocusInfoNotFoundException() throw() {}

const string LocusInfoNotFoundException::getLocusInfoName() const {
	return _locus_name;
}
