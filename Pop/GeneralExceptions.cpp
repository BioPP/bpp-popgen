/*
 * File GeneralExceptions.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 01 2004
 */

#include "GeneralExceptions.h"

// From Utils
#include <Utils/TextTools.h>

BadIdentifierException::BadIdentifierException(const char *text,
		const unsigned int id): Exception("BadIdentifierException: " +
			string(text) + "(" + TextTools::toString(id) + ")"),
			_id(TextTools::toString(id)) {}
			
BadIdentifierException::BadIdentifierException(const string &text,
		const unsigned int id): Exception("BadIdentifierException: " +
			text + "(" + TextTools::toString(id) + ")"),
			_id(TextTools::toString(id)) {}

BadIdentifierException::BadIdentifierException(const char *text,
		const string &id): Exception("BadIdentifierException: " + string(text) +
			"(" + id + ")"),
			_id(id) {}

BadIdentifierException::BadIdentifierException(const string &text,
		const string &id): Exception("BadIdentifierException: " + text +
			"(" + id + ")"),
			_id(id) {}
			
BadIdentifierException::~BadIdentifierException() throw() {}

const string BadIdentifierException::getIdentifier() const {
	return _id;
}
