/*
 * File GeneralExceptions.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 21 2004
 */

#include "GeneralExceptions.h"

// From Utils
#include <Utils/TextTools.h>

//** BadIdentifierException **************************************************/
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

//** LocusNotFoundException **************************************************/
LocusNotFoundException::LocusNotFoundException(const char *text,
		const unsigned int id): BadIdentifierException("LocusNotFoundException: " +
			string(text) + "(" + TextTools::toString(id) + ")",
			id) {} 
			
LocusNotFoundException::LocusNotFoundException(const string &text,
		const unsigned int id): BadIdentifierException("LocusNotFoundException: " +
			text + "(" + TextTools::toString(id) + ")",
			id) {}

LocusNotFoundException::LocusNotFoundException(const char *text,
		const string &id): BadIdentifierException("LocusNotFoundException: " + string(text) +
			"(" + id + ")",
			id) {}

LocusNotFoundException::LocusNotFoundException(const string &text,
		const string &id): BadIdentifierException("LocusNotFoundException: " + text +
			"(" + id + ")",
			id) {}
			
LocusNotFoundException::~LocusNotFoundException() throw() {}

const string LocusNotFoundException::getIdentifier() const {
	return BadIdentifierException::getIdentifier();
}

//** AlleleNotFoundException **************************************************/
AlleleNotFoundException::AlleleNotFoundException(const char *text,
		const unsigned int id): BadIdentifierException("AlleleNotFoundException: " +
			string(text) + "(" + TextTools::toString(id) + ")",
			id) {}
			
AlleleNotFoundException::AlleleNotFoundException(const string &text,
		const unsigned int id): BadIdentifierException("AlleleNotFoundException: " +
			text + "(" + TextTools::toString(id) + ")",
			id) {}

AlleleNotFoundException::AlleleNotFoundException(const char *text,
		const string &id): BadIdentifierException("AlleleNotFoundException: " + string(text) +
			"(" + id + ")",
			id) {}

AlleleNotFoundException::AlleleNotFoundException(const string &text,
		const string &id): BadIdentifierException("AlleleNotFoundException: " + text +
			"(" + id + ")",
			id) {}
			
AlleleNotFoundException::~AlleleNotFoundException() throw() {}

const string AlleleNotFoundException::getIdentifier() const {
	return BadIdentifierException::getIdentifier();
}

//** LocalityNotFoundException **************************************************/
LocalityNotFoundException::LocalityNotFoundException(const char *text,
		const unsigned int id): BadIdentifierException("LocalityNotFoundException: " +
			string(text) + "(" + TextTools::toString(id) + ")",
			id) {}
			
LocalityNotFoundException::LocalityNotFoundException(const string &text,
		const unsigned int id): BadIdentifierException("LocalityNotFoundException: " +
			text + "(" + TextTools::toString(id) + ")",
			id) {}

LocalityNotFoundException::LocalityNotFoundException(const char *text,
		const string &id): BadIdentifierException("LocalityNotFoundException: " + string(text) +
			"(" + id + ")",
			id) {}

LocalityNotFoundException::LocalityNotFoundException(const string &text,
		const string &id): BadIdentifierException("LocalityNotFoundException: " + text +
			"(" + id + ")",
			id) {}
			
LocalityNotFoundException::~LocalityNotFoundException() throw() {}

const string LocalityNotFoundException::getIdentifier() const {
	return BadIdentifierException::getIdentifier();
}

//** IndividualNotFoundException **************************************************/
IndividualNotFoundException::IndividualNotFoundException(const char *text,
		const unsigned int id): BadIdentifierException("IndividualNotFoundException: " +
			string(text) + "(" + TextTools::toString(id) + ")",
			id) {}
			
IndividualNotFoundException::IndividualNotFoundException(const string &text,
		const unsigned int id): BadIdentifierException("IndividualNotFoundException: " +
			text + "(" + TextTools::toString(id) + ")",
			id) {}

IndividualNotFoundException::IndividualNotFoundException(const char *text,
		const string &id): BadIdentifierException("IndividualNotFoundException: " + string(text) +
			"(" + id + ")",
			id) {}

IndividualNotFoundException::IndividualNotFoundException(const string &text,
		const string &id): BadIdentifierException("IndividualNotFoundException: " + text +
			"(" + id + ")",
			id) {}
			
IndividualNotFoundException::~IndividualNotFoundException() throw() {}

const string IndividualNotFoundException::getIdentifier() const {
	return BadIdentifierException::getIdentifier();
}
