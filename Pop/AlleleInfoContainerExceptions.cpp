/*
 * File AlleleInfoContainerExceptions.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 01 2004
 */

#include "AlleleInfoContainerExceptions.h"

// From Utils
#include <Utils/TextTools.h>

AlleleInfoNotFoundException::AlleleInfoNotFoundException(const char *text,
		const unsigned int id): Exception("AlleleInfoNotFoundException: " +
			string(text) + "(" + TextTools::toString(id) + ")"),
			_id(id) {}
			
AlleleInfoNotFoundException::AlleleInfoNotFoundException(const string &text,
		const unsigned int id): Exception("AlleleInfoNotFoundException: " +
			text + "(" + TextTools::toString(id) + ")"),
			_id(id) {}

AlleleInfoNotFoundException::~AlleleInfoNotFoundException() throw() {}

const unsigned int AlleleInfoNotFoundException::getAlleleInfoId() const {
	return _id;
}
