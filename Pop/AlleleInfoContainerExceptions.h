/*
 * File AlleleInfoContainerExceptions.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 01 2004
 */

// Secured inclusion of header's file
#ifndef _ALLELEINFOCONTAINEREXCEPTIONS_H_
#define _ALLELEINFOCONTAINEREXCEPTIONS_H_

// From Utils
#include <Utils/Exceptions.h>


class AlleleInfoNotFoundException : public Exception {
	public:
		// Class constructor
		AlleleInfoNotFoundException(const char *text, const unsigned int id);
		AlleleInfoNotFoundException(const string &text, const unsigned int id);

		// Class destructor
		~AlleleInfoNotFoundException() throw();

	public:
		virtual const unsigned int getAlleleInfoId() const;

	protected:
		const unsigned int _id;
};

#endif // _ALLELEINFOCONTAINEREXCEPTIONS_H_
