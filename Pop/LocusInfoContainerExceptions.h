/*
 * File LocusInfoContainerExceptions.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 01 2004
 */

// Secured inclusion of header's file
#ifndef _LOCUSINFOCONTAINEREXCEPTIONS_H_
#define _LOCUSINFOCONTAINEREXCEPTIONS_H_

// From Utils
#include <Utils/Exceptions.h>


class LocusInfoNotFoundException : public Exception {
	public:
		// Class constructor
		LocusInfoNotFoundException(const char *text, const string &locus_name);
		LocusInfoNotFoundException(const string &text, const string &locus_name);

		// Class destructor
		~LocusInfoNotFoundException() throw();

	public:
		virtual const string getLocusInfoName() const;

	protected:
		const string _locus_name;
};

#endif // _LOCUSINFOCONTAINEREXCEPTIONS_H_
