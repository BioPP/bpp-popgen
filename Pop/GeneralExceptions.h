/*
 * File GeneralExceptions.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 14 2004
 */

// Secured inclusion of header's file
#ifndef _GENERALEXCEPTIONS_H_
#define _GENERALEXCEPTIONS_H_

// From Utils
#include <Utils/Exceptions.h>

/**
 * @brief The BadIdentifierException class.
 *
 * This exception is used when an identifier is not found.
 * The identifier can be either a string or an integer but its
 * value is stored as a string.
 */
class BadIdentifierException : public Exception {
	public:
		// Class constructor
		/**
		 * @brief Build the exception with a numerical identifier.
		 */
		BadIdentifierException(const char *text, const unsigned int id);
		/**
		 * @brief Build the exception with a numerical identifier.
		 */
		BadIdentifierException(const string &text, const unsigned int id);
		
		/**
		 * @brief Build the exception with a textual identifier.
		 */
		BadIdentifierException(const char *text, const string &id);
		/**
		 * @brief Build the exception with a textual identifier.
		 */
		BadIdentifierException(const string &text, const string &id);

		// Class destructor
		~BadIdentifierException() throw();

	public:
		/**
		 * @brief Return the value of the identifier as a string.
		 */
		virtual const string getIdentifier() const;

	protected:
		const string _id;
};

#endif // _GENERALEXCEPTIONS_H_
