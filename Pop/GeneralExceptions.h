/*
 * File GeneralExceptions.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 21 2004
 */

// Secured inclusion of header's file
#ifndef _GENERALEXCEPTIONS_H_
#define _GENERALEXCEPTIONS_H_

// From Utils
#include <Utils/Exceptions.h>
//*****************************************************************************
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

//*****************************************************************************
/**
 * @brief The LocusNotFoundException class.
 */
class LocusNotFoundException : public BadIdentifierException {
	public:
		// Class constructor
		/**
		 * @brief Build the exception with a numerical identifier.
		 */
		LocusNotFoundException(const char *text, const unsigned int id);
		
		/**
		 * @brief Build the exception with a numerical identifier.
		 */
		LocusNotFoundException(const string &text, const unsigned int id);
		
		/**
		 * @brief Build the exception with a textual identifier.
		 */
		LocusNotFoundException(const char *text, const string &id);
		
		/**
		 * @brief Build the exception with a textual identifier.
		 */
		LocusNotFoundException(const string &text, const string &id);

		// Class destructor
		~LocusNotFoundException() throw();

	public:
		/**
		 * @brief Return the value of the identifier as a string.
		 */
		virtual const string getIdentifier() const;
};

//*****************************************************************************
/**
 * @brief The AlleleNotFoundException class.
 */
class AlleleNotFoundException : public BadIdentifierException {
	public:
		// Class constructor
		/**
		 * @brief Build the exception with a numerical identifier.
		 */
		AlleleNotFoundException(const char *text, const unsigned int id);
		
		/**
		 * @brief Build the exception with a numerical identifier.
		 */
		AlleleNotFoundException(const string &text, const unsigned int id);
		
		/**
		 * @brief Build the exception with a textual identifier.
		 */
		AlleleNotFoundException(const char *text, const string &id);
		
		/**
		 * @brief Build the exception with a textual identifier.
		 */
		AlleleNotFoundException(const string &text, const string &id);

		// Class destructor
		~AlleleNotFoundException() throw();

	public:
		/**
		 * @brief Return the value of the identifier as a string.
		 */
		virtual const string getIdentifier() const;
};

//*****************************************************************************
/**
 * @brief The LocalityNotFoundException class.
 */
class LocalityNotFoundException : public BadIdentifierException {
	public:
		// Class constructor
		/**
		 * @brief Build the exception with a numerical identifier.
		 */
		LocalityNotFoundException(const char *text, const unsigned int id);
		
		/**
		 * @brief Build the exception with a numerical identifier.
		 */
		LocalityNotFoundException(const string &text, const unsigned int id);
		
		/**
		 * @brief Build the exception with a textual identifier.
		 */
		LocalityNotFoundException(const char *text, const string &id);
		
		/**
		 * @brief Build the exception with a textual identifier.
		 */
		LocalityNotFoundException(const string &text, const string &id);

		// Class destructor
		~LocalityNotFoundException() throw();

	public:
		/**
		 * @brief Return the value of the identifier as a string.
		 */
		virtual const string getIdentifier() const;
};

//*****************************************************************************
/**
 * @brief The IndividualNotFoundException class.
 */
class IndividualNotFoundException : public BadIdentifierException {
	public:
		// Class constructor
		/**
		 * @brief Build the exception with a numerical identifier.
		 */
		IndividualNotFoundException(const char *text, const unsigned int id);
		
		/**
		 * @brief Build the exception with a numerical identifier.
		 */
		IndividualNotFoundException(const string &text, const unsigned int id);
		
		/**
		 * @brief Build the exception with a textual identifier.
		 */
		IndividualNotFoundException(const char *text, const string &id);
		
		/**
		 * @brief Build the exception with a textual identifier.
		 */
		IndividualNotFoundException(const string &text, const string &id);

		// Class destructor
		~IndividualNotFoundException() throw();

	public:
		/**
		 * @brief Return the value of the identifier as a string.
		 */
		virtual const string getIdentifier() const;
};

#endif // _GENERALEXCEPTIONS_H_
