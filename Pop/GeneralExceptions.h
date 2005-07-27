/*
 * File GeneralExceptions.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
 *
*/
/*
Copyright or © or Copr. CNRS, (November 17, 2004)


This software is a computer program whose purpose is to provide classes
for sequences analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
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

//*****************************************************************************
/**
 * @brief The GroupNotFoundException class.
 */
class GroupNotFoundException : public BadIdentifierException {
	public:
		// Class constructor
		/**
		 * @brief Build the exception with a numerical identifier.
		 */
		GroupNotFoundException(const char *text, const unsigned int id);
		
		/**
		 * @brief Build the exception with a numerical identifier.
		 */
		GroupNotFoundException(const string &text, const unsigned int id);
		
		/**
		 * @brief Build the exception with a textual identifier.
		 */
		GroupNotFoundException(const char *text, const string &id);
		
		/**
		 * @brief Build the exception with a textual identifier.
		 */
		GroupNotFoundException(const string &text, const string &id);

		// Class destructor
		~GroupNotFoundException() throw();

	public:
		/**
		 * @brief Return the value of the identifier as a string.
		 */
		virtual const string getIdentifier() const;
};

#endif // _GENERALEXCEPTIONS_H_
