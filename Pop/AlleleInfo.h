/*
 * File AlleleInfo.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday June 30 2004
 */

// Secured inclusion of header's file
#ifndef _ALLELEINFO_H_
#define _ALLELEINFO_H_

// From STL
#include <string>
using namespace std;

// From Utils
#include <Utils/Clonable.h>

/**
 * @brief The AlleleInfo interface.
 */
class AlleleInfo : public Clonable {
	public: // Destructor

		virtual ~AlleleInfo();

	public: // Methodes
		/**
		 * @brief Set the identifier of the allele.
		 */
		virtual void setId(const string & allele_id) = 0;
		
		/**
		 * @brief Get the identitier of the allele.
		 */
		virtual string getId() const = 0;
};

#endif // _ALLELEINFO_H_
