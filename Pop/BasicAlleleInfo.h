/*
 * File BasicAlleleInfo.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday June 30 2004
 */

// Secured inclusion of header's file
#ifndef _BASICALLELEINFO_H_
#define _BASICALLELEINFO_H_

// From local Pop
#include "AlleleInfo.h"
#include "GeneralExceptions.h"

/**
 * @brief The BasicAlleleInfo class.
 *
 * This is the simplest allele class implementation which contains just an identitier.
 */
class BasicAlleleInfo : public AlleleInfo {
	public: // Constructors and destructor
		/**
		 * @brief Build a new allele.
		 *
		 * @param id The identity number of the allele.
		 */
		BasicAlleleInfo(const string & id);

		/**
		 * @brief The BasicAlleleInfo copy constructor.
		 */
		BasicAlleleInfo(const BasicAlleleInfo &allele);

		virtual ~BasicAlleleInfo();

	public: // Methodes
		/**
		 * @brief The assignation operator.
		 */
		virtual BasicAlleleInfo & operator= (const BasicAlleleInfo & allele);

		/**
		 * @brief The == operator.
		 */
		virtual bool operator== (const BasicAlleleInfo & allele) const;

		/**
		 * @brief The != operator.
		 */
		virtual bool operator!= (const BasicAlleleInfo & allele) const;

		/**
		 * @name The Clonable interface
		 * @{
		 */
		Clonable * clone() const;
		/** @} */
			
		/**
		 * @name The AlleleInfo interface
		 */
		void setId(const string & allele_id);
		string getId() const;
		/** @} */

	private:
		string _id;
};

#endif // _BASICALLELEINFO_H_
