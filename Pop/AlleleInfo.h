/*
 * File AlleleInfo.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 01 2004
 */

// Secured inclusion of header's file
#ifndef _ALLELEINFO_H_
#define _ALLELEINFO_H_

// From Utils
#include <Utils/Clonable.h>

/**
 * @brief The AlleleInfo class.
 *
 * This is the simplest allele class which contains just an identity number.
 */
class AlleleInfo : public Clonable {
	public: // Constructors and destructor
		/**
		 * @brief Build a new allele.
		 *
		 * @param id The identity number of the allele.
		 */
		AlleleInfo(unsigned int id);

		/**
		 * @brief The AlleleInfo copy constructor.
		 */
		AlleleInfo(const AlleleInfo &allele);

		virtual ~AlleleInfo();

	public: // Methodes
		/**
		 * @brief Implement the Cloanble interface.
		 */
		Clonable * clone() const;
			
		/**
		 * @brief The == operator.
		 */
		virtual bool operator== (const AlleleInfo & allele) const;
		
		/**
		 * @brief The !h operator.
		 */
		bool operator!= (const AlleleInfo & allele) const {
			return !(*this == allele);
		}
		
		/**
		 * @brief Get the identity number of the allele.
		 */
		unsigned int getId() const;

	private:
		unsigned int _id;
};

#endif // _ALLELEINFO_H_
