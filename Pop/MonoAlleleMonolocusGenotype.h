/*
 * File MonoAlleleMonolocusGenotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday July 05 2004
 */

// Secured inclusion of header's file
#ifndef _MONOALLELEMONOLOCUSGENOTYPE_H_
#define _MONOALLELEMONOLOCUSGENOTYPE_H_

// From Utils
#include <Utils/Exceptions.h>

// From local
#include "MonolocusGenotype.h"

/**
 * @brief The MonoAlleleMonolocusGenotype class.
 */
class MonoAlleleMonolocusGenotype : public MonolocusGenotype {
	public: // Constructors and destructor
		/**
		 * @brief Build a monolocus genotype containing one allele.
		 */
		MonoAlleleMonolocusGenotype(unsigned int allele_index);

		/**
		 * @brief Build a monolocus genotype containing one allele.
		 */
		MonoAlleleMonolocusGenotype(vector<unsigned int> allele_index) throw (BadIntegerException);

		/**
		 * @brief Copy constructor.
		 */
		MonoAlleleMonolocusGenotype(const MonoAlleleMonolocusGenotype & mmg);
		
		/**
		 * @brief Destroy the MonoAlleleMonolocusGenotype.
		 */
		~MonoAlleleMonolocusGenotype();

	public: // Other methodes
		/**
		 * @brief The affectation operator.
		 */
		MonoAlleleMonolocusGenotype & operator= (const MonoAlleleMonolocusGenotype & mmg);
		
		/**
		 * @brief The == operator.
		 */
		virtual bool operator== (const MonoAlleleMonolocusGenotype & mmg) const;
		
		/**
		 * @name The MonolocusGenotype interface:
		 *
		 * @{
		 */
		vector<unsigned int> getAlleleIndex() const;
		/** @} */

		/**
		 * @name The Clonable interface:
		 *
		 * @{
		 */
		Clonable * clone() const;
		/** @} */

	protected:
		unsigned int _allele_index;
};
#endif // _MONOALLELEMONOLOCUSGENOTYPE_H_
