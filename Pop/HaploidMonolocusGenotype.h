/*
 * File HaploidMonolocusGenotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 22 2004
 */

// Secured inclusion of header's file
#ifndef _HAPLOIDMONOLOCUSGENOTYPE_H_
#define _HAPLOIDMONOLOCUSGENOTYPE_H_

//From local
#include "MonolocusGenotype.h"
#include "LocusInfo.h"

/**
 * @brief The HaploidMonolocusGenotype class.
 */
class HaploidMonolocusGenotype : public MonolocusGenotype {
	public: // Constructors and destructor
		/**
		 * @brief Build a monolocus genotype containing one allele.
		 */
		HaploidMonolocusGenotype(unsigned int allele_key);

		/**
		 * @brief Copy constructor.
		 */
		HaploidMonolocusGenotype(const HaploidMonolocusGenotype & hmg);
		
		/**
		 * @brief Destroy the HaploidMonolocusGenotype.
		 */
		~HaploidMonolocusGenotype();

	public: // Other methodes
		/**
		 * @name The MonolocusGenotype interface:
		 *
		 * @{
		 */
		
		unsigned int getAlleleKey() const;
		
		unsigned int getPloidy() const;

		unsigned int getSize() const;

		/** @} */

		/**
		 * @name The Clonable interface:
		 *
		 * @{
		 */

		/**
		 * @brief Clone the HaploidMonolocusGenotype.
		 */
		Clonable * clone() const;

		/** @} */

	protected:
		unsigned int _allele_key;
};
#endif // _HAPLOIDMONOLOCUSGENOTYPE_H_
