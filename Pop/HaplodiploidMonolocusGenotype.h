/*
 * File HaplodiploidMonolocusGenotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 07 2004
 */

// Secured inclusion of header's file
#ifndef _HAPLODIPLOIDMONOLOCUSGENOTYPE_H_
#define _HAPLODIPLOIDMONOLOCUSGENOTYPE_H_

//From local
#include "MonolocusGenotype.h"
#include "LocusInfo.h"

/**
 * @brief The HaplodiploidMonolocusGenotype class.
 */
class HaplodiploidMonolocusGenotype : public MonolocusGenotype {
	public: // Constructors and destructor
		/**
		 * @brief Build a monolocus genotype containing one allele.
		 */
		HaplodiploidMonolocusGenotype(unsigned int allele_key);

		/**
		 * @brief Copy constructor.
		 */
		HaplodiploidMonolocusGenotype(const HaplodiploidMonolocusGenotype & hdmg);
		
		/**
		 * @brief Destroy the HaplodiploidMonolocusGenotype.
		 */
		~HaplodiploidMonolocusGenotype();

	public: // Other methodes
		/**
		 * @name The MonolocusGenotype interface:
		 *
		 * @{
		 */
		
		/**
		 * @brief Get the allele key.
		 */
		unsigned int getAlleleKey() const;

		/**
		 * @brief Get the ploidy of the locus.
		 */
		unsigned int getPloidy() const;

		/** @} */

		/**
		 * @name The Clonable interface:
		 *
		 * @{
		 */

		/**
		 * @brief Clone the HaplodiploidMonolocusGenotype.
		 */
		Clonable * clone() const;

		/** @} */

	protected:
		unsigned int _allele_key;
};
#endif // _HAPLODIPLOIDMONOLOCUSGENOTYPE_H_
