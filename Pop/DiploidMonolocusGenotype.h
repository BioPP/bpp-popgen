/*
 * File DiploidMonolocusGenotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 07 2004
 */

// Secured inclusion of header's file
#ifndef _DIPLOIDMONOLOCUSGENOTYPE_H_
#define _DIPLOIDMONOLOCUSGENOTYPE_H_

// From STL
#include <vector>
using namespace std;

// From Utils
#include <Utils/Exceptions.h>

//From local
#include "MonolocusGenotype.h"
#include "LocusInfo.h"

/**
 * @brief The DiploidMonolocusGenotype class.
 */
class DiploidMonolocusGenotype : public MonolocusGenotype {
	public: // Constructors and destructor
		/**
		 * @brief Build a monolocus genotype containing two alleles.
		 */
		DiploidMonolocusGenotype(unsigned int first_allele_key,
				unsigned int second_allele_key);
		
		/**
		 * @brief Copy constructor.
		 */
		DiploidMonolocusGenotype(const DiploidMonolocusGenotype & dmg);

		/**
		 * @brief Destroy the DiploidMonolocusGenotype.
		 */
		~DiploidMonolocusGenotype();

	public: // Other methodes
		/**
		 * @brief Get a specific allele key.
		 */
		unsigned int getAlleleKey(unsigned int allele_index) const
			throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Test the homozygozity of the locus.
		 */
		bool isHomozygous() const;

		/**
		 * @name The MonolocusGenotype interface:
		 *
		 * @{
		 */
		
		/**
		 * @brief Get the first allele key.
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
		 * @brief Clone the DiploidMonolocusGenotype.
		 */
		Clonable * clone() const;

		/** @} **/
	protected:
		vector<unsigned int> _allelekeys;
};
#endif // _DIPLOIDMONOLOCUSGENOTYPE_H_
