/*
 * File BiAlleleMonolocusGenotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday July 05 2004
 */

// Secured inclusion of header's file
#ifndef _BIALLELEMONOLOCUSGENOTYPE_H_
#define _BIALLELEMONOLOCUSGENOTYPE_H_

// From STL
#include <vector>
using namespace std;

// From Utils
#include <Utils/Exceptions.h>

//From local
#include "MonolocusGenotype.h"

/**
 * @brief The BiAlleleMonolocusGenotype class.
 */
class BiAlleleMonolocusGenotype : public MonolocusGenotype {
	public: // Constructors and destructor
		/**
		 * @brief Build a monolocus genotype containing two alleles.
		 */
		BiAlleleMonolocusGenotype(unsigned int first_allele_index,
				unsigned int second_allele_index);
		
		/**
		 * @brief Build a monolocus genotype containing two alleles.
		 */
		BiAlleleMonolocusGenotype(vector<unsigned int> allele_index) throw (BadIntegerException);
		
		/**
		 * @brief Copy constructor.
		 */
		BiAlleleMonolocusGenotype(const BiAlleleMonolocusGenotype & bmg);

		/**
		 * @brief Destroy the BiAlleleMonolocusGenotype.
		 */
		~BiAlleleMonolocusGenotype();

	public: // Other methodes
		/**
		 * @brief The affectation operator.
		 */
		BiAlleleMonolocusGenotype & operator= (const BiAlleleMonolocusGenotype & bmg);

		/**
		 * @brief The == operator.
		 */
		bool operator== (const BiAlleleMonolocusGenotype & bmg) const;

		/**
		 * @brief Get the first allele index.
		 */
		unsigned int getFirstAlleleIndex() const;

		/**
		 * @brief Get the second allele index.
		 */
		unsigned int getSecondAlleleIndex() const;
		
		/**
		 * @brief Test the homozygozity of the locus.
		 */
		bool isHomozygous() const;

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
		/** @} **/
	protected:
		vector<unsigned int> _allele_index;
};
#endif // _BIALLELEMONOLOCUSGENOTYPE_H_
