/*
 * File MultilocusGenotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday July 05 2004
 */

// Secured inclusion of header's file
#ifndef _MULTILOCUSGENOTYPE_H_
#define _MULTILOCUSGENOTYPE_H_

// From STL
#include <vector>
using namespace std;

// From Utils
#include <Utils/Exceptions.h>

// From local
#include "MonolocusGenotype.h"
#include "BiAlleleMonolocusGenotype.h"
#include "MonoAlleleMonolocusGenotype.h"
#include "LocusInfo.h"

/**
 * @brief The MultilocusGenotype class.
 *
 * This is a MonolocusGenotype containor.
 */
class MultilocusGenotype {
	public: // Constructors and Destructor
		/**
		 * @brief Build a MultilocusGenotype linked to an AnalyzedLoci object.
		 *
		 * @throw BadIntegerException if loci_number < 1.
		 */
		MultilocusGenotype(unsigned int loci_number) throw (BadIntegerException);

		/**
		 * @brief Copy constructor.
		 */
		MultilocusGenotype(const MultilocusGenotype & genotype);

		/**
		 * @brief Destroy a MultilocusGenotype.
		 */
		~MultilocusGenotype();
		
	public:
		/**
		 * @brief Set a MonolocusGenotype.
		 */
		void setMonolocusGenotype(unsigned int locus_index,
				const MonolocusGenotype & monogen) throw (IndexOutOfBoundsException);

		/**
		 * @brief Set a MonolocusGenotype by allele keys.
		 *
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of loci.
		 * @throw Exception if there is no key in allele_keys.
		 */
		void setMonolocusGenotypeByAlleleKey(unsigned int locus_index,
				const vector<unsigned int> allele_keys) throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype by allele id.
		 *
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of loci.
		 * @throw AlleleNotFoundException if at least one of the id is not found in the LocusInfo.
		 */
		void setMonolocusGenotypeByAlleleId(unsigned int locus_index,
				const vector<string> allele_id, const LocusInfo & locus_info) throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype as missing data.
		 *
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of loci.
		 */
		void setMonolocusGenotypeAsMissing(unsigned int locus_index) throw (IndexOutOfBoundsException);

		/**
		 * @brief Tell if a MonolocusGenotype is a missing data.
		 *
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of loci.
		 */
		bool isMonolocusGenotypeMissing(unsigned int locus_index) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Get a MonolocusGenotype.
		 */
		const MonolocusGenotype * getMonolocusGenotype(unsigned int locus_index) const
			throw (IndexOutOfBoundsException);

	protected:
		vector<MonolocusGenotype *> _loci;
};
#endif // _MULTILOCUSGENOTYPE_H_
