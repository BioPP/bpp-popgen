/*
 * File Genotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 22 2004
 */

// Secured inclusion of header's file
#ifndef _GENOTYPE_H_
#define _GENOTYPE_H_

// From STL
#include <vector>
using namespace std;

// From Utils
#include <Utils/Exceptions.h>

// From local
#include "AnalyzedLoci.h"
#include "MonolocusGenotype.h"
#include "HaplodiploidMonolocusGenotype.h"
#include "HaploidMonolocusGenotype.h"
#include "DiploidMonolocusGenotype.h"
#include "LocusInfo.h"

/**
 * @brief The Genotype class.
 *
 * This is a MonolocusGenotype containor.
 */
class Genotype {
	public: // Constructors and Destructor
		/**
		 * @brief Build a Genotype linked to an AnalyzedLoci object.
		 *
		 * @throw NullPointerException if analyzed_loci is NULL.
		 */
		Genotype(const AnalyzedLoci * analyzed_loci) throw (NullPointerException);

		/**
		 * @brief Copy constructor.
		 */
		Genotype(const Genotype & genotype);

		/**
		 * @brief Destroy a Genotype.
		 */
		~Genotype();
		
	public:
		/**
		 * @brief Set a MonolocusGenotype.
		 */
		void setMonolocusGenotype(unsigned int locus_index,
				const MonolocusGenotype & monogen) throw (IndexOutOfBoundsException);

		/**
		 * @brief Set a MonolocusGenotype by allele keys.
		 */
		void setMonolocusGenotypeByAlleleKey(unsigned int locus_index,
				const vector<unsigned int> allele_keys) throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype by allele id.
		 */
		void setMonolocusGenotypeByAlleleId(unsigned int locus_index,
				const vector<unsigned int> allele_id) throw (Exception);

		/**
		 * @brief Get the ploidy of a locus.
		 */
		unsigned int getPloidy(unsigned int locus_index) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get a MonolocusGenotype.
		 */
		const MonolocusGenotype * getMonolocusGenotype(unsigned int locus_index) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the AnalyzedLoci linked to this Genotype.
		 */
		const AnalyzedLoci * getAnalyzedLoci();

	protected:
		const AnalyzedLoci * _analyzedLoci;
		vector<MonolocusGenotype *> _loci;
};
#endif // _GENOTYPE_H_
