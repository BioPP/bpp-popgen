/*
 * File Genotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 07 2004
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
		 */
		Genotype(AnalyzedLoci * analyzed_loci);

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

	protected:
		AnalyzedLoci * _analyzedLoci;
		vector<MonolocusGenotype *> _loci;
};
#endif // _GENOTYPE_H_
