/*
 * File MultilocusGenotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
 *
 * Copyright (C) 2004 Sylvain Gaillard and the
 *                    PopGenLib Development Core Team
 *
 * This file is part of PopGenLib.
 *
 * PopGenLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopGenLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopGenLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
		void setMonolocusGenotype(unsigned int locus_position,
				const MonolocusGenotype & monogen) throw (IndexOutOfBoundsException);

		/**
		 * @brief Set a MonolocusGenotype by allele keys.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
		 * @throw Exception if there is no key in allele_keys.
		 */
		void setMonolocusGenotypeByAlleleKey(unsigned int locus_position,
				const vector<unsigned int> allele_keys) throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype by allele id.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
		 * @throw AlleleNotFoundException if at least one of the id is not found in the LocusInfo.
		 */
		void setMonolocusGenotypeByAlleleId(unsigned int locus_position,
				const vector<string> allele_id, const LocusInfo & locus_info) throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype as missing data.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
		 */
		void setMonolocusGenotypeAsMissing(unsigned int locus_position) throw (IndexOutOfBoundsException);

		/**
		 * @brief Tell if a MonolocusGenotype is a missing data.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
		 */
		bool isMonolocusGenotypeMissing(unsigned int locus_position) const throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Get a MonolocusGenotype.
		 */
		const MonolocusGenotype * getMonolocusGenotype(unsigned int locus_position) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Count the number of loci.
		 *
		 * Return the size of _loci.
		 */
		unsigned int size() const;

		/**
		 * @brief Count the number of non missing MonolocusGenotype.
		 */
		unsigned int countNonMissingLoci() const;

		/**
		 * @brief Count the number of homozygous MonolocusGenotype.
		 */
		unsigned int countHomozygousLoci() const;

		/**
		 * @brief Count the number of heterozygous MonolocusGenotype.
		 */
		unsigned int countHeterozygousLoci() const;

	protected:
		vector<MonolocusGenotype *> _loci;
};
#endif // _MULTILOCUSGENOTYPE_H_
