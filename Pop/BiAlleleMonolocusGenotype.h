/*
 * File BiAlleleMonolocusGenotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
 *
 * Copyright (C) 2004 Sylvain Gaillard and the
 *                    PopLib Development Core Team
 *
 * This file is part of PopLib.
 *
 * PopLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
