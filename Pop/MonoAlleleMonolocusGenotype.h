/*
 * File MonoAlleleMonolocusGenotype.h
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
