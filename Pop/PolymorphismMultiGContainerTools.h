/*
 * File PolymorphismMultiGContainerTools.h
 * Author : Sylvain Gailard <yragael2001@yahoo.fr>
 * Last modification : Thursday September 30 2004
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

//Secured inclusion of header's file
#ifndef _POLYMORPHISMMULTIGCONTAINERTOOLS_H_
#define _POLYMORPHISMMULTIGCONTAINERTOOLS_H_

// From the STL
#include <vector>
#include <set>
using namespace std;

//From the PolyLib library
#include "PolymorphismMultiGContainer.h"

//From the NumCalc library
#include <NumCalc/RandomTools.h>

/**
 * @brief Tools for PolymorphismMultiGContainer.
 *
 * Provides static methods for permutations.
 */
class PolymorphismMultiGContainerTools {
	public:
		virtual ~PolymorphismMultiGContainerTools();
		
	public:
		/**
		 * @brief Permut the MultilocusGenotype in the whole PolymorphismMultiGContainer.
		 *
		 * @param pmgc The PolymorphismMultiGContainer to permut.
		 * @return A permuted PolymorphismMultiGContainer.
		 */
		static PolymorphismMultiGContainer permutMultiG(const PolymorphismMultiGContainer & pmgc);

		/**
		 * @brief Permut the MonolocusGenotype.
		 *
		 * Permut the MonolocusGenotypes in one or several groups breaking
		 * the links between them.
		 *
		 * @param pmgc The PolymorphismMultiGContainer to permut.
		 * @param groups The groups ids between which the MonolocusGenotypes will be permuted.
		 * @return A permuted PolymorphismMultiGContainer.
		 */
		static PolymorphismMultiGContainer permutMonoG(const PolymorphismMultiGContainer & pmgc, const set<unsigned int> & groups);

		/**
		 * @brief Permut the Alleles.
		 *
		 * Permut the alleles in one or several groups breaking
		 * the links between them.
		 *
		 * @param pmgc The PolymorphismMultiGContainer to permut.
		 * @param groups The groups ids between which the MonolocusGenotypes will be permuted.
		 * @return A permuted PolymorphismMultiGContainer.
		 */
		static PolymorphismMultiGContainer permutAlleles(const PolymorphismMultiGContainer & pmgc, const set<unsigned int> & groups);
};

#endif // _POLYMORPHISMMULTIGCONTAINERTOOLS_H_
