/*
 * File MonolocusGenotype.h
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
#ifndef _MONOLOCUSGENOTYPE_H_
#define _MONOLOCUSGENOTYPE_H_

// From STL
#include <vector>
using namespace std;

// From Utils
#include <Utils/Clonable.h>

/**
 * @brief The MonolocusGenotype virtual class.
 *
 * A MonolocusGenotype containes the Alleles' keys defined in a Locus object.
 * This keys are returned as unsigned integers.
 * This class is an interface for all monolocus genotypes.
 */
class MonolocusGenotype : public Clonable {
	public: // Constructors and Destructor
		/**
		 * @brief Destroy a MonolocusGenotype.
		 */
		virtual ~MonolocusGenotype();
		
	public: // Methodes
		/**
		 * @brief Get the alleles' index.
		 *
		 * The alleles' index are the position of the AlleleInfo in a LocusInfo object.
		 * If no LocusInfo is used, the index are just numbers to identify the alleles.
		 *
		 * @return A vector of unsigned int.
		 *
		 * The size of the vector corresponds to the number of alleles at this locus.
		 */
		virtual vector<unsigned int> getAlleleIndex() const = 0;
};
#endif // _MONOLOCUSGENOTYPE_H_
