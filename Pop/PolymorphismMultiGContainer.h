/*
 * File PolymorphismMultiGContainer.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday September 28 2004
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
#ifndef _POLYMORPHYSMMULTIGCONTAINER_H_
#define _POLYMORPHYSMMULTIGCONTAINER_H_

// From C library
#include <cmath>

// From STL
#include <vector>
#include <map>
#include <set>
#include <algorithm>
using namespace std;

// From Utils
#include <Utils/Clonable.h>
#include <Utils/Exceptions.h>
#include <Utils/MapTools.h>

// From popgenlib
#include "MultilocusGenotype.h"
#include "GeneralExceptions.h"

/**
 * @brief The PolymorphismMultiGContainer class
 *
 * This class is a container of MultilocusGenotype.
 */
class PolymorphismMultiGContainer {
	public: // Constructors and destructor

		/**
		 * @brief Build a new PolymorphismMultilocusGenotypeContainer.
		 */
		PolymorphismMultiGContainer();

		/**
		 * @brief The copy constructor.
		 */
		PolymorphismMultiGContainer(const PolymorphismMultiGContainer & pmgc);

		/**
		 * @brief Destroy a PolymorphismMultilocusGenotypeContainer.
		 */
		~PolymorphismMultiGContainer();

	public:
		/**
		 * @brief The assignation operator=.
		 */
		PolymorphismMultiGContainer & operator= (const PolymorphismMultiGContainer & pmgc);

		/**
		 * @brief Add a MultilocusGenotype to the container.
		 */
		void addMultilocusGenotype(const MultilocusGenotype & mg, unsigned int group);

		/**
		 * @brief Get a MultilocusGenotype at a position.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		const MultilocusGenotype * getMultilocusGenotype(unsigned int position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Remove a MultilocusGenotype.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		MultilocusGenotype * removeMultilocusGenotype(unsigned int position) throw (IndexOutOfBoundsException);

		/**
		 * @brief Delete a MultilocusGenotype.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		void deleteMultilocusGenotype(unsigned int position) throw (IndexOutOfBoundsException);

		/**
		 * @brief Tell if the MultilocusGenotypes are aligned (i.e. same size).
		 */
		bool isAligned() const;

		/**
		 * @brief Get the number of loci if the MultilocusGenotypes are aligned.
		 *
		 * @throw Exception if MultilocusGenotypes are not aligned.
		 */
		unsigned int getNumberOfLoci() const throw (Exception);

		/**
		 * @brief Get the Group of a MultilocusGenotype.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		unsigned int getGroupId(unsigned int position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Set the Group of a MultilocusGenotype.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		void setGroupId(unsigned int position, unsigned int group_id) throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the groups' ids.
		 */
		set<unsigned int> getAllGroupsIds() const;

		/**
		 * @brief Tell if a group exists.
		 */
		bool groupExists(unsigned int group) const;

		/**
		 * @brief Get the number of groups.
		 */
		unsigned int getNumberOfGroups() const;

		/**
		 * @brief Get group size.
		 */
		unsigned int getGroupSize(unsigned int group) const;

		/**
		 * @brief Get the size of a group for a given locus.
		 */
		unsigned int getLocusGroupSize(unsigned int group, unsigned int locus_position) const;
		
		/**
		 * @brief Get the number of MultilocusGenotype.
		 */
		unsigned int size() const;
		
		/**
		 * @brief Clear the container.
		 */
		void clear();

	protected:
		vector<MultilocusGenotype *> _multilocusGenotypes;
		vector<unsigned int> _groups;
};

#endif // _POLYMORPHYSMMULTIGCONTAINER_H_
