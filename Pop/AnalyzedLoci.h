/*
 * File AnalyzedLoci.h
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
#ifndef _ANALYZEDLOCI_H_
#define _ANALYZEDLOCI_H_

// From STL
#include <vector>
using namespace std;

//From Utils
#include <Utils/Exceptions.h>

// From local
#include "LocusInfo.h"
#include "GeneralExceptions.h"

/**
 * @brief The AnalyzedLoci class.
 *
 * This is a LocusInfo container.
 * Its instanciation requires a number of locus wich is fixed
 * and can't be modified.
 */
class AnalyzedLoci {
	public: // Constructors and Destructor
		/**
		 * @brief Build a void AnalyzedLoci with a specific number of loci.
		 */
		AnalyzedLoci(unsigned int number_of_loci);

		/**
		 * @brief Copy constructor.
		 */
		AnalyzedLoci(const AnalyzedLoci & analyzed_loci);

		/**
		 * @brief Destroy the AnalyzedLoci.
		 */
		~AnalyzedLoci();

	public: // Other methodes
		/**
		 * @brief Set a LocusInfo.
		 *
		 * @throw IndexOutOfBoundsException if locus_position is out of bounds.
		 */
		void setLocusInfo(unsigned int locus_position, const LocusInfo & locus)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the position of a LocusInfo.
		 *
		 * @throw BadIdentifierException if locus_name is not found.
		 */
		unsigned int getLocusInfoPosition(const string & locus_name) const
			throw (BadIdentifierException);

		/**
		 * @brief Get a LocusInfo by name.
		 *
		 * @throw BadIdentifierException if locus_name is not found.
		 */
		const LocusInfo * getLocusInfoByName(const string & locus_name) const
			throw (BadIdentifierException);

		/**
		 * @brief Get a LocusInfo by its position.
		 *
		 * @throw NullPointerException if the LocusInfo is not difined.
		 * @throw IndexOutOfBoundsException if locus_position is out of bounds.
		 */
		const LocusInfo * getLocusInfoAtPosition(unsigned int locus_position) const
			throw (Exception);

		/**
		 * @brief Add an AlleleInfo to a LocusInfo by LocusInfo name.
		 *
		 * @throw BadIdentifierException if the allele's id is already in use.
		 * @throw LocusNotFoundException if locus_name is not found.
		 */
		void addAlleleInfoByLocusName(const string & locus_name,
				const AlleleInfo & allele)
			throw (Exception);

		/**
		 * @brief Add an AlleleInfo to a LocusInfo by its position.
		 *
		 * @throw BadIdentifierException if the allele's id is already in use.
		 * @throw IndexOutOfBoundsException if locus_position is out of bounds.
		 */
		void addAlleleInfoByLocusPosition(unsigned int locus_position,
				const AlleleInfo & allele)
			throw (Exception);

		/**
		 * @brief Get the number of loci.
		 */
		unsigned int getNumberOfLoci() const;
		
		/**
		 * @brief Get the number of alleles at each locus.
		 */
		vector<unsigned int> getNumberOfAlleles() const;

		/**
		 * @brief Get the ploidy of a locus by name.
		 *
		 * @throw LocusNotFoundException if locus_name is not found.
		 */
		unsigned int getPloidyByLocusName(const string & locus_name) const
			throw (LocusNotFoundException);

		/**
		 * @brief Get the ploidy of a locus by its position.
		 *
		 * @throw IndexOutOfBoundsException if locus_position is out of bounds.
		 */
		unsigned int getPloidyByLocusPosition(unsigned int locus_position) const
			throw (IndexOutOfBoundsException);
		
	protected:
		vector<LocusInfo *> _loci;
};
#endif // _ANALYZEDLOCI_H_
