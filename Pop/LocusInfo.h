/*
 * File LocusInfo.h
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
#ifndef _LOCUSINFO_H_
#define _LOCUSINFO_H_

//From STL
#include <vector>
using namespace std;

// From local Poplib
#include "AlleleInfo.h"
#include "GeneralExceptions.h"

// From Utils
#include <Utils/Exceptions.h>

/**
 * @brief The LocusInfo class.
 *
 * This is an AlleleInfo container with additionnal data like a name,
 * the ploidy and some comments.
 */
class LocusInfo {
	public: // Constantes
		static unsigned int HAPLODIPLOID;
		static unsigned int HAPLOID;
		static unsigned int DIPLOID;
		
	public: // Constructors and destructor
		/**
		 * @brief Build a new LocusInfo object.
		 *
		 * @param name The name of the locus.
		 * @param ploidy The ploidy of the locus.
		 */
		LocusInfo(const string &name, const unsigned int ploidy = DIPLOID);

		/**
		 * @brief Copy constructor.
		 */
		LocusInfo(const LocusInfo & locus_info);

		/**
		 * @brief Destroy the LocusInfo.
		 */
		virtual ~LocusInfo();
		
	public: // Methodes
		/**
		 * @brief Get the name of the locus.
		 */
		string getName() const;

		/**
		 * @brief Get the ploidy of the locus.
		 *
		 * @return The ploidy as an unsigned integer.
		 */
		unsigned int getPloidy() const;
		
		/**
		 * @brief Add an AlleleInfo to the LocusInfo.
		 *
		 * @throw BadIdentifierException if the AlleleInfo's id already exists.
		 */
		void addAlleleInfo(const AlleleInfo &allele)
			throw (BadIdentifierException);

		/**
		 * @brief Retrieve an AlleleInfo object of the LocusInfo.
		 *
		 * @throw AlleleNotFoundException if the id is not found.
		 */
		AlleleInfo * getAlleleInfoById(const string & id) const
			throw (AlleleNotFoundException);

		/**
		 * @brief Retrieve an AlleleInfo object of the LocusInfo.
		 *
		 * @throw IndexOutOfBoundsException if key excedes the number of alleles.
		 */
		AlleleInfo * getAlleleInfoByKey(unsigned int key) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the position of an AlleleInfo.
		 *
		 * @throw AlleleNotFoundException if the AlleleInfo's id is not found.
		 */
		unsigned int getAlleleInfoKey(const string & id) const
			throw (AlleleNotFoundException);
		
		/**
		 * @brief Get the number of alleles at this locus.
		 */
		unsigned int getNumberOfAlleles() const;

		/**
		 * @brief Delete all alleles from the locus.
		 */
		void clear();

	protected:
		string _name;
		unsigned int _ploidy;
		vector <AlleleInfo *> _alleles;
};

#endif // _LOCUSINFO_H_
