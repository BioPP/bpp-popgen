/*
 * File BasicAlleleInfo.h
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
#ifndef _BASICALLELEINFO_H_
#define _BASICALLELEINFO_H_

// From local Pop
#include "AlleleInfo.h"
#include "GeneralExceptions.h"

/**
 * @brief The BasicAlleleInfo class.
 *
 * This is the simplest allele class implementation which contains just an identitier.
 */
class BasicAlleleInfo : public AlleleInfo {
	public: // Constructors and destructor
		/**
		 * @brief Build a new allele.
		 *
		 * @param id The identity number of the allele.
		 */
		BasicAlleleInfo(const string & id);

		/**
		 * @brief The BasicAlleleInfo copy constructor.
		 */
		BasicAlleleInfo(const BasicAlleleInfo &allele);

		virtual ~BasicAlleleInfo();

	public: // Methodes
		/**
		 * @brief The assignation operator.
		 */
		virtual BasicAlleleInfo & operator= (const BasicAlleleInfo & allele);

		/**
		 * @brief The == operator.
		 */
		virtual bool operator== (const BasicAlleleInfo & allele) const;

		/**
		 * @brief The != operator.
		 */
		virtual bool operator!= (const BasicAlleleInfo & allele) const;

		/**
		 * @name The Clonable interface
		 * @{
		 */
		Clonable * clone() const;
		/** @} */
			
		/**
		 * @name The AlleleInfo interface
		 */
		void setId(const string & allele_id);
		string getId() const;
		/** @} */

	private:
		string _id;
};

#endif // _BASICALLELEINFO_H_
