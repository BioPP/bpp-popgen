/*
 * File AlleleInfo.h
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
#ifndef _ALLELEINFO_H_
#define _ALLELEINFO_H_

// From STL
#include <string>
using namespace std;

// From Utils
#include <Utils/Clonable.h>

/**
 * @brief The AlleleInfo interface.
 */
class AlleleInfo : public Clonable {
	public: // Destructor

		virtual ~AlleleInfo();

	public: // Methodes
		/**
		 * @brief Set the identifier of the allele.
		 */
		virtual void setId(const string & allele_id) = 0;
		
		/**
		 * @brief Get the identitier of the allele.
		 */
		virtual string getId() const = 0;
};

#endif // _ALLELEINFO_H_
