/*
 * File AnalyzedSequences.h
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
#ifndef _ANALYZEDSEQUENCES_H_
#define _ANALYZEDSEQUENCES_H_

// From Seq
#include <Seq/Alphabet.h>
#include <Seq/DNA.h>
#include <Seq/RNA.h>
#include <Seq/ProteicAlphabet.h>

/**
 * @brief The AnalyzedSequences class.
 *
 * This is a class to store info about the sequences.
 */
class AnalyzedSequences {
	public: // Constructor and destructor
		AnalyzedSequences();
		~AnalyzedSequences();

	public:
		/**
		 * @brief Set the alphabet used for the sequences.
		 */
		void setAlphabet(const Alphabet * alpha);

		/**
		 * @brief Set the alphabet used for the sequences by alphabet type.
		 */
		void setAlphabet(const string & alpha_type) throw (Exception);

		/**
		 * @brief Get the alphabet.
		 */
		const Alphabet * getAlphabet() const;

		/**
		 * @brief Get the alphabet type as a string.
		 */
		string getAlphabetType() const;

	protected:
		const Alphabet * _alphabet;
};

#endif // _ANALYZEDSEQUENCES_H_
