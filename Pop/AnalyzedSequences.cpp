/*
 * File AnalyzedSequences.cpp
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

#include "AnalyzedSequences.h"

AnalyzedSequences::AnalyzedSequences() {
	_alphabet = NULL;
}

AnalyzedSequences::~AnalyzedSequences() {
	if (_alphabet != NULL)
		delete _alphabet;
}

void AnalyzedSequences::setAlphabet(const Alphabet * alpha) {
	_alphabet = alpha;
}

void AnalyzedSequences::setAlphabet(const string & alpha_type) throw (Exception) {
	if (alpha_type != string("DNA") && alpha_type != string("RNA") && alpha_type != string("PROTEIN"))
		throw Exception(string("AnalyzedSequences::setAlphabet: bad alphabet type. (") + alpha_type + string(")."));
	Alphabet * alpha = NULL;
	if (alpha_type == string("DNA"))
		alpha = new DNA();
	if (alpha_type == string("RNA"))
		alpha = new RNA();
	if (alpha_type == string("PROTEIN"))
		alpha = new ProteicAlphabet();
	setAlphabet(alpha);
}

const Alphabet * AnalyzedSequences::getAlphabet() const {
	return _alphabet;
}

string AnalyzedSequences::getAlphabetType() const {
	if (_alphabet == NULL)
		return string("---");
	string alpha_type = _alphabet->getAlphabetType();
	int bs = alpha_type.find(" ",0);
	alpha_type = string(alpha_type.begin(), alpha_type.begin() + bs);
	if (alpha_type == "Proteic")
		alpha_type = "PROTEIN";
	return alpha_type;
}
