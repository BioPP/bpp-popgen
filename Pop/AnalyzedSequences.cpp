/*
 * File AnalyzedSequences.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 08 2004
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
