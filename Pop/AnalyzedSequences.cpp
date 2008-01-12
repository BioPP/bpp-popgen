//
// File AnalyzedSequences.cpp
// Author : Sylvain Gaillard
// Last modification : Thursday July 29 2004
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for population genetics analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "AnalyzedSequences.h"

using namespace bpp;

AnalyzedSequences::AnalyzedSequences()
{
	_alphabet = NULL;
}

AnalyzedSequences::~AnalyzedSequences()
{
	if (_alphabet != NULL)
		delete _alphabet;
}

void AnalyzedSequences::setAlphabet(const Alphabet * alpha)
{
	_alphabet = alpha;
}

void AnalyzedSequences::setAlphabet(const string & alpha_type) throw (Exception)
{
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

const Alphabet * AnalyzedSequences::getAlphabet() const
{
	return _alphabet;
}

string AnalyzedSequences::getAlphabetType() const
{
	if (_alphabet == NULL)
		return string("---");
	string alpha_type = _alphabet->getAlphabetType();
	int bs = alpha_type.find(" ",0);
	alpha_type = string(alpha_type.begin(), alpha_type.begin() + bs);
	if (alpha_type == "Proteic")
		alpha_type = "PROTEIN";
	return alpha_type;
}

