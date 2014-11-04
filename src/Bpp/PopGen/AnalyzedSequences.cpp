//
// File AnalyzedSequences.cpp
// Created by: Sylvain Gaillard
// Created on: Thursday July 29 2004
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
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/RNA.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>

using namespace bpp;
using namespace std;

AnalyzedSequences::AnalyzedSequences() : alphabet_(0),
  autoset_(false) {}

AnalyzedSequences::AnalyzedSequences(const Alphabet* alpha) : alphabet_(alpha),
  autoset_(false) {}

AnalyzedSequences::~AnalyzedSequences()
{
  clear_();
}

AnalyzedSequences::AnalyzedSequences(const AnalyzedSequences& as) : alphabet_(0),
  autoset_(false)
{
  if (as.autoset_)
  {
    setAlphabet(as.getAlphabetType());
  }
  else
  {
    alphabet_ = as.alphabet_;
  }
  autoset_ = as.autoset_;
}

AnalyzedSequences& AnalyzedSequences::operator=(const AnalyzedSequences& as)
{
  if (as.autoset_)
  {
    setAlphabet(as.getAlphabetType());
  }
  else
  {
    alphabet_ = as.alphabet_;
  }
  autoset_ = as.autoset_;
  return *this;
}

void AnalyzedSequences::setAlphabet(const Alphabet* alpha)
{
  alphabet_ = alpha;
  autoset_ = false;
}

void AnalyzedSequences::setAlphabet(const std::string& alpha_type) throw (Exception)
{
  if (alpha_type != string("DNA") && alpha_type != string("RNA") && alpha_type != string("PROTEIN"))
    throw Exception(string("AnalyzedSequences::setAlphabet: bad alphabet type. (") + alpha_type + string(")."));
  Alphabet* alpha = 0;
  if (alpha_type == string("DNA"))
    alpha = new DNA();
  if (alpha_type == string("RNA"))
    alpha = new RNA();
  if (alpha_type == string("PROTEIN"))
    alpha = new ProteicAlphabet();
  alphabet_ = alpha;
  autoset_ = true;
}

std::string AnalyzedSequences::getAlphabetType() const
{
  if (alphabet_ == 0)
    return string("---");
  string alpha_type = alphabet_->getAlphabetType();
  size_t bs = alpha_type.find(" ", 0);
  alpha_type = string(alpha_type.begin(), alpha_type.begin() + static_cast<ptrdiff_t>(bs));
  if (alpha_type == "Proteic")
    alpha_type = "PROTEIN";
  return alpha_type;
}

void AnalyzedSequences::clear_()
{
  if (alphabet_ != 0 && autoset_)
  {
    delete alphabet_;
    alphabet_ = 0;
    autoset_ = false;
  }
}
