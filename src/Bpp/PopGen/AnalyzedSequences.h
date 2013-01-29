//
// File AnalyzedSequences.h
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

#ifndef _ANALYZEDSEQUENCES_H_
#define _ANALYZEDSEQUENCES_H_

// From Seq
#include <Bpp/Seq/Alphabet/Alphabet.h>

namespace bpp
{
/**
 * @brief The AnalyzedSequences class.
 *
 * This is a class to store info about the sequences.
 *
 * The object stores a pointer toward a const Alphabet.
 * The way the pointer is managed depend on the method used to set it.
 *
 * If one use a method using a const Alphabet* to set the Alphabet, then he
 * has to take care of the memory management (i.e. freeing the Alphabet
 * object).
 *
 * If one use a method that create an Alphabet object like those using a
 * string description of the Alphabet then the AnalyzedSequences object will
 * delete himself the Alphabet object on destruction.
 *
 * Be carefull when copying an AnalyzedSequences object, the way that the
 * Alphabet object is managed is also copyed then if the initial
 * AnalyzedSequences takes care of its Alphabet member then the copy will hold
 * copy af the Alphabet an manage it else the new AnalyzedSequences will just
 * copy the pointer and it's up to the user to take care of its deletion.
 *
 * @author Sylvain Gaillard
 */
class AnalyzedSequences
{
private:
  const Alphabet* alphabet_;
  bool autoset_;

public:
  // Constructor and destructor
  AnalyzedSequences();
  AnalyzedSequences(const Alphabet* alpha);
  ~AnalyzedSequences();

  // Copie constructor
  AnalyzedSequences(const AnalyzedSequences& as);
  AnalyzedSequences& operator=(const AnalyzedSequences& as);

public:
  /**
   * @brief Set the alphabet used for the sequences.
   */
  void setAlphabet(const Alphabet* alpha);

  /**
   * @brief Set the alphabet used for the sequences by alphabet type.
   */
  void setAlphabet(const std::string& alpha_type) throw (Exception);

  /**
   * @brief Get the alphabet.
   */
  const Alphabet* getAlphabet() const
  {
    return alphabet_;
  }

  /**
   * @brief Get the alphabet type as a string.
   */
  std::string getAlphabetType() const;

private:
  void clear_();
};
} // end of namespace bpp;

#endif // _ANALYZEDSEQUENCES_H_

