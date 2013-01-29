//
// File MonolocusGenotype.h
// Author : Sylvain Gaillard
// Last modification : Thursday July 29 2004
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#ifndef _MONOLOCUSGENOTYPE_H_
#define _MONOLOCUSGENOTYPE_H_

// From STL
#include <vector>

#include <Bpp/Clonable.h>

namespace bpp
{
/**
 * @brief The MonolocusGenotype virtual class.
 *
 * A MonolocusGenotype containes the Alleles' keys defined in a Locus object.
 * This keys are returned as size_tegers.
 * This class is an interface for all monolocus genotypes.
 *
 * @author Sylvain Gaillard
 */
class MonolocusGenotype :
  public Clonable
{
public:
  // Constructors and Destructor
  /**
   * @brief Destroy a MonolocusGenotype.
   */
  virtual ~MonolocusGenotype() {}

public:
  // Methodes
  /**
   * @brief Get the alleles' index.
   *
   * The alleles' index are the position of the AlleleInfo in a LocusInfo object.
   * If no LocusInfo is used, the index are just numbers to identify the alleles.
   *
   * @return A vector of size_t.
   *
   * The size of the vector corresponds to the number of alleles at this locus.
   */
  virtual std::vector<size_t> getAlleleIndex() const = 0;
};
} // end of namespace bpp;

#endif // _MONOLOCUSGENOTYPE_H_

