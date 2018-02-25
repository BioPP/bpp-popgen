//
// File MonoAlleleMonolocusGenotype.h
// Author : Sylvain Gaillard <yragael2001@yahoo.fr>
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

#ifndef _MONOALLELEMONOLOCUSGENOTYPE_H_
#define _MONOALLELEMONOLOCUSGENOTYPE_H_

#include <Bpp/Exceptions.h>

// From local
#include "MonolocusGenotype.h"

namespace bpp
{
/**
 * @brief The MonoAlleleMonolocusGenotype class.
 *
 * @author Sylvain Gaillard
 */
class MonoAlleleMonolocusGenotype :
  public MonolocusGenotype
{
private:
  size_t allele_index_;

public:
  // Constructors and destructor
  /**
   * @brief Build a monolocus genotype containing one allele.
   */
  MonoAlleleMonolocusGenotype(size_t allele_index);

  /**
   * @brief Build a monolocus genotype containing one allele.
   */
  MonoAlleleMonolocusGenotype(std::vector<size_t> allele_index);

  /**
   * @brief Copy constructor.
   */
  MonoAlleleMonolocusGenotype(const MonoAlleleMonolocusGenotype& mmg);

  /**
   * @brief Destroy the MonoAlleleMonolocusGenotype.
   */
  ~MonoAlleleMonolocusGenotype();

public:
  // Other methodes
  /**
   * @brief The affectation operator.
   */
  MonoAlleleMonolocusGenotype& operator=(const MonoAlleleMonolocusGenotype& mmg);

  /**
   * @brief The == operator.
   */
  virtual bool operator==(const MonoAlleleMonolocusGenotype& mmg) const;

  /**
   * @name The MonolocusGenotype interface:
   *
   * @{
   */
  std::vector<size_t> getAlleleIndex() const;
  /** @} */

  /**
   * @name The Clonable interface:
   *
   * @{
   */
  MonoAlleleMonolocusGenotype* clone() const;
  /** @} */
};
} // end of namespace bpp;

#endif // _MONOALLELEMONOLOCUSGENOTYPE_H_
