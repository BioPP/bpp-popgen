//
// File BiAlleleMonolocusGenotype.h
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
// Secured inclusion of header's file
#ifndef _BIALLELEMONOLOCUSGENOTYPE_H_
#define _BIALLELEMONOLOCUSGENOTYPE_H_

// From STL
#include <vector>

#include <Bpp/Exceptions.h>

// From local
#include "MonolocusGenotype.h"

namespace bpp
{
/**
 * @brief The BiAlleleMonolocusGenotype class.
 *
 * @author Sylvain Gaillard
 */
class BiAlleleMonolocusGenotype :
  public virtual MonolocusGenotypeInterface
{
private:
  std::vector<size_t> alleleIndex_;

public:
  // Constructors and destructor
  /**
   * @brief Build a monolocus genotype containing two alleles.
   */
  BiAlleleMonolocusGenotype(size_t firstAlleleIndex,
                            size_t secondAlleleIndex);

  /**
   * @brief Build a monolocus genotype containing two alleles.
   */
  BiAlleleMonolocusGenotype(std::vector<size_t> alleleIndex);

  /**
   * @brief Copy constructor.
   */
  BiAlleleMonolocusGenotype(const BiAlleleMonolocusGenotype& bmg);

  /**
   * @brief Destroy the BiAlleleMonolocusGenotype.
   */
  virtual ~BiAlleleMonolocusGenotype();

public:
  /**
   * @brief The affectation operator.
   */
  BiAlleleMonolocusGenotype& operator=(const BiAlleleMonolocusGenotype& bmg);

  /**
   * @brief The == operator.
   */
  bool operator==(const BiAlleleMonolocusGenotype& bmg) const;

  /**
   * @brief Get the first allele index.
   */
  size_t getFirstAlleleIndex() const
  {
    return alleleIndex_[0];
  }

  /**
   * @brief Get the second allele index.
   */
  size_t getSecondAlleleIndex() const
  {
    return alleleIndex_[1];
  }

  /**
   * @brief Test the homozygozity of the locus.
   */
  bool isHomozygous() const
  {
    return alleleIndex_[0] == alleleIndex_[1];
  }

  /**
   * @name The MonolocusGenotype interface:
   *
   * @{
   */
  std::vector<size_t> getAlleleIndex() const override
  {
    return alleleIndex_;
  }

  /** @} */

  /**
   * @name The Clonable interface:
   *
   * @{
   */
  BiAlleleMonolocusGenotype* clone() const override
  {
    return new BiAlleleMonolocusGenotype(*this);
  }
  /** @} */
};
} // end of namespace bpp;

#endif// _BIALLELEMONOLOCUSGENOTYPE_H_
