// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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

#endif // _BIALLELEMONOLOCUSGENOTYPE_H_
