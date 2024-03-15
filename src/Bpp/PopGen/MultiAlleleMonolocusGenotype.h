// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// Secured inclusion of header's file
#ifndef _MULTIALLELEMONOLOCUSGENOTYPE_H_
#define _MULTIALLELEMONOLOCUSGENOTYPE_H_

// From STL
#include <vector>

#include <Bpp/Exceptions.h>

// From local
#include "MonolocusGenotype.h"

namespace bpp
{
/**
 * @brief The MultiAlleleMonolocusGenotype class.
 *
 * This class is intended to handle monolocus genotype with many alleles
 * like polyploid loci or loci obtained by trace file without cutoff on
 * peaks or other filter.
 *
 * @author Sylvain Gaillard
 */
class MultiAlleleMonolocusGenotype :
  public virtual MonolocusGenotypeInterface
{
private:
  std::vector<size_t> alleleIndex_;

public:
  /**
   * @brief Build a monolocus genotype containing many alleles.
   */
  MultiAlleleMonolocusGenotype(std::vector<size_t> alleleIndex);

  /**
   * @brief Copy constructor.
   */
  MultiAlleleMonolocusGenotype(const MultiAlleleMonolocusGenotype& mmg);

  /**
   * @brief Destroy the MultiAlleleMonolocusGenotype.
   */
  virtual ~MultiAlleleMonolocusGenotype();

public:
  // Other methodes
  /**
   * @brief The affectation operator.
   */
  MultiAlleleMonolocusGenotype& operator=(const MultiAlleleMonolocusGenotype& mmg);

  /**
   * @brief The == operator.
   */
  bool operator==(const MultiAlleleMonolocusGenotype& mmg) const;

  /**
   * @brief Test the homozygozity of the locus (i.e. all allele are identical).
   */
  bool isHomozygous() const;

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
  MultiAlleleMonolocusGenotype* clone() const override
  {
    return new MultiAlleleMonolocusGenotype(*this);
  }
  /** @} */
};
} // end of namespace bpp;

#endif// _MULTIALLELEMONOLOCUSGENOTYPE_H_
