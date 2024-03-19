// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
  public virtual MonolocusGenotypeInterface
{
private:
  size_t alleleIndex_;

public:
  /**
   * @brief Build a monolocus genotype containing one allele.
   */
  MonoAlleleMonolocusGenotype(size_t alleleIndex);

  /**
   * @brief Build a monolocus genotype containing one allele.
   */
  MonoAlleleMonolocusGenotype(std::vector<size_t> alleleIndex);

  /**
   * @brief Copy constructor.
   */
  MonoAlleleMonolocusGenotype(const MonoAlleleMonolocusGenotype& mmg);

  /**
   * @brief Destroy the MonoAlleleMonolocusGenotype.
   */
  virtual ~MonoAlleleMonolocusGenotype() = default;

public:
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
  std::vector<size_t> getAlleleIndex() const override;
  /** @} */

  /**
   * @name The Clonable interface:
   *
   * @{
   */
  MonoAlleleMonolocusGenotype* clone() const override
  {
    return new MonoAlleleMonolocusGenotype(*this);
  }
  /** @} */
};
} // end of namespace bpp;

#endif // _MONOALLELEMONOLOCUSGENOTYPE_H_
