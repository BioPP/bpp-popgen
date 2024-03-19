// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _MULTILOCUSGENOTYPE_H_
#define _MULTILOCUSGENOTYPE_H_

// From STL
#include <vector>
#include <string>
#include <memory>

#include <Bpp/Exceptions.h>

// From Pop
#include "MonolocusGenotype.h"
#include "MonolocusGenotypeTools.h"
#include "BiAlleleMonolocusGenotype.h"
#include "MonoAlleleMonolocusGenotype.h"
#include "LocusInfo.h"

namespace bpp
{
/**
 * @brief The MultilocusGenotype class.
 *
 * This is a MonolocusGenotype container.
 *
 * @author Sylvain Gaillard
 */
class MultilocusGenotype :
  public virtual Clonable
{
private:
  std::vector<std::unique_ptr<MonolocusGenotypeInterface>> loci_;

public:
  /**
   * @brief Build a MultilocusGenotype linked to an AnalyzedLoci object.
   *
   * @throw BadIntegerException if lociNumber < 1.
   */
  MultilocusGenotype(size_t lociNumber);

  /**
   * @brief Copy constructor.
   */
  MultilocusGenotype(const MultilocusGenotype& genotype);

  /**
   * @brief Destroy a MultilocusGenotype.
   */
  virtual ~MultilocusGenotype();

  MultilocusGenotype* clone() const override { return new MultilocusGenotype(*this); }

public:
  /**
   * @brief Set a MonolocusGenotype.
   */
  void setMonolocusGenotype(
      size_t locusPosition,
      const MonolocusGenotypeInterface& monogen);

  /**
   * @brief Set a MonolocusGenotype by allele keys.
   *
   * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
   * @throw Exception if there is no key in allele_keys.
   */
  void setMonolocusGenotypeByAlleleKey(
      size_t locusPosition,
      const std::vector<size_t>& alleleKeys);

  /**
   * @brief Set a MonolocusGenotype by allele id.
   *
   * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
   * @throw AlleleNotFoundException if at least one of the id is not found in the LocusInfo.
   */
  void setMonolocusGenotypeByAlleleId(
      size_t locusPosition,
      const std::vector<std::string>& alleleId,
      const LocusInfo& locusInfo);

  /**
   * @brief Set a MonolocusGenotype as missing data.
   *
   * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
   */
  void setMonolocusGenotypeAsMissing(size_t locusPosition);

  /**
   * @brief Tell if a MonolocusGenotype is a missing data.
   *
   * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
   */
  bool isMonolocusGenotypeMissing(size_t locusPosition) const;

  /**
   * @brief Get a MonolocusGenotype.
   */
  const MonolocusGenotypeInterface& monolocusGenotype(size_t locusPosition) const;

  /**
   * @brief Count the number of loci.
   *
   * Return the size of _loci.
   */
  size_t size() const;

  /**
   * @brief Count the number of non missing MonolocusGenotype.
   */
  size_t countNonMissingLoci() const;

  /**
   * @brief Count the number of homozygous MonolocusGenotype.
   */
  size_t countHomozygousLoci() const;

  /**
   * @brief Count the number of heterozygous MonolocusGenotype.
   */
  size_t countHeterozygousLoci() const;
};
} // end of namespace bpp;

#endif // _MULTILOCUSGENOTYPE_H_
