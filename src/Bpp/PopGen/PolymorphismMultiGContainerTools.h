// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _POLYMORPHISMMULTIGCONTAINERTOOLS_H_
#define _POLYMORPHISMMULTIGCONTAINERTOOLS_H_

// From the STL
#include <set>

// From the PolGenLib library
#include "PolymorphismMultiGContainer.h"

#include <Bpp/Numeric/Random/RandomTools.h>

namespace bpp
{
/**
 * @brief Tools for PolymorphismMultiGContainer.
 *
 * Provides static methods for permutations.
 *
 * @author Sylvain Gaillard
 */
class PolymorphismMultiGContainerTools
{
public:
  /**
   * @brief Permut the MultilocusGenotype in the whole PolymorphismMultiGContainer.
   *
   * @param pmgc The PolymorphismMultiGContainer to permut.
   * @return A permuted PolymorphismMultiGContainer.
   */
  static std::unique_ptr<PolymorphismMultiGContainer> permuteMultiG(const PolymorphismMultiGContainer& pmgc);

  /**
   * @brief Permut the MonolocusGenotype.
   *
   * Permut the MonolocusGenotypes in one or several groups breaking
   * the links between them.
   *
   * @param pmgc The PolymorphismMultiGContainer to permut.
   * @param groups The groups ids between which the MonolocusGenotypes will be permuted.
   * @return A permuted PolymorphismMultiGContainer.
   */
  static std::unique_ptr<PolymorphismMultiGContainer> permuteMonoG(const PolymorphismMultiGContainer& pmgc, const std::set<size_t>& groups);

  /**
   * @brief Permut the MonolocusGenotype between individuals in the same group.
   *
   * Permut the MonolocusGenotypes for a set of groups. The idiv for the other groups
   * are kept intact
   *
   * @param pmgc The PolymorphismMultiGContainer to permut.
   * @param groups The groups ids for which the MonolocusGenotypes will be permuted.
   * @return A permuted PolymorphismMultiGContainer.
   */
  static std::unique_ptr<PolymorphismMultiGContainer> permuteIntraGroupMonoG(const PolymorphismMultiGContainer& pmgc, const std::set<size_t>& groups);

  /**
   * @brief Permut the Alleles.
   *
   * Permut the alleles in one or several groups breaking
   * the links between them.
   *
   * @param pmgc The PolymorphismMultiGContainer to permut.
   * @param groups The groups ids between which the MonolocusGenotypes will be permuted.
   * @return A permuted PolymorphismMultiGContainer.
   */
  static std::unique_ptr<PolymorphismMultiGContainer> permuteAlleles(const PolymorphismMultiGContainer& pmgc, const std::set<size_t>& groups);

  /**
   * @brief Permut the Alleles between individuals in the same group.
   *
   * Permut the alleles in one or several groups
   *
   * @param pmgc The PolymorphismMultiGContainer to permut.
   * @param groups The groups ids between which the MonolocusGenotypes will be permuted.
   * @return A permuted PolymorphismMultiGContainer.
   */
  static std::unique_ptr<PolymorphismMultiGContainer> permuteIntraGroupAlleles(const PolymorphismMultiGContainer& pmgc, const std::set<size_t>& groups);

  static std::unique_ptr<PolymorphismMultiGContainer> extractGroups(const PolymorphismMultiGContainer& pmgc, const std::set<size_t>& groups);
};
} // end of namespace bpp;

#endif // _POLYMORPHISMMULTIGCONTAINERTOOLS_H_
