// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// Secured inclusion of header's file
#ifndef _MonolocusGenotypeTools_h_
#define _MonolocusGenotypeTools_h_

// From STL
#include <vector>
#include <memory>

#include <Bpp/Exceptions.h>

// From Pop
#include "MonolocusGenotype.h"

namespace bpp
{
/**
 * @brief The MonolocusGenotypeTools static class.
 *
 * This class provides tools for MonolocusGenotype manipulation or creation.
 *
 * @author Sylvain Gaillard
 */
class MonolocusGenotypeTools
{
public:
  /**
   * @brief Build a proper MonolocusGenotype accordig to the number of alleles.
   *
   * Return a MonolocusGenotype build according to the number of allels.
   * If one allele key, send a MonoAlleleMonolocusGenotype,
   * if two allele keys, send a BiAlleleMonolocusGenotype,
   * if more allele keys, send a MultiAlleleMonolocusGenotype.
   *
   * @param allele_keys A vector containing thes allele keys to put in the MonolocusGenotype.
   * @return A MonolocusGenotype according to the number of alleles
   */
  static std::unique_ptr<MonolocusGenotypeInterface> buildMonolocusGenotypeByAlleleKey(const std::vector<size_t> alleleKeys);
};
} // end of namespace bpp;

#endif// _MonolocusGenotypeTools_h_
