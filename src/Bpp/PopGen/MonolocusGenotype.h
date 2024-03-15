// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
class MonolocusGenotypeInterface :
  public virtual Clonable
{
public:

  MonolocusGenotypeInterface* clone() const override = 0;

public:
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

#endif// _MONOLOCUSGENOTYPE_H_
