//
// File MultilocusGenotype.h
// Author : Sylvain Gaillard
// Last modification : April 4, 2008
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

#ifndef _MULTILOCUSGENOTYPE_H_
#define _MULTILOCUSGENOTYPE_H_

// From STL
#include <vector>
#include <string>

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
class MultilocusGenotype
{
private:
  std::vector<MonolocusGenotype*> loci_;

public:
  // Constructors and Destructor
  /**
   * @brief Build a MultilocusGenotype linked to an AnalyzedLoci object.
   *
   * @throw BadIntegerException if loci_number < 1.
   */
  MultilocusGenotype(size_t loci_number);

  /**
   * @brief Copy constructor.
   */
  MultilocusGenotype(const MultilocusGenotype& genotype);

  /**
   * @brief Destroy a MultilocusGenotype.
   */
  ~MultilocusGenotype();

public:
  /**
   * @brief Set a MonolocusGenotype.
   */
  void setMonolocusGenotype(size_t locus_position,
                            const MonolocusGenotype& monogen);

  /**
   * @brief Set a MonolocusGenotype by allele keys.
   *
   * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
   * @throw Exception if there is no key in allele_keys.
   */
  void setMonolocusGenotypeByAlleleKey(size_t locus_position,
                                       const std::vector<size_t>& allele_keys);

  /**
   * @brief Set a MonolocusGenotype by allele id.
   *
   * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
   * @throw AlleleNotFoundException if at least one of the id is not found in the LocusInfo.
   */
  void setMonolocusGenotypeByAlleleId(size_t locus_position,
                                      const std::vector<std::string>& allele_id,
                                      const LocusInfo& locus_info);

  /**
   * @brief Set a MonolocusGenotype as missing data.
   *
   * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
   */
  void setMonolocusGenotypeAsMissing(size_t locus_position);

  /**
   * @brief Tell if a MonolocusGenotype is a missing data.
   *
   * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
   */
  bool isMonolocusGenotypeMissing(size_t locus_position) const;

  /**
   * @brief Get a MonolocusGenotype.
   */
  const MonolocusGenotype& getMonolocusGenotype(size_t locus_position) const;

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

#endif// _MULTILOCUSGENOTYPE_H_
