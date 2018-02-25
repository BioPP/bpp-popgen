//
// File AnalyzedLoci.h
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

#ifndef _ANALYZEDLOCI_H_
#define _ANALYZEDLOCI_H_

// From STL
#include <vector>
#include <string>

#include <Bpp/Exceptions.h>

// From local
#include "../LocusInfo.h"
#include "../GeneralExceptions.h"

namespace bpp
{
/**
 * @brief The AnalyzedLoci class.
 *
 * This is a LocusInfo container.
 * Its instanciation requires a number of locus wich is fixed
 * and can't be modified.
 *
 * @author Sylvain Gaillard
 */
class AnalyzedLoci
{
private:
  std::vector<LocusInfo*> loci_;

public:
  // Constructors and Destructor
  /**
   * @brief Build a void AnalyzedLoci with a specific number of loci.
   */
  AnalyzedLoci(size_t number_of_loci);

  /**
   * @brief Copy constructor.
   */
  AnalyzedLoci(const AnalyzedLoci& analyzed_loci);

  /**
   * @brief Destroy the AnalyzedLoci.
   */
  ~AnalyzedLoci();

public:
  // Other methodes
  /**
   * @brief Set a LocusInfo.
   *
   * @throw IndexOutOfBoundsException if locus_position is out of bounds.
   */
  void setLocusInfo(size_t locus_position, const LocusInfo& locus);

  /**
   * @brief Get the position of a LocusInfo.
   *
   * @throw BadIdentifierException if locus_name is not found.
   */
  size_t getLocusInfoPosition(const std::string& locus_name) const;

  /**
   * @brief Get a LocusInfo by name.
   *
   * @throw BadIdentifierException if locus_name is not found.
   */
  const LocusInfo& getLocusInfoByName(const std::string& locus_name) const;

  /**
   * @brief Get a LocusInfo by its position.
   *
   * @throw NullPointerException if the LocusInfo is not difined.
   * @throw IndexOutOfBoundsException if locus_position is out of bounds.
   */
  const LocusInfo& getLocusInfoAtPosition(size_t locus_position) const;

  /**
   * @brief Add an AlleleInfo to a LocusInfo by LocusInfo name.
   *
   * @throw BadIdentifierException if the allele's id is already in use.
   * @throw LocusNotFoundException if locus_name is not found.
   */
  void addAlleleInfoByLocusName(const std::string& locus_name,
                                const AlleleInfo& allele);

  /**
   * @brief Add an AlleleInfo to a LocusInfo by its position.
   *
   * @throw BadIdentifierException if the allele's id is already in use.
   * @throw IndexOutOfBoundsException if locus_position is out of bounds.
   */
  void addAlleleInfoByLocusPosition(size_t locus_position,
                                    const AlleleInfo& allele);

  /**
   * @brief Get the number of loci.
   */
  size_t getNumberOfLoci() const { return loci_.size(); }

  /**
   * @brief Get the number of alleles at each locus.
   */
  std::vector<size_t> getNumberOfAlleles() const;

  /**
   * @brief Get the ploidy of a locus by name.
   *
   * @throw LocusNotFoundException if locus_name is not found.
   */
  unsigned int getPloidyByLocusName(const std::string& locus_name) const;

  /**
   * @brief Get the ploidy of a locus by its position.
   *
   * @throw IndexOutOfBoundsException if locus_position is out of bounds.
   */
  unsigned int getPloidyByLocusPosition(size_t locus_position) const;
};
} // end of namespace bpp;

#endif // _ANALYZEDLOCI_H_

