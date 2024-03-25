// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
 * Its instantiation requires a number of locus which is fixed
 * and can't be modified.
 *
 * @author Sylvain Gaillard
 */
class AnalyzedLoci :
  public virtual Clonable
{
private:
  std::vector<std::unique_ptr<LocusInfo>> loci_;

public:
  // Constructors and Destructor
  /**
   * @brief Build a void AnalyzedLoci with a specific number of loci.
   */
  AnalyzedLoci(size_t numberOfLoci) : loci_(numberOfLoci) {}

  /**
   * @brief Copy constructor.
   */
  AnalyzedLoci(const AnalyzedLoci& analyzedLoci) :
    loci_(analyzedLoci.loci_.size())
  {
    size_t i = 0;
    for (const auto& locus : analyzedLoci.loci_)
    {
      loci_[i++].reset(locus->clone());
    }
  }

  AnalyzedLoci& operator=(const AnalyzedLoci& analyzedLoci)
  {
    loci_.resize(analyzedLoci.loci_.size());
    size_t i = 0;
    for (const auto& locus : analyzedLoci.loci_)
    {
      loci_[i++].reset(locus->clone());
    }
    return *this;
  }

  /**
   * @brief Destroy the AnalyzedLoci.
   */
  virtual ~AnalyzedLoci() = default;

  AnalyzedLoci* clone() const { return new AnalyzedLoci(*this); }

public:
  // Other methods
  /**
   * @brief Set a LocusInfo.
   *
   * @throw IndexOutOfBoundsException if locus_position is out of bounds.
   */
  void setLocusInfo(size_t locusPosition, const LocusInfo& locus);

  /**
   * @brief Get the position of a LocusInfo.
   *
   * @throw BadIdentifierException if locus_name is not found.
   */
  size_t getLocusInfoPosition(const std::string& locusName) const;

  /**
   * @brief Get a LocusInfo by name.
   *
   * @throw BadIdentifierException if locus_name is not found.
   */
  const LocusInfo& getLocusInfoByName(const std::string& locusName) const;

  /**
   * @brief Get a LocusInfo by its position.
   *
   * @throw NullPointerException if the LocusInfo is not defined.
   * @throw IndexOutOfBoundsException if locus_position is out of bounds.
   */
  const LocusInfo& getLocusInfoAtPosition(size_t locusPosition) const;

  /**
   * @brief Add an AlleleInfo to a LocusInfo by LocusInfo name.
   *
   * @throw BadIdentifierException if the allele's id is already in use.
   * @throw LocusNotFoundException if locus_name is not found.
   */
  void addAlleleInfoByLocusName(const std::string& locusName,
      const AlleleInfo& allele);

  /**
   * @brief Add an AlleleInfo to a LocusInfo by its position.
   *
   * @throw BadIdentifierException if the allele's id is already in use.
   * @throw IndexOutOfBoundsException if locus_position is out of bounds.
   */
  void addAlleleInfoByLocusPosition(size_t locusosition,
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
  unsigned int getPloidyByLocusName(const std::string& locusName) const;

  /**
   * @brief Get the ploidy of a locus by its position.
   *
   * @throw IndexOutOfBoundsException if locus_position is out of bounds.
   */
  unsigned int getPloidyByLocusPosition(size_t locusPosition) const;
};
} // end of namespace bpp;

#endif // _ANALYZEDLOCI_H_
