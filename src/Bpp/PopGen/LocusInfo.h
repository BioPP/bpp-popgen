// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _LOCUSINFO_H_
#define _LOCUSINFO_H_

// From STL
#include <string>
#include <vector>
#include <memory>

// From local bpp-popgen
#include "AlleleInfo.h"
#include "GeneralExceptions.h"

#include <Bpp/Exceptions.h>

namespace bpp
{
/**
 * @brief The LocusInfo class.
 *
 * This is an AlleleInfo container with additional data like a name,
 * the ploidy and some comments.
 *
 * @author Sylvain Gaillard
 */
class LocusInfo :
  public virtual Clonable
{
private:
  std::string name_;
  unsigned int ploidy_;
  std::vector<std::unique_ptr<AlleleInfo>> alleles_;

public:
  static unsigned int HAPLODIPLOID;
  static unsigned int HAPLOID;
  static unsigned int DIPLOID;
  static unsigned int UNKNOWN;

public:
  // Constructors and destructor
  /**
   * @brief Build a new LocusInfo object.
   *
   * @param name The name of the locus.
   * @param ploidy The ploidy of the locus.
   */
  LocusInfo(const std::string& name, const unsigned int ploidy = DIPLOID) :
    name_(name),
    ploidy_(ploidy),
    alleles_()
  {}

  /**
   * @brief Copy constructor.
   */
  LocusInfo(const LocusInfo& locusInfo) :
    name_(locusInfo.name_),
    ploidy_(locusInfo.ploidy_),
    alleles_(locusInfo.getNumberOfAlleles())
  {
    for (unsigned int i = 0; i < locusInfo.getNumberOfAlleles(); ++i)
    {
      alleles_[i].reset(locusInfo.getAlleleInfoByKey(i).clone());
    }
  }

  LocusInfo& operator=(const LocusInfo& locusInfo)
  {
    name_ = locusInfo.name_;
    ploidy_ = locusInfo.ploidy_;
    alleles_.resize(locusInfo.getNumberOfAlleles());
    for (unsigned int i = 0; i < locusInfo.getNumberOfAlleles(); ++i)
    {
      alleles_[i].reset(locusInfo.getAlleleInfoByKey(i).clone());
    }
    return *this;
  }

  /**
   * @brief Destroy the LocusInfo.
   */
  virtual ~LocusInfo() = default;

  LocusInfo* clone() const override { return new LocusInfo(*this); }

public:
  /**
   * @brief Get the name of the locus.
   */
  const std::string& getName() const { return name_; }

  /**
   * @brief Get the ploidy of the locus.
   *
   * @return The ploidy as an unsigned integer.
   */
  unsigned int getPloidy() const { return ploidy_; }

  /**
   * @brief Add an AlleleInfo to the LocusInfo.
   *
   * @throw BadIdentifierException if the AlleleInfo's id already exists.
   */
  void addAlleleInfo(const AlleleInfo& allele);

  /**
   * @brief Retrieve an AlleleInfo object of the LocusInfo.
   *
   * @throw AlleleNotFoundException if the id is not found.
   */
  const AlleleInfo& getAlleleInfoById(const std::string& id) const;

  /**
   * @brief Retrieve an AlleleInfo object of the LocusInfo.
   *
   * @throw IndexOutOfBoundsException if key exceeds the number of alleles.
   */
  const AlleleInfo& getAlleleInfoByKey(size_t key) const;

  /**
   * @brief Get the position of an AlleleInfo.
   *
   * @throw AlleleNotFoundException if the AlleleInfo's id is not found.
   */
  unsigned int getAlleleInfoKey(const std::string& id) const;

  /**
   * @brief Get the number of alleles at this locus.
   */
  size_t getNumberOfAlleles() const { return alleles_.size(); }

  /**
   * @brief Delete all alleles from the locus.
   */
  void clear() { alleles_.clear(); }
};
} // end of namespace bpp;

#endif // _LOCUSINFO_H_
