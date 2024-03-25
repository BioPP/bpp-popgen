// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _POLYMORPHYSMMULTIGCONTAINER_H_
#define _POLYMORPHYSMMULTIGCONTAINER_H_

// From Utils
#include <Bpp/Clonable.h>
#include <Bpp/Exceptions.h>
#include <Bpp/Utils/MapTools.h>
#include <Bpp/Text/TextTools.h>

// From popgenlib
#include "MultilocusGenotype.h"
#include "GeneralExceptions.h"

// From STL
#include <string>
#include <vector>
#include <map>
#include <set>

namespace bpp
{
/**
 * @brief The PolymorphismMultiGContainer class
 *
 * This class is a container of MultilocusGenotype.
 *
 * @author Sylvain Gaillard
 */
class PolymorphismMultiGContainer :
  public virtual Clonable
{
private:
  std::vector<std::unique_ptr<MultilocusGenotype>> multilocusGenotypes_;
  std::vector<size_t> groups_; // group id for each multilocusgenotype
  std::map<size_t, std::string> groupsNames_;

public:
  // Constructors and destructor
  /**
   * @brief Build a new PolymorphismMultilocusGenotypeContainer.
   */
  PolymorphismMultiGContainer() :
    multilocusGenotypes_(),
    groups_(std::vector<size_t>()),
    groupsNames_(std::map<size_t, std::string>())
  {}

  /**
   * @brief The copy constructor.
   */
  PolymorphismMultiGContainer(const PolymorphismMultiGContainer& pmgc);

  /**
   * @brief Destroy a PolymorphismMultilocusGenotypeContainer.
   */
  virtual ~PolymorphismMultiGContainer()
  {
    clear();
  }

  PolymorphismMultiGContainer* clone() const override { return new PolymorphismMultiGContainer(*this); }

public:
  /**
   * @brief The assignation operator=.
   */
  PolymorphismMultiGContainer& operator=(const PolymorphismMultiGContainer& pmgc);

  /**
   * @brief Add a MultilocusGenotype to the container.
   */
  void addMultilocusGenotype(std::unique_ptr<MultilocusGenotype>& mg, size_t group);

  /**
   * @brief Get a MultilocusGenotype at a position.
   *
   * @throw IndexOutOfBoundsException if position exceeds the size of the container.
   */
  const MultilocusGenotype& multilocusGenotype(size_t position) const;

  /**
   * @brief Remove a MultilocusGenotype.
   *
   * @throw IndexOutOfBoundsException if position exceeds the size of the container.
   */
  std::unique_ptr<MultilocusGenotype> removeMultilocusGenotype(size_t position);

  /**
   * @brief Delete a MultilocusGenotype.
   *
   * @throw IndexOutOfBoundsException if position exceeds the size of the container.
   */
  void deleteMultilocusGenotype(size_t position);

  /**
   * @brief Tell if the MultilocusGenotypes are aligned (i.e. same size).
   */
  bool isAligned() const;

  /**
   * @brief Get the number of loci if the MultilocusGenotypes are aligned.
   *
   * @throw Exception if MultilocusGenotypes are not aligned.
   */
  size_t getNumberOfLoci() const;

  /**
   * @brief Get the Group id of a MultilocusGenotype.
   *
   * @throw IndexOutOfBoundsException if position exceeds the size of the container.
   */
  size_t getGroupId(size_t position) const;

  /**
   * @brief Set the Group id of a MultilocusGenotype.
   *
   * @throw IndexOutOfBoundsException if position exceeds the size of the container.
   */
  void setGroupId(size_t position, size_t groupId);

  /**
   * @brief Get the groups' ids.
   */
  std::set<size_t> getAllGroupsIds() const;

  /**
   * @brief Get the groups names or ids if not available
   */
  std::vector<std::string> getAllGroupsNames() const;

  /**
   * @brief Tell if a group exists.
   */
  bool groupExists(size_t group) const;

  /**
   * @brief Get the number of groups.
   */
  size_t getNumberOfGroups() const
  {
    return getAllGroupsIds().size();
  }

  /**
   * @brief Get group size.
   */
  size_t getGroupSize(size_t group) const;

  /**
   * @brief Get the group name for a given group id or just the id if not available juste return it's id
   */
  std::string getGroupName(size_t groupId) const;

  /**
   * @brief Set the name for the given group id.
   */
  void setGroupName(size_t groupId, const std::string& name);

  /**
   * @brief Inserts a name for the given group id.
   */
  void addGroupName(size_t groupId, const std::string& name);

  /**
   * @brief Get the size of a group for a given locus.
   */
  size_t getLocusGroupSize(size_t group, size_t locusPosition) const;

  /**
   * @brief Get the number of MultilocusGenotype.
   */
  size_t size() const;

  /**
   * @brief Clear the container.
   */
  void clear();
};
} // end of namespace bpp;

#endif // _POLYMORPHYSMMULTIGCONTAINER_H_
