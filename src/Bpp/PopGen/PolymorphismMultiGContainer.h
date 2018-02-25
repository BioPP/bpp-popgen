//
// File PolymorphismMultiGContainer.h
// Author : Sylvain Gaillard
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
class PolymorphismMultiGContainer
{
private:
  std::vector<MultilocusGenotype*> multilocusGenotypes_;
  std::vector<size_t> groups_; // group id for each multilocusgenotype
  std::map<size_t, std::string> groups_names_;

public:
  // Constructors and destructor
  /**
   * @brief Build a new PolymorphismMultilocusGenotypeContainer.
   */
  PolymorphismMultiGContainer();

  /**
   * @brief The copy constructor.
   */
  PolymorphismMultiGContainer(const PolymorphismMultiGContainer& pmgc);

  /**
   * @brief Destroy a PolymorphismMultilocusGenotypeContainer.
   */
  ~PolymorphismMultiGContainer();

public:
  /**
   * @brief The assignation operator=.
   */
  PolymorphismMultiGContainer& operator=(const PolymorphismMultiGContainer& pmgc);

  /**
   * @brief Add a MultilocusGenotype to the container.
   */
  void addMultilocusGenotype(const MultilocusGenotype& mg, size_t group);

  /**
   * @brief Get a MultilocusGenotype at a position.
   *
   * @throw IndexOutOfBoundsException if position excedes the size of the container.
   */
  const MultilocusGenotype* getMultilocusGenotype(size_t position) const;

  /**
   * @brief Remove a MultilocusGenotype.
   *
   * @throw IndexOutOfBoundsException if position excedes the size of the container.
   */
  MultilocusGenotype* removeMultilocusGenotype(size_t position);

  /**
   * @brief Delete a MultilocusGenotype.
   *
   * @throw IndexOutOfBoundsException if position excedes the size of the container.
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
   * @throw IndexOutOfBoundsException if position excedes the size of the container.
   */
  size_t getGroupId(size_t position) const;

  /**
   * @brief Set the Group id of a MultilocusGenotype.
   *
   * @throw IndexOutOfBoundsException if position excedes the size of the container.
   */
  void setGroupId(size_t position, size_t group_id);

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
  size_t getNumberOfGroups() const;

  /**
   * @brief Get group size.
   */
  size_t getGroupSize(size_t group) const;

  /**
   * @brief Get the group name for a given group id or just the id if not available juste return it's id
   */
  std::string getGroupName(size_t group_id) const;

  /**
   * @brief Set the name for the given group id.
   */
  void setGroupName(size_t group_id, std::string name);

  /**
   * @brief Inserts a name for the given group id.
   */
  void addGroupName(size_t group_id, std::string name);

  /**
   * @brief Get the size of a group for a given locus.
   */
  size_t getLocusGroupSize(size_t group, size_t locus_position) const;

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
