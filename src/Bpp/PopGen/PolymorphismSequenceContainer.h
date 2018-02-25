//
// File: PolymorphismSequenceContainer.h
// Authors: Eric Bazin
//          Sylvain Gaillard
// Created on: Wednesday August 04 2004
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

#ifndef _POLYMORPHISMSEQUENCECONTAINER_H_
#define _POLYMORPHISMSEQUENCECONTAINER_H_

#include <set>
#include <string>

#include <Bpp/Clonable.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>

#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

/**
 * @mainpage
 *
 * @par
 * The bpp-popgen library provides classes for population genetics analysis.
 * It makes intensive use of the bpp-seq library, and adds a dedicated container
 * named bpp::PolymorphismSequenceContainer, which associates frequencies to the
 * sequences in the set. The bpp::PolymorphismSequenceContainerTools and
 * bpp::SequenceStatistics static classes provide several tools for data analysis,
 * including diversity indices and positive selection tests.
 *
 * @section dataset Population and sample data storage and manipulation
 *
 * @par
 * PopGenLib library provides data structure for handling sample and data sets
 * for population genetics.
 * These objects are embedded in the bpp::DataSet object which is a container of bpp::Group
 * of bpp::Individual.
 * Each bpp::Individual can store bpp::Sequence data or allelic data with the dedicated
 * classes bpp::MultilocusGenotype.
 *
 * @section genetics Population genetics data and statistics
 *
 * @par
 * To compute statistics on data, two containers families are provided, one for sequences
 * (bpp::PolymorphismSequenceContainer) and the other for allelic data (bpp::PolymorphismMultiGContainer).
 * Static tools class for both families are provided to compute several common or less
 * common statistics.
 *
 * @section statistics Statistics overview
 *
 * @par heterozygosity
 * @par watterson75 Diversity estimator Theta of Watterson
 * @par tajima83 Diversity estimator Theta of Tajima
 * @par DVH Haplotype diversity of Depaulis and Veuille
 * @par D Tajima's D test
 */

namespace bpp
{
/**
 * @brief The PolymorphismSequenceContainer class.
 *
 * This is a VectorSiteContainer with effectif for each sequence.
 * It also has flag for ingroup and outgroup.
 *
 * @author Sylvain Gaillard
 */
class PolymorphismSequenceContainer :
  public VectorSiteContainer
{
private:
  std::vector<bool> ingroup_;
  std::vector<unsigned int> count_;
  std::vector<size_t> group_;

public:
  // Constructors and destructor
  /**
   * @brief Build a new empty PolymorphismSequenceContainer.
   */
  PolymorphismSequenceContainer(const Alphabet* alpha);

  /**
   * @brief Build a new empty PolymorphismSequenceContainer of given size.
   */
  PolymorphismSequenceContainer(size_t size, const Alphabet* alpha);

  /**
   * @brief Build a new empty PolymorphismSequenceContainer with given sequence names.
   */
  PolymorphismSequenceContainer(const std::vector<std::string>& names, const Alphabet* alpha);

  /**
   * @brief Build a PolymorphismSequenceContainer by copying data from an OrderedSequenceContainer.
   *
   * @param sc Sequence container to convert.
   */
  PolymorphismSequenceContainer(const OrderedSequenceContainer& sc);

  /**
   * @brief Build a PolymorphismSequenceContainer by copying data from an OrderedSequenceContainer.
   *
   * @note In case of count = false, the constructor with additional argument will be more efficient.
   *
   * @param sc Sequence container to convert.
   * @param count Tell if identical sequences should be collapsed and counted.
   *              If not (the historical behavior), sequences are duplicated and stored with a frequency of 1.
   */
  PolymorphismSequenceContainer(const OrderedSequenceContainer& sc, bool count);

  /**
   * @brief Build a PolymorphismSequenceContainer by copying data from a SiteContainer.
   *
   * @param sc Sequence container to convert.
   */
  PolymorphismSequenceContainer(const SiteContainer& sc);

  /**
   * @brief Build a PolymorphismSequenceContainer by copying data from a SiteContainer.
   *
   * @note In case of count = false, the constructor with additional argument will be more efficient.
   *
   * @param sc Sequence container to convert.
   * @param count Tell if identical sequences should be collapsed and counted.
   *              If not (the historical behavior), sequences are duplicated and stored with a frequency of 1.
   */
  PolymorphismSequenceContainer(const SiteContainer& sc, bool count);

  /**
   * @brief Copy constructor.
   */
  PolymorphismSequenceContainer(const PolymorphismSequenceContainer& psc);

  /**
   * @brief Operator= : copy operator.
   */
  PolymorphismSequenceContainer& operator=(const PolymorphismSequenceContainer& psc);

  /**
   * @brief Destroy a PolymorphismSequenceContainer.
   */
  virtual ~PolymorphismSequenceContainer();

  /**
   * @brief Clone a PolymorphismSequenceContainer.
   */
  PolymorphismSequenceContainer* clone() const
  {
    return new PolymorphismSequenceContainer(*this);
  }

public:
  // Other methods
  /**
   * @brief Remove a sequence by index and return a pointer to this removed sequence.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences.
   */
  Sequence* removeSequence(size_t index);

  /**
   * @brief Remove a sequence by name and return a pointer to this removed sequence.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   */
  Sequence* removeSequence(const std::string& name);

  /**
   * @brief Delete a sequence by index.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences.
   */
  void deleteSequence(size_t index);

  /**
   * @brief Delete a sequence by name.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   */
  void deleteSequence(const std::string& name);

  /**
   * @brief Add a sequence to the container.
   *
   * @throw AlphabetMismatchException if the sequence's alphabet doesn't match the container's alphabet.
   * @throw SequenceException if the sequence's size doesn't match the sequence's size of the container.
   * @throw SequenceException if the sequence's name already exists in the container.
   */
  void addSequenceWithFrequency(const Sequence& sequence, unsigned int frequency, bool checkName = true);
  void addSequenceWithFrequency(const Sequence& sequence, size_t sequenceIndex, unsigned int frequency, bool checkName = true);
  void addSequence(const Sequence& sequence, bool checkName = true) {
    addSequenceWithFrequency(sequence, 1, checkName);
  }
  void addSequence(const Sequence& sequence, size_t sequenceIndex, bool checkName = true) {
    addSequenceWithFrequency(sequence, sequenceIndex, 1, checkName);
  }

  /**
   * @brief Clear the container of all its sequences.
   */
  void clear();

  /**
   * @brief Get the group identifier of the sequence.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
   */
  size_t getGroupId(size_t index) const;

  /**
   * @brief Get the group identifier of a sequence.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   */
  size_t getGroupId(const std::string& name) const;

  /**
   * @brief Get all the groups identifiers.
   */
  std::set<size_t> getAllGroupsIds() const;

  /**
   * @brief Set the group identifier of a sequence.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
   */
  void setGroupId(size_t index, size_t group_id);

  /**
   * @brief Set the group identifier of a sequence.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   */
  void setGroupId(const std::string& name, size_t group_id);

  /**
   * @brief Get the number of groups.
   */
  size_t getNumberOfGroups() const;

  /**
   * @return True is the container contains at least one outgroup sequence.
   */
  bool hasOutgroup() const;

  /**
   * @brief Tell if the sequence is ingroup by index.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
   */
  bool isIngroupMember(size_t index) const;

  /**
   * @brief Tell if a sequence is ingroup by name.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   */
  bool isIngroupMember(const std::string& name) const;

  /**
   * @brief Set a sequence as ingroup member by index.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
   */
  void setAsIngroupMember(size_t index);

  /**
   * @brief Set a sequence as ingroup member by name.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   */
  void setAsIngroupMember(const std::string& name);

  /**
   * @brief Set a sequence as outgroup member by index.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
   */
  void setAsOutgroupMember(size_t index);

  /**
   * @brief Set a sequence as outgroup member by name.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   */
  void setAsOutgroupMember(const std::string& name);

  /**
   * @brief Set the count of a sequence by index.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
   * @throw BadIntegerException if count < 1 ... use deleteSequence instead of setting the count to 0.
   */
  void setSequenceCount(size_t index, unsigned int count);

  /**
   * @brief Set the count of a sequence by name.
   *
   * @throw throw SequenceNotFoundException if name is not found among the sequences' names.
   * @throw BadIntegerException if count < 1 ... use deleteSequence instead of setting the count to 0.
   */
  void setSequenceCount(const std::string& name, unsigned int count);

  /**
   * @brief Add 1 to the sequence count.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
   */
  void incrementSequenceCount(size_t index);

  /**
   * @brief Add 1 to the sequence count.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   */
  void incrementSequenceCount(const std::string& name);

   /**
   * @brief Removz 1 to the sequence count.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
   * @throw BadIntegerException if count < 1 ... use deleteSequence instead of setting the count to 0.
   */
  void decrementSequenceCount(size_t index);

  /**
   * @brief Remove 1 to the sequence count.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   * @throw BadIntegerException if count < 1 ... use deleteSequence instead of setting the count to 0.
   */
  void decrementSequenceCount(const std::string& name);

  /**
   * @brief Get the count of a sequence by index.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
   */
  unsigned int getSequenceCount(size_t index) const;

  /**
   * @brief Get the count of a sequence by name.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   */
  unsigned int getSequenceCount(const std::string& name) const;

  /**
   * @brief convert the container to a site container, with sequences dulicated according to their respective frequencies.
   *
   * @return A SiteContainer object, eventually with duplicated sequences. Names of duplicated sequences are happended with _1, _2, etc.
   */
  SiteContainer* toSiteContainer() const;
  
};
} // end of namespace bpp;

#endif  // _POLYMORPHISMSEQUENCECONTAINER_H_

