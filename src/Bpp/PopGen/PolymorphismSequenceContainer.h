// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
  PolymorphismSequenceContainer(std::shared_ptr<const Alphabet> alpha) :
    VectorSiteContainer(alpha),
    ingroup_(std::vector<bool>()),
    count_(0),
    group_(0)
  {}

  /**
   * @brief Build a new empty PolymorphismSequenceContainer of given size.
   */
  PolymorphismSequenceContainer(size_t size, std::shared_ptr<const Alphabet> alpha) :
    VectorSiteContainer(size, alpha),
    ingroup_(size),
    count_(size),
    group_(size)
  {}

  /**
   * @brief Build a new empty PolymorphismSequenceContainer with given sequence names.
   */
  PolymorphismSequenceContainer(const std::vector<std::string>& names, std::shared_ptr<const Alphabet> alpha) :
    VectorSiteContainer(names, alpha),
    ingroup_(names.size()),
    count_(names.size()),
    group_(names.size())
  {}

  /**
   * @brief Build a PolymorphismSequenceContainer by copying data from a SequenceContainer.
   *
   * @param sc Sequence container to convert.
   */
  PolymorphismSequenceContainer(const SequenceContainerInterface& sc) :
    VectorSiteContainer(sc),
    ingroup_(sc.getNumberOfSequences(), true),
    count_(sc.getNumberOfSequences(), 1),
    group_(sc.getNumberOfSequences(), 1)
  {}

  /**
   * @brief Build a PolymorphismSequenceContainer by copying data from a SequenceContainer.
   *
   * @note In case of count = false, the constructor without additional argument will be more efficient.
   *
   * @param sc Sequence container to convert.
   * @param count Tell if identical sequences should be collapsed and counted.
   *              If not (the historical behavior), sequences are duplicated and stored with a frequency of 1.
   */
  PolymorphismSequenceContainer(const SequenceContainerInterface& sc, bool count);

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
  PolymorphismSequenceContainer* clone() const override
  {
    return new PolymorphismSequenceContainer(*this);
  }

public:
  // Other methods
  std::unique_ptr<Sequence> removeSequence(size_t sequencePosition) override;

  std::unique_ptr<Sequence> removeSequence(const std::string& sequenceKey)override;

  void deleteSequence(size_t sequencePosition) override;

  void deleteSequence(const std::string& sequenceKey) override;

  /**
   * @brief Add a sequence to the container.
   *
   * @throw AlphabetMismatchException if the sequence's alphabet doesn't match the container's alphabet.
   * @throw SequenceException if the sequence's size doesn't match the sequence's size of the container.
   * @throw SequenceException if the sequence's name already exists in the container.
   */
  void addSequenceWithFrequency(
      const std::string& sequenceKey,
      std::unique_ptr<Sequence>& sequence,
      unsigned int frequency)
  {
    VectorSiteContainer::addSequence(sequenceKey, sequence);
    count_.push_back(frequency);
    ingroup_.push_back(true);
    group_.push_back(0);
  }


  void insertSequenceWithFrequency(
      size_t sequencePosition,
      std::unique_ptr<Sequence>& sequence,
      const std::string& sequenceKey,
      unsigned int frequency)
  {
    VectorSiteContainer::insertSequence(sequencePosition, sequence, sequenceKey);
    count_.insert(count_.begin() + static_cast<ptrdiff_t>(sequencePosition), frequency);
    ingroup_.insert(ingroup_.begin() + static_cast<ptrdiff_t>(sequencePosition), true);
    group_.insert(group_.begin() + static_cast<ptrdiff_t>(sequencePosition), 0);
  }

  void addSequence(
      const std::string& sequenceKey,
      std::unique_ptr<Sequence>& sequence) override
  {
    addSequenceWithFrequency(sequenceKey, sequence, 1);
  }

  void insertSequence(size_t sequencePosition, std::unique_ptr<Sequence>& sequence, const std::string& sequenceKey) override
  {
    insertSequenceWithFrequency(sequencePosition, sequence, sequenceKey, 1);
  }

  /**
   * @brief Clear the container of all its sequences.
   */
  void clear() override
  {
    VectorSiteContainer::clear();
    count_.clear();
    ingroup_.clear();
    group_.clear();
  }

  /**
   * @brief Get the group identifier of the sequence.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
   */
  size_t getGroupId(size_t index) const
  {
    if (index >= getNumberOfSequences())
      throw IndexOutOfBoundsException("PolymorphismSequenceContainer::getGroupId: index out of bounds.", index, 0, getNumberOfSequences());
    return group_[index];
  }

  /**
   * @brief Get the group identifier of a sequence.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   */
  size_t getGroupId(const std::string& name) const
  {
    try
    {
      return group_[getSequencePosition(name)];
    }
    catch (SequenceNotFoundException& snfe)
    {
      throw SequenceNotFoundException("PolymorphismSequenceContainer::getGroupId.", name);
    }
  }


  /**
   * @brief Get all the groups identifiers.
   */
  std::set<size_t> getAllGroupsIds() const;

  /**
   * @brief Set the group identifier of a sequence.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
   */
  void setGroupId(size_t index, size_t group_id)
  {
    if (index >= getNumberOfSequences())
      throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setGroupId: index out of bounds.", index, 0, getNumberOfSequences());
    group_[index] = group_id;
  }

  /**
   * @brief Set the group identifier of a sequence.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   */
  void setGroupId(const std::string& name, size_t group_id)
  {
    try
    {
      group_[getSequencePosition(name)] = group_id;
    }
    catch (SequenceNotFoundException& snfe)
    {
      throw SequenceNotFoundException("PolymorphismSequenceContainer::setGroupId.", name);
    }
  }

  /**
   * @brief Get the number of groups.
   */
  size_t getNumberOfGroups() const
  {
    return getAllGroupsIds().size();
  }

  /**
   * @return True is the container contains at least one outgroup sequence.
   */
  bool hasOutgroup() const;

  /**
   * @brief Tell if the sequence is ingroup by index.
   *
   * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
   */
  bool isIngroupMember(size_t index) const
  {
    if (index >= getNumberOfSequences())
      throw IndexOutOfBoundsException("PolymorphismSequenceContainer::isIngroupMember: index out of bounds.", index, 0, getNumberOfSequences());
    return ingroup_[index];
  }

  /**
   * @brief Tell if a sequence is ingroup by name.
   *
   * @throw SequenceNotFoundException if name is not found among the sequences' names.
   */
  bool isIngroupMember(const std::string& name) const
  {
    try
    {
      return ingroup_[getSequencePosition(name)];
    }
    catch (SequenceNotFoundException& snfe)
    {
      throw SequenceNotFoundException("PolymorphismSequenceContainer::isIngroupMember.", name);
    }
  }

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
  std::unique_ptr<SiteContainerInterface> toSiteContainer() const;
};
} // end of namespace bpp;

#endif // _POLYMORPHISMSEQUENCECONTAINER_H_
