// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _GROUP_H_
#define _GROUP_H_

// From STL
#include <vector>
#include <memory>

#include <Bpp/Exceptions.h>
#include <Bpp/Graphics/Point2D.h>

// From bpp-seq
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

// From bpp-popgen
#include "Individual.h"
#include "../GeneralExceptions.h"

namespace bpp
{
/**
 * @brief The Group class.
 *
 * A Group is an ensembl of Individuals with some statistics like the average
 * allele number.
 *
 * @author Sylvain Gaillard
 */
class Group :
  public virtual Clonable
{
protected:
  size_t id_;
  std::string name_;
  std::vector<std::unique_ptr<Individual>> individuals_;

public:
  // Constructors and destructor :
  /**
   * @brief Build a void new Group.
   */
  Group(size_t groupId) :
    id_(groupId),
    name_(""),
    individuals_()
  {}

  /**
   * @brief Copy constructor.
   *
   * If you need to use a copy constructor in a DataSet context, use the one
   * which specify a new Group Id.
   */
  Group(const Group& group) :
    id_(group.id_),
    name_(group.name_),
    individuals_()
  {
    for (size_t i = 0; i < group.getNumberOfIndividuals(); ++i)
    {
      addIndividual(group.getIndividualAtPosition(i));
    }
  }

  /**
   * @brief The assignation operator =.
   */
  Group& operator=(const Group& group)
  {
    id_ = group.id_;
    name_ = group.name_;
    individuals_.clear();
    for (size_t i = 0; i < group.getNumberOfIndividuals(); ++i)
    {
      addIndividual(group.getIndividualAtPosition(i));
    }
    return *this;
  }

  /**
   * @brief A duplication constructor with new Group Id.
   */
  Group(const Group& group, size_t groupId) :
    id_(groupId),
    name_(group.name_),
    individuals_()
  {
    for (size_t i = 0; i < group.getNumberOfIndividuals(); ++i)
    {
      addIndividual(group.getIndividualAtPosition(i));
    }
  }


  /**
   * @brief Destroy a Group.
   */
  virtual ~Group() = default;

  Group* clone() const override { return new Group(*this); }

public:
  /**
   * @brief Set the id of the Group.
   *
   * @param group_id The id of the Group as an size_t.
   */
  void setGroupId(size_t groupId){ id_ = groupId; }

  /**
   * @brief Get the name of the Group.
   *
   * @return The name of the Group as a string.
   */
  const std::string& getGroupName() const { return name_; }

  /**
   * @brief Set the name of the Group.
   *
   * @param group_name Name of the Group as string.
   */
  void setGroupName(const std::string& groupName) { name_ = groupName; }

  /**
   * @brief Get the id of the Group.
   *
   * @return The id of the Group as an size_t.
   */
  size_t getGroupId() const { return id_; }

  /**
   * @brief Add an Individual.
   *
   * Add an Individual to the group.
   *
   * @param ind The Individual to add to the Group.
   * @throw BadIdentifierException if individual's identifier is already in use.
   */
  void addIndividual(const Individual& ind);

  /**
   * @brief Add an empty Individual to the Group.
   *
   * @throw BadIdentifierException if individualId is already in use.
   */
  void addEmptyIndividual(const std::string& individualId);

  /**
   * @brief Get the number of Individual in the Group.
   *
   * @return An integer as the number of Individual.
   */
  size_t getNumberOfIndividuals() const { return individuals_.size(); }

  /**
   * @brief Get the maximum number of sequence.
   *
   * Give the value of the highest sequence key. This value is useful to
   * discover the missing sequences data for each individual.
   */
  size_t getMaxNumberOfSequences() const;

  /**
   * @brief Get the position of an Individual.
   *
   * @throw IndividualNotFoundException if individualId is not found.
   */
  size_t getIndividualPosition(const std::string& individualId) const;

  /**
   * @brief Get a reference to an Individual.
   *
   * @param individualId The id of the Individual to find.
   *
   * @return A pointer to the Individual or NULL if the Individual is not found.
   */
  const Individual& getIndividualById(const std::string& individualId) const;

  /**
   * @brief Get a reference to an Individual by its position.
   *
   * @param individualPosition The position of the Individual in the group.
   *
   * @return A pointer to the Individual.
   * @throw IndividualNotFoundException if individualId is not found.
   */
  const Individual& getIndividualAtPosition(size_t individualPosition) const;

  /**
   * @brief Remove an Individual from the Group.
   *
   * @param individualId The id of the Individual to remove.
   *
   * @return An std::auto_ptr to the removed Individual.
   * @throw IndividualNotFoundException if individualId is not found.
   *
   * Search an Individual in the Group by checking the id and remove it
   * if it is found then return a pointer to this Individual.
   */
  std::unique_ptr<Individual> removeIndividualById(const std::string& individualId);

  /**
   * @brief Remove an Individual from the Group.
   *
   * @param individualPosition The position in the Group of the Individual to remove.
   *
   * @return An std::auto_ptr to the removed Individual.
   *
   * Remove the individual at the specified position and return a pointer
   * to this Individual.
   */
  std::unique_ptr<Individual> removeIndividualAtPosition(size_t individualPosition);

  /**
   * @brief Delete an Individual from the Group.
   *
   * @param individualId The id of the Individual to delete.
   * @throw IndividualNotFoundException if individualId is not found.
   *
   * Search an Individual in the Group by checking the id and delete it
   * if it is foundi and free the memory by calling the destructor of the
   * Individual.
   */
  void deleteIndividualById(const std::string& individualId);

  /**
   * @brief Delete an Individual from the Group.
   *
   * @param individualPosition The position in the Group of the Individual to delete.
   *
   * Free the memory by calling the destructor of the Individual.
   */
  void deleteIndividualAtPosition(size_t individualPosition);

  /**
   * @brief Clear the Group.
   *
   * Delete all the Individuals of the group.
   */
  void clear() { individuals_.clear(); }

  // -- Dealing with Individuals -----------------------------
  /**
   * @brief Set the sex of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   */
  void setIndividualSexAtPosition(size_t individualPosition, const unsigned short sex);

  /**
   * @brief Get the sex of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   */
  unsigned short getIndividualSexAtPosition(size_t individualPosition) const;

  /**
   * @brief Set the date of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   */
  void setIndividualDateAtPosition(size_t individualPosition, const Date& date);

  /**
   * @brief Get the date of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if the Individual has no date.
   */
  const Date& getIndividualDateAtPosition(size_t individualPosition) const;

  /**
   * @brief Set the coordinates of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   */
  void setIndividualCoordAtPosition(size_t individualPosition, const Point2D<double>& coord);

  /**
   * @brief Get the coordinates of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if the individual has no coordinate.
   */
  const Point2D<double>& getIndividualCoordAtPosition(size_t individualPosition) const;

  /**
   * @brief Set the locality of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   */
  void setIndividualLocalityAtPosition(
      size_t individualPosition,
      std::shared_ptr<const Locality<double>> locality);

  /**
   * @brief Get the locality of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if the individual has no locality.
   */
  std::shared_ptr<const Locality<double>> getIndividualLocalityAtPosition(size_t individualPosition) const;

  /**
   * @brief Add a sequence to an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw AlphabetMismatchException if the sequence's alphabet doesn't match the container's alphabet.
   * @throw BadIdentifierException if the sequence's name is already in use.
   * @throw BadIntegerException if sequencePosition is already in use.
   */
  void addIndividualSequenceAtPosition(size_t individualPosition,
      size_t sequencePosition,
      std::unique_ptr<Sequence>& sequence);

  /**
   * @brief Get a sequence of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if there is no sequence container defined in the individual.
   * @throw SequenceNotFoundException if sequence_name is not found.
   */
  const Sequence& getIndividualSequenceByName(
      size_t individualPosition,
      const std::string& sequence_name) const;

  /**
   * @brief Get a sequence of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if there is no sequence container defined in the individual.
   * @throw SequenceNotFoundException if sequencePosition is not found.
   */
  const Sequence& getIndividualSequenceAtPosition(
      size_t individualPosition,
      size_t sequencePosition) const;

  /**
   * @brief Delete a sequence of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if there is no sequence container defined in the individual.
   * @throw SequenceNotFoundException if sequence_name is not found.
   */
  void deleteIndividualSequenceByName(size_t individualPosition, const std::string& sequence_name);

  /**
   * @brief Delete a sequence of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if there is no sequence container defined in the individual.
   * @throw SequenceNotFoundException if sequencePosition is not found.
   */
  void deleteIndividualSequenceAtPosition(size_t individualPosition, size_t sequencePosition);

  /**
   * @brief Tell if the Individual has some sequences.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   */
  bool hasIndividualSequences(size_t individualPosition) const;

  /**
   * @brief Get the sequences' names from an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if there is no sequence container defined in the individual.
   */
  std::vector<std::string> getIndividualSequencesNames(size_t individualPosition) const;

  /**
   * @brief Get the position of a sequence in an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if there is no sequence container defined in the individual.
   * @throw SequenceNotFoundException if sequence_name is not found.
   */
  size_t getIndividualSequencePosition(
      size_t individualPosition,
      const std::string& sequence_name) const;

  /**
   * @brief Get the number of sequences in an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if there is no sequence container defined in the individual.
   */
  size_t getIndividualNumberOfSequences(size_t individualPosition) const;

  /**
   * @brief Set all the sequences by copying an OrderedSequenceContainer.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   */
  void setIndividualSequences(size_t individualPosition, const SequenceContainerInterface& sc);

  /**
   * @brief Set the genotype of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   */
  void setIndividualGenotype(size_t individualPosition, const MultilocusGenotype& genotype);

  /**
   * @brief Initialize the genotype of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw BadIntegerException if loci_number < 1.
   * @throw Exception if the individual already has a genotype.
   */
  void initIndividualGenotype(size_t individualPosition, size_t loci_number);

  /**
   * @brief Delete the genotype of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   */
  void deleteIndividualGenotype(size_t individualPosition);

  /**
   * @brief Tell if an Individual has a genotype.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   */
  bool hasIndividualGenotype(size_t individualPosition) const;

  /**
   * @brief Set a MonolocusGenotype of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if the individual has no genotype.
   * @throw IndexOutOfBoundsException if locusPosition exceeds the number of locus.
   */
  void setIndividualMonolocusGenotype(
      size_t individualPosition,
      size_t locusPosition,
      const MonolocusGenotypeInterface& monogen);

  /**
   * @brief Set a MonolocusGenotype of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if the individual has no genotype.
   * @throw IndexOutOfBoundsException if locusPosition exceeds the number of locus.
   * @throw Exception if there is no key in allele_keys.
   */
  void setIndividualMonolocusGenotypeByAlleleKey(
      size_t individualPosition,
      size_t locusPosition,
      const std::vector<size_t>& alleleKeys);

  /**
   * @brief Set a MonolocusGenotype of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if the individual has no genotype.
   * @throw IndexOutOfBoundsException if locusPosition exceeds the number of locus.
   * @throw AlleleNotFoundException if at least one id is not found in locus_info.
   */
  void setIndividualMonolocusGenotypeByAlleleId(
      size_t individualPosition,
      size_t locusPosition,
      const std::vector<std::string>& alleleId,
      const LocusInfo& locusInfo);

  /**
   * @brief Get a MonolocusGenotype of an Individual.
   *
   * @throw IndexOutOfBoundsException if individualPosition exceeds the number of individuals.
   * @throw NullPointerException if the individual has no genotype.
   * @throw IndexOutOfBoundsException if locusPosition exceeds the number of locus.
   */
  const MonolocusGenotypeInterface& getIndividualMonolocusGenotype(
      size_t individualPosition,
      size_t locusPosition) const;

  /**
   * @brief Tell if at least one individual has at least one sequence.
   */
  bool hasSequenceData() const;

  /**
   * @brief Get the alphabet used for the sequences.
   */
  std::shared_ptr<const Alphabet> getAlphabet() const;

  /**
   * @brief Get the number of individual that have a data at the specified locus.
   */
  size_t getGroupSizeForLocus(size_t locusPosition) const;

  /**
   * @brief Get the number of individual that have a sequence at the specified position.
   */
  size_t getGroupSizeForSequence(size_t sequencePosition) const;
};
} // end of namespace bpp;

#endif // _GROUP_H_
