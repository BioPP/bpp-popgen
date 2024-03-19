// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _INDIVIDUAL_H_
#define _INDIVIDUAL_H_

// From STL
#include <vector>
#include <memory>

#include <Bpp/Graphics/Point2D.h>
#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>

// From bpp-seq
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/SequenceExceptions.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>

// From bpp-popgen
#include "Locality.h"
#include "Date.h"
#include "../MultilocusGenotype.h"
#include "../GeneralExceptions.h"

namespace bpp
{
/**
 * @brief The Individual class.
 *
 * This class is designed to store data on a single individual.
 * This individual has only one sequence for each locus ... no information
 * about diploid sequence data.
 * See the no more in use MultiSeqIndividual documentation for an alternative.
 *
 * @author Sylvain Gaillard
 */
class Individual
{
protected:
  std::string id_;
  unsigned short sex_;
  std::unique_ptr<Date> date_;
  std::unique_ptr<Point2D<double>> coord_;
  std::shared_ptr<const Locality<double>> locality_;
  std::unique_ptr<VectorSequenceContainer> sequences_;
  std::unique_ptr<MultilocusGenotype> genotype_;

public:
  // Constructors and destructor :
  /**
   * @brief Build a void new Individual.
   */
  Individual();

  /**
   * @brief Build a new Individual with an identifier.
   */
  Individual(const std::string& id);

  /**
   * @brief Build a new Individual with parameters.
   *
   * @param id The id of the Individual as a string.
   * @param date The date of the Individual as a Date object.
   * @param coord The coordinates of the Individual as a Point2D object.
   * @param locality The locality of the Individual as a pointer to a Locality
   * object.
   * @param sex The sex of the Individual as an unsigned short.
   */
  Individual(const std::string& id,
      const Date& date,
      const Point2D<double>& coord,
      std::shared_ptr<Locality<double>> locality,
      const unsigned short sex);

  /**
   * @brief The Individual copy constructor.
   */
  Individual(const Individual& ind);

  /**
   * @brief Destroy an Individual.
   */
  virtual ~Individual();

public:
  // Methods
  /**
   * @brief The Individual copy operator.
   *
   * @return A ref toward the assigned Individual.
   * Make a copy of each atribute of the Individual.
   */
  Individual& operator=(const Individual& ind);

  /**
   * @brief Set the id of the Individual.
   *
   * @param id The id of the Individual as a string.
   */
  void setId(const std::string& id);

  /**
   * @brief Get the id of the Individual.
   *
   * @return The id of the Individual as a string.
   */
  const std::string& getId() const { return id_; }

  /**
   * @brief Set the sex of the Individual.
   *
   * @param sex An unsigned short coding for the sex.
   */
  void setSex(const unsigned short sex);

  /**
   * @brief Get the sex of the Individual.
   *
   * @return The sex of the Individual as an unsigned short.
   */
  unsigned short getSex() const { return sex_; }

  /**
   * @brief Set the date of the Individual.
   *
   * @param date The date as a Date object.
   */
  void setDate(const Date& date);

  /**
   * @brief Get the date of the Individual.
   *
   * @return A reference toward a Date object if the Individual has a date.
   * Otherwise throw a NullPointerException.
   */
  const Date& date() const;

  /**
   * @brief Tell if this Individual has a date.
   */
  bool hasDate() const;

  /**
   * @brief Set the coodinates of the Individual.
   *
   * @param coord A Point2D object.
   */
  void setCoord(const Point2D<double>& coord);

  /**
   * @brief Set the coordinates of the Individual.
   *
   * @param x The X coordinate as a double.
   * @param y The Y coordinate as a double.
   */
  void setCoord(const double x, const double y);

  /**
   * @brief Get the coordinates of the Induvidual.
   *
   * @return A pointer toward a Point2D object if the Individual has
   * coordinates. Otherwise throw a NullPointerException.
   */
  const Point2D<double>& coord() const;

  /**
   * @brief Tell if this Individual has coordinates.
   */
  bool hasCoord() const;

  /**
   * @brief Set the X coordinate of the Individual.
   *
   * @param x The X coordinate as a double.
   *
   * Set the X coordinate if the Individual has coordinates.
   * Otherwise throw a NullPointerException.
   */
  void setX(const double x);

  /**
   * @brief Set the Y coordinate of th Individual.
   *
   * @param y The Y coordinate as a double.
   *
   * Set the Y coordinate if the Individual has coordinates.
   * Otherwise throw a NullPointerException.
   */
  void setY(const double y);

  /**
   * @brief Get the X coordinate of the Individual.
   *
   * @return The X coordinate as a double if the Individual has coordinates.
   * Otherwise throw a NullPointerException.
   */
  double getX() const;

  /**
   * @brief Get the Y coordinate of the Individual.
   *
   * @return The Y coordinate as a double if the Individual has coordinates.
   * Otherwise throw a NullPointerException.
   */
  double getY() const;

  /**
   * @brief Set the locality of the Individual.
   *
   * @param locality A pointer to a Locality object.
   */
  void setLocality(std::shared_ptr<const Locality<double>> locality)
  {
    locality_ = locality;
  }

  /**
   * @brief Get the locality of the Individual.
   *
   * @return A pointer to the Locality of the Individual.
   */
  std::shared_ptr<const Locality<double>> getLocality() const;

  /**
   * @brief Get the locality of the Individual.
   *
   * @return A pointer to the Locality of the Individual.
   */
  const Locality<double>& locality() const;

  /**
   * @brief Tell if this Individual has a locality.
   */
  bool hasLocality() const
  {
    return locality_ != nullptr;
  }

  /**
   * @brief Add a sequence to the Individual.
   *
   * Creates the sequence container when adding the first sequence.
   * Otherwize add the sequence to the end of the sequence container.
   *
   * @param sequenceKey the place where the sequence will be put.
   * @param sequence The sequence to add.
   * @throw AlphabetMismatchException if the sequence's alphabet doesn't match the container's alphabet.
   * @throw BadIdentifierException if sequence's name is already in use.
   * @throw BadIntegerException if sequence_position is already in use.
   */
  void addSequence(size_t sequenceKey, std::unique_ptr<Sequence>& sequence);

  /**
   * @brief Get a sequence from its name.
   *
   * @param sequenceName The name of the sequence.
   * @return A reference to the sequence.
   * @throw NullPointerException if there is no sequence container defined.
   * @throw SequenceNotFoundException if sequence_name is not found.
   */
  const Sequence& sequenceByName(const std::string& sequenceName) const;

  /**
   * @brief Get a sequence by its position.
   *
   * @param sequencePosition The position of the sequence in the sequence set.
   * @return A reference to the sequence.
   * @throw NullPointerException if there is no sequence container defined.
   * @throw SequenceNotFoundException if sequence_position is not found (i.e. missing data or not used).
   */
  const Sequence& sequenceAtPosition(size_t sequencePosition) const;

  /**
   * @brief Delete a sequence.
   *
   * @param sequenceName The name of the sequence.
   * @throw NullPointerException if there is no sequence container defined.
   * @throw SequenceNotFoundException if sequence_name is not found.
   */
  void deleteSequenceByName(const std::string& sequenceName);

  /**
   * @brief Delete a sequence.
   *
   * @param sequencePosition The position of the sequence.
   * @throw NullPointerException if there is no sequence container defined.
   * @throw SequenceNotFoundException if sequence_postion is not found.
   */
  void deleteSequenceAtPosition(size_t sequencePosition);

  /**
   * @brief Tell if the Individual has some sequences.
   *
   * @return TRUE if the individual has at least one sequence.
   * @return FALSE if the container is empty or undifined.
   */
  bool hasSequences() const;

  /**
   * @brief Tell if the Individual has a sequence at a given position.
   */
  bool hasSequenceAtPosition(size_t position) const;

  /**
   * @brief Return the alphabet of the sequences.
   *
   * @throw NullPointerException if there is no sequence container defined.
   */
  std::shared_ptr<const Alphabet> getSequenceAlphabet() const;

  /**
   * @brief Return the alphabet of the sequences.
   *
   * @throw NullPointerException if there is no sequence container defined.
   */
  const Alphabet& sequenceAlphabet() const;

  /**
   * @brief Get the sequence names.
   *
   * @return All the sequences' names of the individual in a vector of string.
   * @throw NullPointerException if there is no sequence container defined.
   */
  std::vector<std::string> getSequenceNames() const
  {
    if (!sequences_)
      throw NullPointerException("Individual::getSequencesNames: no sequence data.");
    return sequences_->getSequenceNames();
  }

  /**
   * @brief Get the sequences' positions.
   *
   * @return All the positions where a sequence is found.
   * @throw NullPointerException if there is no sequence container defined.
   */
  std::vector<size_t> getSequencePositions() const;

  /**
   * @brief Get the position of a sequence.
   *
   * @throw NullPointerException if there is no sequence container defined.
   * @throw SequenceNotFoundException if sequence_name is not found.
   */
  size_t getSequencePosition(const std::string& sequenceKey) const;

  /**
   * @brief Get the number of sequences.
   */
  size_t getNumberOfSequences() const
  {
    return sequences_ ? sequences_->getNumberOfSequences() : 0;
  }


  /**
   * @brief Set all the sequences with a SequenceContainer. The container will be copied.
   */
  void setSequences(const SequenceContainerInterface& sc)
  {
    sequences_.reset(new VectorSequenceContainer(sc));
  }

  /**
   * @brief Get a reference to the sequence container.
   *
   * @throw NullPointerException if there is no sequence container defined.
   */
  const SequenceContainerInterface& sequences() const
  {
    if (!sequences_)
      throw NullPointerException("Individual::getSequences: no sequence data.");
    return *sequences_;
  }


  /**
   * @brief Set a genotype.
   *
   * @param genotype The MultilocusGenotype which will be copied.
   */
  void setGenotype(const MultilocusGenotype& genotype);

  /**
   * @brief Init the genotype.
   *
   * @throw Exception if the Individual already has a Genotype.
   * @throw BadIntegerException if loci_number < 1.
   */
  void initGenotype(size_t lociNumber);

  /**
   * @brief Get the genotype.
   */
  const MultilocusGenotype& getGenotype() const;

  /**
   * @brief Delete the genotype of the individual.
   */
  void deleteGenotype();

  /**
   * @brief Tell if the Individual has a MultilocusGenotype.
   */
  bool hasGenotype() const;

  /**
   * @brief Set a MonolocusGenotype.
   *
   * @throw NullPointerException if there is no genotype defined.
   * @throw IndexOutOfBoundsException if locusPosition excedes the number of loci.
   */
  void setMonolocusGenotype(size_t locusPosition, const MonolocusGenotypeInterface& monogen);

  /**
   * @brief Set a MonolocusGenotype.
   *
   * @throw NullPointerException if there is no genotype defined.
   * @throw IndexOutOfBoundsException if locusPosition excedes the number of loci.
   * @throw Exception if there is no key in alleleKeys.
   */
  void setMonolocusGenotypeByAlleleKey(
      size_t locusPosition,
      const std::vector<size_t> alleleKeys);

  /**
   * @brief Set a MonolocusGenotype.
   *
   * @throw NullPointerException if there is no genotype defined.
   * @throw IndexOutOfBoundsException if locusPosition excedes the number of loci.
   * @throw AlleleNotFoundException if at least one the id is not found in the LocusInfo.
   */
  void setMonolocusGenotypeByAlleleId(
      size_t locusPosition,
      const std::vector<std::string> alleleId,
      const LocusInfo& locusInfo);

  /**
   * @brief Get a MonolocusGenotype.
   *
   * @throw NullPointerException if there is no genotype defined.
   * @throw IndexOutOfBoundsException if locusPosition excedes the number of loci.
   */
  const MonolocusGenotypeInterface& getMonolocusGenotype(size_t locusPosition);

  /**
   * @brief Count the number of non missing MonolocusGenotype.
   *
   * @throw NullPointerException if there is no genotype defined.
   */
  size_t countNonMissingLoci() const;

  /**
   * @brief Count the number of homozygous MonolocusGenotype.
   *
   * @throw NullPointerException if there is no genotype defined.
   */
  size_t countHomozygousLoci() const;

  /**
   * @brief Count the number of heterozygous MonolocusGenotype.
   *
   * @throw NullPointerException if there is no genotype defined.
   */
  size_t countHeterozygousLoci() const;
};
} // end of namespace bpp;

#endif // _INDIVIDUAL_H_
