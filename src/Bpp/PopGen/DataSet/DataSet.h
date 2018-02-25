//
// File DataSet.h
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

#ifndef _DATASET_H_
#define _DATASET_H_

// From the STL
#include <algorithm>
#include <vector>
#include <map>
#include <string>

#include <Bpp/Exceptions.h>
#include <Bpp/Graphics/Point2D.h>
#include <Bpp/Utils/MapTools.h>

#include "Group.h"
#include "Individual.h"
#include "Locality.h"
#include "../GeneralExceptions.h"
#include "AnalyzedLoci.h"
#include "AnalyzedSequences.h"
#include "../PolymorphismMultiGContainer.h"
#include "../PolymorphismSequenceContainer.h"

namespace bpp
{
/**
 * @brief The DataSet class.
 *
 * A DataSet the object that manage every data on which one can compute
 * some statistics.
 *
 * @author Sylvain Gaillard
 */
class DataSet
{
private:
  AnalyzedLoci* analyzedLoci_;
  AnalyzedSequences* analyzedSequences_;
  std::vector<Locality<double>*> localities_;
  std::vector<Group*> groups_;

public:
  // Constructor and destructor
  /**
   * @brief Build a new void DataSet.
   */
  DataSet();

  /**
   * @brief Destroy a DataSet.
   */
  ~DataSet();

  /**
   * @brief Copy constructor.
   */
  DataSet(const DataSet& ds);

  DataSet& operator=(const DataSet& ds);

public:
  // Methodes
// ** Locality manipulation ***************************************************/
  /**
   * @brief Add a locality to the DataSet.
   *
   * @param locality A Locality object.
   * @throw BadIdentifierException if the locality's name already exists.
   */
  void addLocality(Locality<double>& locality);

  /**
   * @brief Get the position of a locality in the container.
   *
   * @return The locality_position (position) of the Locality.
   * @param name The locality's name to find.
   * @throw LocalityNotFoundException if the locality's name doesn't match any name in the DataSet.
   */
  size_t getLocalityPosition(const std::string& name) const;

  /**
   * @brief Get a Locality by locality_position.
   *
   * @return A const pointer to the locality matching the locality_position.
   * @param locality_position The position of the Locality in the DataSet.
   * @throw IndexOutOfBoundsException if locality_position excedes the number of locality of the DataSet.
   */
  const Locality<double>& getLocalityAtPosition(size_t locality_position) const;

  /**
   * @brief Get a Locality by name.
   *
   * @throw LocalityNotFoundException if the locality's name is not found.
   */
  const Locality<double>& getLocalityByName(const std::string& name) const;

  /**
   * @brief Delete a Locality from the DataSet.
   *
   * @throw IndexOutOfBoundsException if locality_position excedes the number of Locality.
   */
  void deleteLocalityAtPosition(size_t locality_position);

  /**
   * @brief Delete a Locality from the DataSet.
   *
   * @throw LocalityNotFoundException if the locality's name is not found.
   */
  void deleteLocalityByName(const std::string& name);

  /**
   * @brief Get the number of Localities.
   */
  size_t getNumberOfLocalities() const;

  /**
   * @brief Tell if there is at least one locality.
   */
  bool hasLocality() const;

  // ** Group manipulation ******************************************************/
  /**
   * @brief Add a Group to the DataSet.
   *
   * Add a Group to the DataSet.
   *
   * @param group A pointer to the Group to add.
   */
  void addGroup(const Group& group);

  /**
   * @brief Add an empty Group to the DataSet.
   */
  void addEmptyGroup(size_t group_id);

  /**
   * @brief Get a group by identifier.
   */
  const Group& getGroupById(size_t group_id) const;

  /**
   * @brief Get the position of a Group.
   *
   * @throw GroupNotFoundException if the group_id is not found.
   */
  size_t getGroupPosition(size_t group_id) const;

  /**
   * @brief Get the name of a Group. If the name is an empty string it just returns the group_id
   *
   * @throw GroupNotFoundException if the group_id is not found.
   */
  std::string getGroupName(size_t group_id) const;
  /**
   * @brief set the name of a Group.
   *
   * @throw GroupNotFoundException if the group_id is not found.
   */
  void setGroupName(size_t group_id, const std::string& group_name) const;

  /**
   * @brief Get a group by position.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   */
  const Group& getGroupAtPosition(size_t group_position) const;

  /**
   * @brief Delete a Group from the DataSet.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   */
  void deleteGroupAtPosition(size_t group_position);

  /**
   * @brief Get the number of Groups.
   */
  size_t getNumberOfGroups() const;

  /**
   * @brief Merge two groups.
   *
   * This methode merge two groups. The source group is emptied into the target
   * and then is deleted.
   */
  void mergeTwoGroups(size_t source_id, size_t target_id);

  /**
   * @brief Merge some Groups in one.
   *
   * Merge all the groups which are specified in the first one (smallest
   * identifier). When a group is merged to the first, it is deleted from the
   * DataSet.
   *
   * @param group_ids A vector size_t listing the id of groups to merge.
   * @throw IndexOutOfBoundsException if one of the int in groups excedes the number of groups.
   */
  void mergeGroups(std::vector<size_t>& group_ids);

  /**
   * @brief Split a group in two.
   *
   * @param group_id The identifier of the source group.
   * @param individuals_selection The positions of the Individuals to extract from the group to make the new group.
   * @throw GroupNotFoundException if the group_id is not found.
   * @throw IndexOutOfBoundsException if one position of the selection excedes the number of individuals of the group.
   */
  void splitGroup(size_t group_id, std::vector<size_t> individuals_selection);

  // ** Individuals manipulation ************************************************/
  /**
   * @brief Add an Individual to a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw BadIdentifierException if the individual's id is already in use.
   */
  void addIndividualToGroup(size_t group_position, const Individual& individual);

  /**
   * @brief Add an empty Individual to a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw BadIdentifierException if the individual's id is already in use.
   */
  void addEmptyIndividualToGroup(size_t group_position, const std::string& individual_id);

  /**
   * @brief Get the number of Individuals in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   */
  size_t getNumberOfIndividualsInGroup(size_t group_position) const;

  /**
   * @brief Get the position of an Individual in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndividualNotFoundException if individual_id is not found.
   */
  size_t getIndividualPositionInGroup(size_t group_position, const std::string& individual_id) const;
  
  /**
   * @brief Get an Individual from a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   */
  const Individual* getIndividualAtPositionFromGroup(size_t group_position, size_t individual_position) const;

  /**
   * @brief Get an Individual from a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndividualNotFoundException if individual_id is not found.
   */
  const Individual* getIndividualByIdFromGroup(size_t group_position, const std::string& individual_id) const;

  /**
   * @brief Delete an Individual from a group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   */
  void deleteIndividualAtPositionFromGroup(size_t group_position, size_t individual_position);

  /**
   * @brief Delete an Individual from a group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndividualNotFoundException if individual_id is not found.
   */
  void deleteIndividualByIdFromGroup(size_t group_position, const std::string& individual_id);

  /**
   * @brief Set the sex of an Individual in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   */
  void setIndividualSexInGroup(size_t group_position, size_t individual_position, const unsigned short sex);

  /**
   * @brief Get the sex of an Individual in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   */
  unsigned short getIndividualSexInGroup(size_t group_position, size_t individual_position) const;

  /**
   * @brief Set the Date of an Individual in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   */
  void setIndividualDateInGroup(size_t group_position, size_t individual_position, const Date& date);

  /**
   * @brief Get the Date of an Individual in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no date.
   */
  const Date* getIndividualDateInGroup(size_t group_position, size_t individual_position) const;

  /**
   * @brief Set the coordinates of an Individual in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   */
  void setIndividualCoordInGroup(size_t group_position, size_t individual_position, const Point2D<double>& coord);

  /**
   * @brief Get the coordinate of an Individual in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no coordinate.
   */
  const Point2D<double>* getIndividualCoordInGroup(size_t group_position, size_t individual_position) const;

  /**
   * @brief Set the Locality of an Individual in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw LocalityNotFoundException if locality_name is not found.
   */
  void setIndividualLocalityInGroupByName(size_t group_position, size_t individual_position, const std::string& locality_name);

  /**
   * @brief Get the Locality of an Individual in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no locality.
   */
  const Locality<double>* getIndividualLocalityInGroup(size_t group_position, size_t individual_position) const;

  /**
   * @brief Add a Sequence to an Individual in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw AlphabetMismatchException if the sequence's alphabet doesn't match the container's alphabet.
   * @throw BadIdentifierException if the sequence's name is already in use.
   */
  void addIndividualSequenceInGroup(size_t group_position, size_t individual_position,
                                    size_t sequence_position, const Sequence& sequence);

  /**
   * @brief Get a Sequence from an Individual of a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no sequences.
   * @throw SequenceNotFoundException if sequence_name is not found.
   * @throw BadIntegerException if sequence_position is already in use.
   */
  const Sequence& getIndividualSequenceByNameInGroup(size_t group_position, size_t individual_position, const std::string& sequence_name) const;

  /**
   * @brief Get a Sequence from an Individual of a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no sequences.
   * @throw SequenceNotFoundException if sequence_position is not found.
   */
  const Sequence& getIndividualSequenceAtPositionInGroup(size_t group_position, size_t individual_position, size_t sequence_position) const;

  /**
   * @brief Delete a Sequence of an Individual of a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no sequences.
   * @throw SequenceNotFoundException if sequence_name is not found.
   */
  void deleteIndividualSequenceByNameInGroup(size_t group_position, size_t individual_position, const std::string& sequence_name);

  /**
   * @brief Delete a Sequence of an Individual of a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no sequences.
   * @throw SequenceNotFoundException if sequence_position is not found.
   */
  void deleteIndividualSequenceAtPositionInGroup(size_t group_position, size_t individual_position, size_t sequence_position);

  /**
   * @brief Get the Sequences' names from an Individual of a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no sequences.
   */
  std::vector<std::string> getIndividualSequencesNamesInGroup(size_t group_position, size_t individual_position) const;

  /**
   * @brief Get the position of a Sequence in an Individual of a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no sequences.
   * @throw SequenceNotFoundException if sequence_name is not found.
   */
  size_t getIndividualSequencePositionInGroup(size_t group_position, size_t individual_position, const std::string& sequence_name) const;

  /**
   * @brief Get the number of Sequences in an Individual of a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no sequences.
   */
  size_t getIndividualNumberOfSequencesInGroup(size_t group_position, size_t individual_position) const;

  /**
   * @brief Set the MultilocusGenotype of an Individual in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   */
  void setIndividualGenotypeInGroup(size_t group_position, size_t individual_position, const MultilocusGenotype& genotype);

  /**
   * @brief Initialyze the genotype of an Individual in a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw BadIntegerException if the number of loci is < 1;
   * @throw NullPointerException if analyzed_loci is NULL.
   * @throw Exception if the individual already has a genotype.
   */
  void initIndividualGenotypeInGroup(size_t group_position, size_t individual_position);

  /**
   * @brief Delete the MultilocusGenotype of an Individual from a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   */
  void deleteIndividualGenotypeInGroup(size_t group_position, size_t individual_position);

  /**
   * @brief Set a MonolocusGenotype of an Individual from a group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no genotype.
   * @throw IndexOutOfBoundsException if locus_position excedes the number of locus.
   */
  void setIndividualMonolocusGenotypeInGroup(size_t group_position, size_t individual_position, size_t locus_position, const MonolocusGenotype& monogen);

  /**
   * @brief Set a MonolocusGenotype of an Individual from a group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no genotype.
   * @throw IndexOutOfBoundsException if locus_position excedes the number of locus.
   * @throw Exception if the ploidy doesn't match.
   */
  void setIndividualMonolocusGenotypeByAlleleKeyInGroup(size_t group_position, size_t individual_position, size_t locus_position, const std::vector<size_t> allele_keys);

  /**
   * @brief Set a MonolocusGenotype of an Individual from a group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no genotype.
   * @throw IndexOutOfBoundsException if locus_position excedes the number of locus.
   * @throw Exception if there is no key in allele_keys.
   */
  void setIndividualMonolocusGenotypeByAlleleIdInGroup(size_t group_position, size_t individual_position, size_t locus_position, const std::vector<std::string> allele_id);

  /**
   * @brief Get a MonolocusGenotype from an Individual of a Group.
   *
   * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
   * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
   * @throw NullPointerException if the individual has no genotype.
   * @throw IndexOutOfBoundsException if locus_position excedes the number of locus.
   * @throw AlleleNotFoundException if at least one of the id is not found.
   */
  const MonolocusGenotype* getIndividualMonolocusGenotypeInGroup(size_t group_position, size_t individual_position, size_t locus_position) const;

  // ** AnalyzedSequences manipulation ******************************************/
  /**
   * @brief Set the alphabet of the AnalyzedSequences.
   */
  void setAlphabet(const Alphabet* alpha);

  /**
   * @brief Set the alphabet of the AnalyzedSequences by its type..
   */
  void setAlphabet(const std::string& alpha_type);

  /**
   * @brief Get the alphabet if there is sequence data.
   *
   * @throw NullPointerException if there is no sequence data.
   */
  const Alphabet* getAlphabet() const;

  /**
   * @brief Get the alphabet type as a string.
   *
   * @throw NullPointerException if there is no sequence data.
   */
  std::string getAlphabetType() const;

  // ** AnalyzedLoci manipulation ***********************************************/
  /**
   * @brief Set the AnalyzedLoci to the DataSet.
   *
   * @throw Exception if at least one Individual has a genotype refering to the actual AnalyzedLoci.
   */
  void setAnalyzedLoci(const AnalyzedLoci& analyzedLoci);

  /**
   * @brief Initialize the AnalyzedLoci for number of loci.
   *
   * @throw Exception if the AnalyzedLoci has already been initialyzed.
   */
  void initAnalyzedLoci(size_t number_of_loci);

  /**
   * @brief Get the AnalyzedLoci if there is one.
   *
   * @throw NullPointerException if there is no AnalyzedLoci.
   */
  const AnalyzedLoci* getAnalyzedLoci() const;

  /**
   * @brief Delete the AnalyzedLoci.
   */
  void deleteAnalyzedLoci();

  /**
   * @brief Set a LocusInfo.
   *
   * @throw NullPointerException if there is no AnalyzedLoci to setup.
   * @throw IndexOutOfBoundsException if locus_position excedes the total of LocusInfo of the DataSet.
   */
  void setLocusInfo(size_t locus_position, const LocusInfo& locus);

  /**
   * @brief Get a LocusInfo by its name.
   */
  const LocusInfo& getLocusInfoByName(const std::string& locus_name) const;

  /**
   * @brief Get a LocusInfo by its position.
   */
  const LocusInfo& getLocusInfoAtPosition(size_t locus_position) const;

  /**
   * @brief Add an AlleleInfo to a LocusInfo.
   */
  void addAlleleInfoByLocusName(const std::string& locus_name, const AlleleInfo& allele);

  /**
   * @brief Add an AlleleInfo to a LocusInfo.
   */
  void addAlleleInfoByLocusPosition(size_t locus_position, const AlleleInfo& allele);

  /**
   * @brief Get the number of loci.
   */
  size_t getNumberOfLoci() const;

  /**
   * @brief Get the ploidy of a locus.
   */
  size_t getPloidyByLocusName(const std::string& locus_name) const;

  /**
   * @brief Get the ploidy of a locus.
   */
  size_t getPloidyByLocusPosition(size_t locus_position) const;

  // ** Container extraction ***************************************************/
  /**
   * @brief Get a PolymorphismMultiGContainer with all allelic data of the DataSet.
   */
  PolymorphismMultiGContainer* getPolymorphismMultiGContainer() const;

  /**
   * @brief Get a PolymorphismMultiGContainer from a selection of groups and individuals.
   *
   * @param selection A map with groups id as keys and vector of individuals position in each group as values.
   */
  PolymorphismMultiGContainer* getPolymorphismMultiGContainer(const std::map<size_t, std::vector<size_t> >& selection) const;

  /**
   * @brief Get a PolymorphismSequenceContainer from a selection of groups and individuals.
   *
   * All the sequences are ingroup. You may change their state after created the container.
   * @param selection A map with groups id as keys and vector of individuals position in each group as values.
   * @param sequence_position The position of the sequence in the individuals;
   */
  PolymorphismSequenceContainer* getPolymorphismSequenceContainer(
		  const std::map<size_t,
		  std::vector<size_t> >& selection,
		  size_t sequence_position) const;

  // ** General tests **********************************************************/
  /**
   * @brief Tell if at least one individual has at least one sequence.
   */
  bool hasSequenceData() const;

  /**
   * @brief Tell if there is alelelic data.
   */
  bool hasAlleleicData() const;
};
} // end of namespace bpp;

#endif // _DATASET_H_

