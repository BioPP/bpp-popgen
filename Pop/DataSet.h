/*
 * File DataSet.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday July 06 2004
 */

// Secured inclusion of header's file
#ifndef _DATASET_H_
#define _DATASET_H_

// From the STL
#include <algorithm>
#include <vector>
#include <map>
using namespace std;

// From Utils
#include <Utils/Exceptions.h>

// From PopLib (local)
#include "Group.h"
#include "Individual.h"
#include "Locality.h"
#include "GeneralExceptions.h"
#include "AnalyzedLoci.h"

/**
 * @brief The DataSet class.
 *
 * A DataSet the object that manage every data on which one can compute
 * some statistics.
 */
class DataSet {
	public: // Constructor and destructor
		/**
		 * @brief Build a new void DataSet.
		 */
		DataSet();

		/**
		 * @brief Destroy a DataSet.
		 */
		~DataSet();
		
	public: // Methodes

//** Locality manipulation ***************************************************/
		/**
		 * @brief Add a locality to the DataSet.
		 *
		 * @param locality A Locality object.
		 * @throw BadIdentifierException if the locality's name already exists.
		 */
		void addLocality(Locality<double> & locality) throw (BadIdentifierException);

		/**
		 * @brief Get the position of a locality in the container.
		 *
		 * @return The locality_position (position) of the Locality.
		 * @param name The locality's name to find.
		 * @throw LocalityNotFoundException if the locality's name doesn't match any name in the DataSet.
		 */
		unsigned int getLocalityPosition(const string & name) const throw (LocalityNotFoundException);
		
		/**
		 * @brief Get a Locality by locality_position.
		 *
		 * @return A const pointer to the locality matching the locality_position.
		 * @param locality_position The position of the Locality in the DataSet.
		 * @throw IndexOutOfBoundsException if locality_position excedes the number of locality of the DataSet.
		 */
		const Locality<double> * getLocalityAtPosition(unsigned int locality_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Get a Locality by name.
		 *
		 * @throw LocalityNotFoundException if the locality's name is not found.
		 */
		const Locality<double> * getLocalityByName(const string & name) const throw (LocalityNotFoundException);

		/**
		 * @brief Delete a Locality from the DataSet.
		 *
		 * @throw IndexOutOfBoundsException if locality_position excedes the number of Locality.
		 */
		void deleteLocalityAtPosition(unsigned int locality_position) throw (IndexOutOfBoundsException);

		/**
		 * @brief Delete a Locality from the DataSet.
		 *
		 * @throw LocalityNotFoundException if the locality's name is not found.
		 */
		void deleteLocalityByName(const string & name) throw (LocalityNotFoundException);

		/**
		 * @brief Get the number of Localities.
		 */
		unsigned int getNumberOfLocalities() const;

		/**
		 * @brief Tell if there is at least one locality.
		 */
		bool hasLocality() const;
		
//** Group manipulation ******************************************************/
		/**
		 * @brief Add a Group to the DataSet.
		 *
		 * Add a Group to the DataSet.
		 *
		 * @param group A pointer to the Group to add.
		 */
		void addGroup(const Group & group) throw (BadIdentifierException);

		/**
		 * @brief Add an empty Group to the DataSet.
		 */
		void addEmptyGroup(unsigned int group_id) throw (BadIdentifierException);
		
		/**
		 * @brief Get a group by identifier.
		 */
		const Group * getGroupById(unsigned int group_id);
		
		/**
		 * @brief Get the position of a Group.
		 *
		 * @throw GroupNotFoundException if the group_id is not found.
		 */
		unsigned int getGroupPosition(unsigned int group_id) throw (GroupNotFoundException);
		
		/**
		 * @brief Get a group by position.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 */
		const Group * getGroupAtPosition(unsigned int group_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Delete a Group from the DataSet.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 */
		void deleteGroupAtPosition(unsigned int group_position) throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the number of Groups.
		 */
		unsigned int getNumberOfGroups() const;

		/**
		 * @brief Merge two groups.
		 *
		 * This methode merge two groups. The source group is emptied into the target
		 * and then is deleted.
		 */
		void mergeTwoGroups(unsigned int source_id, unsigned int target_id) throw (GroupNotFoundException);
		
		/**
		 * @brief Merge some Groups in one.
		 *
		 * Merge all the groups which are specified in the first one (smallest
		 * identifier). When a group is merged to the first, it is deleted from the
		 * DataSet.
		 *
		 * @param group_ids A vector unsigned int listing the id of groups to merge.
		 * @throw IndexOutOfBoundsException if one of the int in groups excedes the number of groups.
		 */
		void mergeGroups(vector<unsigned int> & group_ids) throw (GroupNotFoundException);

		/**
		 * @brief Split a group in two.
		 *
		 * @param group_id The identifier of the source group.
		 * @param individuals_selection The positions of the Individuals to extract from the group to make the new group.
		 * @throw GroupNotFoundException if the group_id is not found.
		 * @throw IndexOutOfBoundsException if one position of the selection excedes the number of individuals of the group.
		 */
		void splitGroup(unsigned int group_id, vector<unsigned int> individuals_selection) throw (Exception);

//** Individuals manipulation ************************************************/
		/**
		 * @brief Add an Individual to a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw BadIdentifierException if the individual's id is already in use.
		 */
		void addIndividualToGroup(unsigned int group_position, const Individual & individual) throw (Exception);

		/**
		 * @brief Add an empty Individual to a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw BadIdentifierException if the individual's id is already in use.
		 */
		void addEmptyIndividualToGroup(unsigned int group_position, const string & individual_id) throw (Exception);

		/**
		 * @brief Get the number of Individuals in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 */
		unsigned int getNumberOfIndividualsInGroup(unsigned int group_position) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the position of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 */
		unsigned int getIndividualPositionInGroup(unsigned int group_position, const string & individual_id) const
			throw (Exception);

		/**
		 * @brief Get an Individual from a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 */
		const Individual * getIndividualAtPositionFromGroup(unsigned int group_position, unsigned int individual_position) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get an Individual from a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 */
		const Individual * getIndividualByIdFromGroup(unsigned int group_position, const string & individual_id) const
			throw (Exception);

		/**
		 * @brief Delete an Individual from a group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 */
		void deleteIndividualAtPositionFromGroup(unsigned int group_position, unsigned int individual_position)
			throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Delete an Individual from a group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 */
		void deleteIndividualByIdFromGroup(unsigned int group_position, const string & individual_id)
			throw (Exception);

		/**
		 * @brief Set the sex of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 */
		void setIndividualSexInGroup(unsigned int group_position, unsigned int individual_position, const unsigned short sex)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the sex of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 */
		unsigned short getIndividualSexInGroup(unsigned int group_position, unsigned int individual_position) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Set the Date of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 */
		void setIndividualDateInGroup(unsigned int group_position, unsigned int individual_position, const Date & date)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the Date of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no date.
		 */
		const Date * getIndividualDateInGroup(unsigned int group_position, unsigned int individual_position) const
			throw (Exception);

		/**
		 * @brief Set the coordinates of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 */
		void setIndividualCoordInGroup(unsigned int group_position, unsigned int individual_position, const Coord<double> & coord)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the coordinate of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no coordinate.
		 */
		const Coord<double> * getIndividualCoordInGroup(unsigned int group_position, unsigned int individual_position) const
			throw (Exception);

		/**
		 * @brief Set the Locality of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw LocalityNotFoundException if locality_name is not found.
		 */
		void setIndividualLocalityInGroupByName(unsigned int group_position, unsigned int individual_position, const string & locality_name)
			throw (Exception);

		/**
		 * @brief Get the Locality of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no locality.
		 */
		const Locality<double> * getIndividualLocalityInGroup(unsigned int group_position, unsigned int individual_position) const
			throw (Exception);

		/**
		 * @brief Add a Sequence to an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw AlphabetMismatchException if the sequence's alphabet doesn't match the container's alphabet.
		 * @throw BadIdentifierException if the sequence's name is already in use.
		 */
		void addIndividualSequenceInGroup(unsigned int group_position, unsigned int individual_position,
				unsigned int sequence_position, const Sequence & sequence)
			throw (Exception);

		/**
		 * @brief Get a Sequence from an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 * @throw BadIntegerException if sequence_position is already in use.
		 */
		const Sequence * getIndividualSequenceByNameInGroup(unsigned int group_position, unsigned int individual_position, const string & sequence_name) const
			throw (Exception);

		/**
		 * @brief Get a Sequence from an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 * @throw SequenceNotFoundException if sequence_position is not found.
		 */
		const Sequence * getIndividualSequenceAtPositionInGroup(unsigned int group_position, unsigned int individual_position, unsigned int sequence_position) const
			throw (Exception);

		/**
		 * @brief Delete a Sequence of an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 */
		void deleteIndividualSequenceByNameInGroup(unsigned int group_position, unsigned int individual_position, const string & sequence_name)
			throw (Exception);

		/**
		 * @brief Delete a Sequence of an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 * @throw SequenceNotFoundException if sequence_position is not found.
		 */
		void deleteIndividualSequenceAtPositionInGroup(unsigned int group_position, unsigned int individual_position, unsigned int sequence_position)
			throw (Exception);

		/**
		 * @brief Get the Sequences' names from an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 */
		vector<string> getIndividualSequencesNamesInGroup(unsigned int group_position, unsigned int individual_position) const
			throw (Exception);

		/**
		 * @brief Get the position of a Sequence in an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 */
		unsigned int getIndividualSequencePositionInGroup(unsigned int group_position, unsigned int individual_position, const string & sequence_name) const
			throw (Exception);

		/**
		 * @brief Get the number of Sequences in an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 */
		unsigned int getIndividualNumberOfSequencesInGroup(unsigned int group_position, unsigned int individual_position) const
			throw (Exception);

		/**
		 * @brief Set the MultilocusGenotype of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 */
		void setIndividualGenotypeInGroup(unsigned int group_position, unsigned int individual_position, const MultilocusGenotype & genotype)
			throw (Exception);

		/**
		 * @brief Initialyze the genotype of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw BadIntegerException if the number of loci is < 1;
		 * @throw NullPointerException if analyzed_loci is NULL.
		 * @throw Exception if the individual already has a genotype.
		 */
		void initIndividualGenotypeInGroup(unsigned int group_position, unsigned int individual_position)
			throw (Exception);

		/**
		 * @brief Delete the MultilocusGenotype of an Individual from a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 */
		void deleteIndividualGenotypeInGroup(unsigned int group_position, unsigned int individual_position)
			throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Set a MonolocusGenotype of an Individual from a group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no genotype.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of locus.
		 */
		void setIndividualMonolocusGenotypeInGroup(unsigned int group_position, unsigned int individual_position, unsigned int locus_position, const MonolocusGenotype & monogen)
			throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype of an Individual from a group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no genotype.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of locus.
		 * @throw Exception if the ploidy doesn't match.
		 */
		void setIndividualMonolocusGenotypeByAlleleKeyInGroup(unsigned int group_position, unsigned int individual_position, unsigned int locus_position, const vector<unsigned int> allele_keys)
			throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype of an Individual from a group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no genotype.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of locus.
		 * @throw Exception if there is no key in allele_keys.
		 */
		void setIndividualMonolocusGenotypeByAlleleIdInGroup(unsigned int group_position, unsigned int individual_position, unsigned int locus_position, const vector<string> allele_id)
			throw (Exception);

		/**
		 * @brief Get a MonolocusGenotype from an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_position excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no genotype.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of locus.
		 * @throw AlleleNotFoundException if at least one of the id is not found.
		 */
		const MonolocusGenotype * getIndividualMonolocusGenotypeInGroup(unsigned int group_position, unsigned int individual_position, unsigned int locus_position) const
			throw (Exception);

//** AnalyzedLoci manipulation ***********************************************/
		/**
		 * @brief Set the AnalyzedLoci to the DataSet.
		 *
		 * @throw Exception if at least one Individual has a genotype refering to the actual AnalyzedLoci.
		 */
		void setAnalyzedLoci(AnalyzedLoci & analyzedLoci) throw (Exception);

		/**
		 * @brief Initialize the AnalyzedLoci for number of loci.
		 *
		 * @throw Exception if the AnalyzedLoci has already been initialyzed.
		 */
		void initAnalyzedLoci(unsigned int number_of_loci) throw (Exception);

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
		void setLocusInfo(unsigned int locus_position, const LocusInfo & locus)
			throw (Exception);

		/**
		 * @brief Get a LocusInfo by its name.
		 */
		const LocusInfo * getLocusInfoByName(const string & locus_name) const
			throw (Exception);

		/**
		 * @brief Get a LocusInfo by its position.
		 */
		const LocusInfo * getLocusInfoAtPosition(unsigned int locus_position) const
			throw (Exception);

		/**
		 * @brief Add an AlleleInfo to a LocusInfo.
		 */
		void addAlleleInfoByLocusName(const string & locus_name, const AlleleInfo & allele)
			throw (Exception);
		
		/**
		 * @brief Add an AlleleInfo to a LocusInfo.
		 */
		void addAlleleInfoByLocusPosition(unsigned int locus_position, const AlleleInfo & allele)
			throw (Exception);

		/**
		 * @brief Get the number of loci.
		 */
		unsigned int getNumberOfLoci() const throw (NullPointerException);

		/**
		 * @brief Get the ploidy of a locus.
		 */
		unsigned int getPloidyByLocusName(const string & locus_name) const throw (Exception);

		/**
		 * @brief Get the ploidy of a locus.
		 */
		unsigned int getPloidyByLocusPosition(unsigned int locus_position) const throw (Exception);

	//** General tests **********************************************************/
		/**
		 * @brief Tell if at least one individual has at least one sequence.
		 */
		bool hasSequenceData() const;

		/**
		 * @brief Get the alphabet if there is sequence data.
		 *
		 * @throw NullPointerException if there is no sequence data.
		 */
		const Alphabet * getAlphabet() const throw (NullPointerException);
		
		/**
		 * @brief Tell if there is alelelic data.
		 */
		bool hasAlleleicData() const;
		
	protected:
		AnalyzedLoci * _analyzedLoci;
		vector<Locality<double> *> _localities;
		vector<Group *> _groups;
};

#endif // _DATASET_H_
