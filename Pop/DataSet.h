/*
 * File DataSet.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 22 2004
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

//** File manipulation *******************************************************/
		/**
		 * @brief Read the DataSet from a file.
		 *
		 * @param path The path to the file which contains the data
		 * to store in the DataSet.
		 */
		void readFile(const string & path);

		/**
		 * @brief Write the DataSet in a flat file.
		 *
		 * @param path The path to the file.
		 */
		void writeFile(const string & path);
		
//** Locality manipulation ***************************************************/
		/**
		 * @brief Add a locality to the DataSet.
		 *
		 * @param locality A Locality object.
		 * @throw BadIdentifierException if the locality's name already exists.
		 */
		void addLocality(Locality<double> & locality) throw (BadIdentifierException);

		/**
		 * @brief Get the position (index) of a locality in the container.
		 *
		 * @return The index (position) of the Locality.
		 * @param name The locality's name to find.
		 * @throw LocalityNotFoundException if the locality's name doesn't match any name in the DataSet.
		 */
		unsigned int getLocalityPosition(const string & name) const throw (LocalityNotFoundException);
		
		/**
		 * @brief Get a Locality by index.
		 *
		 * @return A const pointer to the locality matching the index.
		 * @param index The position (index) of the Locality in the DataSet.
		 * @throw IndexOutOfBoundsException if index excedes the number of locality of the DataSet.
		 */
		const Locality<double> * getLocalityByIndex(unsigned int index) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Get a Locality by name.
		 *
		 * @throw LocalityNotFoundException if the locality's name is not found.
		 */
		const Locality<double> * getLocalityByName(const string & name) const throw (LocalityNotFoundException);

		/**
		 * @brief Delete a Locality from the DataSet.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of Locality.
		 */
		void deleteLocalityByIndex(unsigned int index) throw (IndexOutOfBoundsException);

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
		
//** Group manipulation ******************************************************/
		/**
		 * @brief Add a Group to the DataSet.
		 *
		 * Add a Group to the DataSet.
		 *
		 * @param group A pointer to the Group to add.
		 */
		void addGroup(const Group & group);

		/**
		 * @brief Get a group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 */
		const Group * getGroup(unsigned int group_index) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Delete a Group from the DataSet.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 */
		void deleteGroup(unsigned int group_index) throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the number of Groups.
		 */
		unsigned int getNumberOfGroups() const;

		/**
		 * @brief Merge some Groups in one.
		 *
		 * Merge all the groups which are specified in the first one (smallest
		 * index). When a group is merged to the first, it is deleted from the
		 * DataSet.
		 *
		 * @param groups A vector unsigned int listing the indices of groups to merge.
		 * @throw IndexOutOfBoundsException if one of the int in groups excedes the number of groups.
		 */
		void mergeGroups(vector<unsigned int> & groups) throw (IndexOutOfBoundsException);

//** Individuals manipulation ************************************************/
		/**
		 * @brief Add an Individual to a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw BadIdentifierException if the individual's id is already in use.
		 */
		void addIndividualToGroup(unsigned int group_index, const Individual & individual) throw (Exception);

		/**
		 * @brief Add an empty Individual to a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw BadIdentifierException if the individual's id is already in use.
		 */
		void addEmptyIndividualToGroup(unsigned int group_index, const string & individual_id) throw (Exception);

		/**
		 * @brief Get the number of Individuals in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 */
		unsigned int getNumberOfIndividualsInGroup(unsigned int group_index) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the position (index) of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 */
		unsigned int getIndividualPositionInGroup(unsigned int group_index, const string & individual_id) const
			throw (Exception);

		/**
		 * @brief Get an Individual from a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 */
		const Individual * getIndividualByIndexFromGroup(unsigned int group_index, unsigned int individual_index) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get an Individual from a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 */
		const Individual * getIndividualByIdFromGroup(unsigned int group_index, const string & individual_id) const
			throw (Exception);

		/**
		 * @brief Delete an Individual from a group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 */
		void deleteIndividualByIndexFromGroup(unsigned int group_index, unsigned int individual_index)
			throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Delete an Individual from a group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 */
		void deleteIndividualByIdFromGroup(unsigned int group_index, const string & individual_id)
			throw (Exception);

		/**
		 * @brief Set the sex of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 */
		void setIndividualSexInGroup(unsigned int group_index, unsigned int individual_index, const unsigned short sex)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the sex of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 */
		unsigned short getIndividualSexInGroup(unsigned int group_index, unsigned int individual_index) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Set the Date of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 */
		void setIndividualDateInGroup(unsigned int group_index, unsigned int individual_index, const Date & date)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the Date of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no date.
		 */
		const Date * getIndividualDateInGroup(unsigned int group_index, unsigned int individual_index) const
			throw (Exception);

		/**
		 * @brief Set the coordinates of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 */
		void setIndividualCoordInGroup(unsigned int group_index, unsigned int individual_index, const Coord<double> & coord)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the coordinate of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no coordinate.
		 */
		const Coord<double> * getIndividualCoordInGroup(unsigned int group_index, unsigned int individual_index) const
			throw (Exception);

		/**
		 * @brief Set the Locality of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 */
		void setIndividualLocalityInGroup(unsigned int group_index, unsigned int individual_index, const Locality<double> * locality)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the Locality of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no locality.
		 */
		const Locality<double> * getIndividualLocalityInGroup(unsigned int group_index, unsigned int individual_index) const
			throw (Exception);

		/**
		 * @brief Add a Sequence to an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw AlphabetMismatchException if the sequence's alphabet doesn't match the container's alphabet.
		 * @throw BadIdentifierException if the sequence's name is already in use.
		 */
		void addIndividualSequenceInGroup(unsigned int group_index, unsigned int individual_index, const Sequence & sequence)
			throw (Exception);

		/**
		 * @brief Get a Sequence from an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 */
		const Sequence * getIndividualSequenceByNameInGroup(unsigned int group_index, unsigned int individual_index, const string & sequence_name) const
			throw (Exception);

		/**
		 * @brief Get a Sequence from an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 * @throw IndexOutOfBoundsException if sequence_index excedes the number of sequences.
		 */
		const Sequence * getIndividualSequenceByIndexInGroup(unsigned int group_index, unsigned int individual_index, unsigned int sequence_index) const
			throw (Exception);

		/**
		 * @brief Delete a Sequence of an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 */
		void deleteIndividualSequenceByNameInGroup(unsigned int group_index, unsigned int individual_index, const string & sequence_name)
			throw (Exception);

		/**
		 * @brief Delete a Sequence of an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 * @throw IndexOutOfBoundsException if sequence_index excedes the number of sequences.
		 */
		void deleteIndividualSequenceByIndexInGroup(unsigned int group_index, unsigned int individual_index, unsigned int sequence_index)
			throw (Exception);

		/**
		 * @brief Get the Sequences' names from an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 */
		vector<string> getIndividualSequencesNamesInGroup(unsigned int group_index, unsigned int individual_index) const
			throw (Exception);

		/**
		 * @brief Get the position (index) of a Sequence in an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 */
		unsigned int getIndividualSequencePositionInGroup(unsigned int group_index, unsigned int individual_index, const string & sequence_name) const
			throw (Exception);

		/**
		 * @brief Get the number of Sequences in an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no sequences.
		 */
		unsigned int getIndividualNumberOfSequencesInGroup(unsigned int group_index, unsigned int individual_index) const
			throw (Exception);

		/**
		 * @brief Add a Genotype to an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw Exception if the individual already has a genotype.
		 */
		void addIndividualGenotypeInGroup(unsigned int group_index, unsigned int individual_index, const Genotype & genotype)
			throw (Exception);

		/**
		 * @brief Initialyze the genotype of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if analyzed_loci is NULL.
		 * @throw Exception if the individual already has a genotype.
		 */
		void initIndividualGenotypeInGroup(unsigned int group_index, unsigned int individual_index, const AnalyzedLoci * analyzed_loci)
			throw (Exception);

		/**
		 * @brief Delete the Genotype of an Individual from a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 */
		void deleteIndividualGenotypeInGroup(unsigned int group_index, unsigned int individual_index)
			throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Set a MonolocusGenotype of an Individual from a group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no genotype.
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of locus.
		 */
		void setIndividualMonolocusGenotypeInGroup(unsigned int group_index, unsigned int individual_index, unsigned int locus_index, const MonolocusGenotype & monogen)
			throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype of an Individual from a group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no genotype.
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of locus.
		 * @throw Exception if the ploidy doesn't match.
		 */
		void setIndividualMonolocusGenotypeByAlleleKeyInGroup(unsigned int group_index, unsigned int individual_index, unsigned int locus_index, const vector<unsigned int> allele_keys)
			throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype of an Individual from a group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no genotype.
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of locus.
		 * @throw Exception if the ploidy doesn't match.
		 */
		void setIndividualMonolocusGenotypeByAlleleIdInGroup(unsigned int group_index, unsigned int individual_index, unsigned int locus_index, const vector<unsigned int> allele_id)
			throw (Exception);

		/**
		 * @brief Get a MonolocusGenotype from an Individual of a Group.
		 *
		 * @throw IndexOutOfBoundsException if group_index excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 * @throw NullPointerException if the individual has no genotype.
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of locus.
		 */
		const MonolocusGenotype * getIndividualMonolocusGenotypeInGroup(unsigned int group_index, unsigned int individual_index, unsigned int locus_index) const
			throw (Exception);

//** AnalyzedLoci manipulation ***********************************************/
		/**
		 * @brief Set the AnalyzedLoci to the DataSet.
		 *
		 * @throw Exception if at least one Individual has a genotype refering to the actual AnalyzedLoci.
		 */
		void setAnalyzedLoci(AnalyzedLoci & analyzeLoci) throw (Exception);

		/**
		 * @brief Initialize the AnalyzedLoci for number of loci.
		 *
		 * @throw Exception if the AnalyzedLoci has already been initialyzed.
		 */
		void initAnalyzedLoci(unsigned int number_of_loci) throw (Exception);

		/**
		 * @brief Delete the AnalyzedLoci.
		 *
		 * @throw Exception if at least one Individual has a genotype refering to the actual AnalyzedLoci.
		 */
		void deleteAnalyzedLoci() throw (Exception);
		
		/**
		 * @brief Set a LocusInfo.
		 *
		 * @throw NullPointerException if there is no AnalyzedLoci to setup.
		 * @throw Exception if at least ont Individual has a genotype refering to the actual AnalyzedLoci.
		 * @throw IndexOutOfBoundsException if locus_index excedes the total of LocusInfo of the DataSet.
		 */
		void setLocusInfo(unsigned int locus_index, const LocusInfo & locus)
			throw (Exception);

		/**
		 * @brief Get a LocusInfo by its name.
		 */
		const LocusInfo * getLocusInfoByName(const string & locus_name) const
			throw (Exception);

		/**
		 * @brief Get a LocusInfo by its index.
		 */
		const LocusInfo * getLocusInfoByIndex(unsigned int locus_index) const
			throw (Exception);

		/**
		 * @brief Add an AlleleInfo to a LocusInfo.
		 */
		void addAlleleInfoByLocusName(const string & locus_name, const AlleleInfo & allele)
			throw (Exception);
		
		/**
		 * @brief Add an AlleleInfo to a LocusInfo.
		 */
		void addAlleleInfoByLocusIndex(unsigned int locus_index, const AlleleInfo & allele)
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
		unsigned int getPloidyByLocusIndex(unsigned int locus_index) const throw (Exception);
		
	protected:
		AnalyzedLoci * _analyzedLoci;
		vector<Locality<double> *> _localities;
		vector<Group *> _groups;
};

#endif // _DATASET_H_
