/*
 * File DataSet.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday June 18 2004
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
		 * @throw IndexOutOfBoundsException if index excedes the number of groups.
		 */
		const Group * getGroup(unsigned int index) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Delete a Group from the DataSet.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of groups.
		 */
		void deleteGroup(unsigned int index) throw (IndexOutOfBoundsException);

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

		/**
		 * @brief Add an Individual to a Group.
		 *
		 * @throw IndexOutOfBoundsException if group excedes the number of groups.
		 * @throw BadIdentifierException if the individual's id is already in use.
		 */
		void addIndividualToGroup(unsigned int group, const Individual & individual) throw (Exception);

		/**
		 * @brief Get the position (index) of an Individual in a Group.
		 *
		 * @throw IndexOutOfBoundsException if group excedes the number of groups.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 */
		unsigned int getIndividualPositionFromGroup(const string & individual_id, unsigned int group) const
			throw (Exception);

		/**
		 * @brief Get an Individual from a Group.
		 *
		 * @throw IndexOutOfBoundsException if group excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 */
		const Individual * getIndividualByIndexFromGroup(unsigned int individual_index, unsigned int group) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get an Individual from a Group.
		 *
		 * @throw IndexOutOfBoundsException if group excedes the number of groups.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 */
		const Individual * getIndividualByIdFromGroup(const string & individual_id, unsigned int group) const
			throw (Exception);

		/**
		 * @brief Delete an Individual from a group.
		 *
		 * @throw IndexOutOfBoundsException if group excedes the number of groups.
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individual in the group.
		 */
		void deleteIndividualByIndexFromGroup(unsigned int individual_index, unsigned int group)
			throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Delete an Individual from a group.
		 *
		 * @throw IndexOutOfBoundsException if group excedes the number of groups.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 */
		void deleteIndividualByIdFromGroup(const string & individual_id, unsigned int group)
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
