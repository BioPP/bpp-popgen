/*
 * File DataSet.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 01 2004
 */

// Secured inclusion of header's file
#ifndef _DATASET_H_
#define _DATASET_H_

// From the STL
#include <vector>
#include <map>
using namespace std;

// From PopLib (local)
#include "Group.h"
#include "Individual.h"
#include "Locality.h"

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
		
		/**
		 * @brief Add a locality to the DataSet.
		 *
		 * @param locality A Locality object.
		 */
		void addLocality(Locality<double> & locality);

		/**
		 * @brief Remove a Locality from the DataSet.
		 *
		 * Remove a Locality from the DataSet and return a pointer to this Locality.
		 */
		Locality<double> * removeLocality(const Locality<double> * locality) const;

		/**
		 * @brief Delete a Locality from the DataSet.
		 */
		void deleteLocality(const Locality<double> * locality);
		
		/**
		 * @brief Add a Group to the DataSet.
		 *
		 * Add a Group to the DataSet.
		 *
		 * @param group A pointer to the Group to add.
		 */
		void addGroup(const Group * group);

		/**
		 * @brief Delete a Group from the DataSet.
		 */
		void deleteGroup(const Group * group);

		/**
		 * @brief Merge some Groups in one.
		 *
		 * @param groups A vector of pointers to groups to merge.
		 */
		Group * mergeGroups(vector<Group *> & groups) const;
		
	protected:
		map <string, Locality *> _localities;
		vector <Group *> _groups;
};

#endif // _DATASET_H_
