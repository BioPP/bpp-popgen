/*
 * File DataSet.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday April 06 2004
 */

// Secured inclusion of header's file
#ifndef _DATASET_H_
#define _DATASET_H_

// From the STL
#include <vector>
using namespace std;

// From PopLib (local)
#include "Group.h"

/**
 * @brief The DataSet class.
 *
 * A DataSet is a group of Group with the possibility to set
 * which one is ingroup and which ones are outgroup.
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
		 * @brief Add a Group to the DataSet.
		 *
		 * Add a Group to the DataSet.
		 *
		 * @param group A pointer to the Group to add.
		 */
		void addGroup(const Group * group);
		
		/**
		 * @brief Remove all Groups from the DataSet.
		 *
		 * Remove all the Groups from the DataSet.
		 */
		void clear();

		/**
		 * @brief Test if the group is set as ingroup.
		 *
		 * @param group A pointer to the Group to test.
		 * @return The test result as a Boolean.
		 */
		bool isInGroup(const Group * group);

	protected:
		vector <Group *> _groups;
		vector <unsigned short> _status;
};

#endif // _DATASET_H_
