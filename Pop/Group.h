/*
 * File Group.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday April 05 2004
 */

// Secured inclusion of header's file
#ifndef _GROUP_H_
#define _GROUP_H_

// From STL
#include <vector>
using namespace std;

/**
 * @brief The Group class.
 *
 * A Group is a set of Individuals with some statistics like  the average allele number.
 */
class Group {
	public: // Constructors and destructor :
		
		/**
		 * @brief Build a void new Group.
		 */
		Group();

		/**
		 * @brief Destroy an Group.
		 */
		~Group();
		
	public:
		
		/**
		 * @brief Add an Individual.
		 *
		 * Add an Individual to the group.
		 *
		 * @param ind The Individual to add to the Group.
		 */
		void addIndividual(const Individual & ind);

		/**
		 * @brief Clear the Group.
		 *
		 * Remove all the Individuals of the group.
		 */
		void clear();
	protected:
		vector <Individual *> _individuals;
		unsigned int _averageAlleleNumber;
};
#endif // _GROUP_H_
