/*
 * File Group.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 01 2004
 */

// Secured inclusion of header's file
#ifndef _GROUP_H_
#define _GROUP_H_

// From STL
#include <vector>
#include <map>
using namespace std;

// From SeqLib
#include <VectorSequenceContainer>

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
		 * @brief Set the id of the Group.
		 *
		 * @param group_id The id of the Group as a string.
		 */
		void setGroupId(const string group_id);

		/**
		 * @brief Get the id of the Group.
		 *
		 * @return The id of the Group as a String.
		 */
		string getGroupId() const;
		
		/**
		 * @brief Add an Individual.
		 *
		 * Add an Individual to the group.
		 *
		 * @param ind The Individual to add to the Group.
		 */
		void addIndividual(const Individual & ind);

		/**
		 * @brief Get a pointer to an Individual.
		 *
		 * @param individual_id The id of the Individual to find.
		 *
		 * @return A pointer to the removed Individual.
		 */
		Individual * getIndividual(const string individual_id) const;

		/**
		 * @brief Delete an Individual from the Group.
		 *
		 * @param individual_id The id of the Individual to delete.
		 *
		 * Search an Individual in the Group by cheking the id and delete it
		 * if it is found.
		 */
		void deleteIndividualById(const string individual_id);

		/**
		 * @brief Clear the Group.
		 *
		 * Remove all the Individuals of the group.
		 */
		void clear();

		/**
		 * @brief Get a VectorSequenceContainer from a named sequence set.
		 *
		 * @return A pointer to a VectorSequenceContainer.
		 */
		VectorSequenceContainer * getVectorSequenceContainer(const string & id);

		/**
		 * @brief Get a VectorSiteContainer from a named sequence set.
		 *
		 * @return A pointer to a VectorSiteContainer.
		 */
		VectorSiteContainer * getVectorSiteContainer(const string & id);

		/**
		 * @brief Get the number of Individual in the Group.
		 *
		 * @return An integer as the number of Individual.
		 */
		int getNumberOfIndividuals();

	protected:
		string _id;
		vector <Individual *> _individuals;
};
#endif // _GROUP_H_
