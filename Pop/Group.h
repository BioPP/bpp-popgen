/*
 * File Group.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 21 2004
 */

// Secured inclusion of header's file
#ifndef _GROUP_H_
#define _GROUP_H_

// From STL
#include <vector>
using namespace std;

// From Utils
#include <Utils/Exceptions.h>

// From SeqLib
#include <Seq/VectorSequenceContainer.h>
#include <Seq/VectorSiteContainer.h>
#include <Seq/SequenceContainerTools.h>

// From local
#include "Individual.h"
#include "GeneralExceptions.h"

/**
 * @brief The Group class.
 *
 * A Group is an ensemble of Individuals with some statistics like the average
 * allele number.
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
		 * @throw BadIdentifierException if individual's identifier is already in use.
		 */
		void addIndividual(const Individual & ind) throw (BadIdentifierException);

		/**
		 * @brief Add an empty Individual to the Group.
		 *
		 * @throw BadIdentifierException if individual_id is already in use.
		 */
		void addEmptyIndividual(const string & individual_id) throw (BadIdentifierException);
		
		/**
		 * @brief Get the position (index) of an Individual.
		 *
		 * @throw IndividualNotFoundException if individual_id is not found.
		 */
		unsigned int getIndividualPosition(const string individual_id) const
			throw (IndividualNotFoundException);

		/**
		 * @brief Get a pointer to an Individual.
		 *
		 * @param individual_id The id of the Individual to find.
		 *
		 * @return A pointer to the Individual.
		 */
		const Individual * getIndividualById(const string individual_id) const
			throw (IndividualNotFoundException);

		/**
		 * @brief Get a pointer to an Individual by index.
		 *
		 * @param index The index of the Individual.
		 *
		 * @return A pointer to the Individual.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 */
		const Individual * getIndividualByIndex(unsigned int index) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Remove an Individual from the Group.
		 *
		 * @param individual_id The id of the Individual to remove.
		 *
		 * @return A pointer to the removed Individual.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 *
		 * Search an Individual in the Group by cheking the id and remove it
		 * if it is found then return a pointer to this Individual.
		 */
		Individual * removeIndividualById(const string individual_id) throw (IndividualNotFoundException);

		/**
		 * @brief Remove an Individual from the Group.
		 *
		 * @param index The index of the Individual to remove.
		 *
		 * @return A pointer to the removed Individual.
		 *
		 * Remove the individual at the index position and return a pointer
		 * to this Individual.
		 */
		Individual * removeIndividualByIndex(unsigned int index) throw (IndexOutOfBoundsException);

		/**
		 * @brief Delete an Individual from the Group.
		 *
		 * @param individual_id The id of the Individual to delete.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 *
		 * Search an Individual in the Group by cheking the id and delete it
		 * if it is foundi and free the memory by calling the destructor of the
		 * Individual.
		 */
		void deleteIndividualById(const string individual_id) throw (IndividualNotFoundException);

		/**
		 * @brief Delete an Individual from the Group.
		 *
		 * @param index The index of the Individual to delete.
		 *
		 * Free the memory by calling the destructor of the Individual.
		 */
		void deleteIndividualByIndex(unsigned int index) throw (IndexOutOfBoundsException);

		/**
		 * @brief Clear the Group.
		 *
		 * Delete all the Individuals of the group.
		 */
		void clear();

		/**
		 * @brief Append a group at the end of another.
		 */
		void append(const Group & group);

//-- Dealing with Individuals -----------------------------
		/**
		 * @brief Set the sex of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individuals.
		 */
		void setIndividualSexByIndex(unsigned int individual_index, const unsigned short sex)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the sex of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individuals.
		 */
		unsigned short getIndividualSexByIndex(unsigned int individual_index) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Set the date of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individuals.
		 */
		void setIndividualDateByIndex(unsigned int individual_index, const Date & date)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the date of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individuals.
		 * @throw NullPointerException if the Individual has no date.
		 */
		const Date * getIndividualDateByIndex(unsigned int individual_index) const
			throw (Exception);

		/**
		 * @brief Set the coordinates of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individuals.
		 */
		void setIndividualCoordByIndex(unsigned int individual_index, const Coord<double> & coord)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the coordinates of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individuals.
		 * @throw NullPointerException if the individual has no coordinate.
		 */
		const Coord<double> * getIndividualCoordByIndex(unsigned int individual_index) const
			throw (Exception);
		
		/**
		 * @brief Set the locality of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individuals.
		 */
		void setIndividualLocalityByIndex(unsigned int individual_index, const Locality<double> * locality)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the locality of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individuals.
		 * @throw NullPointerException if the individual has no locality.
		 */
		const Locality<double> * getIndividualLocalityByIndex(unsigned int individual_index) const
			throw (Exception);

		/**
		 * @brief Add a sequence to an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individuals.
		 * @throw AlphabetMismatchException if the sequence's alphabet doesn't match the container's alphabet.
		 * @throw BadIdentifierException if the sequence's name is already in use.
		 */
		void addSequenceToIndividualByIndex(unsigned int individual_index,
				const string & seq_set,
				const Sequence & sequence)
			throw (Exception);
		/**
		 * @brief get the keys of the sequence set.
		 *
		 * @throw IndexOutOfBoundsException if individual_index excedes the number of individuals.
		 */
		vector<string> getIndividualSequencesKeys(unsigned int individual_index) const
			throw (IndexOutOfBoundsException);
	
//-- Dealing with sequence containers ---------------------
		/**
		 * @brief Get a VectorSequenceContainer from a named sequence set.
		 *
		 * @return A pointer to a VectorSequenceContainer.
		 */
		VectorSequenceContainer * getVectorSequenceContainer(const string & seqset_id) const throw (Exception);

		/**
		 * @brief Get a VectorSiteContainer from a named sequence set.
		 *
		 * @return A pointer to a VectorSiteContainer.
		 */
		VectorSiteContainer * getVectorSiteContainer(const string & seqset_id) const throw (Exception);

		/**
		 * @brief Get the number of Individual in the Group.
		 *
		 * @return An integer as the number of Individual.
		 */
		unsigned int getNumberOfIndividuals() const;

	protected:
		string _id;
		vector <Individual *> _individuals;
};
#endif // _GROUP_H_
