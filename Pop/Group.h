/*
 * File Group.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 *        : Khalid Belkhir
 * Last modification : Thursday July 29 2004
 *
*/
/*
Copyright or © or Copr. CNRS, (November 17, 2004)


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
		Group(unsigned int group_id);

		/**
		 * @brief Copy constructor.
		 *
		 * If you need to use a copy constructor in a DataSet context, use the one
		 * which specify a new Group Id.
		 */
		Group(const Group & group);

		/**
		 * @brief A duplication constructor with new Group Id.
		 */
		Group(const Group & group, unsigned int group_id);

		/**
		 * @brief Destroy an Group.
		 */
		~Group();
		
	public:
		/**
		 * @brief The assignation operator =.
		 */
		Group & operator= (const Group & group);

		/**
		 * @brief Set the id of the Group.
		 *
		 * @param group_id The id of the Group as an unsigned int.
		 */
		void setGroupId(unsigned int group_id);
		
		/**
		 * @brief Get the name of the Group.
		 *
		 * @return The name of the Group as a string.
		 */
		string getGroupName() const;
		
		/**
		 * @brief Set the name of the Group.
		 *
		 * @param group_name Name of the Group as string.
		 */
		void setGroupName(string group_name);
		
		/**
		 * @brief Get the id of the Group.
		 *
		 * @return The id of the Group as an unsigned int.
		 */
		unsigned int getGroupId() const;
		
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
		 * @brief Get the number of Individual in the Group.
		 *
		 * @return An integer as the number of Individual.
		 */
		unsigned int getNumberOfIndividuals() const;

		/**
		 * @brief Get the maximum number of sequence.
		 *
		 * Give the value of the highest sequence key. This value is usefull to
		 * discover the missing sequences data for each individual.
		 */
		unsigned int getMaxNumberOfSequences() const;
		
		/**
		 * @brief Get the position of an Individual.
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
		 * @return A pointer to the Individual or NULL if the Individual is not found.
		 */
		const Individual * getIndividualById(const string individual_id) const;

		/**
		 * @brief Get a pointer to an Individual by its position.
		 *
		 * @param individual_position The position of the Individual in the group.
		 *
		 * @return A pointer to the Individual.
		 * @throw IndividualNotFoundException if individual_id is not found.
		 */
		const Individual * getIndividualAtPosition(unsigned int individual_position) const
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
		 * @param individual_position The position in the Group of the Individual to remove.
		 *
		 * @return A pointer to the removed Individual.
		 *
		 * Remove the individual at the specified position and return a pointer
		 * to this Individual.
		 */
		Individual * removeIndividualAtPosition(unsigned int individual_position) throw (IndexOutOfBoundsException);

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
		 * @param individual_position The position in the Group of the Individual to delete.
		 *
		 * Free the memory by calling the destructor of the Individual.
		 */
		void deleteIndividualAtPosition(unsigned int individual_position) throw (IndexOutOfBoundsException);

		/**
		 * @brief Clear the Group.
		 *
		 * Delete all the Individuals of the group.
		 */
		void clear();

//-- Dealing with Individuals -----------------------------
		/**
		 * @brief Set the sex of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 */
		void setIndividualSexAtPosition(unsigned int individual_position, const unsigned short sex)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the sex of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 */
		unsigned short getIndividualSexAtPosition(unsigned int individual_position) const
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Set the date of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 */
		void setIndividualDateAtPosition(unsigned int individual_position, const Date & date)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the date of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if the Individual has no date.
		 */
		const Date * getIndividualDateAtPosition(unsigned int individual_position) const
			throw (Exception);

		/**
		 * @brief Set the coordinates of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 */
		void setIndividualCoordAtPosition(unsigned int individual_position, const Coord<double> & coord)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the coordinates of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if the individual has no coordinate.
		 */
		const Coord<double> * getIndividualCoordAtPosition(unsigned int individual_position) const
			throw (Exception);
		
		/**
		 * @brief Set the locality of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 */
		void setIndividualLocalityAtPosition(unsigned int individual_position, const Locality<double> * locality)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the locality of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if the individual has no locality.
		 */
		const Locality<double> * getIndividualLocalityAtPosition(unsigned int individual_position) const
			throw (Exception);

		/**
		 * @brief Add a sequence to an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw AlphabetMismatchException if the sequence's alphabet doesn't match the container's alphabet.
		 * @throw BadIdentifierException if the sequence's name is already in use.
		 * @throw BadIntegerException if sequence_position is already in use.
		 */
		void addIndividualSequenceAtPosition(unsigned int individual_position,
				unsigned int sequence_position, const Sequence & sequence)
			throw (Exception);

		/**
		 * @brief Get a sequence of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if there is no sequence container defined in the individual.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 */
		const Sequence * getIndividualSequenceByName(unsigned int individual_position, const string & sequence_name) const
			throw (Exception);

		/**
		 * @brief Get a sequence of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if there is no sequence container defined in the individual.
		 * @throw SequenceNotFoundException if sequence_position is not found.
		 */
		const Sequence * getIndividualSequenceAtPosition(unsigned int individual_position, unsigned int sequence_position) const
			throw (Exception);

		/**
		 * @brief Delete a sequence of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if there is no sequence container defined in the individual.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 */
		void deleteIndividualSequenceByName(unsigned int individual_position, const string & sequence_name)
			throw (Exception);

		/**
		 * @brief Delete a sequence of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if there is no sequence container defined in the individual.
		 * @throw SequenceNotFoundException if sequence_position is not found.
		 */
		void deleteIndividualSequenceAtPosition(unsigned int individual_position, unsigned int sequence_position)
			throw (Exception);

		/**
		 * @brief Tell if the Individual has some sequences.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 */
		bool hasIndividualSequences(unsigned int individual_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the sequences' names from an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if there is no sequence container defined in the individual.
		 */
		vector<string> getIndividualSequencesNames(unsigned int individual_position) const throw (Exception);

		/**
		 * @brief Get the position of a sequence in an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if there is no sequence container defined in the individual.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 */
		unsigned int getIndividualSequencePosition(unsigned int individual_position, const string & sequence_name) const
			throw (Exception);
		
		/**
		 * @brief Get the number of sequences in an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if there is no sequence container defined in the individual.
		 */
		unsigned int getIndividualNumberOfSequences(unsigned int individual_position) const throw (Exception);

		/**
		 * @brief Set all the sequences by copying an OrderedSequenceContainer.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 */
		void setIndividualSequences(unsigned int individual_position, const MapSequenceContainer & msc)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Set the genotype of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 */
		void setIndividualGenotype(unsigned int individual_position, const MultilocusGenotype & genotype) throw (IndexOutOfBoundsException);

		/**
		 * @brief Initialyze the genotype of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw BadIntegerException if loci_number < 1.
		 * @throw Exception if the individual already has a genotype.
		 */
		void initIndividualGenotype(unsigned int individual_position, unsigned int loci_number)
			throw (Exception);

		/**
		 * @brief Delete the genotype of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 */
		void deleteIndividualGenotype(unsigned int individual_position) throw (IndexOutOfBoundsException);

		/**
		 * @brief Tell if an Individual has a genotype.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 */
		bool hasIndividualGenotype(unsigned int individual_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Set a MonolocusGenotype of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if the individual has no genotype.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of locus.
		 */
		void setIndividualMonolocusGenotype(unsigned int individual_position, unsigned int locus_position,
				const MonolocusGenotype & monogen) throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if the individual has no genotype.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of locus.
		 * @throw Exception if there is no key in allele_keys.
		 */
		void setIndividualMonolocusGenotypeByAlleleKey(unsigned int individual_position, unsigned int locus_position,
				const vector<unsigned int> allele_keys) throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if the individual has no genotype.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of locus.
		 * @throw AlleleNotFoundException if at least one id is not found in locus_info.
		 */
		void setIndividualMonolocusGenotypeByAlleleId(unsigned int individual_position, unsigned int locus_position,
				const vector<string> allele_id, const LocusInfo & locus_info) throw (Exception);

		/**
		 * @brief Get a MonolocusGenotype of an Individual.
		 *
		 * @throw IndexOutOfBoundsException if individual_position excedes the number of individuals.
		 * @throw NullPointerException if the individual has no genotype.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of locus.
		 */
		const MonolocusGenotype * getIndividualMonolocusGenotype(unsigned int individual_position,
				unsigned int locus_position) const throw (Exception);

		/**
		 * @brief Tell if at least one individual has at least one sequence.
		 */
		bool hasSequenceData() const;

		/**
		 * @brief Get the alphabet used for the sequences.
		 */
		const Alphabet * getAlphabet() const throw (NullPointerException);

		/**
		 * @brief Get the number of individual that have a data at the specified locus.
		 */
		unsigned int getGroupSizeForLocus(unsigned int locus_position) const;
		
		/**
		 * @brief Get the number of individual that have a sequence at the specified position.
		 */
		unsigned int getGroupSizeForSequence(unsigned int sequence_position) const;

	protected:
		unsigned int _id;
		string _name;
		vector <Individual *> _individuals;
};
#endif // _GROUP_H_
