//
// File: PolymorphismSequenceContainer.h
// Authors: bazin <bazin@univ-montp2.fr>
//          Sylvain Gaillard <yragael2001@yahoo.fr>
// Last modification : Wednesday June 16 2004
//

#ifndef _POLYMORPHISMSEQUENCECONTAINER_H_
#define _POLYMORPHISMSEQUENCECONTAINER_H_

#include <Utils/Clonable.h>
#include <Seq/VectorSiteContainer.h>
#include <Seq/SequenceContainerTools.h>

/**
 * @brief The PolymorphismSequenceContainer class.
 *
 * This is a VectorSiteContainer with effectif for each sequence.
 * It also has flag for ingroup and outgroup.
 */
class PolymorphismSequenceContainer : public VectorSiteContainer
{
	public: // Constructors and destructor
		/**
		 * @brief Build a new empty PolymorphismSequenceContainer.
		 */
		PolymorphismSequenceContainer(const Alphabet *alpha);

		/**
		 * @brief Build a new empty PolymorphismSequenceContainer of given size.
		 */	
		PolymorphismSequenceContainer(unsigned int size, const Alphabet *alpha);

		/**
		 * @brief Build a PolymorphismSequenceContainer by copying data from an OrderedSequenceContainer.
		 */
		PolymorphismSequenceContainer(const OrderedSequenceContainer & sc);
		
		/**
		 * @brief Build a PolymorphismSequenceContainer by copying data from a SiteContainer.
		 */
		PolymorphismSequenceContainer(const SiteContainer & sc);

		/**
		 * @brief Copy constructor.
		 */
		PolymorphismSequenceContainer(const PolymorphismSequenceContainer & psc);
		
		/**
		 * @brief Operator= : copy operator.
		 */
		PolymorphismSequenceContainer & operator= (const PolymorphismSequenceContainer & psc);
		
		/**
		 * @brief Destroy a PolymorphismSequenceContainer.
		 */
		virtual ~PolymorphismSequenceContainer();
		
		/**
		 * @brief Clone a PolymorphismSequenceContainer.
		 */
		Clonable * clone() const;
		
	public: // Other methodes
		/**
		 * @brief Remove a sequence by index and return a pointer to this removed sequence.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences.
		 */
		Sequence * removeSequence(unsigned int index) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Remove a sequence by name and return a pointer to this removed sequence.
		 *
		 * @brief throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		Sequence * removeSequence(const string &name) throw (SequenceNotFoundException);
		
		/**
		 * @brief Delete a sequence by index.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences.
		 */
		void deleteSequence(unsigned int index) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Delete a sequence by name.
		 *
		 * @brief throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		void deleteSequence(const string &name) throw (SequenceNotFoundException);
		
		/**
		 * @brief Add a sequence to the container.
		 *
		 * @throw AlphabetMismatchException if the sequence's alphabet doesn't match the container's alphabet.
		 * @throw SequenceException if the sequence's size doesn't match the sequence's size of the container.
		 * @throw SequenceException if the sequence's name already exists in the container.
		 */
		void addSequence(const Sequence &sequence,unsigned int effectif=1,  bool checkNames=true) throw (Exception);
		
		/**
		 * @brief Clear the container of all its sequences.
		 */
		void clear();
		
		/**
		 * @brief Tell if the sequence is ingroup by index.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 */
		bool isIngroupMember(unsigned int index) const throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Tell if a sequence is ingroup by name.
		 *
		 * @brief throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		bool isIngroupMember(const string &name) const throw (SequenceNotFoundException);
		
		/**
		 * @brief Set a sequence as ingroup member by index.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 */
		bool setAsIngroupMember(unsigned int index) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Set a sequence as ingroup member by name.
		 *
		 * @brief throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		bool setAsIngroupMember(const string &name) throw (SequenceNotFoundException);
		
		/**
		 * @brief Set a sequence as outgroup member by index.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 */
		bool setAsOutgroupMember(unsigned int index) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Set a sequence as outgroup member by name.
		 *
		 * @brief throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		bool setAsOutgroupMember(const string &name) throw (SequenceNotFoundException);
		
		/**
		 * @brief Set the strength of a sequence by index.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 * @throw BadIntegerException if stregth < 1 ... use deleteSequence instead of setting the strength to 0.
		 */
		void setSequenceStrength(unsigned int index, unsigned int strength) throw (Exception);
		
		/**
		 * @brief Set the strength of a sequence by name.
		 *
		 * @throw throw SequenceNotFoundException if name is not found among the sequences' names.
		 * @throw BadIntegerException if stregth < 1 ... use deleteSequence instead of setting the strength to 0.
		 */
		void setSequenceStrength(const string &name, unsigned int strength) throw (Exception);

		/**
		 * @brief Add 1 to the sequence strength.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 */
		void incrementSequenceStrength(unsigned int index) throw (IndexOutOfBoundsException);

		/**
		 * @brief Add 1 to the sequence strength.
		 *
		 * @throw throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		void incrementSequenceStrength(const string &name) throw (SequenceNotFoundException);

		/**
		 * @brief Remove 1 to the sequence strength.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 * @throw BadIntegerException if stregth < 1 ... use deleteSequence instead of setting the strength to 0.
		 */
		void decrementSequenceStrength(unsigned int index) throw (Exception);

		/**
		 * @brief Remove 1 to the sequence strength.
		 *
		 * @throw throw SequenceNotFoundException if name is not found among the sequences' names.
		 * @throw BadIntegerException if stregth < 1 ... use deleteSequence instead of setting the strength to 0.
		 */
		void decrementSequenceStrength(const string & name) throw (Exception);

		/**
		 * @brief Get the strength of a sequence by index.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 */
		unsigned int getSequenceStrength(unsigned int index) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the strength of a sequence by name.
		 *
		 * @brief throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		unsigned int getSequenceStrength(const string &name) const throw (SequenceNotFoundException);
		
	protected:
		vector<bool> _ingroup;
		vector<unsigned int> _strength;
};

#endif	//_POLYMORPHISMSEQUENCECONTAINER_H_
