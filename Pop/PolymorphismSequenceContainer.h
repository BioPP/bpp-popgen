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
		 * @brief Build a PolymorphismSequenceContainer by copying data from an OrderedSequenceContainer.
		 */
		PolymorphismSequenceContainer(const OrderedSequenceContainer & sc);
		
		/**
		 * @brief Build a PolymorphismSequenceContainer by copying data from a SiteContainer.
		 */
		PolymorphismSequenceContainer(const SiteContainer & sc);
		
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
		Clonable* clone() const;
		
	public: // Other methodes
		/**
		 * @brief Remove a sequence by index and return a pointer to this removed sequence.
		 */
		Sequence * removeSequence(unsigned int index) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Remove a sequence by name and return a pointer to this removed sequence.
		 */
		Sequence * removeSequence(const string &name) throw (SequenceNotFoundException);
		
		/**
		 * @brief Delete a sequence by index.
		 */
		void deleteSequence(unsigned int index) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Delete a sequence by name.
		 */
		void deleteSequence(const string &name) throw (SequenceNotFoundException);
		
		/**
		 * @brief Add a sequence to the container.
		 */
		void addSequence(const Sequence &sequence,unsigned int effectif=1,  bool checkNames=true) throw (Exception);
		
		/**
		 * @brief Clear the container of all its sequences.
		 */
		void clear();
		
		/**
		 * @brief Tell if the sequence is ingroup by index.
		 */
		bool isIngroup(unsigned int index) const throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Tell if a sequence is ingroup by name.
		 */
		bool isIngroup(const string &name) const throw (SequenceNotFoundException);
		
		/**
		 * @brief Toggle the ingroup state of a sequence by index.
		 */
		bool toggleIngroup(unsigned int index) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Toggle the ingroup state of a sequence by name.
		 */
		bool toggleIngroup(const string &name) throw (SequenceNotFoundException);
		
		/**
		 * @brief Set the effectif of a sequence by index.
		 */
		void setEffectif(unsigned int sequence, unsigned int effectif) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Set the effectif of a sequence by name.
		 */
		void setEffectif(const string &name, unsigned int effectif) throw (SequenceNotFoundException);

		/**
		 * @brief Get the effectif of a sequence by index.
		 */
		unsigned int getEffectif(unsigned int index) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Get theeffectif of a sequence by name.
		 */
		unsigned int getEffectif(const string &name) const throw (SequenceNotFoundException);
		
	protected:
		vector<bool> _ingroup;
		vector<unsigned int> _effectif;
};

#endif	//_POLYMORPHISMSEQUENCECONTAINER_H_
