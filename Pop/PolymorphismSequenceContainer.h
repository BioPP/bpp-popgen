/*
 * File: PolymorphismSequenceContainer.h
 * Authors: Eric Bazin <bazin@univ-montp2.fr>
 *          Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday August 04 2004
 *
 * Copyright (C) 2004 Eric Bazin, Sylvain Gaillard and the
 *                    PopLib Development Core Team
 *
 * This file is part of PopLib.
 *
 * PopLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef _POLYMORPHISMSEQUENCECONTAINER_H_
#define _POLYMORPHISMSEQUENCECONTAINER_H_

#include <set>
using namespace std;

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
		 * @throw SequenceNotFoundException if name is not found among the sequences' names.
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
		 * @throw SequenceNotFoundException if name is not found among the sequences' names.
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
		 * @brief Get the group identifier of the sequence.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 */
		unsigned int getGroupId(unsigned int index) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the group identifier of a sequence.
		 *
		 * @throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		unsigned int getGroupId(const string & name) const throw (SequenceNotFoundException);

		/**
		 * @brief Get all the groups identifiers.
		 */
		set<unsigned int> getAllGroupsIds() const;
		
		/**
		 * @brief Set the group identifier of a sequence.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 */
		void setGroupId(unsigned int index, unsigned int group_id) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Set the group identifier of a sequence.
		 *
		 * @throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		void setGroupId(const string & name, unsigned int group_id) throw (SequenceNotFoundException);
		
		/**
		 * @brief Get the number of groups.
		 */
		unsigned int getNumberOfGroups() const;

		/**
		 * @brief Tell if the sequence is ingroup by index.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 */
		bool isIngroupMember(unsigned int index) const throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Tell if a sequence is ingroup by name.
		 *
		 * @throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		bool isIngroupMember(const string &name) const throw (SequenceNotFoundException);
		
		/**
		 * @brief Set a sequence as ingroup member by index.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 */
		void setAsIngroupMember(unsigned int index) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Set a sequence as ingroup member by name.
		 *
		 * @throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		void setAsIngroupMember(const string &name) throw (SequenceNotFoundException);
		
		/**
		 * @brief Set a sequence as outgroup member by index.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 */
		void setAsOutgroupMember(unsigned int index) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Set a sequence as outgroup member by name.
		 *
		 * @throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		void setAsOutgroupMember(const string &name) throw (SequenceNotFoundException);
		
		/**
		 * @brief Set the count of a sequence by index.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 * @throw BadIntegerException if count < 1 ... use deleteSequence instead of setting the count to 0.
		 */
		void setSequenceCount(unsigned int index, unsigned int count) throw (Exception);
		
		/**
		 * @brief Set the count of a sequence by name.
		 *
		 * @throw throw SequenceNotFoundException if name is not found among the sequences' names.
		 * @throw BadIntegerException if count < 1 ... use deleteSequence instead of setting the count to 0.
		 */
		void setSequenceCount(const string &name, unsigned int count) throw (Exception);

		/**
		 * @brief Add 1 to the sequence count.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 */
		void incrementSequenceCount(unsigned int index) throw (IndexOutOfBoundsException);

		/**
		 * @brief Add 1 to the sequence count.
		 *
		 * @throw throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		void incrementSequenceCount(const string &name) throw (SequenceNotFoundException);

		/**
		 * @brief Remove 1 to the sequence count.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 * @throw BadIntegerException if count < 1 ... use deleteSequence instead of setting the count to 0.
		 */
		void decrementSequenceCount(unsigned int index) throw (Exception);

		/**
		 * @brief Remove 1 to the sequence count.
		 *
		 * @throw throw SequenceNotFoundException if name is not found among the sequences' names.
		 * @throw BadIntegerException if count < 1 ... use deleteSequence instead of setting the count to 0.
		 */
		void decrementSequenceCount(const string & name) throw (Exception);

		/**
		 * @brief Get the count of a sequence by index.
		 *
		 * @throw IndexOutOfBoundsException if index excedes the number of sequences in the container.
		 */
		unsigned int getSequenceCount(unsigned int index) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the count of a sequence by name.
		 *
		 * @throw SequenceNotFoundException if name is not found among the sequences' names.
		 */
		unsigned int getSequenceCount(const string &name) const throw (SequenceNotFoundException);
		
	protected:
		vector<bool> _ingroup;
		vector<unsigned int> _count;
		vector<unsigned int> _group;
};

#endif	//_POLYMORPHISMSEQUENCECONTAINER_H_
