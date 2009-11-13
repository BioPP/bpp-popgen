//
// File: PolymorphismSequenceContainer.h
// Created by: Eric Bazin
//             Sylvain Gaillard
// Created on: Wednesday August 04 2004
//

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

#include "PolymorphismSequenceContainer.h"

using namespace bpp;

/******************************************************************************/

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const Alphabet *alpha):
  VectorSiteContainer(alpha), ingroup_(vector<bool>()), count_(vector<unsigned int>()), group_(vector<unsigned int>()) {}

/******************************************************************************/

  PolymorphismSequenceContainer::PolymorphismSequenceContainer(unsigned int size, const Alphabet *alpha):
    VectorSiteContainer(size, alpha), ingroup_(vector<bool>(size)), count_(vector<unsigned int>(size)), group_(vector<unsigned int>(size)) {}

/******************************************************************************/

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const OrderedSequenceContainer & sc):
  VectorSiteContainer(sc), ingroup_(vector<bool>(sc.getNumberOfSequences(), true)), count_(vector<unsigned int>(sc.getNumberOfSequences(), 1)), group_(vector<unsigned int>(sc.getNumberOfSequences(), 1)) {}

/******************************************************************************/

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const SiteContainer & sc):
  VectorSiteContainer(sc), ingroup_(vector<bool>(sc.getNumberOfSequences(), true)), count_(vector<unsigned int>(sc.getNumberOfSequences(), 1)), group_(vector<unsigned int>(sc.getNumberOfSequences(), 1)) {}

/******************************************************************************/

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const PolymorphismSequenceContainer & psc):
  VectorSiteContainer(psc), ingroup_(vector<bool>(psc.getNumberOfSequences())), count_(vector<unsigned int>(psc.getNumberOfSequences())), group_(vector<unsigned int>(psc.getNumberOfSequences()))
{
  for (unsigned int i = 0 ; i < psc.getNumberOfSequences() ; i++)
  {
    count_[i] = psc.getSequenceCount(i);
    ingroup_[i] = psc.isIngroupMember(i);
    group_[i] = psc.getGroupId(i);
  }
}

/******************************************************************************/

PolymorphismSequenceContainer& PolymorphismSequenceContainer::operator= (const PolymorphismSequenceContainer & psc)
{
  VectorSiteContainer::operator=(psc);
  // Setting up the sequences comments, numbers and ingroup state
  unsigned int nbSeq = psc.getNumberOfSequences();
  count_.resize(nbSeq);
  ingroup_.resize(nbSeq);
  group_.resize(nbSeq);
  for(unsigned int i = 0; i < nbSeq; i++)
  {
    count_[i] = psc.getSequenceCount(i);
    ingroup_[i] = psc.isIngroupMember(i);
    group_[i] = psc.getGroupId(i);
  }
  return *this;
}

/******************************************************************************/

//** Class destructor: *******************************************************/

PolymorphismSequenceContainer::~PolymorphismSequenceContainer()
{
  clear();
}

/******************************************************************************/

//** Other methodes: *********************************************************/

Sequence* PolymorphismSequenceContainer::removeSequence(unsigned int index) throw (IndexOutOfBoundsException)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::removeSequence: index out of bounds.", index, 0, getNumberOfSequences());
  count_.erase(count_.begin() + index);
  ingroup_.erase(ingroup_.begin() + index);
  return VectorSiteContainer::removeSequence(index);
}

/******************************************************************************/

Sequence* PolymorphismSequenceContainer::removeSequence(const std::string& name) throw (SequenceNotFoundException)
{
  try {
    return removeSequence(getSequencePosition(name));
  }
  catch (SequenceNotFoundException & snfe) {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::removeSequence.", name);
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::deleteSequence(unsigned int index) throw (IndexOutOfBoundsException)
{
  try {
    delete removeSequence(index);
  }
  catch (IndexOutOfBoundsException & ioobe) {
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::deleteSequence.", index, 0, getNumberOfSequences());
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::deleteSequence(const std::string& name) throw (SequenceNotFoundException)
{
  try {
    delete removeSequence(name);
  }
  catch (SequenceNotFoundException & snfe) {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::deleteSequence.", name);
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::addSequence(const Sequence& sequence, unsigned int effectif, bool checkNames) throw (Exception)
{
  try {
    VectorSiteContainer::addSequence(sequence, checkNames);
  }
  catch (Exception & e) {
    throw e;
  }
  count_.push_back(effectif);
  ingroup_.push_back(true);
  group_.push_back(0);
}

/******************************************************************************/

void PolymorphismSequenceContainer::clear()
{
  VectorSiteContainer::clear();
  count_.clear();
  ingroup_.clear();
  group_.clear();
}

/******************************************************************************/

unsigned int PolymorphismSequenceContainer::getGroupId(unsigned int index) const throw (IndexOutOfBoundsException)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::getGroupId: index out of bounds.", index, 0, getNumberOfSequences());
  return group_[index];
}

/******************************************************************************/

unsigned int PolymorphismSequenceContainer::getGroupId(const std::string& name) const throw (SequenceNotFoundException)
{
  try {
    return group_[getSequencePosition(name)];
  }
  catch (SequenceNotFoundException & snfe) {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::getGroupId.", name);
  }
}

/******************************************************************************/

std::set<unsigned int> PolymorphismSequenceContainer::getAllGroupsIds() const
{
  set<unsigned int> grp_ids;
  for (unsigned int i = 0 ; i < group_.size() ; i++)
    grp_ids.insert(group_[i]);
  return grp_ids;
}

/******************************************************************************/

void PolymorphismSequenceContainer::setGroupId(unsigned int index, unsigned int group_id) throw (IndexOutOfBoundsException)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setGroupId: index out of bounds.", index, 0, getNumberOfSequences());
  group_[index] = group_id;
}

/******************************************************************************/

void PolymorphismSequenceContainer::setGroupId(const std::string& name, unsigned int group_id) throw (SequenceNotFoundException)
{
  try {
    group_[getSequencePosition(name)] = group_id;
  }
  catch (SequenceNotFoundException & snfe) {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::setGroupId.", name);
  }
}

/******************************************************************************/

unsigned int PolymorphismSequenceContainer::getNumberOfGroups() const
{
  return getAllGroupsIds().size();
}

/******************************************************************************/

bool PolymorphismSequenceContainer::isIngroupMember(unsigned int index) const throw (IndexOutOfBoundsException)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::isIngroupMember: index out of bounds.", index, 0, getNumberOfSequences());
  return ingroup_[index];
}

/******************************************************************************/

bool PolymorphismSequenceContainer::isIngroupMember(const std::string& name) const throw (SequenceNotFoundException)
{
  try {
    return ingroup_[getSequencePosition(name)];
  }
  catch (SequenceNotFoundException & snfe) {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::isIngroupMember.", name);
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::setAsIngroupMember(unsigned int index) throw (IndexOutOfBoundsException)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setAsIngroupMember.", index, 0, getNumberOfSequences());
  ingroup_[index] = true;
}

/******************************************************************************/

void PolymorphismSequenceContainer::setAsIngroupMember(const std::string& name) throw (SequenceNotFoundException)
{
  try {
    unsigned int seqPos = getSequencePosition(name);
    ingroup_[seqPos] = true;
  }
  catch (SequenceNotFoundException & snfe) {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::setAsIngroupMember.", name);
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::setAsOutgroupMember(unsigned int index) throw (IndexOutOfBoundsException)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setAsOutgroupMember.", index, 0, getNumberOfSequences());
  ingroup_[index] = false;
}

/******************************************************************************/

void PolymorphismSequenceContainer::setAsOutgroupMember(const std::string& name) throw (SequenceNotFoundException)
{
  try {
    unsigned int seqPos = getSequencePosition(name);
    ingroup_[seqPos] = false;
  }
  catch (SequenceNotFoundException & snfe) {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::setAsOutgroupMember.", name);
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::setSequenceCount(unsigned int index, unsigned int count) throw (Exception)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setSequenceCount.", index, 0, getNumberOfSequences());
  if (count < 1)
    throw BadIntegerException("PolymorphismSequenceContainer::setSequenceCount: count can't be < 1.", count);
  count_[index] = count;
}

/******************************************************************************/

void PolymorphismSequenceContainer::setSequenceCount(const std::string& name, unsigned int count) throw (Exception)
{
  try {
    setSequenceCount(getSequencePosition(name), count);
  }
  catch (BadIntegerException & bie) {
    throw bie;
  }
  catch (SequenceNotFoundException & snfe) {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::setSequenceCount.", name);
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::incrementSequenceCount(unsigned int index) throw (IndexOutOfBoundsException)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::incrementSequenceCount.", index, 0, getNumberOfSequences());
  count_[index]++;
}

/******************************************************************************/

void PolymorphismSequenceContainer::incrementSequenceCount(const std::string& name) throw (SequenceNotFoundException)
{
  try {
    incrementSequenceCount(getSequencePosition(name));
  }
  catch (SequenceNotFoundException & snfe) {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::incrementSequenceCount.", name);
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::decrementSequenceCount(unsigned int index) throw (Exception)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::decrementSequenceCount.", index, 0, getNumberOfSequences());
  if (count_[index]-1 < 1)
    throw BadIntegerException("PolymorphismSequenceContainer::decrementSequenceCount: count can't be < 1.", count_[index]-1);
  count_[index]--;
}

/******************************************************************************/

void PolymorphismSequenceContainer::decrementSequenceCount(const std::string& name) throw (Exception)
{
  try {
    decrementSequenceCount(getSequencePosition(name));
  }
  catch (BadIntegerException & bie) {
    throw bie;
  }
  catch (SequenceNotFoundException & snfe) {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::decrementSequenceCount.", name);
  }
}

/******************************************************************************/

unsigned int PolymorphismSequenceContainer::getSequenceCount(unsigned int index) const throw (IndexOutOfBoundsException)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::getSequenceCount.", index, 0, getNumberOfSequences());
  return count_[index];
}

/******************************************************************************/

unsigned int PolymorphismSequenceContainer::getSequenceCount(const std::string& name) const throw (SequenceNotFoundException)
{
  try {
    return getSequenceCount(getSequencePosition(name));
  }
  catch (SequenceNotFoundException & snfe) {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::getSequenceCount.", name);
  }
}

/******************************************************************************/

