//
// File: PolymorphismSequenceContainer.h
// Created by: Eric Bazin
//             Sylvain Gaillard
// Created on: Wednesday August 04 2004
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)


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

#include <Bpp/Seq/SequenceTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const SequenceContainerInterface& sc, bool count) :
  VectorSiteContainer(sc.getAlphabet()),
  ingroup_(),
  count_(),
  group_()
{
  if (sc.getNumberOfSequences() == 0)
    return; // done.

  // Add first sequence:
  auto seqTmp0 = unique_ptr<Sequence>(sc.sequence(0).clone());
  addSequenceWithFrequency(seqTmp0->getName(), seqTmp0, 1);
  for (size_t i = 1; i < sc.getNumberOfSequences(); ++i)
  {
    const Sequence& seq = sc.sequence(i);
    // Check if this sequence already exists in this container:
    bool exists = false;
    for (size_t j = 0; !exists && j < getNumberOfSequences(); ++j)
    {
      if (SequenceTools::areSequencesIdentical(sequence(j), seq))
      {
        incrementSequenceCount(j); // We increase frequency, meaning that we discard this sequence name.
        exists = true;
      }
    }
    if (!exists)
    {
      auto seqTmp = unique_ptr<Sequence>(seq.clone());
      addSequenceWithFrequency(seqTmp->getName(), seqTmp, 1);
    }
  }
  ingroup_.resize(getNumberOfSequences(), true);
  group_.resize(getNumberOfSequences());
}

/******************************************************************************/

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const PolymorphismSequenceContainer& psc) :
  VectorSiteContainer(psc),
  ingroup_(psc.getNumberOfSequences()),
  count_(psc.getNumberOfSequences()),
  group_(psc.getNumberOfSequences())
{
  for (size_t i = 0; i < psc.getNumberOfSequences(); i++)
  {
    count_[i] = psc.getSequenceCount(i);
    ingroup_[i] = psc.isIngroupMember(i);
    group_[i] = psc.getGroupId(i);
  }
}

/******************************************************************************/

PolymorphismSequenceContainer& PolymorphismSequenceContainer::operator=(const PolymorphismSequenceContainer& psc)
{
  VectorSiteContainer::operator=(psc);
  // Setting up the sequences comments, numbers and ingroup state
  size_t nbSeq = psc.getNumberOfSequences();
  count_.resize(nbSeq);
  ingroup_.resize(nbSeq);
  group_.resize(nbSeq);
  for (size_t i = 0; i < nbSeq; i++)
  {
    count_[i] = psc.getSequenceCount(i);
    ingroup_[i] = psc.isIngroupMember(i);
    group_[i] = psc.getGroupId(i);
  }
  return *this;
}

/******************************************************************************/

/** Class destructor: *********************************************************/

PolymorphismSequenceContainer::~PolymorphismSequenceContainer()
{
  clear();
}

/*****************************************************************************/

/** Other methods: ***********************************************************/

std::unique_ptr<Sequence> PolymorphismSequenceContainer::removeSequence(size_t sequencePosition)
{
  if (sequencePosition >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::removeSequence: index out of bounds.", sequencePosition, 0, getNumberOfSequences());
  count_.erase(count_.begin() + static_cast<ptrdiff_t>(sequencePosition));
  ingroup_.erase(ingroup_.begin() + static_cast<ptrdiff_t>(sequencePosition));
  group_.erase(group_.begin() + static_cast<ptrdiff_t>(sequencePosition));
  return VectorSiteContainer::removeSequence(sequencePosition);
}

/******************************************************************************/

std::unique_ptr<Sequence> PolymorphismSequenceContainer::removeSequence(const std::string& sequenceName)
{
  try
  {
    return removeSequence(getSequencePosition(sequenceName));
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::removeSequence.", sequenceName);
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::deleteSequence(size_t sequencePosition)
{
  try
  {
    removeSequence(sequencePosition); //This returns a smart pointer, which automatically deletes the object.
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::deleteSequence.", sequencePosition, 0, getNumberOfSequences());
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::deleteSequence(const std::string& sequenceName)
{
  try
  {
    removeSequence(sequenceName); //This returns a smart pointer, which automatically deletes the object.
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::deleteSequence.", sequenceName);
  }
}

/******************************************************************************/

std::set<size_t> PolymorphismSequenceContainer::getAllGroupsIds() const
{
  set<size_t> grp_ids;
  for (size_t i = 0; i < group_.size(); i++)
  {
    grp_ids.insert(group_[i]);
  }
  return grp_ids;
}

/******************************************************************************/

bool PolymorphismSequenceContainer::hasOutgroup() const
{
  for (auto i : ingroup_)
  {
    if (!i)
      return true;
  }
  return false;
}

/******************************************************************************/

void PolymorphismSequenceContainer::setAsIngroupMember(size_t index)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setAsIngroupMember.", index, 0, getNumberOfSequences());
  ingroup_[index] = true;
}

/******************************************************************************/

void PolymorphismSequenceContainer::setAsIngroupMember(const std::string& name)
{
  try
  {
    size_t seqPos = getSequencePosition(name);
    ingroup_[seqPos] = true;
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::setAsIngroupMember.", name);
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::setAsOutgroupMember(size_t index)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setAsOutgroupMember.", index, 0, getNumberOfSequences());
  ingroup_[index] = false;
}

/******************************************************************************/

void PolymorphismSequenceContainer::setAsOutgroupMember(const std::string& name)
{
  try
  {
    size_t seqPos = getSequencePosition(name);
    ingroup_[seqPos] = false;
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::setAsOutgroupMember.", name);
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::setSequenceCount(size_t index, unsigned int count)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setSequenceCount.", index, 0, getNumberOfSequences());
  if (count < 1)
    throw BadIntegerException("PolymorphismSequenceContainer::setSequenceCount: count can't be < 1.", static_cast<int>(count));
  count_[index] = count;
}

/******************************************************************************/

void PolymorphismSequenceContainer::setSequenceCount(const std::string& name, unsigned int count)
{
  try
  {
    setSequenceCount(getSequencePosition(name), count);
  }
  catch (BadIntegerException& bie)
  {
    throw bie;
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::setSequenceCount.", name);
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::incrementSequenceCount(size_t index)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::incrementSequenceCount.", index, 0, getNumberOfSequences());
  count_[index]++;
}

/******************************************************************************/

void PolymorphismSequenceContainer::incrementSequenceCount(const std::string& name)
{
  try
  {
    incrementSequenceCount(getSequencePosition(name));
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::incrementSequenceCount.", name);
  }
}

/******************************************************************************/

void PolymorphismSequenceContainer::decrementSequenceCount(size_t index)
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::decrementSequenceCount.", index, 0, getNumberOfSequences());
  if (count_[index] - 1 < 1)
    throw BadIntegerException("PolymorphismSequenceContainer::decrementSequenceCount: count can't be < 1.", static_cast<int>(count_[index] - 1));
  count_[index]--;
}

/******************************************************************************/

void PolymorphismSequenceContainer::decrementSequenceCount(const std::string& name)
{
  try
  {
    decrementSequenceCount(getSequencePosition(name));
  }
  catch (BadIntegerException& bie)
  {
    throw bie;
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::decrementSequenceCount.", name);
  }
}

/******************************************************************************/

unsigned int PolymorphismSequenceContainer::getSequenceCount(size_t index) const
{
  if (index >= getNumberOfSequences())
    throw IndexOutOfBoundsException("PolymorphismSequenceContainer::getSequenceCount.", index, 0, getNumberOfSequences());
  return count_[index];
}

/******************************************************************************/

unsigned int PolymorphismSequenceContainer::getSequenceCount(const std::string& name) const
{
  try
  {
    return getSequenceCount(getSequencePosition(name));
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("PolymorphismSequenceContainer::getSequenceCount.", name);
  }
}

/******************************************************************************/

unique_ptr<SiteContainerInterface> PolymorphismSequenceContainer::toSiteContainer() const
{
  unique_ptr<SiteContainerInterface> sites(new VectorSiteContainer(getAlphabet()));
  for (size_t i = 0; i < getNumberOfSequences(); ++i)
  {
    const auto& seq = sequence(i);
    unsigned int freq = getSequenceCount(i);
    if (freq > 1)
    {
      for (unsigned int j = 0; j < freq; ++j)
      {
        unique_ptr<Sequence> seqdup(seq.clone());
        seqdup->setName(seq.getName() + "_" + TextTools::toString(j + 1));
        sites->addSequence(seqdup->getName(), seqdup);
      }
    }
    else
    {
      unique_ptr<Sequence> seqdup(seq.clone());
      sites->addSequence(seqdup->getName(), seqdup);
    }
  }
  sites->setSiteCoordinates(getSiteCoordinates());
  return sites;
}

/******************************************************************************/
