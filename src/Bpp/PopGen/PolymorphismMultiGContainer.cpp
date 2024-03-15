// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "PolymorphismMultiGContainer.h"

using namespace bpp;
using namespace std;

// ** Constructors : **********************************************************/

PolymorphismMultiGContainer::PolymorphismMultiGContainer(const PolymorphismMultiGContainer& pmgc) :
  multilocusGenotypes_(pmgc.size()),
  groups_(pmgc.size()),
  groupsNames_()
{
  for (size_t i = 0; i < pmgc.size(); ++i)
  {
    multilocusGenotypes_[i].reset(pmgc.multilocusGenotype(i).clone());
    groups_[i] = pmgc.getGroupId(i);
  }
  for (auto& id : pmgc.getAllGroupsIds())
  {
    string name = pmgc.getGroupName(id);
    groupsNames_[id] = name;
  }
}

// ** Other methods: ********************************************************/

PolymorphismMultiGContainer& PolymorphismMultiGContainer::operator=(const PolymorphismMultiGContainer& pmgc)
{
  clear();
  for (size_t i = 0; i < pmgc.size(); ++i)
  {
    multilocusGenotypes_.push_back(make_unique<MultilocusGenotype>(pmgc.multilocusGenotype(i)));
    groups_.push_back(pmgc.getGroupId(i));
  }
  for (auto& id : pmgc.getAllGroupsIds())
  {
    string name = pmgc.getGroupName(id);
    groupsNames_[id] = name;
  }

  return *this;
}

/******************************************************************************/

void PolymorphismMultiGContainer::addMultilocusGenotype(unique_ptr<MultilocusGenotype>& mg, size_t group)
{
  multilocusGenotypes_.push_back(std::move(mg));
  groups_.push_back(group);
  auto it = groupsNames_.find(group);
  if (!(it != groupsNames_.end()) )
  {
    //Add group with empty name
    groupsNames_[group] = "";
  }
}

/******************************************************************************/

const MultilocusGenotype& PolymorphismMultiGContainer::multilocusGenotype(size_t position) const
{
  if (position >= size())
    throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getMultilocusGenotype: position out of bounds.", position, 0, size() - 1);
  return *multilocusGenotypes_[position];
}

/******************************************************************************/

unique_ptr<MultilocusGenotype> PolymorphismMultiGContainer::removeMultilocusGenotype(size_t position)
{
  if (position >= size())
    throw IndexOutOfBoundsException("PolymorphismMultiGContainer::removeMultilocusGenotype: position out of bounds.", position, 0, size() - 1);
  unique_ptr<MultilocusGenotype> tmpMg = std::move(multilocusGenotypes_[position]);
  multilocusGenotypes_.erase(multilocusGenotypes_.begin() + static_cast<ptrdiff_t>(position));
  groups_.erase(groups_.begin() + static_cast<ptrdiff_t>(position));
  return tmpMg;
}

/******************************************************************************/

void PolymorphismMultiGContainer::deleteMultilocusGenotype(size_t position)
{
  if (position >= size())
    throw IndexOutOfBoundsException("PolymorphismMultiGContainer::deleteMultilocusGenotype: position out of bounds.", position, 0, size() - 1);
  multilocusGenotypes_.erase(multilocusGenotypes_.begin() + static_cast<ptrdiff_t>(position));
  groups_.erase(groups_.begin() + static_cast<ptrdiff_t>(position));
}

/******************************************************************************/

bool PolymorphismMultiGContainer::isAligned() const
{
  size_t value = 0;
  for (size_t i = 0; i < size(); i++)
  {
    if (i == 0)
      value = multilocusGenotypes_[i]->size();
    else if (multilocusGenotypes_[i]->size() != value)
      return false;
  }
  return true;
}

/******************************************************************************/

size_t PolymorphismMultiGContainer::getNumberOfLoci() const
{
  if (!isAligned())
    throw Exception("MultilocusGenotypes are not aligned.");
  if (size() < 1)
    return 0;
  return multilocusGenotypes_[0]->size();
}

/******************************************************************************/

size_t PolymorphismMultiGContainer::getGroupId(size_t position) const
{
  if (position >= size())
    throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getGroupId: position out of bounds.", position, 0, size() - 1);
  return groups_[position];
}

/******************************************************************************/

void PolymorphismMultiGContainer::setGroupId(size_t position, size_t group_id)
{
  if (position >= size())
    throw IndexOutOfBoundsException("PolymorphismMultiGContainer::setGroupId: position out of bounds.", position, 0, size() - 1);
  groups_[position] = group_id;
}

/******************************************************************************/

std::set<size_t> PolymorphismMultiGContainer::getAllGroupsIds() const
{
  set<size_t> groups_ids;
  for (size_t i = 0; i < size(); i++)
  {
    groups_ids.insert(groups_[i]);
  }
  return groups_ids;
}

/******************************************************************************/

std::vector<std::string> PolymorphismMultiGContainer::getAllGroupsNames() const
{
  vector<string> grpsNames;
  for (auto& it : groupsNames_)
  {
    string name = it.second;
    if (!name.empty())
      grpsNames.push_back(name);
    else
      grpsNames.push_back(TextTools::toString(it.first) );
  }

  return grpsNames;
}

/******************************************************************************/

bool PolymorphismMultiGContainer::groupExists(size_t group) const
{
  for (auto& g : groups_)
  {
    if (g == group)
      return true;
  }
  return false;
}

/******************************************************************************/

size_t PolymorphismMultiGContainer::getGroupSize(size_t group) const
{
  size_t counter = 0;
  for (auto& g : groups_)
  {
    if (g == group)
      counter++;
  }
  return counter;
}

/******************************************************************************/

std::string PolymorphismMultiGContainer::getGroupName(size_t groupId) const
{
  string name = TextTools::toString(groupId); //return group id per default.
  auto it = groupsNames_.find(groupId);
  if (it != groupsNames_.end() )
    name = it->second;
  else
    throw GroupNotFoundException("PolymorphismMultiGContainer::getGroupName: group not found.", groupId);
  return name;
}

/******************************************************************************/

void PolymorphismMultiGContainer::setGroupName(size_t groupId, const string& name)
{
  auto it = groupsNames_.find(groupId);
  if (it != groupsNames_.end() )
    it->second = name;
  else
    throw GroupNotFoundException("PolymorphismMultiGContainer::getGroupName: group not found.", groupId);
}

/******************************************************************************/

void PolymorphismMultiGContainer::addGroupName(size_t groupId, const string& name)
{
  groupsNames_[groupId] = name;
}

/******************************************************************************/

size_t PolymorphismMultiGContainer::getLocusGroupSize(size_t group, size_t locusPosition) const
{
  size_t counter = 0;
  for (size_t i = 0; i < size(); ++i)
  {
    try
    {
      if (groups_[i] == group && !multilocusGenotypes_[i]->isMonolocusGenotypeMissing(locusPosition))
        counter++;
    }
    catch (IndexOutOfBoundsException& ioobe)
    {
      throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getGroupSize: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
    }
  }
  return counter;
}

/******************************************************************************/

size_t PolymorphismMultiGContainer::size() const
{
  return multilocusGenotypes_.size();
}

/******************************************************************************/

void PolymorphismMultiGContainer::clear()
{
  multilocusGenotypes_.clear();
  groups_.clear();
  groupsNames_.clear();
}

/******************************************************************************/
