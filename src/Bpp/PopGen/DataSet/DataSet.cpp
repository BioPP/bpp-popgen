// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include "DataSet.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

DataSet::DataSet(const DataSet& ds) :
  analyzedLoci_(nullptr),
  sequenceAlphabet_(ds.sequenceAlphabet_),
  localities_(),
  groups_()
{
  if (ds.analyzedLoci_)
    analyzedLoci_.reset(ds.analyzedLoci_->clone());
  for (const auto& locality : ds.localities_)
    localities_.push_back(unique_ptr<Locality<double>>(locality->clone()));
  for (const auto& group : ds.groups_)
    groups_.push_back(unique_ptr<Group>(group->clone()));
}

/******************************************************************************/

DataSet& DataSet::operator=(const DataSet& ds)
{
  if (ds.analyzedLoci_)
    analyzedLoci_.reset(ds.analyzedLoci_->clone());
  else
    analyzedLoci_.reset(nullptr);

  sequenceAlphabet_ = ds.sequenceAlphabet_;
 
  localities_.clear(); 
  for (const auto& locality : ds.localities_)
    localities_.push_back(unique_ptr<Locality<double>>(locality->clone()));

  groups_.clear();
  for (const auto& group : ds.groups_)
    groups_.push_back(unique_ptr<Group>(group->clone()));
  
  return *this;
}

// ** Other methods: *********************************************************/

// Dealing with Localities ---------------------------------
void DataSet::addLocality(const Locality<double>& locality)
{
  for (const auto& existingLocality : localities_)
  {
    if (existingLocality->getName() == locality.getName())
      throw BadIdentifierException("DataSet::addLocality: locality name already in use.", locality.getName());
  }
  localities_.push_back(make_unique<Locality<double>>(locality));
}

/******************************************************************************/

size_t DataSet::getLocalityPosition(const std::string& name) const
{
  for (size_t i = 0; i < localities_.size(); ++i)
  {
    if (localities_[i]->getName() == name)
      return i;
  }
  throw LocalityNotFoundException("DataSet::getLocalityPosition: Locality not found.", name);
}

/******************************************************************************/

shared_ptr<const Locality<double>> DataSet::getLocalityAtPosition(size_t localityPosition) const
{
  if (localityPosition >= localities_.size())
    throw IndexOutOfBoundsException("DataSet::getLocalityAtPosition: localityPosition out of bounds.", localityPosition, 0, localities_.size());
  return localities_[localityPosition];
}

/******************************************************************************/

const Locality<double>& DataSet::localityAtPosition(size_t localityPosition) const
{
  if (localityPosition >= localities_.size())
    throw IndexOutOfBoundsException("DataSet::getLocalityAtPosition: localityPosition out of bounds.", localityPosition, 0, localities_.size());
  return *localities_[localityPosition];
}

/******************************************************************************/

shared_ptr<const Locality<double>> DataSet::getLocalityByName(const std::string& name) const
{
  try
  {
    return getLocalityAtPosition(getLocalityPosition(name));
  }
  catch (LocalityNotFoundException& lnfe)
  {
    throw LocalityNotFoundException("DataSet::getLocalityByName: Locality not found.", name);
  }
}

/******************************************************************************/

const Locality<double>& DataSet::localityByName(const std::string& name) const
{
  try
  {
    return localityAtPosition(getLocalityPosition(name));
  }
  catch (LocalityNotFoundException& lnfe)
  {
    throw LocalityNotFoundException("DataSet::getLocalityByName: Locality not found.", name);
  }
}

/******************************************************************************/

void DataSet::deleteLocalityAtPosition(size_t localityPosition)
{
  if (localityPosition >= localities_.size())
    throw IndexOutOfBoundsException("DataSet::deleteLocalityAtPosition: localityPosition out of bounds.", localityPosition, 0, localities_.size());
  localities_.erase(localities_.begin() + static_cast<ptrdiff_t>(localityPosition));
}

/******************************************************************************/

void DataSet::deleteLocalityByName(const std::string& name)
{
  try
  {
    deleteLocalityAtPosition(getLocalityPosition(name));
  }
  catch (LocalityNotFoundException& lnfe)
  {
    throw LocalityNotFoundException("DataSet::deleteLocalityByName: Locality not found.", name);
  }
}

/******************************************************************************/

// Dealing with groups -------------------------------------
void DataSet::addGroup(const Group& group)
{
  for (const auto& existingGroup : groups_)
  {
    if (group.getGroupId() == existingGroup->getGroupId())
      throw BadIdentifierException("DataSet::addGroup: group id already in use.", group.getGroupId());
  }
  groups_.push_back(make_unique<Group>(group));
}

/******************************************************************************/

void DataSet::addEmptyGroup(size_t groupId)
{
  for (const auto& existingGroup : groups_)
  {
    if (groupId == existingGroup->getGroupId())
      throw BadIdentifierException("DataSet::addEmptyGroup: groupId already in use.", groupId);
  }
  groups_.push_back(make_unique<Group>(groupId));
}

/******************************************************************************/

const Group& DataSet::getGroupById(size_t groupId) const
{
  for (size_t i = 0; i < groups_.size(); i++)
  for (const auto& group : groups_)
  {
    if (groupId == group->getGroupId())
      return *group;
  }
  throw GroupNotFoundException("DataSet::getGroupById: groupId not found.", groupId);
}

/******************************************************************************/

string DataSet::getGroupName(size_t groupId) const
{
  string name;
  name = getGroupById(groupId).getGroupName();
  if (!name.empty() )
    return name;
  else
    return TextTools::toString(groupId);
  throw GroupNotFoundException("DataSet::getGroupName: groupId not found.", groupId);
}

/******************************************************************************/

void DataSet::setGroupName(size_t groupId, const std::string& group_name) const
{
  for (size_t i = 0; i < groups_.size(); i++)
  {
    if (groupId == groups_[i]->getGroupId())
    {
      groups_[i]->setGroupName(group_name);
      return;
    }
  }
  throw GroupNotFoundException("DataSet::setGroupName: groupId not found.", groupId);
}

/******************************************************************************/

size_t DataSet::getGroupPosition(size_t groupId) const
{
  for (size_t i = 0; i < groups_.size(); i++)
  {
    if (groupId == groups_[i]->getGroupId())
      return i;
  }
  throw GroupNotFoundException("DataSet::getGroupPosition: groupId not found.", groupId);
}

/******************************************************************************/

const Group& DataSet::getGroupAtPosition(size_t groupPosition) const
{
  if (groupPosition >= groups_.size())
    throw IndexOutOfBoundsException("DataSet::getGroup.", groupPosition, 0, groups_.size());
  return *(groups_[groupPosition]);
}

/******************************************************************************/

void DataSet::deleteGroupAtPosition(size_t groupPosition)
{
  if (groupPosition >= groups_.size())
    throw IndexOutOfBoundsException("DataSet::deleteGroup.", groupPosition, 0, groups_.size());
  groups_.erase(groups_.begin() + static_cast<ptrdiff_t>(groupPosition));
}

/******************************************************************************/

size_t DataSet::getNumberOfGroups() const
{
  return groups_.size();
}

/******************************************************************************/

void DataSet::mergeTwoGroups(size_t source_id, size_t target_id)
{
  // Test the existance of the two groups.
  try
  {
    getGroupById(source_id);
  }
  catch (GroupNotFoundException& e)
  {
    throw GroupNotFoundException("DataSet::mergeTwoGroups: source_id not found.", source_id);
  }
  try
  {
    getGroupById(target_id);
  }
  catch (GroupNotFoundException& e)
  {
    throw GroupNotFoundException("DataSet::mergeTwoGroups: target_id not found.", target_id);
  }
  // Emptie the source into the target
  size_t source_pos = getGroupPosition(source_id);
  size_t target_pos = getGroupPosition(target_id);
  for (size_t i = 0; i < groups_[source_pos]->getNumberOfIndividuals(); i++)
  {
    groups_[target_pos]->addIndividual(groups_[source_pos]->getIndividualAtPosition(i));
  }
  deleteGroupAtPosition(source_pos);
}

/******************************************************************************/

void DataSet::mergeGroups(std::vector<size_t>& groupIds)
{
  // Test if all group id exists in the DataSet
  for (size_t i = 0; i < groupIds.size(); i++)
  {
    try
    {
      getGroupById(groupIds[i]);
    }
    catch (GroupNotFoundException& e)
    {
      throw GroupNotFoundException("DataSet::mergeGroups: group not found.", groupIds[i]);
    }
  }
  // Sort the group id
  sort(groupIds.begin(), groupIds.end());
  // Merge all the groups in the first
  size_t pos_first = getGroupPosition(groupIds[0]);
  for (size_t i = 1; i < groupIds.size(); i++)
  {
    size_t pos_current = getGroupPosition(groupIds[i]);
    for (size_t j = 0; j < getGroupAtPosition(pos_current).getNumberOfIndividuals(); j++)
    {
      groups_[pos_first]->addIndividual(getGroupAtPosition(pos_current).getIndividualAtPosition(j));
    }
    deleteGroupAtPosition(pos_current);
  }
}

/******************************************************************************/

void DataSet::splitGroup(size_t groupId, vector<size_t> individualSelection)
{
  size_t sourcePos;
  try
  {
    sourcePos = getGroupPosition(groupId);
  }
  catch (GroupNotFoundException& gnfe)
  {
    throw GroupNotFoundException("DataSet::splitGroup: groupId not found.", gnfe.getIdentifier());
  }
  size_t newGroupId = 0;
  for (size_t i = 0; i < groups_.size(); i++)
  {
    if (groups_[i]->getGroupId() > newGroupId)
      newGroupId = groups_[i]->getGroupId();
  }
  newGroupId++;
  Group newGroup(newGroupId);
  for (size_t i = 0; i < individualSelection.size(); i++)
  {
    if (individualSelection[i] >= groups_[sourcePos]->getNumberOfIndividuals())
      throw IndexOutOfBoundsException("DataSet::splitGroup: individuals_selection excedes the number of individual in the group.", individualSelection[i], 0, groups_[sourcePos]->getNumberOfIndividuals());
  }
  for (size_t i = 0; i < individualSelection.size(); i++)
  {
    newGroup.addIndividual(*groups_[sourcePos]->removeIndividualAtPosition(individualSelection[i]));
    groups_[sourcePos]->deleteIndividualAtPosition(individualSelection[i]);
  }
  addGroup(newGroup);
}

/******************************************************************************/

// Dealing with individuals -------------------------------

void DataSet::addIndividualToGroup(size_t group, const Individual& individual)
{
  if (group >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::addIndividualToGroup: group out of bounds.", group, 0, getNumberOfGroups());
  try
  {
    groups_[group]->addIndividual(individual);
    if (individual.hasSequences())
      setAlphabet(individual.getSequenceAlphabet());
  }
  catch (BadIdentifierException& bie)
  {
    throw BadIdentifierException("DataSet::addIndividualToGroup: individual's id already in use in this group.", bie.getIdentifier());
  }
}

/******************************************************************************/

void DataSet::addEmptyIndividualToGroup(size_t group, const std::string& individual_id)
{
  if (group >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::addEmptyIndividual: group out of bounds.", group, 0, getNumberOfGroups());
  try
  {
    groups_[group]->addEmptyIndividual(individual_id);
  }
  catch (BadIdentifierException& bie)
  {
    throw BadIdentifierException("DataSet::addEmptyIndividual: individual_id already in use.", bie.getIdentifier());
  }
}

/******************************************************************************/

size_t DataSet::getNumberOfIndividualsInGroup(size_t groupPosition) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getNumberOfIndividualsInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  return groups_[groupPosition]->getNumberOfIndividuals();
}

/******************************************************************************/

size_t DataSet::getIndividualPositionInGroup(size_t groupPosition, const std::string& individual_id) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualPositionFromGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualPosition(individual_id);
  }
  catch (IndividualNotFoundException& infe)
  {
    throw IndividualNotFoundException("DataSet::getIndividualPositionFromGroup: individual_id not found.", infe.getIdentifier());
  }
}

/******************************************************************************/

const Individual& DataSet::getIndividualAtPositionFromGroup(size_t groupPosition, size_t individualPosition) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualAtPositionFromGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualAtPosition(individualPosition);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualAtPositionFromGroup: individualPosition out of bouds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

const Individual& DataSet::getIndividualByIdFromGroup(size_t groupPosition, const std::string& individualId) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualByIdFromGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualById(individualId);
  }
  catch (IndividualNotFoundException& infe)
  {
    throw IndividualNotFoundException("DataSet::getIndividualByIdFromGroup: individualId not found.", infe.getIdentifier());
  }
}

/******************************************************************************/

void DataSet::deleteIndividualAtPositionFromGroup(size_t groupPosition, size_t individualPosition)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::deleteIndividualAtPositionFromGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->deleteIndividualAtPosition(individualPosition);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::deleteIndividualAtPositionFromGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

void DataSet::deleteIndividualByIdFromGroup(size_t groupPosition, const std::string& individual_id)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::deleteIndividualByIdFromGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->deleteIndividualById(individual_id);
  }
  catch (IndividualNotFoundException& infe)
  {
    throw IndividualNotFoundException("DataSet::deleteIndividualByIdFromGroup: individual_id not found.", infe.getIdentifier());
  }
}

/******************************************************************************/

void DataSet::setIndividualSexInGroup(size_t groupPosition, size_t individualPosition, const unsigned short sex)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualSexInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->setIndividualSexAtPosition(individualPosition, sex);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::setIndividualSexInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

unsigned short DataSet::getIndividualSexInGroup(size_t groupPosition, size_t individualPosition) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualSexInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualSexAtPosition(individualPosition);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualSexInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

void DataSet::setIndividualDateInGroup(size_t groupPosition, size_t individualPosition, const Date& date)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualDateInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->setIndividualDateAtPosition(individualPosition, date);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::setIndividualDateInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

const Date& DataSet::individualDateInGroup(size_t groupPosition, size_t individualPosition) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualDateInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualDateAtPosition(individualPosition);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualDateInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualDateInGroup: individual has no date.");
  }
}

/******************************************************************************/

void DataSet::setIndividualCoordInGroup(size_t groupPosition, size_t individualPosition, const Point2D<double>& coord)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualCoordInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->setIndividualCoordAtPosition(individualPosition, coord);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::setIndividualCoordInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

const Point2D<double>& DataSet::individualCoordInGroup(size_t groupPosition, size_t individualPosition) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualCoordInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualCoordAtPosition(individualPosition);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualCoordAtPosition: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualCoordInGroup: individual has no coordinate.");
  }
}

/******************************************************************************/

void DataSet::setIndividualLocalityInGroupByName(size_t groupPosition, size_t individualPosition, const std::string& locality_name)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualLocalityInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->setIndividualLocalityAtPosition(individualPosition, getLocalityByName(locality_name));
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::setIndividualLocalityInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (LocalityNotFoundException& lnfe)
  {
    throw LocalityNotFoundException("DataSet::setIndividualLocalityInGroup: locality_name not found.", lnfe.getIdentifier());
  }
}

/******************************************************************************/

shared_ptr<const Locality<double>> DataSet::getIndividualLocalityInGroup(size_t groupPosition, size_t individualPosition) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualLocalityInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualLocalityAtPosition(individualPosition);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualLocalityInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualLocalityInGroup: individual has no locality.");
  }
}

/******************************************************************************/

void DataSet::addIndividualSequenceInGroup(size_t groupPosition, size_t individualPosition, size_t sequence_position, unique_ptr<Sequence>& sequence)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::addIndividualSequenceInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->addIndividualSequenceAtPosition(individualPosition, sequence_position, sequence);
    setAlphabet(sequence->getAlphabet());
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::addIndividualSequenceInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (AlphabetMismatchException& ame)
  {
    throw AlphabetMismatchException("DataSet::addIndividualSequenceInGroup: sequence's alphabet doesn't match.", ame.getFirstAlphabet(), ame.getSecondAlphabet());
  }
  catch (BadIdentifierException& bie)
  {
    throw BadIdentifierException("DataSet::addIndividualSequenceInGroup: sequence's name already in use.", bie.getIdentifier());
  }
  catch (BadIntegerException& bie)
  {
    throw BadIntegerException("DataSet::addIndividualSequenceInGroup: sequence_position already in use.", bie.getBadInteger());
  }
}

/******************************************************************************/

const Sequence& DataSet::getIndividualSequenceByNameInGroup(size_t groupPosition, size_t individualPosition, const std::string& sequence_name) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualSequenceByNameInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualSequenceByName(individualPosition, sequence_name);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualSequenceByNameInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualSequenceByNameInGroup: individual has no sequences.");
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("DataSet::getIndividualSequenceByNameInGroup: sequence_name not found.", snfe.getSequenceId());
  }
}

/******************************************************************************/

const Sequence& DataSet::getIndividualSequenceAtPositionInGroup(size_t groupPosition, size_t individualPosition, size_t sequence_position) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualSequenceAtPositionInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualSequenceAtPosition(individualPosition, sequence_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    if (string(ioobe.what()).find("individualPosition") < string(ioobe.what()).size())
      throw IndexOutOfBoundsException("DataSet::getIndividualSequenceAtPositionInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
    // if (string(ioobe.what()).find("sequence_position") < string(ioobe.what()).size())
    else
      throw IndexOutOfBoundsException("DataSet::getIndividualSequenceAtPositionInGroup: sequence_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualSequenceAtPositionInGroup: individual has no sequences.");
  }
}

/******************************************************************************/

void DataSet::deleteIndividualSequenceByNameInGroup(size_t groupPosition, size_t individualPosition, const std::string& sequence_name)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceByNameInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->deleteIndividualSequenceByName(individualPosition, sequence_name);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceByNameInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::deleteIndividualSequenceByNameInGroup: individual has no sequences.");
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("DataSet::deleteIndividualSequenceByNameInGroup: sequence_name not found.", snfe.getSequenceId());
  }
}

/******************************************************************************/

void DataSet::deleteIndividualSequenceAtPositionInGroup(size_t groupPosition, size_t individualPosition, size_t sequence_position)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceAtPositionInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->deleteIndividualSequenceAtPosition(individualPosition, sequence_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    if (string(ioobe.what()).find("individualPosition") < string(ioobe.what()).size())
      throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceAtPositionInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
    // if (string(ioobe.what()).find("sequence_position") < string(ioobe.what()).size())
    else
      throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceAtPositionInGroup: sequence_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::deleteIndividualSequenceAtPositionInGroup: individual has no sequences.");
  }
}

/******************************************************************************/

std::vector<std::string> DataSet::getIndividualSequencesNamesInGroup(size_t groupPosition, size_t individualPosition) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualSequencesNamesInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualSequencesNames(individualPosition);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualSequencesNamesInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualSequencesNamesInGroup: individual has no sequences.");
  }
}

/******************************************************************************/

size_t DataSet::getIndividualSequencePositionInGroup(size_t groupPosition, size_t individualPosition, const std::string& sequence_name) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualSequencePositionInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualSequencePosition(individualPosition, sequence_name);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualSequencePositionInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualSequencePositionInGroup: individual has no sequences.");
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("DataSet::getIndividualSequencePositionInGroup: sequence_name not found.", snfe.getSequenceId());
  }
}

/******************************************************************************/

size_t DataSet::getIndividualNumberOfSequencesInGroup(size_t groupPosition, size_t individualPosition) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualNumberOfSequencesInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualNumberOfSequences(individualPosition);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualNumberOfSequencesInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualNumberOfSequencesInGroup: individual has no sequences.");
  }
}

/******************************************************************************/

void DataSet::setIndividualGenotypeInGroup(size_t groupPosition, size_t individualPosition, const MultilocusGenotype& genotype)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualGenotypeInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->setIndividualGenotype(individualPosition, genotype);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::setIndividualGenotypeInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

void DataSet::initIndividualGenotypeInGroup(size_t groupPosition, size_t individualPosition)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::initIndividualGenotypeInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->initIndividualGenotype(individualPosition, getNumberOfLoci());
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::initIndividualGenotypeInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (BadIntegerException& bie)
  {
    throw BadIntegerException("DataSet::initIndividualGenotypeInGroup: number of loci must be > 0.", bie.getBadInteger());
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::initIndividualGenotypeInGroup: analyzed_loci is NULL.");
  }
  catch (Exception&)
  {
    throw Exception("DataSet::initIndividualGenotypeInGroup: individual already has a genotype.");
  }
}

/******************************************************************************/

void DataSet::deleteIndividualGenotypeInGroup(size_t groupPosition, size_t individualPosition)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::deleteIndividualGenotypeInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->deleteIndividualGenotype(individualPosition);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::deleteIndividualGenotypeInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

void DataSet::setIndividualMonolocusGenotypeInGroup(
    size_t groupPosition,
    size_t individualPosition,
    size_t locusPosition,
    const MonolocusGenotypeInterface& monogen)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->setIndividualMonolocusGenotype(individualPosition, locusPosition, monogen);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    if (string(ioobe.what()).find("individualPosition") < string(ioobe.what()).size())
      throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
    // if (string(ioobe.what()).find("locus_position") < string(ioobe.what()).size())
    else
      throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeInGroup: locus_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::setIndividualMonolocusGenotypeInGroup: individual has no genotype.");
  }
}

/******************************************************************************/

void DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup(size_t groupPosition, size_t individualPosition, size_t locus_position, const std::vector<size_t> allele_keys)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    groups_[groupPosition]->setIndividualMonolocusGenotypeByAlleleKey(individualPosition, locus_position, allele_keys);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    if (string(ioobe.what()).find("individualPosition") < string(ioobe.what()).size())
      throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
    // if (string(ioobe.what()).find("locus_position") < string(ioobe.what()).size())
    else
      throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup: locus_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup: individual has no genotype.");
  }
  catch (Exception&)
  {
    throw Exception("DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup: no key in allele_keys.");
  }
}

/******************************************************************************/

void DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup(size_t groupPosition, size_t individualPosition, size_t locus_position, const std::vector<std::string> allele_id)
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  const LocusInfo& locus_info = getLocusInfoAtPosition(locus_position);
  try
  {
    groups_[groupPosition]->setIndividualMonolocusGenotypeByAlleleId(individualPosition, locus_position, allele_id, locus_info);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    if (string(ioobe.what()).find("individualPosition") < string(ioobe.what()).size())
      throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
    // if (string(ioobe.what()).find("locus_position") < string(ioobe.what()).size())
    else
      throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup: locus_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup: individual has no genotype.");
  }
  catch (AlleleNotFoundException& anfe)
  {
    throw AlleleNotFoundException("DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup: id not found.", anfe.getIdentifier());
  }
}

/******************************************************************************/

const MonolocusGenotypeInterface& DataSet::getIndividualMonolocusGenotypeInGroup(
    size_t groupPosition,
    size_t individualPosition,
    size_t locusPosition) const
{
  if (groupPosition >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualMonolocusGenotypeInGroup: groupPosition out of bounds.", groupPosition, 0, getNumberOfGroups());
  try
  {
    return groups_[groupPosition]->getIndividualMonolocusGenotype(individualPosition, locusPosition);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    if (string(ioobe.what()).find("individualPosition") < string(ioobe.what()).size())
      throw IndexOutOfBoundsException("DataSet::getIndividualMonolocusGenotypeInGroup: individualPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
    // if (string(ioobe.what()).find("locus_position") < string(ioobe.what()).size())
    else
      throw IndexOutOfBoundsException("DataSet::getIndividualMonolocusGenotypeInGroup: locus_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualMonolocusGenotypeInGroup: individual has no genotype.");
  }
}

/******************************************************************************/

void DataSet::setAlphabet(const std::string& alphaType)
{
  if (alphaType != string("DNA") && alphaType != string("RNA") && alphaType != string("PROTEIN"))
    throw Exception(string("DataSet::setAlphabet: bad alphabet type. (") + alphaType + string(")."));
  if (alphaType == string("DNA"))
    sequenceAlphabet_ = AlphabetTools::DNA_ALPHABET;
  if (alphaType == string("RNA"))
    sequenceAlphabet_ = AlphabetTools::RNA_ALPHABET;
  if (alphaType == string("PROTEIN"))
    sequenceAlphabet_ = AlphabetTools::PROTEIN_ALPHABET;
}

/******************************************************************************/

void DataSet::setLocusInfo(size_t locus_position, const LocusInfo& locus)
{
  if (analyzedLoci_ == 0)
    throw NullPointerException("DataSet::setLocusInfo: there's no AnalyzedLoci to setup.");
  try
  {
    analyzedLoci_->setLocusInfo(locus_position, locus);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::setLocusInfo: locus_position out of bounds.", locus_position, 0, analyzedLoci_->getNumberOfLoci());
  }
}

/******************************************************************************/

const LocusInfo& DataSet::getLocusInfoByName(const std::string& locus_name) const
{
  if (analyzedLoci_ == 0)
    throw NullPointerException("DataSet::getLocusInfoByName: there's no AnalyzedLoci.");
  try
  {
    return analyzedLoci_->getLocusInfoByName(locus_name);
  }
  catch (LocusNotFoundException& lnfe)
  {
    throw LocusNotFoundException("DataSet::getLocusInfoByName: locus_name not found", locus_name);
  }
}

/******************************************************************************/

const LocusInfo& DataSet::getLocusInfoAtPosition(size_t locus_position) const
{
  if (analyzedLoci_ == 0)
    throw NullPointerException("DataSet::getLocusInfoAtPosition: there's no AnalyzedLoci.");
  try
  {
    return analyzedLoci_->getLocusInfoAtPosition(locus_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getLocusInfoAtPosition: locus_position out of bounds.", locus_position, ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("DataSet::getLocusInfoAtPosition: no locus defined here");
  }
}

/******************************************************************************/

void DataSet::addAlleleInfoByLocusName(const std::string& locus_name, const AlleleInfo& allele)
{
  if (analyzedLoci_ == 0)
    throw NullPointerException("DataSet::addAlleleInfoByLocusName: there's no AnalyzedLoci.");
  try
  {
    analyzedLoci_->addAlleleInfoByLocusName(locus_name, allele);
  }
  catch (LocusNotFoundException& lnfe)
  {
    throw LocusNotFoundException("DataSet::addAlleleInfoByLocusName: locus_name not found.", lnfe.getIdentifier());
  }
  catch (BadIdentifierException& bie)
  {
    throw BadIdentifierException("DataSet::addAlleleInfoByLocusName: allele's id already in use.", bie.getIdentifier());
  }
}

/******************************************************************************/

void DataSet::addAlleleInfoByLocusPosition(size_t locus_position, const AlleleInfo& allele)
{
  if (analyzedLoci_ == 0)
    throw NullPointerException("DataSet::addAlleleInfoByLocusPosition: there's no AnalyzedLoci.");
  try
  {
    analyzedLoci_->addAlleleInfoByLocusPosition(locus_position, allele);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::addAlleleInfoByLocusPosition: locus_position out of bounds.", locus_position, ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (BadIdentifierException& bie)
  {
    throw BadIdentifierException("DataSet::addAlleleInfoByLocusPosition: allele'e id already in use.", bie.getIdentifier());
  }
}

/******************************************************************************/

size_t DataSet::getNumberOfLoci() const
{
  if (analyzedLoci_ == 0)
    throw NullPointerException("DataSet::getNumberOfLoci: there's no AnalyzedLoci.");
  return analyzedLoci_->getNumberOfLoci();
}

/******************************************************************************/

size_t DataSet::getPloidyByLocusName(const std::string& locus_name) const
{
  if (analyzedLoci_ == 0)
    throw NullPointerException("DataSet::getPloidyByLocusName: there's no AnalyzedLoci.");
  try
  {
    return analyzedLoci_->getPloidyByLocusName(locus_name);
  }
  catch (LocusNotFoundException& lnfe)
  {
    throw LocusNotFoundException("DataSet::getPloidyByLocusName: locus_name not found.", lnfe.getIdentifier());
  }
}

/******************************************************************************/

size_t DataSet::getPloidyByLocusPosition(size_t locus_position) const
{
  if (analyzedLoci_ == 0)
    throw NullPointerException("DataSet::getPloidyByLocusPosition: there's no AnalyzedLoci.");
  try
  {
    return analyzedLoci_->getPloidyByLocusPosition(locus_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getPloidyByLocusPosition: locus_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

unique_ptr<PolymorphismMultiGContainer> DataSet::getPolymorphismMultiGContainer() const
{
  auto pmgc = make_unique<PolymorphismMultiGContainer>();
  for (size_t i = 0; i < getNumberOfGroups(); i++)
  {
    string name = groups_[i]->getGroupName();
    pmgc->addGroupName(i, name);
    for (size_t j = 0; j < getNumberOfIndividualsInGroup(i); ++j)
    {
      const auto& tmpInd = getIndividualAtPositionFromGroup(i, j);
      if (tmpInd.hasGenotype())
      {
        auto tmpMg = unique_ptr<MultilocusGenotype>(tmpInd.getGenotype().clone());
        pmgc->addMultilocusGenotype(tmpMg, i);
      }
    }
  }
  return pmgc;
}

/******************************************************************************/

unique_ptr<PolymorphismMultiGContainer> DataSet::getPolymorphismMultiGContainer(const std::map<size_t, std::vector<size_t> >& selection) const
{
  auto pmgc = make_unique<PolymorphismMultiGContainer>();
  for (const auto& it : selection)
  {
    size_t i;
    try
    {
      i = getGroupPosition(it.first);
    }
    catch (GroupNotFoundException& gnfe)
    {
      throw gnfe;
    }
    string name = groups_[i]->getGroupName();
    pmgc->addGroupName(i, name);
    for (size_t j = 0; j < it.second.size(); j++)
    {
      try
      {
        const auto& tmpInd = getIndividualAtPositionFromGroup(i, j);
        if (tmpInd.hasGenotype())
        {
          auto tmpMg = unique_ptr<MultilocusGenotype>(tmpInd.getGenotype().clone());
          pmgc->addMultilocusGenotype(tmpMg, i);
        }
      }
      catch (IndexOutOfBoundsException& ioobe)
      {
        throw ioobe;
      }
    }
  }
  return pmgc;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> DataSet::getPolymorphismSequenceContainer(
    const std::map<size_t, std::vector<size_t> >& selection,
    size_t sequencePosition) const
{
  auto psc = make_unique<PolymorphismSequenceContainer>(getAlphabet());
  for (auto& it : selection)
  {
    size_t i;
    try
    {
      i = getGroupPosition(it.first);
    }
    catch (GroupNotFoundException& gnfe)
    {
      throw gnfe;
    }
    for (size_t j = 0; j < it.second.size(); ++j)
    {
      try
      {
        const auto& tmpInd = getIndividualAtPositionFromGroup(i, j);
        if (tmpInd.hasSequenceAtPosition(sequencePosition))
        {
	  auto tmpSeq = unique_ptr<Sequence>(tmpInd.sequenceAtPosition(sequencePosition).clone());
          psc->addSequenceWithFrequency(tmpSeq->getName(), tmpSeq, 1);
          psc->setGroupId(tmpSeq->getName(), it.first);
        }
      }
      catch (IndexOutOfBoundsException& ioobe)
      {
        throw ioobe;
      }
    }
  }
  return psc;
}

/******************************************************************************/

