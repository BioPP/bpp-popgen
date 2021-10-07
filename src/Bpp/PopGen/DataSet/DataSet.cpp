//
// File DataSet.cpp
// Author : Sylvain Gaillard
//          Khalid Belkhir
// Last modification : November 10, 2008
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

#include "DataSet.h"

using namespace bpp;
using namespace std;

// ** Class constructor: *******************************************************/

DataSet::DataSet() : analyzedLoci_(0),
  analyzedSequences_(0),
  localities_(vector<Locality<double>*>()),
  groups_(vector<Group*>()) {}

/******************************************************************************/

DataSet::DataSet(const DataSet& ds) : analyzedLoci_(0),
  analyzedSequences_(0),
  localities_(vector<Locality<double>*>()),
  groups_(vector<Group*>())
{
  if (ds.analyzedLoci_ != 0)
    analyzedLoci_ = new AnalyzedLoci(*(ds.analyzedLoci_));
  if (ds.analyzedSequences_ != 0)
    analyzedSequences_ = new AnalyzedSequences(*(ds.analyzedSequences_));
  if (ds.localities_.size() != 0)
    for (size_t i = 0; i < ds.localities_.size(); i++)
    {
      localities_.push_back(new Locality<double>(*(ds.localities_[i])));
    }
  if (ds.groups_.size() != 0)
    for (size_t i = 0; i < ds.groups_.size(); i++)
    {
      groups_.push_back(new Group(*(ds.groups_[i])));
    }
}

/******************************************************************************/

DataSet& DataSet::operator=(const DataSet& ds)
{
  if (ds.analyzedLoci_ != 0)
    analyzedLoci_ = new AnalyzedLoci(*(ds.analyzedLoci_));
  if (ds.analyzedSequences_ != 0)
    analyzedSequences_ = new AnalyzedSequences(*(ds.analyzedSequences_));
  if (ds.localities_.size() != 0)
    for (size_t i = 0; i < ds.localities_.size(); i++)
    {
      localities_.push_back(new Locality<double>(*(ds.localities_[i])));
    }
  if (ds.groups_.size() != 0)
    for (size_t i = 0; i < ds.groups_.size(); i++)
    {
      groups_.push_back(new Group(*(ds.groups_[i])));
    }
  return *this;
}

// ** Class destructor: *******************************************************/
DataSet::~DataSet()
{
  if (getNumberOfGroups() > 0)
    for (size_t i = 0; i < getNumberOfGroups(); i++)
    {
      delete groups_[i];
    }
  if (analyzedLoci_ != 0)
    delete analyzedLoci_;
  if (getNumberOfLocalities() > 0)
    for (size_t i = 0; i < getNumberOfLocalities(); i++)
    {
      delete localities_[i];
    }
  if (analyzedSequences_ != 0)
    delete analyzedSequences_;
}

// ** Other methodes: *********************************************************/

// Dealing with Localities ---------------------------------
void DataSet::addLocality(Locality<double>& locality)
{
  for (size_t i = 0; i < localities_.size(); i++)
  {
    if (localities_[i]->getName() == locality.getName())
      throw BadIdentifierException("DataSet::addLocality: locality name already in use.", locality.getName());
  }
  localities_.push_back(new Locality<double>(locality));
}

/******************************************************************************/

size_t DataSet::getLocalityPosition(const std::string& name) const
{
  for (size_t i = 0; i < localities_.size(); i++)
  {
    if (localities_[i]->getName() == name)
      return i;
  }
  throw LocalityNotFoundException("DataSet::getLocalityPosition: Locality not found.", name);
}

/******************************************************************************/

const Locality<double>& DataSet::getLocalityAtPosition(size_t locality_position) const
{
  if (locality_position >= localities_.size())
    throw IndexOutOfBoundsException("DataSet::getLocalityAtPosition: locality_position out of bounds.", locality_position, 0, localities_.size());
  return *(localities_[locality_position]);
}

/******************************************************************************/

const Locality<double>& DataSet::getLocalityByName(const std::string& name) const
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

void DataSet::deleteLocalityAtPosition(size_t locality_position)
{
  if (locality_position >= localities_.size())
    throw IndexOutOfBoundsException("DataSet::deleteLocalityAtPosition: locality_position out of bounds.", locality_position, 0, localities_.size());
  delete localities_[locality_position];
  localities_.erase(localities_.begin() + static_cast<ptrdiff_t>(locality_position));
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

size_t DataSet::getNumberOfLocalities() const
{
  return localities_.size();
}

/******************************************************************************/

bool DataSet::hasLocality() const
{
  return getNumberOfLocalities() > 0;
}

/******************************************************************************/

// Dealing with groups -------------------------------------
void DataSet::addGroup(const Group& group)
{
  for (size_t i = 0; i < groups_.size(); i++)
  {
    if (group.getGroupId() == groups_[i]->getGroupId())
      throw BadIdentifierException("DataSet::addGroup: group id already in use.", group.getGroupId());
  }
  groups_.push_back(new Group(group));
}

/******************************************************************************/

void DataSet::addEmptyGroup(size_t group_id)
{
  for (size_t i = 0; i < groups_.size(); i++)
  {
    if (group_id == groups_[i]->getGroupId())
      throw BadIdentifierException("DataSet::addEmptyGroup: group_id already in use.", group_id);
  }
  groups_.push_back(new Group(group_id));
}

/******************************************************************************/

const Group& DataSet::getGroupById(size_t group_id) const
{
  for (size_t i = 0; i < groups_.size(); i++)
  {
    if (group_id == groups_[i]->getGroupId())
      return *(groups_[i]);
  }
  throw GroupNotFoundException("DataSet::getGroupById: group_id not found.", group_id);
}

/******************************************************************************/

string DataSet::getGroupName(size_t group_id) const
{
  string name;
  name = getGroupById(group_id).getGroupName();
  if (!name.empty() )
    return name;
  else
    return TextTools::toString(group_id);
  throw GroupNotFoundException("DataSet::getGroupName: group_id not found.", group_id);
}

/******************************************************************************/

void DataSet::setGroupName(size_t group_id, const std::string& group_name) const
{
  for (size_t i = 0; i < groups_.size(); i++)
  {
    if (group_id == groups_[i]->getGroupId())
    {
      groups_[i]->setGroupName(group_name);
      return;
    }
  }
  throw GroupNotFoundException("DataSet::setGroupName: group_id not found.", group_id);
}

/******************************************************************************/

size_t DataSet::getGroupPosition(size_t group_id) const
{
  for (size_t i = 0; i < groups_.size(); i++)
  {
    if (group_id == groups_[i]->getGroupId())
      return i;
  }
  throw GroupNotFoundException("DataSet::getGroupPosition: group_id not found.", group_id);
}

/******************************************************************************/

const Group& DataSet::getGroupAtPosition(size_t group_position) const
{
  if (group_position >= groups_.size())
    throw IndexOutOfBoundsException("DataSet::getGroup.", group_position, 0, groups_.size());
  return *(groups_[group_position]);
}

/******************************************************************************/

void DataSet::deleteGroupAtPosition(size_t group_position)
{
  if (group_position >= groups_.size())
    throw IndexOutOfBoundsException("DataSet::deleteGroup.", group_position, 0, groups_.size());
  delete groups_[group_position];
  groups_.erase(groups_.begin() + static_cast<ptrdiff_t>(group_position));
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

void DataSet::mergeGroups(std::vector<size_t>& group_ids)
{
  // Test if all group id exists in the DataSet
  for (size_t i = 0; i < group_ids.size(); i++)
  {
    try
    {
      getGroupById(group_ids[i]);
    }
    catch (GroupNotFoundException& e)
    {
      throw GroupNotFoundException("DataSet::mergeGroups: group not found.", group_ids[i]);
    }
  }
  // Sort the group id
  sort(group_ids.begin(), group_ids.end());
  // Merge all the groups in the first
  size_t pos_first = getGroupPosition(group_ids[0]);
  for (size_t i = 1; i < group_ids.size(); i++)
  {
    size_t pos_current = getGroupPosition(group_ids[i]);
    for (size_t j = 0; j < getGroupAtPosition(pos_current).getNumberOfIndividuals(); j++)
    {
      groups_[pos_first]->addIndividual(getGroupAtPosition(pos_current).getIndividualAtPosition(j));
    }
    deleteGroupAtPosition(pos_current);
  }
}

/******************************************************************************/

void DataSet::splitGroup(size_t group_id, std::vector<size_t> individuals_selection)
{
  size_t source_pos;
  try
  {
    source_pos = getGroupPosition(group_id);
  }
  catch (GroupNotFoundException& gnfe)
  {
    throw GroupNotFoundException("DataSet::splitGroup: group_id not found.", gnfe.getIdentifier());
  }
  size_t new_group_id = 0;
  for (size_t i = 0; i < groups_.size(); i++)
  {
    if (groups_[i]->getGroupId() > new_group_id)
      new_group_id = groups_[i]->getGroupId();
  }
  new_group_id++;
  Group new_group(new_group_id);
  for (size_t i = 0; i < individuals_selection.size(); i++)
  {
    if (individuals_selection[i] >= groups_[source_pos]->getNumberOfIndividuals())
      throw IndexOutOfBoundsException("DataSet::splitGroup: individuals_selection excedes the number of individual in the group.", individuals_selection[i], 0, groups_[source_pos]->getNumberOfIndividuals());
  }
  for (size_t i = 0; i < individuals_selection.size(); i++)
  {
    new_group.addIndividual(*groups_[source_pos]->removeIndividualAtPosition(individuals_selection[i]));
    groups_[source_pos]->deleteIndividualAtPosition(individuals_selection[i]);
  }
  addGroup(new_group);
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

size_t DataSet::getNumberOfIndividualsInGroup(size_t group_position) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getNumberOfIndividualsInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  return groups_[group_position]->getNumberOfIndividuals();
}

/******************************************************************************/

size_t DataSet::getIndividualPositionInGroup(size_t group_position, const std::string& individual_id) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualPositionFromGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return groups_[group_position]->getIndividualPosition(individual_id);
  }
  catch (IndividualNotFoundException& infe)
  {
    throw IndividualNotFoundException("DataSet::getIndividualPositionFromGroup: individual_id not found.", infe.getIdentifier());
  }
}

/******************************************************************************/

const Individual* DataSet::getIndividualAtPositionFromGroup(size_t group_position, size_t individual_position) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualAtPositionFromGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return &groups_[group_position]->getIndividualAtPosition(individual_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualAtPositionFromGroup: individual_position out of bouds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

const Individual* DataSet::getIndividualByIdFromGroup(size_t group_position, const std::string& individual_id) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualByIdFromGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return &groups_[group_position]->getIndividualById(individual_id);
  }
  catch (IndividualNotFoundException& infe)
  {
    throw IndividualNotFoundException("DataSet::getIndividualByIdFromGroup: individual_id not found.", infe.getIdentifier());
  }
}

/******************************************************************************/

void DataSet::deleteIndividualAtPositionFromGroup(size_t group_position, size_t individual_position)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::deleteIndividualAtPositionFromGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->deleteIndividualAtPosition(individual_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::deleteIndividualAtPositionFromGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

void DataSet::deleteIndividualByIdFromGroup(size_t group_position, const std::string& individual_id)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::deleteIndividualByIdFromGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->deleteIndividualById(individual_id);
  }
  catch (IndividualNotFoundException& infe)
  {
    throw IndividualNotFoundException("DataSet::deleteIndividualByIdFromGroup: individual_id not found.", infe.getIdentifier());
  }
}

/******************************************************************************/

void DataSet::setIndividualSexInGroup(size_t group_position, size_t individual_position, const unsigned short sex)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualSexInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->setIndividualSexAtPosition(individual_position, sex);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::setIndividualSexInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

unsigned short DataSet::getIndividualSexInGroup(size_t group_position, size_t individual_position) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualSexInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return groups_[group_position]->getIndividualSexAtPosition(individual_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualSexInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

void DataSet::setIndividualDateInGroup(size_t group_position, size_t individual_position, const Date& date)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualDateInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->setIndividualDateAtPosition(individual_position, date);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::setIndividualDateInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

const Date* DataSet::getIndividualDateInGroup(size_t group_position, size_t individual_position) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualDateInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return &groups_[group_position]->getIndividualDateAtPosition(individual_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualDateInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualDateInGroup: individual has no date.");
  }
}

/******************************************************************************/

void DataSet::setIndividualCoordInGroup(size_t group_position, size_t individual_position, const Point2D<double>& coord)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualCoordInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->setIndividualCoordAtPosition(individual_position, coord);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::setIndividualCoordInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

const Point2D<double>* DataSet::getIndividualCoordInGroup(size_t group_position, size_t individual_position) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualCoordInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return &groups_[group_position]->getIndividualCoordAtPosition(individual_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualCoordAtPosition: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualCoordInGroup: individual has no coordinate.");
  }
}

/******************************************************************************/

void DataSet::setIndividualLocalityInGroupByName(size_t group_position, size_t individual_position, const std::string& locality_name)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualLocalityInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->setIndividualLocalityAtPosition(individual_position, &getLocalityByName(locality_name));
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::setIndividualLocalityInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (LocalityNotFoundException& lnfe)
  {
    throw LocalityNotFoundException("DataSet::setIndividualLocalityInGroup: locality_name not found.", lnfe.getIdentifier());
  }
}

/******************************************************************************/

const Locality<double>* DataSet::getIndividualLocalityInGroup(size_t group_position, size_t individual_position) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualLocalityInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return &groups_[group_position]->getIndividualLocalityAtPosition(individual_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualLocalityInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualLocalityInGroup: individual has no locality.");
  }
}

/******************************************************************************/

void DataSet::addIndividualSequenceInGroup(size_t group_position, size_t individual_position, size_t sequence_position, const Sequence& sequence)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::addIndividualSequenceInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->addIndividualSequenceAtPosition(individual_position, sequence_position, sequence);
    setAlphabet(sequence.getAlphabet());
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::addIndividualSequenceInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (AlphabetMismatchException& ame)
  {
    throw AlphabetMismatchException("DataSet::addIndividualSequenceInGroup: sequence's alphabet doesn't match.", ame.getAlphabets()[0], ame.getAlphabets()[1]);
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

const Sequence& DataSet::getIndividualSequenceByNameInGroup(size_t group_position, size_t individual_position, const std::string& sequence_name) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualSequenceByNameInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return groups_[group_position]->getIndividualSequenceByName(individual_position, sequence_name);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualSequenceByNameInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
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

const Sequence& DataSet::getIndividualSequenceAtPositionInGroup(size_t group_position, size_t individual_position, size_t sequence_position) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualSequenceAtPositionInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return groups_[group_position]->getIndividualSequenceAtPosition(individual_position, sequence_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    if (string(ioobe.what()).find("individual_position") < string(ioobe.what()).size())
      throw IndexOutOfBoundsException("DataSet::getIndividualSequenceAtPositionInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
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

void DataSet::deleteIndividualSequenceByNameInGroup(size_t group_position, size_t individual_position, const std::string& sequence_name)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceByNameInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->deleteIndividualSequenceByName(individual_position, sequence_name);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceByNameInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
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

void DataSet::deleteIndividualSequenceAtPositionInGroup(size_t group_position, size_t individual_position, size_t sequence_position)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceAtPositionInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->deleteIndividualSequenceAtPosition(individual_position, sequence_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    if (string(ioobe.what()).find("individual_position") < string(ioobe.what()).size())
      throw IndexOutOfBoundsException("DataSet::deleteIndividualSequenceAtPositionInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
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

std::vector<std::string> DataSet::getIndividualSequencesNamesInGroup(size_t group_position, size_t individual_position) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualSequencesNamesInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return groups_[group_position]->getIndividualSequencesNames(individual_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualSequencesNamesInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualSequencesNamesInGroup: individual has no sequences.");
  }
}

/******************************************************************************/

size_t DataSet::getIndividualSequencePositionInGroup(size_t group_position, size_t individual_position, const std::string& sequence_name) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualSequencePositionInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return groups_[group_position]->getIndividualSequencePosition(individual_position, sequence_name);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualSequencePositionInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
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

size_t DataSet::getIndividualNumberOfSequencesInGroup(size_t group_position, size_t individual_position) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualNumberOfSequencesInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return groups_[group_position]->getIndividualNumberOfSequences(individual_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::getIndividualNumberOfSequencesInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (NullPointerException&)
  {
    throw NullPointerException("DataSet::getIndividualNumberOfSequencesInGroup: individual has no sequences.");
  }
}

/******************************************************************************/

void DataSet::setIndividualGenotypeInGroup(size_t group_position, size_t individual_position, const MultilocusGenotype& genotype)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualGenotypeInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->setIndividualGenotype(individual_position, genotype);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::setIndividualGenotypeInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

void DataSet::initIndividualGenotypeInGroup(size_t group_position, size_t individual_position)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::initIndividualGenotypeInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->initIndividualGenotype(individual_position, getNumberOfLoci());
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::initIndividualGenotypeInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
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

void DataSet::deleteIndividualGenotypeInGroup(size_t group_position, size_t individual_position)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::deleteIndividualGenotypeInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->deleteIndividualGenotype(individual_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("DataSet::deleteIndividualGenotypeInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

void DataSet::setIndividualMonolocusGenotypeInGroup(size_t group_position, size_t individual_position, size_t locus_position, const MonolocusGenotype& monogen)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->setIndividualMonolocusGenotype(individual_position, locus_position, monogen);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    if (string(ioobe.what()).find("individual_position") < string(ioobe.what()).size())
      throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
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

void DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup(size_t group_position, size_t individual_position, size_t locus_position, const std::vector<size_t> allele_keys)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    groups_[group_position]->setIndividualMonolocusGenotypeByAlleleKey(individual_position, locus_position, allele_keys);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    if (string(ioobe.what()).find("individual_position") < string(ioobe.what()).size())
      throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleKeyInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
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

void DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup(size_t group_position, size_t individual_position, size_t locus_position, const std::vector<std::string> allele_id)
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  const LocusInfo& locus_info = getLocusInfoAtPosition(locus_position);
  try
  {
    groups_[group_position]->setIndividualMonolocusGenotypeByAlleleId(individual_position, locus_position, allele_id, locus_info);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    if (string(ioobe.what()).find("individual_position") < string(ioobe.what()).size())
      throw IndexOutOfBoundsException("DataSet::setIndividualMonolocusGenotypeByAlleleIdInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
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

const MonolocusGenotype* DataSet::getIndividualMonolocusGenotypeInGroup(size_t group_position, size_t individual_position, size_t locus_position) const
{
  if (group_position >= getNumberOfGroups())
    throw IndexOutOfBoundsException("DataSet::getIndividualMonolocusGenotypeInGroup: group_position out of bounds.", group_position, 0, getNumberOfGroups());
  try
  {
    return &groups_[group_position]->getIndividualMonolocusGenotype(individual_position, locus_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    if (string(ioobe.what()).find("individual_position") < string(ioobe.what()).size())
      throw IndexOutOfBoundsException("DataSet::getIndividualMonolocusGenotypeInGroup: individual_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
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

// Dealing with AnalyzedSequences --------------------------

void DataSet::setAlphabet(const Alphabet* alpha)
{
  if (analyzedSequences_ == 0)
    analyzedSequences_ = new AnalyzedSequences();
  analyzedSequences_->setAlphabet(alpha);
}

/******************************************************************************/

void DataSet::setAlphabet(const std::string& alpha_type)
{
  if (analyzedSequences_ == 0)
    analyzedSequences_ = new AnalyzedSequences();
  analyzedSequences_->setAlphabet(alpha_type);
}

/******************************************************************************/

const Alphabet* DataSet::getAlphabet() const
{
  if (analyzedSequences_ != 0)
    return analyzedSequences_->getAlphabet();
  throw NullPointerException("DataSet::getAlphabet: no sequence data.");
}

/******************************************************************************/

std::string DataSet::getAlphabetType() const
{
  if (analyzedSequences_ != 0)
    return analyzedSequences_->getAlphabetType();
  throw NullPointerException("DataSet::getAlphabetType: no sequence data.");
}

/******************************************************************************/

// Dealing with AnalyzedLoci -------------------------------

void DataSet::setAnalyzedLoci(const AnalyzedLoci& analyzedLoci)
{
  if (analyzedLoci_ != 0)
  {
    try
    {
      deleteAnalyzedLoci();
    }
    catch (Exception& e)
    {
      throw Exception ("DataSet::setAnalyzedLoci: at least one individual has a genotype of the actual AnalyzedLoci.");
    }
  }
  analyzedLoci_ = new AnalyzedLoci(analyzedLoci);
}

/******************************************************************************/

void DataSet::initAnalyzedLoci(size_t number_of_loci)
{
  if (analyzedLoci_ != 0)
    throw Exception("DataSet::initAnalyzedLoci: analyzedLoci_ already initialyzed.");
  analyzedLoci_ = new AnalyzedLoci(number_of_loci);
}

/******************************************************************************/

const AnalyzedLoci* DataSet::getAnalyzedLoci() const
{
  if (analyzedLoci_ != 0)
    return analyzedLoci_;
  throw NullPointerException("DataSet::getAnalyzedLoci: no loci initialized.");
}

/******************************************************************************/

void DataSet::deleteAnalyzedLoci()
{
  if (analyzedLoci_ != 0)
    delete analyzedLoci_;
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

// Container extraction -----------------------------------

PolymorphismMultiGContainer* DataSet::getPolymorphismMultiGContainer() const
{
  PolymorphismMultiGContainer* pmgc = new PolymorphismMultiGContainer();
  for (size_t i = 0; i < getNumberOfGroups(); i++)
  {
    // nommer les groupes khalid
    string name = groups_[i]->getGroupName();
    pmgc->addGroupName(i, name);
    for (size_t j = 0; j < getNumberOfIndividualsInGroup(i); j++)
    {
      const Individual* tmp_ind = getIndividualAtPositionFromGroup(i, j);
      if (tmp_ind->hasGenotype())
      {
        const MultilocusGenotype& tmp_mg = tmp_ind->getGenotype();
        pmgc->addMultilocusGenotype(tmp_mg, i);
      }
    }
  }
  return pmgc;
}

/******************************************************************************/

PolymorphismMultiGContainer* DataSet::getPolymorphismMultiGContainer(const std::map<size_t, std::vector<size_t> >& selection) const
{
  PolymorphismMultiGContainer* pmgc = new PolymorphismMultiGContainer();
  for (map<size_t, vector<size_t> >::const_iterator it = selection.begin(); it != selection.end(); it++)
  {
    size_t i;
    try
    {
      i = getGroupPosition(it->first);
    }
    catch (GroupNotFoundException& gnfe)
    {
      throw gnfe;
    }
    string name = groups_[i]->getGroupName();
    pmgc->addGroupName(i, name);
    for (size_t j = 0; j < it->second.size(); j++)
    {
      const Individual* tmp_ind = 0;
      try
      {
        tmp_ind = getIndividualAtPositionFromGroup(i, j);
      }
      catch (IndexOutOfBoundsException& ioobe)
      {
        throw ioobe;
      }
      if (tmp_ind->hasGenotype())
      {
        const MultilocusGenotype& tmp_mg = tmp_ind->getGenotype();
        pmgc->addMultilocusGenotype(tmp_mg, i);
      }
    }
  }
  return pmgc;
}

/******************************************************************************/

PolymorphismSequenceContainer* DataSet::getPolymorphismSequenceContainer(const std::map<size_t, std::vector<size_t> >& selection, size_t sequence_position) const
{
  PolymorphismSequenceContainer* psc = new PolymorphismSequenceContainer(getAlphabet());
  for (map<size_t, vector<size_t> >::const_iterator it = selection.begin(); it != selection.end(); it++)
  {
    size_t i;
    try
    {
      i = getGroupPosition(it->first);
    }
    catch (GroupNotFoundException& gnfe)
    {
      delete psc;
      throw gnfe;
    }
    for (size_t j = 0; j < it->second.size(); j++)
    {
      const Individual* tmp_ind = 0;
      try
      {
        tmp_ind = getIndividualAtPositionFromGroup(i, j);
      }
      catch (IndexOutOfBoundsException& ioobe)
      {
        delete psc;
        throw ioobe;
      }
      if (tmp_ind->hasSequenceAtPosition(sequence_position))
      {
        const Sequence& tmp_seq = tmp_ind->getSequenceAtPosition(sequence_position);
        psc->addSequenceWithFrequency(tmp_seq, 1, false);
        psc->setGroupId(tmp_seq.getName(), it->first);
      }
    }
  }
  return psc;
}

/******************************************************************************/

// General tests ------------------------------------------

bool DataSet::hasSequenceData() const
{
  return analyzedSequences_ != 0;
}

/******************************************************************************/

bool DataSet::hasAlleleicData() const
{
  return analyzedLoci_ != 0;
}

/******************************************************************************/
