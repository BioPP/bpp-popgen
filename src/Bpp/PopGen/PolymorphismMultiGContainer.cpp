//
// File PolymorphismMultiGContainer.cpp
// Author : Sylvain Gaillard
//          Khalid Belkhir
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

#include "PolymorphismMultiGContainer.h"

using namespace bpp;
using namespace std;

// ** Constructors : **********************************************************/

PolymorphismMultiGContainer::PolymorphismMultiGContainer() : multilocusGenotypes_(std::vector<MultilocusGenotype*>()),
  groups_(std::vector<size_t>()),
  groups_names_(std::map<size_t, std::string>()) {}

PolymorphismMultiGContainer::PolymorphismMultiGContainer(const PolymorphismMultiGContainer& pmgc) : multilocusGenotypes_(std::vector<MultilocusGenotype*>(pmgc.size())),
  groups_(std::vector<size_t>(pmgc.size())),
  groups_names_(std::map<size_t, std::string>())
{
  for (size_t i = 0; i < pmgc.size(); i++)
  {
    multilocusGenotypes_[i] = new MultilocusGenotype(*pmgc.getMultilocusGenotype(i));
    groups_[i] = pmgc.getGroupId(i);
  }
  set<size_t> grp_ids = pmgc.getAllGroupsIds();
  for (set<size_t>::iterator it = grp_ids.begin(); it != grp_ids.end(); it++)
  {
    size_t id = *it;
    string name = pmgc.getGroupName(id);
    groups_names_[id] = name;
  }
}

// ** Destructor : ************************************************************/

PolymorphismMultiGContainer::~PolymorphismMultiGContainer()
{
  clear();
}

// ** Other methodes : ********************************************************/

PolymorphismMultiGContainer& PolymorphismMultiGContainer::operator=(const PolymorphismMultiGContainer& pmgc)
{
  clear();
  for (size_t i = 0; i < pmgc.size(); i++)
  {
    multilocusGenotypes_.push_back(new MultilocusGenotype(*pmgc.getMultilocusGenotype(i)));
    groups_.push_back(pmgc.getGroupId(i));
  }
  set<size_t> grp_ids = pmgc.getAllGroupsIds();
  for (set<size_t>::iterator it = grp_ids.begin(); it != grp_ids.end(); it++)
  {
    size_t id = *it;
    string name = pmgc.getGroupName(id);
    groups_names_[id] = name;
  }

  return *this;
}

/******************************************************************************/

void PolymorphismMultiGContainer::addMultilocusGenotype(const MultilocusGenotype& mg, size_t group)
{
  multilocusGenotypes_.push_back(new MultilocusGenotype(mg));
  groups_.push_back(group);
  map<size_t, string>::const_iterator it = groups_names_.find(group);
  if (!(it != groups_names_.end()) )
  {
    // ajouter ce groupe avec un nom vide
    groups_names_[group] = "";
  }
}

/******************************************************************************/

const MultilocusGenotype* PolymorphismMultiGContainer::getMultilocusGenotype(size_t position) const throw (IndexOutOfBoundsException)
{
  if (position >= size())
    throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getMultilocusGenotype: position out of bounds.", position, 0, size() - 1);
  return multilocusGenotypes_[position];
}

/******************************************************************************/

MultilocusGenotype* PolymorphismMultiGContainer::removeMultilocusGenotype(size_t position) throw (IndexOutOfBoundsException)
{
  if (position >= size())
    throw IndexOutOfBoundsException("PolymorphismMultiGContainer::removeMultilocusGenotype: position out of bounds.", position, 0, size() - 1);
  MultilocusGenotype* tmp_mg = multilocusGenotypes_[position];
  multilocusGenotypes_.erase(multilocusGenotypes_.begin() + position);
  groups_.erase(groups_.begin() + position);
  return tmp_mg;
}

/******************************************************************************/

void PolymorphismMultiGContainer::deleteMultilocusGenotype(size_t position) throw (IndexOutOfBoundsException)
{
  if (position >= size())
    throw IndexOutOfBoundsException("PolymorphismMultiGContainer::deleteMultilocusGenotype: position out of bounds.", position, 0, size() - 1);
  delete multilocusGenotypes_[position];
  multilocusGenotypes_.erase(multilocusGenotypes_.begin() + position);
  groups_.erase(groups_.begin() + position);
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

size_t PolymorphismMultiGContainer::getNumberOfLoci() const throw (Exception)
{
  if (!isAligned())
    throw Exception("MultilocusGenotypes are not aligned.");
  if (size() < 1)
    return 0;
  return multilocusGenotypes_[0]->size();
}

/******************************************************************************/

size_t PolymorphismMultiGContainer::getGroupId(size_t position) const throw (IndexOutOfBoundsException)
{
  if (position >= size())
    throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getGroupId: position out of bounds.", position, 0, size() - 1);
  return groups_[position];
}

/******************************************************************************/

void PolymorphismMultiGContainer::setGroupId(size_t position, size_t group_id) throw (IndexOutOfBoundsException)
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
  vector<string> grps_names;
  map<size_t, string>::const_iterator it;
  for (it = groups_names_.begin(); it != groups_names_.end(); it++)
  {
    string name = it->second;
    if (!name.empty())
      grps_names.push_back(name);
    else
      grps_names.push_back(TextTools::toString(it->first) );
  }

  return grps_names;
}

/******************************************************************************/

bool PolymorphismMultiGContainer::groupExists(size_t group) const
{
  for (size_t i = 0; i < size(); i++)
  {
    if (groups_[i] == group)
      return true;
  }
  return false;
}

/******************************************************************************/

size_t PolymorphismMultiGContainer::getNumberOfGroups() const
{
  return getAllGroupsIds().size();
}

/******************************************************************************/

size_t PolymorphismMultiGContainer::getGroupSize(size_t group) const
{
  size_t counter = 0;
  for (size_t i = 0; i < size(); i++)
  {
    if (groups_[i] == group)
      counter++;
  }
  return counter;
}

/******************************************************************************/

std::string PolymorphismMultiGContainer::getGroupName(size_t group_id) const throw (GroupNotFoundException)
{
  string name = TextTools::toString(group_id); // par defaut on retourne le n° de groupe
  map<size_t, string>::const_iterator it = groups_names_.find(group_id);
  if (it != groups_names_.end() )
    name = it->second;
  else
    throw GroupNotFoundException("PolymorphismMultiGContainer::getGroupName: group not found.", group_id);
  return name;
}

/******************************************************************************/

void PolymorphismMultiGContainer::setGroupName(size_t group_id, std::string name)  throw (GroupNotFoundException)
{
  map<size_t, string>::iterator it = groups_names_.find(group_id);
  if (it != groups_names_.end() )
    it->second = name;
  else
    throw GroupNotFoundException("PolymorphismMultiGContainer::getGroupName: group not found.", group_id);
  return;
}

/******************************************************************************/

void PolymorphismMultiGContainer::addGroupName(size_t group_id, std::string name)
{
  groups_names_[group_id] = name;
  return;
}

/******************************************************************************/

size_t PolymorphismMultiGContainer::getLocusGroupSize(size_t group, size_t locus_position) const
{
  size_t counter = 0;
  for (size_t i = 0; i < size(); i++)
  {
    try
    {
      if (groups_[i] == group && !multilocusGenotypes_[i]->isMonolocusGenotypeMissing(locus_position))
        counter++;
    }
    catch (IndexOutOfBoundsException& ioobe)
    {
      throw IndexOutOfBoundsException("PolymorphismMultiGContainer::getGroupSize: locus_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
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
  for (size_t i = 0; i < multilocusGenotypes_.size(); i++)
  {
    delete multilocusGenotypes_[i];
  }
  multilocusGenotypes_.clear();
  groups_.clear();
  groups_names_.clear();
}

/******************************************************************************/

