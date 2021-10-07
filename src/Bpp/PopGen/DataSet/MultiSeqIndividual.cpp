//
// File MultiSeqIndividual.cpp
// Author : Sylvain Gaillard
// Last modification : Tuesday August 03 2004
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

#include "MultiSeqIndividual.h"

using namespace bpp;
using namespace std;

// ** Class constructor: *******************************************************/

MultiSeqIndividual::MultiSeqIndividual() : id_(""),
  sex_(0),
  date_(0),
  coord_(0),
  locality_(0),
  sequences_(map<string, VectorSequenceContainer*>()),
  genotype_(0) {}

MultiSeqIndividual::MultiSeqIndividual(const std::string& id) : id_(id),
  sex_(0),
  date_(0),
  coord_(0),
  locality_(0),
  sequences_(map<string, VectorSequenceContainer*>()),
  genotype_(0) {}

MultiSeqIndividual::MultiSeqIndividual(
  const std::string& id,
  const Date& date,
  const Point2D<double>& coord,
  Locality<double>* locality,
  const unsigned short sex) :
  id_(id),
  sex_(sex),
  date_(new Date(date)),
  coord_(new Point2D<double>(coord)),
  locality_(locality),
  sequences_(map<string, VectorSequenceContainer*>()),
  genotype_(0) {}

MultiSeqIndividual::MultiSeqIndividual(const MultiSeqIndividual& ind) : id_(ind.getId()),
  sex_(ind.getSex()),
  date_(0),
  coord_(0),
  locality_(0),
  sequences_(map<string, VectorSequenceContainer*>()),
  genotype_(0)
{
  try
  {
    setDate(*ind.getDate());
  }
  catch (NullPointerException&)
  {
    date_ = 0;
  }
  try
  {
    setCoord(*ind.getCoord());
  }
  catch (NullPointerException&)
  {
    coord_ = 0;
  }
  try
  {
    setLocality(ind.getLocality());
  }
  catch (NullPointerException&)
  {
    locality_ = 0;
  }
  if (ind.hasSequences())
  {
    vector<string> keys = ind.getSequencesKeys();
    for (size_t i = 0; i < keys.size(); i++)
    {
      sequences_[keys[i]] = new VectorSequenceContainer(*const_cast<const VectorSequenceContainer*>(ind.getVectorSequenceContainer(keys[i])));
    }
  }
  genotype_ = ind.hasGenotype() ? new MultilocusGenotype(*ind.getGenotype()) : 0;
}

// ** Class destructor: *******************************************************/

MultiSeqIndividual::~MultiSeqIndividual()
{
  delete date_;
  delete coord_;
}

// ** Other methodes: *********************************************************/

MultiSeqIndividual& MultiSeqIndividual::operator=(const MultiSeqIndividual& ind)
{
  setId(ind.getId());
  setSex(ind.getSex());
  try
  {
    setDate(*ind.getDate());
  }
  catch (NullPointerException&)
  {
    date_ = 0;
  }
  try
  {
    setCoord(*ind.getCoord());
  }
  catch (NullPointerException&)
  {
    coord_ = 0;
  }
  try
  {
    setLocality(ind.getLocality());
  }
  catch (NullPointerException&)
  {
    locality_ = 0;
  }
  if (ind.hasSequences())
  {
    vector<string> keys = ind.getSequencesKeys();
    for (size_t i = 0; i < keys.size(); i++)
    {
      sequences_[keys[i]] = new VectorSequenceContainer(*const_cast<const VectorSequenceContainer*>(ind.getVectorSequenceContainer(keys[i])));
    }
  }
  genotype_ = ind.hasGenotype() ? new MultilocusGenotype(*ind.getGenotype()) : 0;
  return *this;
}

/******************************************************************************/

// Id
void MultiSeqIndividual::setId(const std::string id)
{
  id_ = id;
}

/******************************************************************************/

std::string MultiSeqIndividual::getId() const
{
  return id_;
}

/******************************************************************************/

// Sex
void MultiSeqIndividual::setSex(const unsigned short sex)
{
  sex_ = sex;
}

/******************************************************************************/

unsigned short MultiSeqIndividual::getSex() const
{
  return sex_;
}

/******************************************************************************/

// Date
void MultiSeqIndividual::setDate(const Date& date)
{
  if (!hasDate())
  {
    date_ = new Date(date);
  }
  else if (*date_ != date)
  {
    delete date_;
    date_ = new Date(date);
  }
}

/******************************************************************************/

const Date* MultiSeqIndividual::getDate() const
{
  if (hasDate())
    return new Date(*date_);
  else
    throw (NullPointerException("MultiSeqIndividual::getDate: no date associated to this individual."));
}

/******************************************************************************/

bool MultiSeqIndividual::hasDate() const
{
  return date_ != 0;
}

/******************************************************************************/

// Coord
void MultiSeqIndividual::setCoord(const Point2D<double>& coord)
{
  if (!hasCoord())
  {
    coord_ = new Point2D<double>(coord);
  }
  else if  (*coord_ != coord)
  {
    delete coord_;
    coord_ = new Point2D<double>(coord);
  }
}

/******************************************************************************/

void MultiSeqIndividual::setCoord(const double x, const double y)
{
  if (!hasCoord())
  {
    coord_ = new Point2D<double>(x, y);
  }
  else if (this->getX() != x || this->getY() != y)
  {
    delete coord_;
    coord_ = new Point2D<double>(x, y);
  }
}

/******************************************************************************/

const Point2D<double>* MultiSeqIndividual::getCoord() const
{
  if (hasCoord())
    return new Point2D<double>(*coord_);
  else
    throw (NullPointerException("MultiSeqIndividual::getCoord: no coord associated to this individual."));
}

/******************************************************************************/

bool MultiSeqIndividual::hasCoord() const
{
  return coord_ != 0;
}

/******************************************************************************/

void MultiSeqIndividual::setX(const double x)
{
  if (hasCoord())
    coord_->setX(x);
  else
    throw (NullPointerException("MultiSeqIndividual::setX: no coord associated to this individual."));
}

/******************************************************************************/

void MultiSeqIndividual::setY(const double y)
{
  if (hasCoord())
    coord_->setY(y);
  else
    throw (NullPointerException("MultiSeqIndividual::setY: no coord associated to this individual."));
}

/******************************************************************************/

double MultiSeqIndividual::getX() const
{
  if (hasCoord())
    return coord_->getX();
  else
    throw (NullPointerException("MultiSeqIndividual::getX: no coord associated to this individual."));
}

/******************************************************************************/

double MultiSeqIndividual::getY() const
{
  if (hasCoord())
    return coord_->getY();
  else
    throw (NullPointerException("MultiSeqIndividual::getY: no coord associated to this individual."));
}

/******************************************************************************/

// Locality
void MultiSeqIndividual::setLocality(const Locality<double>* locality)
{
  locality_ = locality;
}

/******************************************************************************/

const Locality<double>* MultiSeqIndividual::getLocality() const
{
  if (hasLocality())
    return locality_;
  else
    throw (NullPointerException("MultiSeqIndividual::getLocality: no locality associated to this individual."));
}

/******************************************************************************/

bool MultiSeqIndividual::hasLocality() const
{
  return locality_ != 0;
}

/******************************************************************************/

// Sequences
const VectorSequenceContainer* MultiSeqIndividual::getVectorSequenceContainer(const std::string& id) const
{
  map<string, VectorSequenceContainer*>::const_iterator it;
  it = sequences_.find(id);
  // Test existence of id in the map.
  if (it == sequences_.end())
  {
    string mes = "MultiSeqIndividual::getSequence: sequence set not found (" + id
                 + ").";
    throw (Exception(mes));
  }
  return const_cast<const VectorSequenceContainer*>(it->second);
}

/******************************************************************************/

void MultiSeqIndividual::addSequence(const std::string& id, const Sequence& sequence)
{
  try
  {
    sequences_[id]->addSequence(sequence);
  }
  catch (AlphabetMismatchException& ame)
  {
    throw (AlphabetMismatchException("MultiSeqIndividual::addSequence: alphabets don't match.", ame.getAlphabets()[0], ame.getAlphabets()[1]));
  }
  catch (Exception& e)
  {
    throw (BadIdentifierException("MultiSeqIndividual::addSequence: sequence's name already in use.", sequence.getName()));
  }
}

/******************************************************************************/

const Sequence& MultiSeqIndividual::getSequence(const std::string& id, const std::string& name) const
{
  map<string, VectorSequenceContainer*>::const_iterator it;
  it = sequences_.find(id);
  // Test existence of id in the map.
  if (it == sequences_.end())
  {
    string mes = "MultiSeqIndividual::getSequence: sequence set not found (" + id
                 + ").";
    throw (Exception(mes));
  }
  try
  {
    return const_cast<const VectorSequenceContainer*>(it->second)->getSequence(name);
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw (snfe);
  }
}

/******************************************************************************/

const Sequence& MultiSeqIndividual::getSequence(const std::string& id, size_t i) const
{
  map<string, VectorSequenceContainer*>::const_iterator it;
  it = sequences_.find(id);
  // Test existence of id in the map.
  if (it == sequences_.end())
  {
    string mes = "MultiSeqIndividual::getSequence: sequence set not found (" + id
                 + ").";
    throw (Exception(mes));
  }
  try
  {
    return const_cast<const VectorSequenceContainer*>(it->second)->getSequence(i);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw (ioobe);
  }
}

/******************************************************************************/

std::vector<std::string> MultiSeqIndividual::getSequencesKeys() const
{
  vector<string> keys;
  map<string, VectorSequenceContainer*>::const_iterator it;
  for (it = sequences_.begin(); it != sequences_.end(); it++)
  {
    keys.push_back(it->first);
  }
  return keys;
}

/******************************************************************************/

bool MultiSeqIndividual::hasSequences() const
{
  return sequences_.size() != 0;
}

/******************************************************************************/

size_t MultiSeqIndividual::getNumberOfSequenceSet() const
{
  return sequences_.size();
}

/******************************************************************************/

size_t MultiSeqIndividual::getNumberOfSequences(const std::string& id) const
{
  map<string, VectorSequenceContainer*>::const_iterator it;
  it = sequences_.find(id);
  // Test existence of id in the map.
  if (it == sequences_.end())
  {
    string mes = "MultiSeqIndividual::getSequence: sequence set not found (" + id
                 + ").";
    throw (Exception(mes));
  }

  return const_cast<const VectorSequenceContainer*>(it->second)->getNumberOfSequences();
}

/******************************************************************************/

// MultilocusGenotype

void MultiSeqIndividual::addGenotype(const MultilocusGenotype& genotype)
{
  genotype_ = new MultilocusGenotype(genotype);
}

/******************************************************************************/

const MultilocusGenotype* MultiSeqIndividual::getGenotype() const
{
  return genotype_;
}

/******************************************************************************/

bool MultiSeqIndividual::hasGenotype() const
{
  return genotype_ != 0;
}

/******************************************************************************/
