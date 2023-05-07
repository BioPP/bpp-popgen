//
// File Group.cpp
// Author : Sylvain Gaillard
// Last modification : Thursday July 29 2004
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

#include "Group.h"

using namespace bpp;
using namespace std;

// ** Class constructors: ******************************************************/
Group::Group(size_t group_id) : id_(group_id),
  name_(""),
  individuals_(vector<Individual*>()) {}

Group::Group(const Group& group) : id_(group.getGroupId()),
  name_(group.getGroupName()),
  // individuals_(vector<Individuals*>(group.getNumberOfIndividuals()))
  individuals_(vector<Individual*>())
{
  for (size_t i = 0; i < group.getNumberOfIndividuals(); i++)
  {
    addIndividual(group.getIndividualAtPosition(i));
  }
}

Group::Group(const Group& group, size_t group_id) : id_(group_id),
  name_(group.getGroupName()),
  individuals_(vector<Individual*>())
{
  for (size_t i = 0; i < group.getNumberOfIndividuals(); i++)
  {
    addIndividual(group.getIndividualAtPosition(i));
  }
}

// ** Class destructor: ********************************************************/

Group::~Group () {}

// ** Other methodes: **********************************************************/

Group& Group::operator=(const Group& group)
{
  setGroupId(group.getGroupId());
  for (size_t i = 0; i < group.getNumberOfIndividuals(); i++)
  {
    addIndividual(group.getIndividualAtPosition(i));
  }
  return *this;
}

void Group::setGroupId(size_t group_id)
{
  id_ = group_id;
}

void Group::setGroupName(const std::string& group_name)
{
  name_ = group_name;
}

void Group::addIndividual(const Individual& ind)
{
  try
  {
    getIndividualPosition(ind.getId());
    throw BadIdentifierException("Group::addIndividual: individual id already used.", ind.getId());
  }
  catch (BadIdentifierException& bie)
  {}
  individuals_.push_back(new Individual(ind));
}

void Group::addEmptyIndividual(const std::string& individual_id)
{
  for (size_t i = 0; i < getNumberOfIndividuals(); i++)
  {
    if (individuals_[i]->getId() == individual_id)
      throw BadIdentifierException("Group::addEmptyIndividual: individual_id already in use.", individual_id);
  }
  individuals_.push_back(new Individual(individual_id));
}

size_t Group::getIndividualPosition(const std::string& individual_id) const
{
  for (size_t i = 0; i < getNumberOfIndividuals(); i++)
  {
    if (individuals_[i]->getId() == individual_id)
      return i;
  }
  throw IndividualNotFoundException("Group::getIndividualPosition: individual_id not found.", individual_id);
}

std::unique_ptr<Individual> Group::removeIndividualById(const std::string& individual_id)
{
  try
  {
    size_t indPos = getIndividualPosition(individual_id);
    unique_ptr<Individual> ind(individuals_[indPos]);
    individuals_.erase(individuals_.begin() + static_cast<ptrdiff_t>(indPos));
    return ind;
  }
  catch (IndividualNotFoundException& infe)
  {
    throw IndividualNotFoundException("Group::removeIndividualById: individual_id not found.", individual_id);
  }
}

std::unique_ptr<Individual> Group::removeIndividualAtPosition(size_t individual_position)
{
  if (individual_position >= individuals_.size())
    throw IndexOutOfBoundsException("Group::removeIndividualAtPosition.", individual_position, 0, individuals_.size());
  unique_ptr<Individual> ind(individuals_[individual_position]);
  individuals_.erase(individuals_.begin() + static_cast<ptrdiff_t>(individual_position));
  return ind;
}

void Group::deleteIndividualById(const std::string& individual_id)
{
  try
  {
    removeIndividualById(individual_id);
  }
  catch (IndividualNotFoundException& infe)
  {
    throw IndividualNotFoundException("Group::deleteIndividualById: individual_id not found.", individual_id);
  }
}

void Group::deleteIndividualAtPosition(size_t individual_position)
{
  try
  {
    removeIndividualAtPosition(individual_position);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Group::deleteIndividualAtPosition.", individual_position, 0, getNumberOfIndividuals());
  }
}

void Group::clear()
{
  for (size_t i = 0; i < individuals_.size(); i++)
  {
    delete (individuals_[i]);
  }
  individuals_.clear();
}

const Individual& Group::getIndividualById(const std::string& individual_id) const
{
  for (size_t i = 0; i < individuals_.size(); i++)
  {
    if (individuals_[i]->getId() == individual_id)
      return getIndividualAtPosition(i);
  }
  throw IndividualNotFoundException("Group::getIndividualById: individual_id not found.", individual_id);
}

const Individual& Group::getIndividualAtPosition(size_t individual_position) const
{
  if (individual_position >= individuals_.size())
    throw IndexOutOfBoundsException("Group::getIndividualAtPosition: individual_position out of bounds.", individual_position, 0, individuals_.size());
  return *individuals_[individual_position];
}

size_t Group::getNumberOfIndividuals() const
{
  return individuals_.size();
}

size_t Group::getMaxNumberOfSequences() const
{
  size_t maxnum = 0;
  for (size_t i = 0; i < getNumberOfIndividuals(); i++)
  {
    vector<size_t> seqpos = individuals_[i]->getSequencePositions();
    for (size_t j = 0; j < seqpos.size(); j++)
    {
      if (maxnum < seqpos[j])
        maxnum = seqpos[j];
    }
  }
  return maxnum + 1;
}

// -- Dealing with individual's properties -----------------

void Group::setIndividualSexAtPosition(size_t individual_position, const unsigned short sex)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualSexAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  individuals_[individual_position]->setSex(sex);
}

unsigned short Group::getIndividualSexAtPosition(size_t individual_position) const
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualSexAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  return individuals_[individual_position]->getSex();
}

void Group::setIndividualDateAtPosition(size_t individual_position, const Date& date)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualDateAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  individuals_[individual_position]->setDate(date);
}

const Date& Group::getIndividualDateAtPosition(size_t individual_position) const
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualDateAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individual_position]->getDate();
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualDateAtPosition: individual has no date.");
  }
}

void Group::setIndividualCoordAtPosition(size_t individual_position, const Point2D<double>& coord)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualCoordAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  individuals_[individual_position]->setCoord(coord);
}

const Point2D<double>& Group::getIndividualCoordAtPosition(size_t individual_position) const
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualCoordAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individual_position]->getCoord();
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualCoordAtPosition: individual has no coordinates.");
  }
}

void Group::setIndividualLocalityAtPosition(size_t individual_position, const Locality<double>* locality)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualLocalityAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  individuals_[individual_position]->setLocality(locality);
}

const Locality<double>& Group::getIndividualLocalityAtPosition(size_t individual_position) const
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualLocalityAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    return *individuals_[individual_position]->getLocality();
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualLocalityAtPosition: individuals has no locality.");
  }
}

void Group::addIndividualSequenceAtPosition(size_t individual_position, size_t sequence_position, const Sequence& sequence)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::addIndividualSequenceAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individual_position]->addSequence(sequence_position, sequence);
  }
  catch (AlphabetMismatchException& ame)
  {
    throw AlphabetMismatchException("Group::addIndividualSequenceAtPosition: sequence's alphabet doesn't match.", ame.getAlphabets()[0], ame.getAlphabets()[1]);
  }
  catch (BadIdentifierException& bie)
  {
    throw BadIdentifierException("Group::addIndividualSequenceAtPosition: sequence's name already in use.", bie.getIdentifier());
  }
  catch (BadIntegerException& bie)
  {
    throw BadIntegerException("Group::addIndividualSequenceAtPosition: sequence_position already in use.", bie.getBadInteger());
  }
}

const Sequence& Group::getIndividualSequenceByName(size_t individual_position, const string& sequence_name) const
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualSequenceByName: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individual_position]->getSequenceByName(sequence_name);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualSequenceByName: no sequence data in individual.");
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Group::getIndividualSequenceByName: sequence_name not found.", snfe.getSequenceId());
  }
}

const Sequence& Group::getIndividualSequenceAtPosition(size_t individual_position, size_t sequence_position) const
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individual_position]->getSequenceAtPosition(sequence_position);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualSequenceAtPosition: no sequence data in individual.");
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Group::getIndividualSequenceAtPosition: sequence_position not found.", snfe.getSequenceId());
  }
}

void Group::deleteIndividualSequenceByName(size_t individual_position, const string& sequence_name)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::deleteIndividualSequenceByName: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individual_position]->deleteSequenceByName(sequence_name);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::deleteSequenceByName: no sequence data in individual.");
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Group::deleteSequenceByName: sequence_name not found.", snfe.getSequenceId());
  }
}

void Group::deleteIndividualSequenceAtPosition(size_t individual_position, size_t sequence_position)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::deleteIndividualSequenceAtPosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individual_position]->deleteSequenceAtPosition(sequence_position);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::deleteSequenceAtPosition: no sequence data in individual.");
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Group::deleteSequenceAtPosition: sequence_position not found.", snfe.getSequenceId());
  }
}

bool Group::hasIndividualSequences(size_t individual_position) const
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::hasIndividualSequences: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  return individuals_[individual_position]->hasSequences();
}

vector<string> Group::getIndividualSequencesNames(size_t individual_position) const
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualSequencesNames: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individual_position]->getSequenceNames();
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getSequencesNames: no sequence data in individual.");
  }
}

size_t Group::getIndividualSequencePosition(size_t individual_position, const string& sequence_name) const
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualSequencePosition: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individual_position]->getSequencePosition(sequence_name);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getSequencePosition: no sequence data in individual.");
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Group::getSequencePosition: sequence_name not found.", snfe.getSequenceId());
  }
}

size_t Group::getIndividualNumberOfSequences(size_t individual_position) const
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualNumberOfSequences: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individual_position]->getNumberOfSequences();
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualNumberOfSequences: no sequence data in individual.");
  }
}

void Group::setIndividualSequences(size_t individual_position, const OrderedSequenceContainer& osc)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualSequences: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  individuals_[individual_position]->setSequences(osc);
}

void Group::setIndividualGenotype(size_t individual_position, const MultilocusGenotype& genotype)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualGenotype: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  individuals_[individual_position]->setGenotype(genotype);
}

void Group::initIndividualGenotype(size_t individual_position, size_t loci_number)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::initIndividualGenotype: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individual_position]->initGenotype(loci_number);
  }
  catch (BadIntegerException& bie)
  {
    throw BadIntegerException("Group::initIndividualGenotype: loci_number must be > 0.", bie.getBadInteger());
  }
  catch (Exception&)
  {
    throw Exception("Group::initIndividualGenotype: individual already has a genotype.");
  }
}

void Group::deleteIndividualGenotype(size_t individual_position)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::deleteIndividualGenotype: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  individuals_[individual_position]->deleteGenotype();
}

bool Group::hasIndividualGenotype(size_t individual_position) const
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::hasIndividualGenotype: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  return individuals_[individual_position]->hasGenotype();
}

void Group::setIndividualMonolocusGenotype(size_t individual_position, size_t locus_position, const MonolocusGenotype& monogen)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotype: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individual_position]->setMonolocusGenotype(locus_position, monogen);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::setIndividualMonolocusGenotype: individual has no genotype.");
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotype: locus_position excedes the number of locus.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

void Group::setIndividualMonolocusGenotypeByAlleleKey(size_t individual_position, size_t locus_position, const std::vector<size_t>& allele_keys)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleKey: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individual_position]->setMonolocusGenotypeByAlleleKey(locus_position, allele_keys);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::setIndividualMonolocusGenotypeByAlleleKey: individual has no genotype.");
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleKey: locus_position excedes the number of locus.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (Exception&)
  {
    throw Exception("Group::setIndividualMonolocusGenotypeByAlleleKey: no key in allele_keys.");
  }
}

void Group::setIndividualMonolocusGenotypeByAlleleId(size_t individual_position, size_t locus_position, const std::vector<std::string>& allele_id, const LocusInfo& locus_info)
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleId: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individual_position]->setMonolocusGenotypeByAlleleId(locus_position, allele_id, locus_info);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::setIndividualMonolocusGenotypeByAlleleId: individual has no genotype.");
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleId: locus_position excedes the number of locus.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (AlleleNotFoundException& anfe)
  {
    throw AlleleNotFoundException("Group::setIndividualMonolocusGenotypeByAlleleId: id not found.", anfe.getIdentifier());
  }
}

const MonolocusGenotype&  Group::getIndividualMonolocusGenotype(size_t individual_position, size_t locus_position) const
{
  if (individual_position >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualMonolocusGenotype: individual_position out of bounds.", individual_position, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individual_position]->getMonolocusGenotype(locus_position);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualMonolocusGenotype: individual has no genotype.");
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Group::getIndividualMonolocusGenotype: locus_position excedes the number of locus.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

bool Group::hasSequenceData() const
{
  for (size_t i = 0; i < getNumberOfIndividuals(); i++)
  {
    if (hasIndividualSequences(i))
      return true;
  }
  return false;
}

const Alphabet* Group::getAlphabet() const
{
  for (size_t i = 0; i < getNumberOfIndividuals(); i++)
  {
    if (hasIndividualSequences(i))
      return individuals_[i]->getSequenceAlphabet();
  }
  throw NullPointerException("Group::getAlphabet: individual has no sequence data.");
}

size_t Group::getGroupSizeForLocus(size_t locus_position) const
{
  size_t count = 0;
  for (size_t i = 0; i < individuals_.size(); i++)
  {
    if (individuals_[i]->hasGenotype() && !individuals_[i]->getGenotype().isMonolocusGenotypeMissing(locus_position))
      count++;
  }
  return count;
}

size_t Group::getGroupSizeForSequence(size_t sequence_position) const
{
  size_t count = 0;
  for (size_t i = 0; i < individuals_.size(); i++)
  {
    if (individuals_[i]->hasSequences())
    {
      try
      {
        individuals_[i]->getSequenceAtPosition(sequence_position);
        count++;
      }
      catch (...)
      {}
    }
  }
  return count;
}
