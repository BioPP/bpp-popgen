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

// ** Other methods: **********************************************************/

void Group::addIndividual(const Individual& ind)
{
  try
  {
    getIndividualPosition(ind.getId());
    throw BadIdentifierException("Group::addIndividual: individual id already used.", ind.getId());
  }
  catch (BadIdentifierException& bie)
  {}
  individuals_.push_back(make_unique<Individual>(ind));
}

void Group::addEmptyIndividual(const std::string& individualId)
{
  for (size_t i = 0; i < getNumberOfIndividuals(); ++i)
  {
    if (individuals_[i]->getId() == individualId)
      throw BadIdentifierException("Group::addEmptyIndividual: individualId already in use.", individualId);
  }
  individuals_.push_back(make_unique<Individual>(individualId));
}

size_t Group::getIndividualPosition(const std::string& individualId) const
{
  for (size_t i = 0; i < getNumberOfIndividuals(); ++i)
  {
    if (individuals_[i]->getId() == individualId)
      return i;
  }
  throw IndividualNotFoundException("Group::getIndividualPosition: individualId not found.", individualId);
}

unique_ptr<Individual> Group::removeIndividualById(const std::string& individualId)
{
  try
  {
    size_t indPos = getIndividualPosition(individualId);
    unique_ptr<Individual> ind = move(individuals_[indPos]);
    individuals_.erase(individuals_.begin() + static_cast<ptrdiff_t>(indPos));
    return ind;
  }
  catch (IndividualNotFoundException& infe)
  {
    throw IndividualNotFoundException("Group::removeIndividualById: individualId not found.", individualId);
  }
}

std::unique_ptr<Individual> Group::removeIndividualAtPosition(size_t individualPosition)
{
  if (individualPosition >= individuals_.size())
    throw IndexOutOfBoundsException("Group::removeIndividualAtPosition.", individualPosition, 0, individuals_.size());
  unique_ptr<Individual> ind = move(individuals_[individualPosition]);
  individuals_.erase(individuals_.begin() + static_cast<ptrdiff_t>(individualPosition));
  return ind;
}

void Group::deleteIndividualById(const std::string& individualId)
{
  try
  {
    removeIndividualById(individualId);
  }
  catch (IndividualNotFoundException& infe)
  {
    throw IndividualNotFoundException("Group::deleteIndividualById: individualId not found.", individualId);
  }
}

void Group::deleteIndividualAtPosition(size_t individualPosition)
{
  try
  {
    removeIndividualAtPosition(individualPosition);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Group::deleteIndividualAtPosition.", individualPosition, 0, getNumberOfIndividuals());
  }
}

const Individual& Group::getIndividualById(const std::string& individualId) const
{
  for (size_t i = 0; i < individuals_.size(); ++i)
  {
    if (individuals_[i]->getId() == individualId)
      return getIndividualAtPosition(i);
  }
  throw IndividualNotFoundException("Group::getIndividualById: individualId not found.", individualId);
}

const Individual& Group::getIndividualAtPosition(size_t individualPosition) const
{
  if (individualPosition >= individuals_.size())
    throw IndexOutOfBoundsException("Group::getIndividualAtPosition: individualPosition out of bounds.", individualPosition, 0, individuals_.size());
  return *individuals_[individualPosition];
}

size_t Group::getMaxNumberOfSequences() const
{
  size_t maxnum = 0;
  for (size_t i = 0; i < getNumberOfIndividuals(); ++i)
  {
    vector<size_t> seqpos = individuals_[i]->getSequencePositions();
    for (size_t j = 0; j < seqpos.size(); ++j)
    {
      if (maxnum < seqpos[j])
        maxnum = seqpos[j];
    }
  }
  return maxnum + 1;
}

// -- Dealing with individual's properties -----------------

void Group::setIndividualSexAtPosition(size_t individualPosition, const unsigned short sex)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualSexAtPosition: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  individuals_[individualPosition]->setSex(sex);
}

unsigned short Group::getIndividualSexAtPosition(size_t individualPosition) const
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualSexAtPosition: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  return individuals_[individualPosition]->getSex();
}

void Group::setIndividualDateAtPosition(size_t individualPosition, const Date& date)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualDateAtPosition: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  individuals_[individualPosition]->setDate(date);
}

const Date& Group::getIndividualDateAtPosition(size_t individualPosition) const
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualDateAtPosition: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individualPosition]->date();
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualDateAtPosition: individual has no date.");
  }
}

void Group::setIndividualCoordAtPosition(size_t individualPosition, const Point2D<double>& coord)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualCoordAtPosition: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  individuals_[individualPosition]->setCoord(coord);
}

const Point2D<double>& Group::getIndividualCoordAtPosition(size_t individualPosition) const
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualCoordAtPosition: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individualPosition]->coord();
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualCoordAtPosition: individual has no coordinates.");
  }
}

void Group::setIndividualLocalityAtPosition(size_t individualPosition, shared_ptr<const Locality<double>> locality)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualLocalityAtPosition: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  individuals_[individualPosition]->setLocality(locality);
}

std::shared_ptr<const Locality<double>> Group::getIndividualLocalityAtPosition(size_t individualPosition) const
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualLocalityAtPosition: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individualPosition]->getLocality();
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualLocalityAtPosition: individuals has no locality.");
  }
}

void Group::addIndividualSequenceAtPosition(size_t individualPosition, size_t sequencePosition, unique_ptr<Sequence>& sequence)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::addIndividualSequenceAtPosition: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individualPosition]->addSequence(sequencePosition, sequence);
  }
  catch (AlphabetMismatchException& ame)
  {
    throw AlphabetMismatchException("Group::addIndividualSequenceAtPosition: sequence's alphabet doesn't match.", ame.getFirstAlphabet(), ame.getSecondAlphabet());
  }
  catch (BadIdentifierException& bie)
  {
    throw BadIdentifierException("Group::addIndividualSequenceAtPosition: sequence's name already in use.", bie.getIdentifier());
  }
  catch (BadIntegerException& bie)
  {
    throw BadIntegerException("Group::addIndividualSequenceAtPosition: sequencePosition already in use.", bie.getBadInteger());
  }
}

const Sequence& Group::getIndividualSequenceByName(size_t individualPosition, const string& sequenceName) const
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualSequenceByName: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individualPosition]->sequenceByName(sequenceName);
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

const Sequence& Group::getIndividualSequenceAtPosition(size_t individualPosition, size_t sequencePosition) const
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualAtPosition: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individualPosition]->sequenceAtPosition(sequencePosition);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualSequenceAtPosition: no sequence data in individual.");
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Group::getIndividualSequenceAtPosition: sequencePosition not found.", snfe.getSequenceId());
  }
}

void Group::deleteIndividualSequenceByName(size_t individualPosition, const string& sequence_name)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::deleteIndividualSequenceByName: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individualPosition]->deleteSequenceByName(sequence_name);
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

void Group::deleteIndividualSequenceAtPosition(size_t individualPosition, size_t sequencePosition)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::deleteIndividualSequenceAtPosition: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individualPosition]->deleteSequenceAtPosition(sequencePosition);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::deleteSequenceAtPosition: no sequence data in individual.");
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Group::deleteSequenceAtPosition: sequencePosition not found.", snfe.getSequenceId());
  }
}

bool Group::hasIndividualSequences(size_t individualPosition) const
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::hasIndividualSequences: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  return individuals_[individualPosition]->hasSequences();
}

vector<string> Group::getIndividualSequencesNames(size_t individualPosition) const
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualSequencesNames: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individualPosition]->getSequenceNames();
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getSequencesNames: no sequence data in individual.");
  }
}

size_t Group::getIndividualSequencePosition(size_t individualPosition, const string& sequence_name) const
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualSequencePosition: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individualPosition]->getSequencePosition(sequence_name);
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

size_t Group::getIndividualNumberOfSequences(size_t individualPosition) const
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualNumberOfSequences: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individualPosition]->getNumberOfSequences();
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualNumberOfSequences: no sequence data in individual.");
  }
}

void Group::setIndividualSequences(size_t individualPosition, const SequenceContainerInterface& sc)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualSequences: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  individuals_[individualPosition]->setSequences(sc);
}

void Group::setIndividualGenotype(size_t individualPosition, const MultilocusGenotype& genotype)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualGenotype: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  individuals_[individualPosition]->setGenotype(genotype);
}

void Group::initIndividualGenotype(size_t individualPosition, size_t lociNumber)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::initIndividualGenotype: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individualPosition]->initGenotype(lociNumber);
  }
  catch (BadIntegerException& bie)
  {
    throw BadIntegerException("Group::initIndividualGenotype: lociNumber must be > 0.", bie.getBadInteger());
  }
  catch (Exception&)
  {
    throw Exception("Group::initIndividualGenotype: individual already has a genotype.");
  }
}

void Group::deleteIndividualGenotype(size_t individualPosition)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::deleteIndividualGenotype: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  individuals_[individualPosition]->deleteGenotype();
}

bool Group::hasIndividualGenotype(size_t individualPosition) const
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::hasIndividualGenotype: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  return individuals_[individualPosition]->hasGenotype();
}

void Group::setIndividualMonolocusGenotype(
    size_t individualPosition,
    size_t locusPosition,
    const MonolocusGenotypeInterface& monogen)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotype: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individualPosition]->setMonolocusGenotype(locusPosition, monogen);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::setIndividualMonolocusGenotype: individual has no genotype.");
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotype: locusPosition excedes the number of locus.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

void Group::setIndividualMonolocusGenotypeByAlleleKey(
    size_t individualPosition, 
    size_t locusPosition,
    const vector<size_t>& alleleKeys)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleKey: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individualPosition]->setMonolocusGenotypeByAlleleKey(locusPosition, alleleKeys);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::setIndividualMonolocusGenotypeByAlleleKey: individual has no genotype.");
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleKey: locusPosition excedes the number of locus.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (Exception&)
  {
    throw Exception("Group::setIndividualMonolocusGenotypeByAlleleKey: no key in allele_keys.");
  }
}

void Group::setIndividualMonolocusGenotypeByAlleleId(size_t individualPosition, size_t locusPosition, const std::vector<std::string>& allele_id, const LocusInfo& locus_info)
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleId: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    individuals_[individualPosition]->setMonolocusGenotypeByAlleleId(locusPosition, allele_id, locus_info);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::setIndividualMonolocusGenotypeByAlleleId: individual has no genotype.");
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Group::setIndividualMonolocusGenotypeByAlleleId: locusPosition excedes the number of locus.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (AlleleNotFoundException& anfe)
  {
    throw AlleleNotFoundException("Group::setIndividualMonolocusGenotypeByAlleleId: id not found.", anfe.getIdentifier());
  }
}

const MonolocusGenotypeInterface& Group::getIndividualMonolocusGenotype(
    size_t individualPosition, 
    size_t locusPosition) const
{
  if (individualPosition >= getNumberOfIndividuals())
    throw IndexOutOfBoundsException("Group::getIndividualMonolocusGenotype: individualPosition out of bounds.", individualPosition, 0, getNumberOfIndividuals());
  try
  {
    return individuals_[individualPosition]->getMonolocusGenotype(locusPosition);
  }
  catch (NullPointerException& npe)
  {
    throw NullPointerException("Group::getIndividualMonolocusGenotype: individual has no genotype.");
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Group::getIndividualMonolocusGenotype: locusPosition excedes the number of locus.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

bool Group::hasSequenceData() const
{
  for (size_t i = 0; i < getNumberOfIndividuals(); ++i)
  {
    if (hasIndividualSequences(i))
      return true;
  }
  return false;
}

shared_ptr<const Alphabet> Group::getAlphabet() const
{
  for (size_t i = 0; i < getNumberOfIndividuals(); ++i)
  {
    if (hasIndividualSequences(i))
      return individuals_[i]->getSequenceAlphabet();
  }
  throw NullPointerException("Group::getAlphabet: individual has no sequence data.");
}

size_t Group::getGroupSizeForLocus(size_t locusPosition) const
{
  size_t count = 0;
  for (size_t i = 0; i < individuals_.size(); ++i)
  {
    if (individuals_[i]->hasGenotype() && !individuals_[i]->getGenotype().isMonolocusGenotypeMissing(locusPosition))
      count++;
  }
  return count;
}

size_t Group::getGroupSizeForSequence(size_t sequencePosition) const
{
  size_t count = 0;
  for (size_t i = 0; i < individuals_.size(); ++i)
  {
    if (individuals_[i]->hasSequences())
    {
      try
      {
        individuals_[i]->sequenceAtPosition(sequencePosition);
        count++;
      }
      catch (...)
      {}
    }
  }
  return count;
}
