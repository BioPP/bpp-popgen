// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include "Individual.h"

using namespace bpp;
using namespace std;

// ** Class constructor: *******************************************************/
Individual::Individual() : 
  id_(""),
  sex_(0),
  date_(),
  coord_(),
  locality_(0),
  sequences_(),
  genotype_() {}

Individual::Individual(const std::string& id) : id_(id),
  sex_(0),
  date_(),
  coord_(),
  locality_(0),
  sequences_(),
  genotype_() {}

Individual::Individual(const string& id,
                       const Date& date,
                       const Point2D<double>& coord,
                       shared_ptr<Locality<double>> locality,
                       const unsigned short sex) :
  id_(id),
  sex_(sex),
  date_(new Date(date)),
  coord_(new Point2D<double>(coord)),
  locality_(locality),
  sequences_(),
  genotype_() {}

Individual::Individual(const Individual& ind) :
  id_(ind.id_),
  sex_(ind.sex_),
  date_(),
  coord_(),
  locality_(ind.locality_),
  sequences_(),
  genotype_()
{
  try
  {
    setDate(ind.date());
  }
  catch (...)
  {}
  try
  {
    setCoord(ind.coord());
  }
  catch (...)
  {}
  try
  {
    setSequences(ind.sequences());
  }
  catch (...)
  {}
  if (ind.hasGenotype())
    genotype_.reset(new MultilocusGenotype(ind.getGenotype()));
}

// ** Class destructor: *******************************************************/
Individual::~Individual () {}

// ** Other methodes: *********************************************************/

Individual& Individual::operator=(const Individual& ind)
{
  setId(ind.getId());
  setSex(ind.getSex());
  try
  {
    setDate(ind.date());
  }
  catch (NullPointerException&)
  {
    date_.reset();
  }
  try
  {
    setCoord(ind.coord());
  }
  catch (NullPointerException&)
  {
    coord_.reset();
  }
  try
  {
    setLocality(ind.getLocality());
  }
  catch (NullPointerException&)
  {
    locality_ = 0;
  }
  try
  {
    setSequences(ind.sequences());
  }
  catch (NullPointerException&)
  {
    sequences_.reset();
  }
  genotype_.reset(ind.hasGenotype() ? new MultilocusGenotype(ind.getGenotype()) : 0);
  return *this;
}

/******************************************************************************/

// Id
void Individual::setId(const std::string& id)
{
  id_ = id;
}

/******************************************************************************/

// Sex
void Individual::setSex(const unsigned short sex)
{
  sex_ = sex;
}

/******************************************************************************/

// Date
void Individual::setDate(const Date& date)
{
  date_.reset(new Date(date));
}

/******************************************************************************/

const Date& Individual::date() const
{
  if (hasDate())
    return *date_;
  else
    throw (NullPointerException("Individual::getDate: no date associated to this individual."));
}

/******************************************************************************/

bool Individual::hasDate() const
{
  return date_ != nullptr;
}

/******************************************************************************/

// Coord
void Individual::setCoord(const Point2D<double>& coord)
{
  coord_.reset(new Point2D<double>(coord));
}

/******************************************************************************/

void Individual::setCoord(const double x, const double y)
{
  coord_.reset(new Point2D<double>(x, y));
}

/******************************************************************************/

const Point2D<double>& Individual::coord() const
{
  if (hasCoord())
    return *coord_;
  else
    throw (NullPointerException("Individual::getCoord: no coord associated to this individual."));
}

/******************************************************************************/

bool Individual::hasCoord() const
{
  return coord_ != nullptr;
}

/******************************************************************************/

void Individual::setX(const double x)
{
  if (hasCoord())
    coord_->setX(x);
  else
    throw (NullPointerException("Individual::setX: no coord associated to this individual."));
}

/******************************************************************************/

void Individual::setY(const double y)
{
  if (hasCoord())
    coord_->setY(y);
  else
    throw (NullPointerException("Individual::setY: no coord associated to this individual."));
}

/******************************************************************************/

double Individual::getX() const
{
  if (hasCoord())
    return coord_->getX();
  else
    throw (NullPointerException("Individual::getX: no coord associated to this individual."));
}

/******************************************************************************/

double Individual::getY() const
{
  if (hasCoord())
    return coord_->getY();
  else
    throw (NullPointerException("Individual::getY: no coord associated to this individual."));
}

/******************************************************************************/

shared_ptr<const Locality<double>> Individual::getLocality() const
{
  if (hasLocality())
    return locality_;
  else
    throw (NullPointerException("Individual::getLocality: no locality associated to this individual."));
}

/******************************************************************************/

const Locality<double>& Individual::locality() const
{
  if (hasLocality())
    return *locality_;
  else
    throw (NullPointerException("Individual::locality: no locality associated to this individual."));
}

/******************************************************************************/

// Sequences
void Individual::addSequence(size_t sequenceKey, unique_ptr<Sequence>& sequence)
{
  if (!sequences_)
    sequences_ = make_unique<VectorSequenceContainer>(sequence->getAlphabet());
  try
  {
    sequences_->addSequence(TextTools::toString(sequenceKey), sequence);
  }
  catch (AlphabetMismatchException& ame)
  {
    throw (AlphabetMismatchException("Individual::addSequence: alphabets don't match.", ame.getFirstAlphabet(), ame.getSecondAlphabet()));
  }
  catch (Exception& e)
  {
    if (string(e.what()).find("name") < string(e.what()).size())
      throw (BadIdentifierException("Individual::addSequence: sequence's name already in use.", sequence->getName()));
    // if (string(e.what()).find("key") < string(e.what()).size())
    else
      throw (Exception("Individual::addSequence: sequence_key already in use:" + TextTools::toString(sequenceKey)));
  }
}

/******************************************************************************/

const Sequence& Individual::sequenceByName(const std::string& sequenceName) const
{
  if (!sequences_)
    throw NullPointerException("Individual::getSequenceByName: no sequence data.");
  try
  {
    size_t pos = VectorTools::which(sequences_->getSequenceNames(), sequenceName);
    return sequences_->sequence(pos);
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Individual::getSequenceByName: sequence_name not found.", snfe.getSequenceId());
  }
}

/******************************************************************************/

const Sequence& Individual::sequenceAtPosition(size_t sequencePosition) const
{
  if (!sequences_)
    throw NullPointerException("Individual::getSequenceAtPosition: no sequence data.");
  try
  {
    return sequences_->sequence(TextTools::toString(sequencePosition));
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Individual::getSequenceAtPosition: sequence_position not found", snfe.getSequenceId());
  }
}

/******************************************************************************/

void Individual::deleteSequenceByName(const std::string& sequenceName)
{
  if (!sequences_)
    throw NullPointerException("Individual::deleteSequenceByName: no sequence data.");
  try
  {
    size_t pos = VectorTools::which(sequences_->getSequenceNames(), sequenceName);
    sequences_->deleteSequence(pos);
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Individual::deleteSequenceByName: sequence_name not found.", snfe.getSequenceId());
  }
}

/******************************************************************************/

void Individual::deleteSequenceAtPosition(size_t sequence_position)
{
  if (!sequences_)
    throw NullPointerException("Individual::deleteSequenceAtPosition: no sequence data.");
  try
  {
    sequences_->removeSequence(TextTools::toString(sequence_position));
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Individual::deleteSequenceAtPosition: sequence_position not found.", snfe.getSequenceId());
  }
}

/******************************************************************************/

std::vector<size_t> Individual::getSequencePositions() const
{
  if (!sequences_)
    throw NullPointerException("Individual::getSequencesPositions: no sequence data.");
  vector<size_t> seqpos;
  vector<string> seqkeys = sequences_->getSequenceKeys();
  for (size_t i = 0; i < seqkeys.size(); i++)
  {
    seqpos.push_back((size_t) TextTools::toInt(seqkeys[i]));
  }
  return seqpos;
}

/******************************************************************************/

size_t Individual::getSequencePosition(const std::string& sequenceName) const
{
  if (sequences_.get() == 0)
    throw NullPointerException("Individual::getSequencePosition: no sequence data.");
  try
  {
    return TextTools::to<size_t>(sequences_->sequenceKey(getSequencePosition(sequenceName)));
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Individual::getSequencePosition: sequence_name not found.", snfe.getSequenceId());
  }
}

/******************************************************************************/

bool Individual::hasSequences() const
{
  return !(getNumberOfSequences() == 0);
}

/******************************************************************************/

bool Individual::hasSequenceAtPosition(size_t position) const
{
  if (hasSequences())
  {
    vector<size_t> pos = getSequencePositions();
    for (size_t i = 0; i < pos.size(); ++i)
    {
      if (pos[i] == position)
        return true;
    }
  }
  return false;
}

/******************************************************************************/

shared_ptr<const Alphabet> Individual::getSequenceAlphabet() const
{
  if (!sequences_)
    throw NullPointerException("Individual::getSequenceAlphabet: no sequence data.");
  return sequences_->getAlphabet();
}

/******************************************************************************/

const Alphabet& Individual::sequenceAlphabet() const
{
  if (!sequences_)
    throw NullPointerException("Individual::getSequenceAlphabet: no sequence data.");
  return sequences_->alphabet();
}

/******************************************************************************/

// MultilocusGenotype

void Individual::setGenotype(const MultilocusGenotype& genotype)
{
  genotype_.reset(new MultilocusGenotype(genotype));
}

/******************************************************************************/

void Individual::initGenotype(size_t loci_number)
{
  if (hasGenotype())
    throw Exception("Individual::initGenotype: individual already has a genotype.");
  try
  {
    genotype_.reset(new MultilocusGenotype(loci_number));
  }
  catch (BadIntegerException& bie)
  {
    throw BadIntegerException("Individual::initGenotype: loci_number must be > 0.", bie.getBadInteger());
  }
}

/******************************************************************************/

const MultilocusGenotype& Individual::getGenotype() const
{
  if (!hasGenotype())
    throw NullPointerException("Individual::getGenotype: individual has no genotype.");
  return *genotype_;
}

/******************************************************************************/

void Individual::deleteGenotype()
{
  genotype_.reset();
}

/******************************************************************************/

bool Individual::hasGenotype() const
{
  return genotype_.get() != 0;
}

/******************************************************************************/

void Individual::setMonolocusGenotype(
    size_t locusPosition,
    const MonolocusGenotypeInterface& monogen)
{
  if (!hasGenotype())
    throw NullPointerException("Individual::setMonolocusGenotype: individual has no genotype.");
  try
  {
    genotype_->setMonolocusGenotype(locusPosition, monogen);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Individual::setMonolocusGenotype: locus_position out of boubds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

void Individual::setMonolocusGenotypeByAlleleKey(size_t locus_position, const std::vector<size_t> allele_keys)
{
  if (!hasGenotype())
    throw NullPointerException("Individual::setMonolocusGenotypeByAlleleKey: individual has no genotype.");
  try
  {
    genotype_->setMonolocusGenotypeByAlleleKey(locus_position, allele_keys);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Individual::setMonolocusGenotypeByAlleleKey: locus_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (Exception&)
  {
    throw Exception("Individual::setMonolocusGenotypeByAlleleKey: no key in allele_keys.");
  }
}

/******************************************************************************/

void Individual::setMonolocusGenotypeByAlleleId(size_t locus_position, const std::vector<std::string> allele_id, const LocusInfo& locus_info)
{
  if (!hasGenotype())
    throw NullPointerException("Individual::setMonolocusGenotypeByAlleleId: individual has no genotype.");
  try
  {
    genotype_->setMonolocusGenotypeByAlleleId(locus_position, allele_id, locus_info);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Individual::setMonolocusGenotypeByAlleleId: locus_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (AlleleNotFoundException& anfe)
  {
    throw AlleleNotFoundException("Individual::setMonolocusGenotypeByAlleleId: id not found.", anfe.getIdentifier());
  }
}

/******************************************************************************/

const MonolocusGenotypeInterface& Individual::getMonolocusGenotype(size_t locusPosition)
{
  if (!hasGenotype())
    throw NullPointerException("Individual::getMonolocusGenotype: individual has no genotype.");
  try
  {
    return genotype_->monolocusGenotype(locusPosition);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("Individual::getMonolocusGenotype: locus_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

/******************************************************************************/

size_t Individual::countNonMissingLoci() const
{
  if (!hasGenotype())
    throw NullPointerException("Individual::countNonMissingLoci: individual has no genotype.");
  return genotype_->countNonMissingLoci();
}

/******************************************************************************/

size_t Individual::countHomozygousLoci() const
{
  if (!hasGenotype())
    throw NullPointerException("Individual::countHomozygousLoci: individual has no genotype.");
  return genotype_->countHomozygousLoci();
}

/******************************************************************************/

size_t Individual::countHeterozygousLoci() const
{
  if (!hasGenotype())
    throw NullPointerException("Individual::countHeterozygousLoci: individual has no genotype.");
  return genotype_->countHeterozygousLoci();
}

/******************************************************************************/
