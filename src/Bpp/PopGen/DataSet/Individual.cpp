//
// File Individual.cpp
// Author : Sylvain Gaillard
// Last modification : Tuesday August 03 2004
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

#include "Individual.h"

using namespace bpp;

using namespace std;

// ** Class constructor: *******************************************************/
Individual::Individual() : id_(""),
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
                       Locality<double>* locality,
                       const unsigned short sex) :
  id_(id),
  sex_(sex),
  date_(new Date(date)),
  coord_(new Point2D<double>(coord)),
  locality_(locality),
  sequences_(),
  genotype_() {}

Individual::Individual(const Individual& ind) : id_(ind.getId()),
  sex_(ind.getSex()),
  date_(),
  coord_(),
  locality_(0),
  sequences_(),
  genotype_()
{
  try
  {
    setDate(ind.getDate());
  }
  catch (...)
  {}
  try
  {
    setCoord(ind.getCoord());
  }
  catch (...)
  {}
  try
  {
    setLocality(ind.getLocality());
  }
  catch (...)
  {}
  try
  {
    setSequences(dynamic_cast<const MapSequenceContainer&>(ind.getSequences()));
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
    setDate(ind.getDate());
  }
  catch (NullPointerException)
  {
    date_.reset();
  }
  try
  {
    setCoord(ind.getCoord());
  }
  catch (NullPointerException)
  {
    coord_.reset();
  }
  try
  {
    setLocality(ind.getLocality());
  }
  catch (NullPointerException)
  {
    locality_ = 0;
  }
  try
  {
    setSequences(dynamic_cast<const MapSequenceContainer&>(ind.getSequences()));
  }
  catch (NullPointerException)
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

const Date& Individual::getDate() const
{
  if (hasDate())
    return *date_.get();
  else
    throw (NullPointerException("Individual::getDate: no date associated to this individual."));
}

/******************************************************************************/

bool Individual::hasDate() const
{
  return date_.get() != 0;
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

const Point2D<double>& Individual::getCoord() const
{
  if (hasCoord())
    return *coord_.get();
  else
    throw (NullPointerException("Individual::getCoord: no coord associated to this individual."));
}

/******************************************************************************/

bool Individual::hasCoord() const
{
  return coord_.get() != 0;
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

// Locality
void Individual::setLocality(const Locality<double>* locality)
{
  locality_ = locality;
}

/******************************************************************************/

const Locality<double>* Individual::getLocality() const
{
  if (hasLocality())
    return locality_;
  else
    throw (NullPointerException("Individual::getLocality: no locality associated to this individual."));
}

/******************************************************************************/

bool Individual::hasLocality() const
{
  return locality_ != 0;
}

/******************************************************************************/

// Sequences
void Individual::addSequence(size_t sequence_key, const Sequence& sequence)
{
  if (sequences_.get() == 0)
    sequences_.reset(new MapSequenceContainer(sequence.getAlphabet()));
  try
  {
    sequences_->addSequence(TextTools::toString(sequence_key), sequence);
  }
  catch (AlphabetMismatchException& ame)
  {
    throw (AlphabetMismatchException("Individual::addSequence: alphabets don't match.", ame.getAlphabets()[0], ame.getAlphabets()[1]));
  }
  catch (Exception& e)
  {
    if (string(e.what()).find("name") < string(e.what()).size())
      throw (BadIdentifierException("Individual::addSequence: sequence's name already in use.", sequence.getName()));
    // if (string(e.what()).find("key") < string(e.what()).size())
    else
      throw (Exception("Individual::addSequence: sequence_key already in use:" + TextTools::toString(sequence_key)));
  }
}

/******************************************************************************/

const Sequence& Individual::getSequenceByName(const std::string& sequence_name) const
{
  if (sequences_.get() == 0)
    throw NullPointerException("Individual::getSequenceByName: no sequence data.");
  try
  {
    return sequences_->getSequence(sequence_name);
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Individual::getSequenceByName: sequence_name not found.", snfe.getSequenceId());
  }
}

/******************************************************************************/

const Sequence& Individual::getSequenceAtPosition(size_t sequence_position) const
{
  if (sequences_.get() == 0)
    throw NullPointerException("Individual::getSequenceAtPosition: no sequence data.");
  try
  {
    return sequences_->getSequenceByKey(TextTools::toString(sequence_position));
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Individual::getSequenceAtPosition: sequence_position not found", snfe.getSequenceId());
  }
}

/******************************************************************************/

void Individual::deleteSequenceByName(const std::string& sequence_name)
{
  if (sequences_.get() == 0)
    throw NullPointerException("Individual::deleteSequenceByName: no sequence data.");
  try
  {
    sequences_->deleteSequence(sequence_name);
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Individual::deleteSequenceByName: sequence_name not found.", snfe.getSequenceId());
  }
}

/******************************************************************************/

void Individual::deleteSequenceAtPosition(size_t sequence_position)
{
  if (sequences_.get() == 0)
    throw NullPointerException("Individual::deleteSequenceAtPosition: no sequence data.");
  try
  {
    sequences_->deleteSequenceByKey(TextTools::toString(sequence_position));
  }
  catch (SequenceNotFoundException& snfe)
  {
    throw SequenceNotFoundException("Individual::deleteSequenceAtPosition: sequence_position not found.", snfe.getSequenceId());
  }
}

/******************************************************************************/

std::vector<std::string> Individual::getSequencesNames() const
{
  if (sequences_.get() == 0)
    throw NullPointerException("Individual::getSequencesNames: no sequence data.");
  return sequences_->getSequencesNames();
}

/******************************************************************************/

std::vector<size_t> Individual::getSequencesPositions() const
{
  if (sequences_.get() == 0)
    throw NullPointerException("Individual::getSequencesPositions: no sequence data.");
  vector<size_t> seqpos;
  vector<string> seqkeys = sequences_->getKeys();
  for (size_t i = 0; i < seqkeys.size(); i++)
  {
    seqpos.push_back((size_t) TextTools::toInt(seqkeys[i]));
  }
  return seqpos;
}

/******************************************************************************/

size_t Individual::getSequencePosition(const std::string& sequence_name) const
{
  if (sequences_.get() == 0)
    throw NullPointerException("Individual::getSequencePosition: no sequence data.");
  try
  {
    return (size_t) TextTools::toInt(sequences_->getKey(getSequencePosition(sequence_name)));
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
    vector<size_t> pos = getSequencesPositions();
    for (size_t i = 0; i < pos.size(); i++)
    {
      if (pos[i] == position)
        return true;
    }
  }
  return false;
}

/******************************************************************************/

const Alphabet* Individual::getSequenceAlphabet() const
{
  if (sequences_.get() == 0)
    throw NullPointerException("Individual::getSequenceAlphabet: no sequence data.");
  return sequences_->getAlphabet();
}

/******************************************************************************/

size_t Individual::getNumberOfSequences() const
{
  if (sequences_.get() == 0)
    return 0;
  return sequences_->getNumberOfSequences();
}

/******************************************************************************/

void Individual::setSequences(const MapSequenceContainer& msc)
{
  sequences_.reset(new MapSequenceContainer(msc));
}

/******************************************************************************/

const OrderedSequenceContainer& Individual::getSequences() const
{
  if (sequences_.get() == 0)
    throw NullPointerException("Individual::getSequences: no sequence data.");
  return *sequences_;
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

void Individual::setMonolocusGenotype(size_t locus_position, const MonolocusGenotype& monogen)
{
  if (!hasGenotype())
    throw NullPointerException("Individual::setMonolocusGenotype: individual has no genotype.");
  try
  {
    genotype_->setMonolocusGenotype(locus_position, monogen);
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
  catch (Exception)
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

const MonolocusGenotype& Individual::getMonolocusGenotype(size_t locus_position)
{
  if (!hasGenotype())
    throw NullPointerException("Individual::getMonolocusGenotype: individual has no genotype.");
  try
  {
    return genotype_->getMonolocusGenotype(locus_position);
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

