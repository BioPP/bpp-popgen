// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "MultilocusGenotype.h"

using namespace bpp;
using namespace std;

// ** Class constructor: *******************************************************/

MultilocusGenotype::MultilocusGenotype(size_t loci_number) :
  loci_(loci_number)
{
  if (loci_number < 1)
    throw BadIntegerException("MultilocusGenotype::MultilocusGenotype: loci_number must be > 0.", static_cast<int>(loci_number));

  // Set all the loci_ pointers to nullptr
  for (size_t i = 0; i < loci_number; ++i)
  {
    loci_[i] = nullptr;
  }
}

MultilocusGenotype::MultilocusGenotype(const MultilocusGenotype& genotype) :
  loci_(genotype.size())
{
  for (size_t i = 0; i < genotype.size(); ++i)
  {
    if (!genotype.isMonolocusGenotypeMissing(i))
      loci_[i].reset(genotype.monolocusGenotype(i).clone());
    else
      loci_[i] = nullptr;
  }
}

// ** Class destructor: *******************************************************/

MultilocusGenotype::~MultilocusGenotype()
{
  loci_.clear();
}

// ** Other methodes: *********************************************************/

void MultilocusGenotype::setMonolocusGenotype(
    size_t locusPosition,
    const MonolocusGenotypeInterface& monogen)
{
  if (locusPosition < loci_.size())
    loci_[locusPosition].reset(monogen.clone());
  else
    throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotype: locusPosition out of bounds.", locusPosition, 0, loci_.size());
}

void MultilocusGenotype::setMonolocusGenotypeByAlleleKey(
    size_t locusPosition,
    const std::vector<size_t>& alleleKeys)
{
  if (alleleKeys.size() < 1)
    throw Exception("MultilocusGenotype::setMonolocusGenotypeByAlleleKey: no key in alleleKeys.");

  if (locusPosition < loci_.size())
  {
    setMonolocusGenotype(locusPosition, *MonolocusGenotypeTools::buildMonolocusGenotypeByAlleleKey(alleleKeys));
  }
  else
    throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotype: locusPosition out of bounds.", locusPosition, 0, loci_.size());
}

void MultilocusGenotype::setMonolocusGenotypeByAlleleId(
    size_t locusPosition,
    const vector<string>& alleleId,
    const LocusInfo& locusInfo)
{
  vector<size_t> alleleKeys;
  for (size_t i = 0; i < alleleId.size(); i++)
  {
    try
    {
      alleleKeys.push_back(locusInfo.getAlleleInfoKey(alleleId[i]));
    }
    catch (AlleleNotFoundException& anfe)
    {
      throw AlleleNotFoundException("MultilocusGenotype::setMonolocusGenotypeByAlleleId: id not found.", anfe.getIdentifier());
    }
  }
  try
  {
    setMonolocusGenotypeByAlleleKey(locusPosition, alleleKeys);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotypeByAlleleId: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

void MultilocusGenotype::setMonolocusGenotypeAsMissing(size_t locusPosition)
{
  if (locusPosition >= loci_.size())
    throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotypeAsMissing: locusPosition out of bounds.", locusPosition, 0, loci_.size());
  loci_[locusPosition] = nullptr;
}

bool MultilocusGenotype::isMonolocusGenotypeMissing(size_t locusPosition) const
{
  if (locusPosition >= loci_.size())
    throw IndexOutOfBoundsException("MultilocusGenotype::isMonolocusGenotypeMissing: locusPosition out of bounds.", locusPosition, 0, loci_.size());
  return loci_[locusPosition] == nullptr;
}

const MonolocusGenotypeInterface& MultilocusGenotype::monolocusGenotype(size_t locusPosition) const
{
  if (locusPosition >= loci_.size())
    throw IndexOutOfBoundsException("MultilocusGenotype::getMonolocusGenotype: locusPosition out of bounds", locusPosition, 0, loci_.size());
  return *loci_[locusPosition];
}

size_t MultilocusGenotype::size() const
{
  return loci_.size();
}

size_t MultilocusGenotype::countNonMissingLoci() const
{
  size_t count = 0;
  for (size_t i = 0; i < loci_.size(); ++i)
  {
    if (loci_[i])
      count++;
  }
  return count;
}

size_t MultilocusGenotype::countHomozygousLoci() const
{
  size_t count = 0;
  for (size_t i = 0; i < loci_.size(); ++i)
  {
    try
    {
      if (dynamic_cast<BiAlleleMonolocusGenotype&>(*loci_[i]).isHomozygous())
        count++;
    }
    catch (...)
    {}
  }
  return count;
}

size_t MultilocusGenotype::countHeterozygousLoci() const
{
  size_t count = 0;
  for (size_t i = 0; i < loci_.size(); i++)
  {
    try
    {
      if (!(dynamic_cast<BiAlleleMonolocusGenotype&>(*loci_[i]).isHomozygous()))
        count++;
    }
    catch (...)
    {}
  }
  return count;
}
