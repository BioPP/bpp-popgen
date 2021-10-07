//
// File MultilocusGenotype.cpp
// Author : Sylvain Gaillard <sylvain.gaillard@angers.inra.fr>
// Last modification : April 4, 2008
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

#include "MultilocusGenotype.h"

using namespace bpp;
using namespace std;

// ** Class constructor: *******************************************************/

MultilocusGenotype::MultilocusGenotype(size_t loci_number) : loci_(vector<MonolocusGenotype*>(loci_number))
{
  if (loci_number < 1)
    throw BadIntegerException("MultilocusGenotype::MultilocusGenotype: loci_number must be > 0.", static_cast<int>(loci_number));

  // Set all the loci_ pointers to nullptr
  for (size_t i = 0; i < loci_number; i++)
  {
    loci_[i] = 0;
  }
}

MultilocusGenotype::MultilocusGenotype(const MultilocusGenotype& genotype) : loci_(vector<MonolocusGenotype*>(genotype.size()))
{
  for (size_t i = 0; i < genotype.size(); i++)
  {
    if (!genotype.isMonolocusGenotypeMissing(i))
      loci_[i] = dynamic_cast<MonolocusGenotype*>(genotype.getMonolocusGenotype(i).clone());
    else
      loci_[i] = 0;
  }
}

// ** Class destructor: *******************************************************/

MultilocusGenotype::~MultilocusGenotype()
{
  for (size_t i = 0; i < loci_.size(); i++)
  {
    delete loci_[i];
  }
  loci_.clear();
}

// ** Other methodes: *********************************************************/

void MultilocusGenotype::setMonolocusGenotype(size_t locus_position,
                                              const MonolocusGenotype& monogen)
{
  if (locus_position < loci_.size())
    loci_[locus_position] = dynamic_cast<MonolocusGenotype*>(monogen.clone());
  else
    throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotype: locus_position out of bounds.",
                                    locus_position, 0, loci_.size());
}

void MultilocusGenotype::setMonolocusGenotypeByAlleleKey(size_t locus_position,
                                                         const std::vector<size_t>& allele_keys)
{
  if (allele_keys.size() < 1)
    throw Exception("MultilocusGenotype::setMonolocusGenotypeByAlleleKey: no key in allele_keys.");

  if (locus_position < loci_.size())
  {
    setMonolocusGenotype(locus_position, *MonolocusGenotypeTools::buildMonolocusGenotypeByAlleleKey(allele_keys));
  }
  else
    throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotype: locus_position out of bounds.",
                                    locus_position, 0, loci_.size());
}

void MultilocusGenotype::setMonolocusGenotypeByAlleleId(size_t locus_position,
                                                        const std::vector<std::string>& allele_id, const LocusInfo& locus_info)
{
  vector<size_t> allele_keys;
  for (size_t i = 0; i < allele_id.size(); i++)
  {
    try
    {
      allele_keys.push_back(locus_info.getAlleleInfoKey(allele_id[i]));
    }
    catch (AlleleNotFoundException& anfe)
    {
      throw AlleleNotFoundException("MultilocusGenotype::setMonolocusGenotypeByAlleleId: id not found.", anfe.getIdentifier());
    }
  }
  try
  {
    setMonolocusGenotypeByAlleleKey(locus_position, allele_keys);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotypeByAlleleId: locus_position out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
}

void MultilocusGenotype::setMonolocusGenotypeAsMissing(size_t locus_position)
{
  if (locus_position >= loci_.size())
    throw IndexOutOfBoundsException("MultilocusGenotype::setMonolocusGenotypeAsMissing: locus_position out of bounds.", locus_position, 0, loci_.size());
  if (loci_[locus_position] != NULL)
    delete loci_[locus_position];
  loci_[locus_position] = NULL;
}

bool MultilocusGenotype::isMonolocusGenotypeMissing(size_t locus_position) const
{
  if (locus_position >= loci_.size())
    throw IndexOutOfBoundsException("MultilocusGenotype::isMonolocusGenotypeMissing: locus_position out of bounds.", locus_position, 0, loci_.size());
  return loci_[locus_position] == NULL;
}

const MonolocusGenotype& MultilocusGenotype::getMonolocusGenotype(size_t locus_position) const
{
  if (locus_position >= loci_.size())
    throw IndexOutOfBoundsException("MultilocusGenotype::getMonolocusGenotype: locus_position out of bounds", locus_position, 0, loci_.size());
  return *loci_[locus_position];
}

size_t MultilocusGenotype::size() const
{
  return loci_.size();
}

size_t MultilocusGenotype::countNonMissingLoci() const
{
  size_t count = 0;
  for (size_t i = 0; i < loci_.size(); i++)
  {
    if (loci_[i] != NULL)
      count++;
  }
  return count;
}

size_t MultilocusGenotype::countHomozygousLoci() const
{
  size_t count = 0;
  for (size_t i = 0; i < loci_.size(); i++)
  {
    try
    {
      if (dynamic_cast<BiAlleleMonolocusGenotype*>(loci_[i])->isHomozygous())
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
      if (!(dynamic_cast<BiAlleleMonolocusGenotype*>(loci_[i])->isHomozygous()))
        count++;
    }
    catch (...)
    {}
  }
  return count;
}
