//
// File AnalyzedLoci.cpp
// Author : Sylvain Gaillard
// Last modification : Thursday July 29 2004
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

#include "AnalyzedLoci.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

AnalyzedLoci::AnalyzedLoci(size_t number_of_loci) : loci_(vector<LocusInfo*>(number_of_loci))
{
  for (size_t i = 0; i < loci_.size(); i++)
  {
    loci_[i] = 0;
  }
}

/******************************************************************************/

AnalyzedLoci::AnalyzedLoci(const AnalyzedLoci& analyzed_loci) : loci_(vector<LocusInfo*>(analyzed_loci.loci_.size()))
{
  for (size_t i = 0; i < analyzed_loci.getNumberOfLoci(); i++)
  {
    loci_[i] = new LocusInfo(analyzed_loci.getLocusInfoAtPosition(i));
  }
}

/******************************************************************************/

AnalyzedLoci::~AnalyzedLoci()
{
  for (size_t i = 0; i < loci_.size(); i++)
  {
    delete loci_[i];
  }
}

/******************************************************************************/

void AnalyzedLoci::setLocusInfo(
  size_t locus_position,
  const LocusInfo& locus)
throw (IndexOutOfBoundsException)
{
  if (locus_position < loci_.size())
    loci_[locus_position] = new LocusInfo(locus);
  else
    throw IndexOutOfBoundsException("AnalyzedLoci::setLocusInfo: locus_position out of bounds",
                                    locus_position, 0, loci_.size());
}

/******************************************************************************/

size_t AnalyzedLoci::getLocusInfoPosition(
  const std::string& locus_name) const
throw (BadIdentifierException)
{
  for (size_t i = 0; i < loci_.size(); i++)
  {
    if (loci_[i] != NULL && loci_[i]->getName() == locus_name)
      return i;
  }
  throw BadIdentifierException("AnalyzedLoci::getLocusInfoPosition: locus not found.", locus_name);
}

/******************************************************************************/

const LocusInfo& AnalyzedLoci::getLocusInfoByName(
  const std::string& locus_name) const
throw (BadIdentifierException)
{
  for (size_t i = 0; i < loci_.size(); i++)
  {
    if (loci_[i] != NULL && loci_[i]->getName() == locus_name)
      return *(loci_[i]);
  }
  throw BadIdentifierException("AnalyzedLoci::getLocusInfo: locus not found.",
                               locus_name);
}

/******************************************************************************/

const LocusInfo& AnalyzedLoci::getLocusInfoAtPosition(
  size_t locus_position) const
throw (Exception)
{
  if (locus_position >= loci_.size())
    throw IndexOutOfBoundsException("AnalyzedLoci::getLocusInfoAtPosition: locus_position out of bounds.", locus_position, 0, loci_.size());
  if (loci_[locus_position] != NULL)
    return *(loci_[locus_position]);
  else
    throw NullPointerException("AnalyzedLoci::getLocusInfo: no locus defined here.");
}

/******************************************************************************/

// AlleleInfo
void AnalyzedLoci::addAlleleInfoByLocusName(const std::string& locus_name,
                                            const AlleleInfo& allele)
throw (Exception)
{
  bool locus_found = false;
  for (vector<LocusInfo*>::iterator it = loci_.begin(); it != loci_.end(); it++)
  {
    if ((*it)->getName() == locus_name)
    {
      locus_found = true;
      try
      {
        (*it)->addAlleleInfo(allele);
      }
      catch (BadIdentifierException& bie)
      {
        throw BadIdentifierException("AnalyzedLoci::addAlleleInfoByLocusName: allele id already in use.", bie.getIdentifier());
      }
    }
  }
  if (!locus_found)
    throw LocusNotFoundException("AnalyzedLoci::addAlleleInfoByLocusName: locus_name not found.",
                                 locus_name);
}

/******************************************************************************/

void AnalyzedLoci::addAlleleInfoByLocusPosition(size_t locus_position,
                                                const AlleleInfo& allele)
throw (Exception)
{
  if (locus_position < loci_.size())
  {
    try
    {
      loci_[locus_position]->addAlleleInfo(allele);
    }
    catch (BadIdentifierException& bie)
    {
      throw BadIdentifierException("AnalyzedLoci::addAlleleInfoByLocusPosition: allele id is already in use.", bie.getIdentifier());
    }
  }
  else
    throw IndexOutOfBoundsException("AnalyzedLoci::addAlleleInfoByLocusPosition: locus_position out of bounds.",
                                    locus_position, 0, loci_.size());
}

/******************************************************************************/

std::vector<size_t> AnalyzedLoci::getNumberOfAlleles() const
{
  vector<size_t> allele_count;
  for (size_t i = 0; i < loci_.size(); i++)
  {
    allele_count.push_back(loci_[i]->getNumberOfAlleles());
  }
  return allele_count;
}

/******************************************************************************/

unsigned int AnalyzedLoci::getPloidyByLocusName(const std::string& locus_name) const
throw (LocusNotFoundException)
{
  for (size_t i = 0; i < loci_.size(); i++)
  {
    if (loci_[i] != NULL && loci_[i]->getName() == locus_name)
      return loci_[i]->getPloidy();
  }
  throw LocusNotFoundException("AnalyzedLoci::getLocusInfo: locus_name not found.",
                               locus_name);
}

/******************************************************************************/

unsigned int AnalyzedLoci::getPloidyByLocusPosition(size_t locus_position) const
throw (IndexOutOfBoundsException)
{
  if (locus_position >= loci_.size())
    throw IndexOutOfBoundsException("AnalyzedLoci::getPloidyByLocusPosition: locus_position out of bounds.", locus_position, 0, loci_.size());
  return loci_[locus_position]->getPloidy();
}

/******************************************************************************/

