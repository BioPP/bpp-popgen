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

//** Constructors: ***********************************************************/

AnalyzedLoci::AnalyzedLoci(unsigned int number_of_loci)
{
  _loci.resize(number_of_loci);
  for(unsigned int i = 0 ; i < _loci.size() ; i++)
    _loci[i] = NULL;
}

AnalyzedLoci::AnalyzedLoci(const AnalyzedLoci & analyzed_loci)
{
  for(unsigned int i = 0; i < analyzed_loci.getNumberOfLoci(); i++)
    _loci.push_back(new LocusInfo(* analyzed_loci.getLocusInfoAtPosition(i)));
}

//** Destructor: *************************************************************/

AnalyzedLoci::~AnalyzedLoci()
{
  for(unsigned int i = 0; i < _loci.size(); i++)
    delete _loci[i];
}

//** Other methodes: *********************************************************/
// LocusInfo
  void AnalyzedLoci::setLocusInfo(unsigned int locus_position, const LocusInfo & locus)
throw (IndexOutOfBoundsException)
{
  if (locus_position >= 0 && locus_position < _loci.size())
    _loci[locus_position] = new LocusInfo(locus);
  else
    throw IndexOutOfBoundsException("AnalyzedLoci::setLocusInfo: locus_position out of bounds",
        locus_position, 0, _loci.size());
}

  unsigned int AnalyzedLoci::getLocusInfoPosition(const string & locus_name) const
throw (BadIdentifierException)
{
  for(unsigned int i = 0; i < _loci.size(); i++)
    if (_loci[i] != NULL && _loci[i]->getName() == locus_name)
      return i;
  throw BadIdentifierException("AnalyzedLoci::getLocusInfoPosition: locus not found.", locus_name);
}

  const LocusInfo * AnalyzedLoci::getLocusInfoByName(const string & locus_name) const
throw (BadIdentifierException)
{
  for(unsigned int i = 0; i < _loci.size(); i++)
    if (_loci[i] != NULL && _loci[i]->getName() == locus_name)
      return _loci[i];
  throw BadIdentifierException("AnalyzedLoci::getLocusInfo: locus not found.",
      locus_name);
}

  const LocusInfo * AnalyzedLoci::getLocusInfoAtPosition(unsigned int locus_position) const
throw (Exception)
{
  if(locus_position >= _loci.size())
    throw IndexOutOfBoundsException("AnalyzedLoci::getLocusInfoAtPosition: locus_position out of bounds.", locus_position, 0, _loci.size());
  if(_loci[locus_position] != NULL)
    return _loci[locus_position];
  else
    throw NullPointerException("AnalyzedLoci::getLocusInfo: no locus defined here.");
}

// AlleleInfo
void AnalyzedLoci::addAlleleInfoByLocusName(const string & locus_name,
    const AlleleInfo& allele)
throw (Exception)
{
  bool locus_found = false;
  for (vector<LocusInfo *>::iterator it = _loci.begin(); it != _loci.end() ; it++)
  {
    if((*it)->getName() == locus_name)
    {
      locus_found = true;
      try
      {
        (*it)->addAlleleInfo(allele);
      }
      catch (BadIdentifierException & bie)
      {
        throw BadIdentifierException("AnalyzedLoci::addAlleleInfoByLocusName: allele id already in use.", bie.getIdentifier());
      }
    }
  }
  if(!locus_found)
    throw LocusNotFoundException("AnalyzedLoci::addAlleleInfoByLocusName: locus_name not found.",
        locus_name);
}

void AnalyzedLoci::addAlleleInfoByLocusPosition(unsigned int locus_position,
    const AlleleInfo & allele)
throw (Exception)
{
  if(locus_position >= 0 && locus_position < _loci.size())
  {
    try
    {
      _loci[locus_position]->addAlleleInfo(allele);
    }
    catch (BadIdentifierException & bie)
    {
      throw BadIdentifierException("AnalyzedLoci::addAlleleInfoByLocusPosition: allele id is already in use.", bie.getIdentifier());
    }
  }
  else
    throw IndexOutOfBoundsException("AnalyzedLoci::addAlleleInfoByLocusPosition: locus_position out of bounds.",
        locus_position, 0, _loci.size());
}

// General
unsigned int AnalyzedLoci::getNumberOfLoci() const {
  return _loci.size();
}

vector<unsigned int> AnalyzedLoci::getNumberOfAlleles() const
{
  vector<unsigned int> allele_count;
  for (unsigned int i = 0; i < _loci.size(); i++)
    allele_count.push_back(_loci[i]->getNumberOfAlleles());
  return allele_count;
}

  unsigned int AnalyzedLoci::getPloidyByLocusName(const string & locus_name) const
throw (LocusNotFoundException)
{
  for (unsigned int i = 0; i < _loci.size(); i++)
    if(_loci[i] != NULL && _loci[i]->getName() == locus_name)
      return _loci[i]->getPloidy();
  throw LocusNotFoundException("AnalyzedLoci::getLocusInfo: locus_name not found.",
      locus_name);
}

  unsigned int AnalyzedLoci::getPloidyByLocusPosition(unsigned int locus_position) const
throw (IndexOutOfBoundsException)
{
  if(locus_position >= _loci.size())
    throw IndexOutOfBoundsException("AnalyzedLoci::getPloidyByLocusPosition: locus_position out of bounds.", locus_position, 0, _loci.size());
  return _loci[locus_position]->getPloidy();
}

