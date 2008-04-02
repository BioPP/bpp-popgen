//
// File LocusInfo.cpp
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

// From Utils
#include <Utils/TextTools.h>

#include "LocusInfo.h"
#include "GeneralExceptions.h"

using namespace bpp;

unsigned int LocusInfo::HAPLODIPLOID = 0;
unsigned int LocusInfo::HAPLOID = 1;
unsigned int LocusInfo::DIPLOID = 2;
unsigned int LocusInfo::UNKNOWN = 9999;

//** Class constructor: *******************************************************/

LocusInfo::LocusInfo(const string &name, const unsigned int ploidy)
{
  _name = name;
  _ploidy = ploidy;
}

LocusInfo::LocusInfo(const LocusInfo & locus_info)
{
  _name = locus_info.getName();
  _ploidy = locus_info.getPloidy();
  for (unsigned int i = 0 ; i < locus_info.getNumberOfAlleles() ; i++) {
    Clonable * tmp_allele = locus_info.getAlleleInfoByKey(i)->clone();
    _alleles.push_back(dynamic_cast<AlleleInfo *>(tmp_allele));
  }
}

//** Class destructor: *******************************************************/

LocusInfo::~LocusInfo()
{
  for(unsigned int i = 0; i < _alleles.size(); i++)
    delete _alleles[i];
  _alleles.clear();
}

//** Other methodes: *********************************************************/
// Name
string LocusInfo::getName() const
{
  return _name;
}

// Ploidie
unsigned int LocusInfo::getPloidy() const
{
  return _ploidy;
}

// AlleleInfos
void LocusInfo::addAlleleInfo(const AlleleInfo &allele) throw (BadIdentifierException)
{
  // Check if the allele id is not already in use
  for (unsigned int i = 0 ; i < _alleles.size() ; i++)
    if (_alleles[i]->getId() == allele.getId())
      throw BadIdentifierException("LocusInfo::addAlleleInfo: Id already in use.",allele.getId());
  _alleles.push_back(dynamic_cast<AlleleInfo *>(allele.clone()));
}

  AlleleInfo * LocusInfo::getAlleleInfoById(const string & id) const
throw (AlleleNotFoundException)
{
  for (unsigned int i = 0; i < _alleles.size(); i++)
    if (_alleles[i]->getId() == id)
      return _alleles[i];
  throw AlleleNotFoundException("LocusInfo::getAlleleInfoById: AlleleInfo id unknown.", id);
}

AlleleInfo * LocusInfo::getAlleleInfoByKey(unsigned int key) const throw (IndexOutOfBoundsException)
{
  if (key >= _alleles.size())
    throw IndexOutOfBoundsException("LocusInfo::getAlleleInfoByKey: key out of bounds.", key, 0, _alleles.size());
  return _alleles[key];
}

  unsigned int LocusInfo::getAlleleInfoKey(const string & id) const
throw (AlleleNotFoundException)
{
  for(unsigned int i = 0; i < _alleles.size(); i++)
    if (_alleles[i]->getId() == id)
      return i;
  throw AlleleNotFoundException("LocusInfo::getAlleleInfoKey: AlleleInfo id not found.", id);
}

unsigned int LocusInfo::getNumberOfAlleles() const
{
  return _alleles.size();
}

void LocusInfo::clear()
{
  for(unsigned int i = 0; i < _alleles.size(); i++)
    delete _alleles[i];
  _alleles.clear();
}

