//
// File MultiAlleleMonolocusGenotype.cpp
// Author : Sylvain Gaillard <sylvain.gaillard@angers.inra.fr>
// Last modification : Wednesday March 5 2008
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

#include "MultiAlleleMonolocusGenotype.h"

using namespace bpp;

//** Class constructor: *******************************************************/

MultiAlleleMonolocusGenotype::MultiAlleleMonolocusGenotype(vector<unsigned int> allele_index)
{
  for (unsigned int i=0 ; i<allele_index.size() ; ++i)
    _allele_index.push_back(allele_index[i]);
}

MultiAlleleMonolocusGenotype::MultiAlleleMonolocusGenotype(const MultiAlleleMonolocusGenotype & mmg)
{
  for (unsigned int i=0 ; i<mmg.getAlleleIndex().size() ; ++i)
    _allele_index.push_back(mmg.getAlleleIndex()[i]);
}

//** Class destructor: ********************************************************/

MultiAlleleMonolocusGenotype::~MultiAlleleMonolocusGenotype()
{
  _allele_index.clear();
}

//** Other methodes: **********************************************************/

MultiAlleleMonolocusGenotype & MultiAlleleMonolocusGenotype::operator= (const MultiAlleleMonolocusGenotype & mmg)
{
  for(unsigned int i=0 ; i<mmg.getAlleleIndex().size() ; ++i)
    _allele_index.push_back(mmg.getAlleleIndex()[i]);
  return * this;
}

bool MultiAlleleMonolocusGenotype::operator== (const MultiAlleleMonolocusGenotype & mmg) const
{
  return ((_allele_index[0] == mmg.getAlleleIndex()[0] && _allele_index[1] == mmg.getAlleleIndex()[1])
      || (_allele_index[0] == mmg.getAlleleIndex()[1] && _allele_index[1] == mmg.getAlleleIndex()[0]));
}

bool MultiAlleleMonolocusGenotype::isHomozygous() const
{
  for (unsigned int i=1 ; i<_allele_index.size() ; ++i)
    if (_allele_index[i-1] != _allele_index[i])
      return false;
  return true;
}

vector<unsigned int> MultiAlleleMonolocusGenotype::getAlleleIndex() const
{
  return _allele_index;
}

Clonable * MultiAlleleMonolocusGenotype::clone() const
{
  return new MultiAlleleMonolocusGenotype(* this);
}

