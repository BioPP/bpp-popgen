//
// File MultiAlleleMonolocusGenotype.cpp
// Author : Sylvain Gaillard <sylvain.gaillard@angers.inra.fr>
// Last modification : Wednesday March 5 2008
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

#include "MultiAlleleMonolocusGenotype.h"

using namespace bpp;
using namespace std;

// ** Class constructor: *******************************************************/

MultiAlleleMonolocusGenotype::MultiAlleleMonolocusGenotype(std::vector<size_t> allele_index) : alleleIndex_(vector<size_t>(allele_index.size()))
{
  for (size_t i = 0; i < allele_index.size(); ++i)
  {
    alleleIndex_[i] = allele_index[i];
  }
}

MultiAlleleMonolocusGenotype::MultiAlleleMonolocusGenotype(const MultiAlleleMonolocusGenotype& mmg) : alleleIndex_(vector<size_t>(mmg.alleleIndex_.size()))
{
  for (size_t i = 0; i < mmg.getAlleleIndex().size(); ++i)
  {
    alleleIndex_[i] = mmg.getAlleleIndex()[i];
  }
}

// ** Class destructor: ********************************************************/

MultiAlleleMonolocusGenotype::~MultiAlleleMonolocusGenotype()
{
  alleleIndex_.clear();
}

// ** Other methodes: **********************************************************/

MultiAlleleMonolocusGenotype& MultiAlleleMonolocusGenotype::operator=(const MultiAlleleMonolocusGenotype& mmg)
{
  for (size_t i = 0; i < mmg.getAlleleIndex().size(); ++i)
  {
    alleleIndex_.push_back(mmg.getAlleleIndex()[i]);
  }
  return *this;
}

bool MultiAlleleMonolocusGenotype::operator==(const MultiAlleleMonolocusGenotype& mmg) const
{
  return (alleleIndex_[0] == mmg.getAlleleIndex()[0] && alleleIndex_[1] == mmg.getAlleleIndex()[1])
         || (alleleIndex_[0] == mmg.getAlleleIndex()[1] && alleleIndex_[1] == mmg.getAlleleIndex()[0]);
}

bool MultiAlleleMonolocusGenotype::isHomozygous() const
{
  for (size_t i = 1; i < alleleIndex_.size(); ++i)
  {
    if (alleleIndex_[i - 1] != alleleIndex_[i])
      return false;
  }
  return true;
}

