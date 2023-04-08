//
// File DataSetTools.cpp
// Author : Sylvain Gaillard
// Last modification : Wednesday August 04 2004
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

#include "DataSetTools.h"

using namespace bpp;
using namespace std;

std::unique_ptr<DataSet> DataSetTools::buildDataSet(const SequenceContainerInterface& sc)
{
  auto dataset = make_unique<DataSet>();
  dataset->addEmptyGroup(0);
  for (size_t i = 0; i < sc.getNumberOfSequences(); ++i)
  {
    dataset->addEmptyIndividualToGroup(0, string("Individual_") + TextTools::toString(i + 1));
    try
    {
      auto tmpSeq = unique_ptr<Sequence>(sc.sequence(i).clone());
      dataset->addIndividualSequenceInGroup(0, i, 0, tmpSeq);
    }
    catch (Exception& e)
    {
      throw e;
    }
  }
  return dataset;
}

std::unique_ptr<DataSet> DataSetTools::buildDataSet(const PolymorphismSequenceContainer& psc)
{
  auto dataset = make_unique<DataSet>();
  set<size_t> grpIds = psc.getAllGroupsIds();
  for (auto& it : grpIds)
  {
    dataset->addEmptyGroup(it);
  }
  size_t indCount = 0;
  for (size_t i = 0; i < psc.getNumberOfSequences(); ++i)
  {
    for (size_t j = 0; j < psc.getSequenceCount(i); ++j)
    {
      dataset->addEmptyIndividualToGroup(psc.getGroupId(i), string("Individual_") + TextTools::toString(indCount++));
      try
      {
	auto tmpSeq = unique_ptr<Sequence>(psc.sequence(i).clone());
        dataset->addIndividualSequenceInGroup(psc.getGroupId(i), i, 0, tmpSeq);
      }
      catch (Exception& e)
      {
        throw e;
      }
    }
  }
  return dataset;
}

