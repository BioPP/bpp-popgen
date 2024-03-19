// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
