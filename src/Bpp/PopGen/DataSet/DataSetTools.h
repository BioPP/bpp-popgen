// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _DATASETTOOLS_H_
#define _DATASETTOOLS_H_

// From STL
#include <set>
#include <memory>

#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>

// From bpp-seq
#include <Bpp/Seq/Container/SequenceContainer.h>

// From local bpp-popgen
#include "DataSet.h"
#include "../PolymorphismSequenceContainer.h"

namespace bpp
{
/**
 * @brief A set of tools for DataSet.
 *
 * @author Sylvain Gaillard
 */
class DataSetTools
{
public:
  /**
   * @brief General method to build a DataSet from an OrderedSequenceContainer.
   */
  static std::unique_ptr<DataSet> buildDataSet(const SequenceContainerInterface& sc);

  /**
   * @brief Specific methode to build a DataSet from a PolymorphismSequenceContainer.
   */
  static std::unique_ptr<DataSet> buildDataSet(const PolymorphismSequenceContainer& psc);
};
} // end of namespace bpp;

#endif// _DATASETTOOLS_H_
