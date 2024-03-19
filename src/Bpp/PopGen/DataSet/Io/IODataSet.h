// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// Secured inclusion of header's file
#ifndef _IODATASET_H_
#define _IODATASET_H_

#include "../DataSet.h"

#include <Bpp/Io/IoFormat.h>

// From STL
#include <iostream>
#include <fstream>

namespace bpp
{
/**
 * @brief Interface for input/ouput with DataSet.
 *
 * IODataSet is a virtual class.
 * This is an interface to declare commune methodes for in/out action on DataSet.
 *
 * @author Sylvain Gaillard
 */
class IODataSet : public virtual IOFormat
{
  /**
   * @name The IOFormat interface.
   * @{
   */
  const std::string getDataType() const { return "DataSet for population genetics"; }
  /**
   * @}
   */
};
} // end of namespace bpp;

#endif // _IODATASET_H_
