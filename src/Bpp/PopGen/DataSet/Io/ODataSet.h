// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ODATASET_H_
#define _ODATASET_H_

#include "IODataSet.h"

#include <Bpp/Exceptions.h>

namespace bpp
{
/**
 * @brief The ODataSet interface.
 *
 * @author Sylvain Gaillard
 */
class ODataSet :
  public virtual IODataSet
{
public:
  // Class destructor
  virtual ~ODataSet() {}

public:
  /**
   * @brief Write a DataSet on ostream.
   */
  virtual void write(std::ostream& os, const DataSet& data_set) const = 0;

  /**
   * @brief Write a DataSet in a text file.
   */
  virtual void write(const std::string& path, const DataSet& data_set, bool overwrite) const = 0;
};
} // end of namespace bpp;

#endif// _ODATASET_H_
