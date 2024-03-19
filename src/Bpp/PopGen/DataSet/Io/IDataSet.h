// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _IDATASET_H_
#define _IDATASET_H_

#include "IODataSet.h"

#include <Bpp/Exceptions.h>

namespace bpp
{
/**
 * @brief The IDataSet interface.
 *
 * @author Sylvain Gaillard
 */
class IDataSet :
  public virtual IODataSet
{
public:
  // Class destructor
  virtual ~IDataSet() {}

public:
  /**
   * @brief Read a DataSet on istream.
   */
  virtual void read(std::istream& is, DataSet& data_set) = 0;

  /**
   * @brief Read a DataSet from a text file.
   */
  virtual void read(const std::string& path, DataSet& data_set) = 0;

  /**
   * @brief Read istream and return a DataSet.
   */
  virtual DataSet* read(std::istream& is) = 0;

  /**
   * @brief Read a text file and return a DataSet.
   */
  virtual DataSet* read(const std::string& path) = 0;
};
} // end of namespace bpp;

#endif // _IDATASET_H_
