// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ABSTRACTIDATASET_H_
#define _ABSTRACTIDATASET_H_

#include "IDataSet.h"

#include <Bpp/Exceptions.h>

namespace bpp
{
/**
 * @brief Partial implementation of the DataSet Input interface
 *
 * @author Sylvain Gaillard
 */
class AbstractIDataSet :
  public IDataSet
{
public:
  // Class destructor
  virtual ~AbstractIDataSet();

public:
  /**
   * @name The IDataSet interface.
   * @{
   */
  virtual void read(std::istream& is, DataSet& data_set) = 0;

  virtual void read(const std::string& path, DataSet& data_set);

  virtual DataSet* read(std::istream& is);

  virtual DataSet* read(const std::string& path);
  /**
   * @}
   */
};
} // end of namespace bpp;

#endif// _ABSTRACTIDATASET_H_
