// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ABSTRACTODATASET_H_
#define _ABSTRACTODATASET_H_

#include "ODataSet.h"

namespace bpp
{
/**
 * @brief Partial implementation of the DataSet Output interface.
 *
 * @author Sylvain Gaillard
 */
class AbstractODataSet :
  public ODataSet
{
public:
  virtual ~AbstractODataSet();

public:
  /**
   * @name The ODataSet interface.
   * @{
   */
  virtual void write(std::ostream& os, const DataSet& data_set) const = 0;
  virtual void write(const std::string& path, const DataSet& data_set, bool overwrite) const;
  /**
   * @}
   */
};
} // end of namespace bpp;

#endif // _ABSTRACTODATASET_H_
