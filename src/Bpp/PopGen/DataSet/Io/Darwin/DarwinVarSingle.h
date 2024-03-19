// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _DARWIN_VAR_SINGLE_H_
#define _DARWIN_VAR_SINGLE_H_

#include <Bpp/Exceptions.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>

// From local Pop
#include "../AbstractODataSet.h"

namespace bpp
{
/**
 * @brief The Darwin .don output format for popgenlib.
 *
 * @author Sylvain Gaillard
 */
class DarwinVarSingle :
  public virtual AbstractODataSet
{
private:
  size_t missingData_;

public:
  // Constructor and destructor
  DarwinVarSingle(size_t missingData = 999);
  ~DarwinVarSingle();

public:
  /**
   * @name The ODataSet interface.
   * @{
   */
  void write(std::ostream& os, const DataSet& data_set) const;
  void write(const std::string& path, const DataSet& data_set, bool overwrite) const;
  /**
   * @}
   */

  /**
   * @name The IOFormat interface
   * @{
   */
  virtual const std::string getFormatName() const
  {
    return "Darwin .var single data";
  }
  virtual const std::string getFormatDescription() const
  {
    return "Darwin .var file store data for each marker in each individual (1 variable per allele).";
  }
  /**
   * @}
   */
};
} // end of namespace bpp;

#endif // _DARWIN_VAR_SINGLE_H_
