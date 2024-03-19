// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _DARWIN_DON_H_
#define _DARWIN_DON_H_

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
class DarwinDon :
  public virtual AbstractODataSet
{
public:
  // Constructor and destructor
  DarwinDon();
  ~DarwinDon();

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
  const std::string getFormatName() const
  {
    return "Darwin .don";
  }
  const std::string getFormatDescription() const
  {
    return "Darwin .don file store data identifying individuals.";
  }
  /**
   * @}
   */
};
} // end of namespace bpp;

#endif // _DARWIN_DON_H_
