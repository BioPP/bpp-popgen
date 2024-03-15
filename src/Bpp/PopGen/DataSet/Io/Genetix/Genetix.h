// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _GENETIX_H_
#define _GENETIX_H_

#include <Bpp/Exceptions.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>

// From local Pop
#include "../AbstractIDataSet.h"
#include "../../../BasicAlleleInfo.h"

namespace bpp
{
/**
 * @brief The Genetix input format for popgenlib.
 *
 * @author Sylvain Gaillard
 */
class Genetix :
  public AbstractIDataSet
{
public:
  // Constructor and destructor
  Genetix();
  ~Genetix();

public:
  /**
   * @name The IDataSet interface.
   * @{
   */
  void read(std::istream& is, DataSet& data_set);
  void read(const std::string& path, DataSet& data_set);
  DataSet* read(std::istream& is);
  DataSet* read(const std::string& path);
  /**
   * @}
   */

  /**
   * @name The IOFormat interface
   * @{
   */
  const std::string getFormatName() const
  {
    return "Genetix ver 4.05";
  }
  const std::string getFormatDescription() const
  {
    return "Genetix is a software for populations genetic for Windows(tm)";
  }
  /**
   * @}
   */
};
} // end of namespace bpp;

#endif// _GENETIX_H_
