// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _GENEPOP_H_
#define _GENEPOP_H_

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
 * @brief The Genepop input format for popgenlib.
 *
 * @author Sylvain Gaillard
 */
class Genepop :
  public AbstractIDataSet
{
public:
  // Constructor and destructor
  Genepop();
  ~Genepop();

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
    return "Genepop ver 3.4";
  }

  const std::string getFormatDescription() const
  {
    return "Genepop is a software for populations genetic for DOS operating system";
  }
  /**
   * @}
   */
};
} // end of namespace bpp;

#endif // _GENEPOP_H_
