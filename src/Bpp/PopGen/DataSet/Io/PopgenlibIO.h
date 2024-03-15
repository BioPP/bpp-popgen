// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _POPULIBIO_H_
#define _POPULIBIO_H_

#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Io/FileTools.h>

// From Seq
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>

// From local Pop
#include "AbstractIDataSet.h"
#include "AbstractODataSet.h"
#include "../../BasicAlleleInfo.h"

namespace bpp
{
/**
 * @brief The native I/O format for popgenlib.
 *
 * @author Sylvain Gaillard
 */
class PopgenlibIO :
  public AbstractIDataSet,
  public AbstractODataSet
{
public:
  // Constantes
  static const std::string WHITESPACE;
  static const std::string TAB;
  static const std::string COMA;
  static const std::string SEMICOLON;

  static const std::string DIPLOID;
  static const std::string HAPLOID;
  static const std::string HAPLODIPLOID;
  static const std::string UNKNOWN;

private:
  char data_separator_;
  char missing_data_symbol_;

  std::vector<std::string> getValues_(std::string& param_line, const std::string& delim);
  void parseGeneral_(const std::vector<std::string>& in, DataSet& data_set);
  void parseLocality_(const std::vector<std::string>& in, DataSet& data_set);
  void parseSequence_(const std::vector<std::string>& in, VectorSequenceContainer& vsc);
  void parseLoci_(const std::vector<std::string>& in, std::vector<LocusInfo>& locus_info);
  void parseIndividual_(const std::vector<std::string>& in, DataSet& data_set, const VectorSequenceContainer& vsc);

public:
  // Constructor and destructor
  PopgenlibIO();
  PopgenlibIO(const std::string& missing_data_symbol, const std::string& data_separator);
  ~PopgenlibIO();

public:
  /**
   * @brief Get the code for missing data.
   */
  std::string getMissingDataSymbol() const;

  /**
   * @brief Get the code for data separator.
   */
  std::string getDataSeparator() const;

  /**
   * @brief Get the character for missing data.
   */
  char getMissingDataChar() const;

  /**
   * @brief Get the data separator char.
   */
  char getDataSeparatorChar() const;

  /**
   * @brief Set the code for missing data.
   *
   * The character used to code missing data can be every single non numerical
   * character and can't be the same used as data separator or a white space
   * or a tabulation.
   * The default value is '$'.
   *
   * @throw Excpetion if missing_data_symbol is a not allowed character.
   * @throw Exception if missing_data_symbol contains more than one character.
   */
  void setMissingDataSymbol(const std::string& missing_data_symbol);

  /**
   * @brief Set the code for data separator.
   *
   * The character used to separate data can be every single non numerical
   * character and can't be the same used for coding missing data.
   * Most common characters used are:
   * <ul><li>the white space: "WHITESPACE"</li>
   * <li>the tabulation: "TAB"</li>
   * <li>the coma: "COMA"</li>
   * <li>the semicolon: "SEMICOLON"</li></ul>
   * The default value is "WHITESPACE".
   *
   * @throw Exception if data_separator is a not allowed character.
   * @throw Exception if data_separator containes more than one character other than the two codes defined upper.
   */
  void setDataSeparator(const std::string& data_separator);

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
    return "PopgenlibIO ver 0.1";
  }
  const std::string getFormatDescription() const
  {
    return "IO format used to store DataSets inspired from Arlequin and Fasta";
  }
  /**
   * @}
   */
};
} // end of namespace bpp;

#endif// _POPULIBIO_H_
