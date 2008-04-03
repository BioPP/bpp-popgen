//
// File PopgenlibIO.h
// Author : Sylvain Gaillard
// Last modification : Thursday July 29 2004
//

/*
   Copyright or © or Copr. CNRS, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for population genetics analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use, 
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info". 

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability. 

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or 
   data to be ensured and,  more generally, to use and operate it in the 
   same conditions as regards security. 

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
   */

#ifndef _POPULIBIO_H_
#define _POPULIBIO_H_

// From Utils
#include <Utils/TextTools.h>
#include <Utils/FileTools.h>
#include <Utils/Exceptions.h>

// From Seq
#include <Seq/Fasta.h>
#include <Seq/VectorSequenceContainer.h>

// From local Pop
#include "AbstractIDataSet.h"
#include "AbstractODataSet.h"
#include "BasicAlleleInfo.h"

namespace bpp
{

  /**
   * @brief The natif I/O format for popgenlib.
   */
  class PopgenlibIO:
    public AbstractIDataSet,
    public AbstractODataSet
  {
    public: // Constantes
      static const string WHITESPACE;
      static const string TAB;
      static const string COMA;
      static const string SEMICOLON;

      static const string DIPLOID;
      static const string HAPLOID;
      static const string HAPLODIPLOID;
      static const string UNKNOWN;

    public: // Constructor and destructor
      PopgenlibIO();
      PopgenlibIO(const string & missing_data_symbol, const string & data_separator) throw (Exception);
      ~PopgenlibIO();

    public:
      /**
       * @brief Get the code for missing data.
       */
      string getMissingDataSymbol() const;

      /**
       * @brief Get the code for data separator.
       */
      string getDataSeparator() const;

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
      void setMissingDataSymbol(const string & missing_data_symbol)
        throw (Exception);

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
      void setDataSeparator(const string & data_separator)
        throw (Exception);

      /**
       * @name The IDataSet interface.
       * @{
       */
      void read(istream & is, DataSet & data_set) throw (Exception);
      void read(const string & path, DataSet & data_set) throw (Exception);
      DataSet * read(istream & is) throw (Exception);
      DataSet * read(const string & path) throw (Exception);
      /**
       * @}
       */

      /**
       * @name The ODataSet interface.
       * @{
       */
      void write(ostream & os, const DataSet & data_set) const throw (Exception);

      void write(const string & path, const DataSet & data_set, bool overwrite) const throw (Exception);
      /**
       * @}
       */

      /**
       * @name The IODataSet interface
       * @{
       */
      virtual const string getFormatName();
      virtual const string getFormatDescription();
      /**
       * @}
       */

    protected:
      char _missing_data_symbol;
      char _data_separator;

      vector<string> _getValues(string & param_line, const string & delim);
      void _parseGeneral(const vector<string> & in, DataSet & data_set);
      void _parseLocality(const vector<string> & in, DataSet & data_set);
      void _parseSequence(const vector<string> & in, VectorSequenceContainer & vsc);
      void _parseLoci(const vector<string> & in, vector<LocusInfo> & locus_info);
      void _parseIndividual(const vector<string> & in, DataSet & data_set, const VectorSequenceContainer & vsc);
  };

} //end of namespace bpp;

#endif // _POPULIBIO_H_

