/*
 * File PoplibIO.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 08 2004
 */

// Secured inclusion of header's file
#ifndef _POPLIBIO_H_
#define _POPLIBIO_H_

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

/**
 * @brief The natif I/O format for poplib.
 */
class PoplibIO : public AbstractIDataSet, public AbstractODataSet {
	public: // Constantes
		static const string WHITESPACE;
		static const string TAB;
		static const string COMA;
		static const string SEMICOLON;

	public: // Constructor and destructor
		PoplibIO();
		PoplibIO(const string & missing_data_symbol, const string & data_separator) throw (Exception);
		~PoplibIO();

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
		void _parseLoci(const vector<string> & in, DataSet & data_set);
		void _parseIndividual(const vector<string> & in, DataSet & data_set, const VectorSequenceContainer & vsc);
};

#endif // _POPLIBIO_H_
