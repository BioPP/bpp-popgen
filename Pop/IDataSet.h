/*
 * File IDataSet.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 24 2004
 */

// Secured inclusion of header's file
#ifndef _IDATASET_H_
#define _IDATASET_H_

#include "IODataSet.h"

// From Utils
#include <Utils/Exceptions.h>

/**
 * @brief The IDataSet interface.
 */
class IDataSet : public IODataSet {
	public: // Class destructor
		virtual ~IDataSet();

	public:
		/**
		 * @brief Read a DataSet on istream.
		 */
		virtual void read(istream & is, DataSet & data_set) throw (Exception) = 0;

		/**
		 * @brief Read a DataSet from a text file.
		 */
		virtual void read(const string & path, DataSet & data_set) throw (Exception) = 0;

		/**
		 * @brief Read istream and return a DataSet.
		 */
		virtual DataSet * read(istream & is) throw (Exception) = 0;

		/**
		 * @brief Read a text file and return a DataSet.
		 */
		virtual DataSet * read(const string & path) throw (Exception) = 0;
};

#endif // _IDATASET_H_
