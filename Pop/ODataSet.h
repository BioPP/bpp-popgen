/*
 * File ODataSet.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 24 2004
 */

// Secured inclusion of header's file
#ifndef _ODATASET_H_
#define _ODATASET_H_

#include "IODataSet.h"

// From Utils
#include <Utils/Exceptions.h>

/**
 * @brief The ODataSet interface.
 */
class ODataSet : public IODataSet {
	public: // Class destructor
		virtual ~ODataSet();

	public:
		/**
		 * @brief Write a DataSet on ostream.
		 */
		virtual void write(ostream & os, const DataSet & data_set) const throw (Exception) = 0;

		/**
		 * @brief Write a DataSet in a text filz.
		 */
		virtual void write(const string & path, const DataSet & data_set, bool overwrite) const throw (Exception) = 0;
};

#endif // _ODATASET_H_
