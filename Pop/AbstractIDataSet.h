/*
 * File AbstractIDataSet.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 24 2004
 */

// Secured inclusion of header's file
#ifndef _ABSTRACTIDATASET_H_
#define _ABSTRACTIDATASET_H_

#include "IDataSet.h"

// From Utils
#include <Utils/Exceptions.h>

class AbstractIDataSet : public IDataSet {
	public: // Class destructor
		virtual ~AbstractIDataSet();

	public:
		/**
		 * @name The IDataSet interface.
		 * @{
		 */
		virtual void read(istream & is, DataSet & data_set) throw (Exception) = 0;

		virtual void read(const string & path, DataSet & data_set) throw (Exception);

		virtual DataSet * read(istream & is) throw (Exception);

		virtual DataSet * read(const string & path) throw (Exception);
		/**
		 * @}
		 */
};

#endif // _ABSTRACTIDATASET_H_
