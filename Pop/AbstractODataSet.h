/*
 * File AbstractODataSet.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 24 2004
 */

// Secured inclusion of header's file
#ifndef _ABSTRACTODATASET_H_
#define _ABSTRACTODATASET_H_

#include "ODataSet.h"

class AbstractODataSet : public ODataSet {
	public: // Calss destructor
		virtual ~AbstractODataSet();

	public:
		/**
		 * @name The ODataSet interface.
		 * @{
		 */
		virtual void write(ostream & os, const DataSet & data_set) const throw (Exception) = 0;
		virtual void write(const string & path, const DataSet & data_set, bool overwrite) const throw (Exception);
		/**
		 * @}
		 */
};

#endif // _ABSTRACTODATASET_H_
