/*
 * File IODataSet.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday June 23 2004
 */

// Secured inclusion of header's file
#ifndef _IODATASET_H_
#define _IODATASET_H_

#include "DataSet.h"

// From STL
#include <iostream>
#include <fstream>
using namespace std;

/**
 * @brief Interface for input/ouput with DataSet.
 *
 * IODataSet is a virtual class.
 * This is an interface to declare commune methodes for in/out action on DataSet.
 */
class IODataSet {
	public: // Class destructor
		virtual ~IODataSet();

	public:
		/**
		 * @brief Get the format's name.
		 */
		virtual const string getFormatName() = 0;

		/**
		 * @brief Get a description of the format.
		 */
		virtual const string getFormatDescription() = 0;
};

#endif // _IODATASET_H_
