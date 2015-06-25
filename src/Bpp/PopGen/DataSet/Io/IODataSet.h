//
// File IODataSet.h
// Author : Sylvain Gaillard
// Last modification : Thursday July 29 2004
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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
// Secured inclusion of header's file
#ifndef _IODATASET_H_
#define _IODATASET_H_

#include "../DataSet.h"

#include <Bpp/Io/IoFormat.h>

// From STL
#include <iostream>
#include <fstream>

namespace bpp
{
/**
 * @brief Interface for input/ouput with DataSet.
 *
 * IODataSet is a virtual class.
 * This is an interface to declare commune methodes for in/out action on DataSet.
 *
 * @author Sylvain Gaillard
 */
class IODataSet : public virtual IOFormat
{
  /**
   * @name The IOFormat interface.
   * @{
   */
  const std::string getDataType() const { return "DataSet for population genetics"; }
  /**
   * @}
   */
};
} // end of namespace bpp;

#endif // _IODATASET_H_

