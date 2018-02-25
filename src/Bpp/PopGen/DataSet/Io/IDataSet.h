//
// File IDataSet.h
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

#ifndef _IDATASET_H_
#define _IDATASET_H_

#include "IODataSet.h"

#include <Bpp/Exceptions.h>

namespace bpp
{
/**
 * @brief The IDataSet interface.
 *
 * @author Sylvain Gaillard
 */
class IDataSet :
  public virtual IODataSet
{
public:
  // Class destructor
  virtual ~IDataSet() {}

public:
  /**
   * @brief Read a DataSet on istream.
   */
  virtual void read(std::istream& is, DataSet& data_set) = 0;

  /**
   * @brief Read a DataSet from a text file.
   */
  virtual void read(const std::string& path, DataSet& data_set) = 0;

  /**
   * @brief Read istream and return a DataSet.
   */
  virtual DataSet* read(std::istream& is) = 0;

  /**
   * @brief Read a text file and return a DataSet.
   */
  virtual DataSet* read(const std::string& path) = 0;
};
} // end of namespace bpp;

#endif // _IDATASET_H_

