//
// File DarwinDon.h
// Author : Sylvain Gaillard
// Last modification : April 7, 2008
//

/*
   Copyright or © or Copr. CNRS, (April 7, 2008)

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

#ifndef _DARWIN_DON_H_
#define _DARWIN_DON_H_

// From Utils
#include <Utils/TextTools.h>
#include <Utils/FileTools.h>
#include <Utils/Exceptions.h>
#include <Utils/StringTokenizer.h>

// From local Pop
#include "AbstractODataSet.h"

namespace bpp
{

  /**
   * @brief The Darwin .don output format for popgenlib.
   *
   * @author Sylvain Gaillard
   */
  class DarwinDon:
    public AbstractODataSet
  {

    public: // Constructor and destructor
      DarwinDon();
      ~DarwinDon();

    public:
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
  };

} //end of namespace bpp;

#endif // _DARWIN_DON_H_

