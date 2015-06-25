//
// File DataSetTools.h
// Author : Sylvain Gaillard
// Last modification : Wednesday August 04 2004
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

#ifndef _DATASETTOOLS_H_
#define _DATASETTOOLS_H_

// From STL
#include <set>
#include <memory>

#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>

// From SeqLib
#include <Bpp/Seq/Container/OrderedSequenceContainer.h>

// From local bpp-popgen
#include "DataSet.h"
#include "../PolymorphismSequenceContainer.h"

namespace bpp
{
/**
 * @brief A set of tools for DataSet.
 *
 * @author Sylvain Gaillard
 */
class DataSetTools
{
public:
  /**
   * @brief General method to build a DataSet from an OrderedSequenceContainer.
   */
  static std::auto_ptr<DataSet> buildDataSet(const OrderedSequenceContainer& osc) throw (Exception);

  /**
   * @brief Specific methode to build a DataSet from a PolymorphismSequenceContainer.
   */
  static std::auto_ptr<DataSet> buildDataSet(const PolymorphismSequenceContainer& psc) throw (Exception);
};
} // end of namespace bpp;

#endif // _DATASETTOOLS_H_

