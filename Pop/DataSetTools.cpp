//
// File DataSetTools.cpp
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

#include "DataSetTools.h"

using namespace bpp;

DataSet * DataSetTools::buildDataSet(const OrderedSequenceContainer & osc) throw (Exception)
{
  DataSet * d_s = new DataSet();
  d_s->addEmptyGroup(0);
  for (unsigned int i = 0 ; i < osc.getNumberOfSequences() ; i++) {
    d_s->addEmptyIndividualToGroup(0, string("Individual_") + TextTools::toString(i + 1));
    try {
      d_s->addIndividualSequenceInGroup(0, i, 0, osc.getSequence(i));
    }
    catch (Exception & e) {
      throw e;
    }
  }
  return d_s;
}

DataSet * DataSetTools::buildDataSet(const PolymorphismSequenceContainer & psc) throw (Exception)
{
  DataSet * d_s = new DataSet();
  set<unsigned int> grp_ids = psc.getAllGroupsIds();
  for (set<unsigned int>::iterator it = grp_ids.begin() ; it != grp_ids.end() ; it++)
    d_s->addEmptyGroup(* it);
  unsigned int ind_count = 0;
  for (unsigned int i = 0 ; i < psc.getNumberOfSequences() ; i++) {
    for (unsigned int j = 0 ; j < psc.getSequenceCount(i) ; j++) {
      d_s->addEmptyIndividualToGroup(psc.getGroupId(i), string("Individual_") + TextTools::toString(ind_count++));
      try {
        d_s->addIndividualSequenceInGroup(psc.getGroupId(i), i, 0, psc.getSequence(i));
      }
      catch (Exception & e) {
        throw e;
      }
    }
  }
  return d_s;
}

