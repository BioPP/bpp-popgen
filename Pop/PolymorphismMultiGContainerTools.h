//
// File PolymorphismMultiGContainerTools.h
// Authors : Sylvain Gailard
//           Khalid Belkhir
// Last modification : june 15 2006
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

#ifndef _POLYMORPHISMMULTIGCONTAINERTOOLS_H_
#define _POLYMORPHISMMULTIGCONTAINERTOOLS_H_

// From the STL
#include <vector>
#include <set>

using namespace std;

//From the PolGenLib library
#include "PolymorphismMultiGContainer.h"

//From the NumCalc library
#include <NumCalc/RandomTools.h>

namespace bpp
{

  /**
   * @brief Tools for PolymorphismMultiGContainer.
   *
   * Provides static methods for permutations.
   *
   * @author Sylvain Gaillard
   */
  class PolymorphismMultiGContainerTools
  {
    public:
      virtual ~PolymorphismMultiGContainerTools();

    public:
      /**
       * @brief Permut the MultilocusGenotype in the whole PolymorphismMultiGContainer.
       *
       * @param pmgc The PolymorphismMultiGContainer to permut.
       * @return A permuted PolymorphismMultiGContainer.
       */
      static PolymorphismMultiGContainer permutMultiG(const PolymorphismMultiGContainer & pmgc);

      /**
       * @brief Permut the MonolocusGenotype.
       *
       * Permut the MonolocusGenotypes in one or several groups breaking
       * the links between them.
       *
       * @param pmgc The PolymorphismMultiGContainer to permut.
       * @param groups The groups ids between which the MonolocusGenotypes will be permuted.
       * @return A permuted PolymorphismMultiGContainer.
       */
      static PolymorphismMultiGContainer permutMonoG(const PolymorphismMultiGContainer & pmgc, const set<unsigned int> & groups);

      /**
       * @brief Permut the MonolocusGenotype between individuals in the same group.
       *
       * Permut the MonolocusGenotypes for a set of groups. The idiv for the other groups
       * are kept intact
       *
       * @param pmgc The PolymorphismMultiGContainer to permut.
       * @param groups The groups ids for which the MonolocusGenotypes will be permuted.
       * @return A permuted PolymorphismMultiGContainer.
       */
      static PolymorphismMultiGContainer permutIntraGroupMonoG(const PolymorphismMultiGContainer & pmgc, const set<unsigned int> & groups); 

      /**
       * @brief Permut the Alleles.
       *
       * Permut the alleles in one or several groups breaking
       * the links between them.
       *
       * @param pmgc The PolymorphismMultiGContainer to permut.
       * @param groups The groups ids between which the MonolocusGenotypes will be permuted.
       * @return A permuted PolymorphismMultiGContainer.
       */
      static PolymorphismMultiGContainer permutAlleles(const PolymorphismMultiGContainer & pmgc, const set<unsigned int> & groups);

      /**
       * @brief Permut the Alleles between individuals in the same group.
       *
       * Permut the alleles in one or several groups
       *
       * @param pmgc The PolymorphismMultiGContainer to permut.
       * @param groups The groups ids between which the MonolocusGenotypes will be permuted.
       * @return A permuted PolymorphismMultiGContainer.
       */
      static PolymorphismMultiGContainer permutIntraGroupAlleles(const PolymorphismMultiGContainer & pmgc, const set<unsigned int> & groups) ;
      static PolymorphismMultiGContainer extractGroups(const PolymorphismMultiGContainer & pmgc, const set<unsigned int> & groups) ;
  };

} //end of namespace bpp;

#endif // _POLYMORPHISMMULTIGCONTAINERTOOLS_H_

