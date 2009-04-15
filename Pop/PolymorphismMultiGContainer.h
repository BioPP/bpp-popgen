//
// File PolymorphismMultiGContainer.h
// Author : Sylvain Gaillard <yragael2001@yahoo.fr>
// Last modification : Tuesday September 28 2004
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

#ifndef _POLYMORPHYSMMULTIGCONTAINER_H_
#define _POLYMORPHYSMMULTIGCONTAINER_H_

// From Utils
#include <Utils/Clonable.h>
#include <Utils/Exceptions.h>
#include <Utils/MapTools.h>
#include <Utils/TextTools.h>

// From popgenlib
#include "MultilocusGenotype.h"
#include "GeneralExceptions.h"

// From STL
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

namespace bpp
{

  /**
   * @brief The PolymorphismMultiGContainer class
   *
   * This class is a container of MultilocusGenotype.
   *
   * @author Sylvain Gaillard
   */
  class PolymorphismMultiGContainer
  {
    public: // Constructors and destructor

      /**
       * @brief Build a new PolymorphismMultilocusGenotypeContainer.
       */
      PolymorphismMultiGContainer();

      /**
       * @brief The copy constructor.
       */
      PolymorphismMultiGContainer(const PolymorphismMultiGContainer & pmgc);

      /**
       * @brief Destroy a PolymorphismMultilocusGenotypeContainer.
       */
      ~PolymorphismMultiGContainer();

    public:
      /**
       * @brief The assignation operator=.
       */
      PolymorphismMultiGContainer & operator= (const PolymorphismMultiGContainer & pmgc);

      /**
       * @brief Add a MultilocusGenotype to the container.
       */
      void addMultilocusGenotype(const MultilocusGenotype & mg, unsigned int group);

      /**
       * @brief Get a MultilocusGenotype at a position.
       *
       * @throw IndexOutOfBoundsException if position excedes the size of the container.
       */
      const MultilocusGenotype * getMultilocusGenotype(unsigned int position) const throw (IndexOutOfBoundsException);

      /**
       * @brief Remove a MultilocusGenotype.
       *
       * @throw IndexOutOfBoundsException if position excedes the size of the container.
       */
      MultilocusGenotype * removeMultilocusGenotype(unsigned int position) throw (IndexOutOfBoundsException);

      /**
       * @brief Delete a MultilocusGenotype.
       *
       * @throw IndexOutOfBoundsException if position excedes the size of the container.
       */
      void deleteMultilocusGenotype(unsigned int position) throw (IndexOutOfBoundsException);

      /**
       * @brief Tell if the MultilocusGenotypes are aligned (i.e. same size).
       */
      bool isAligned() const;

      /**
       * @brief Get the number of loci if the MultilocusGenotypes are aligned.
       *
       * @throw Exception if MultilocusGenotypes are not aligned.
       */
      unsigned int getNumberOfLoci() const throw (Exception);

      /**
       * @brief Get the Group id of a MultilocusGenotype.
       *
       * @throw IndexOutOfBoundsException if position excedes the size of the container.
       */
      unsigned int getGroupId(unsigned int position) const throw (IndexOutOfBoundsException);

      /**
       * @brief Set the Group id of a MultilocusGenotype.
       *
       * @throw IndexOutOfBoundsException if position excedes the size of the container.
       */
      void setGroupId(unsigned int position, unsigned int group_id) throw (IndexOutOfBoundsException);

      /**
       * @brief Get the groups' ids.
       */
      set<unsigned int> getAllGroupsIds() const;

      /**
       * @brief Get the groups names or ids if not available
       */
      vector<string> getAllGroupsNames() const;

      /**
       * @brief Tell if a group exists.
       */
      bool groupExists(unsigned int group) const;

      /**
       * @brief Get the number of groups.
       */
      unsigned int getNumberOfGroups() const;

      /**
       * @brief Get group size.
       */
      unsigned int getGroupSize(unsigned int group) const ;

      /**
       * @brief Get the group name for a given group id or just the id if not available juste return it's id
       */
      string getGroupName(unsigned int group_id) const throw (GroupNotFoundException);

      /**
       * @brief Set the name for the given group id.
       */
      void setGroupName(unsigned int group_id, string name)  throw (GroupNotFoundException);

      /**
       * @brief Inserts a name for the given group id.
       */
      void addGroupName(unsigned int group_id, string name)  ;

      /**
       * @brief Get the size of a group for a given locus.
       */
      unsigned int getLocusGroupSize(unsigned int group, unsigned int locus_position) const;

      /**
       * @brief Get the number of MultilocusGenotype.
       */
      unsigned int size() const;

      /**
       * @brief Clear the container.
       */
      void clear();

    protected:
      vector<MultilocusGenotype *> _multilocusGenotypes;
      vector<unsigned int> _groups;//group id for each multilocusgenotype
      map<unsigned int, string> _groups_names;

  };

} //end of namespace bpp;

#endif // _POLYMORPHYSMMULTIGCONTAINER_H_
