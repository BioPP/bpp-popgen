//
// File BasicAlleleInfo.h
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

#ifndef _BASICALLELEINFO_H_
#define _BASICALLELEINFO_H_

// From local Pop
#include "AlleleInfo.h"
#include "GeneralExceptions.h"

namespace bpp
{

  /**
   * @brief The BasicAlleleInfo class.
   *
   * This is the simplest allele class implementation which contains just an identitier.
   *
   * @author Sylvain Gaillard
   */
  class BasicAlleleInfo:
    public AlleleInfo
  {
    public: // Constructors and destructor
      /**
       * @brief Build a new allele.
       *
       * @param id The identity number of the allele.
       */
      BasicAlleleInfo(const std::string& id);

      /**
       * @brief The BasicAlleleInfo copy constructor.
       */
      BasicAlleleInfo(const BasicAlleleInfo& allele);

      virtual ~BasicAlleleInfo();

    public: // Methodes
      /**
       * @brief The assignation operator.
       */
      virtual BasicAlleleInfo& operator= (const BasicAlleleInfo& allele);

      /**
       * @brief The == operator.
       */
      virtual bool operator== (const BasicAlleleInfo& allele) const;

      /**
       * @brief The != operator.
       */
      virtual bool operator!= (const BasicAlleleInfo& allele) const;

      /**
       * @name The Clonable interface
       * @{
       */
      BasicAlleleInfo* clone() const { return new BasicAlleleInfo(*this); }
      /** @} */

      /**
       * @name The AlleleInfo interface
       */
      void setId(const std::string& allele_id);
      std::string getId() const;
      /** @} */

  private:
      std::string id_;
};

} //end of namespace bpp;

#endif // _BASICALLELEINFO_H_

