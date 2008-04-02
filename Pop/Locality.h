//
// File Locality.h
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

#ifndef _LOCALITY_H_
#define _LOCALITY_H_

#include "Coord.h"
#include <string>

using namespace std;

namespace bpp
{

  /**
   * @brief The Locality class.
   *
   * This is a class derivated from the Coord class.
   * It's a Coord with a name.
   */
  template <class T> class Locality:
    public Coord<T>
  {
    protected:
      string _name;

    public: // Constructors and destructor
      /**
       * @brief Build a new locality with name and coordinates.
       *
       * @param name The name of the locality.
       * @param x The longitude.
       * @param y The latitude.
       */
      Locality<T>(const string name, const T x=0, const T y=0):
        Coord<T>(x, y), _name(name) {}

      /**
       * @brief Build a new locality with name and coordinates.
       *
       * @param name The name of the locality.
       * @param coord The coordinates of the locality.
       */
      Locality<T>(const string name, const Coord<T> & coord):
        Coord<T>(coord), _name(name) {}

      /**
       * @brief Destroy a locality.
       */
      virtual ~Locality<T>() {}

    public: // Methodes
      /**
       * @brief Implements the Clonable interface.
       */
      Locality<T> * clone() const { return new Locality<T>(* this); }

      /**
       * @brief The == operator.
       *
       * returns true if both name and coordinates are identical between the two Locality objects.
       */
      virtual bool operator== (const Locality<T> & locality) const
      {
        return this->_x == locality.getX() && this->_y == locality.getY() && _name == locality._name;
      }

      /**
       * @brief The != operator.
       */
      virtual bool operator!= (const Locality<T> & locality) const
      {
        return !(locality == *this);
      }

      /**
       * @brief Set the name of the locality.
       */
      void setName(const string & name) { _name = name; }

      /**
       * @brief Get the name of the locality.
       */
      string getName() const { return _name; }

  };

} //end of namespace bpp;

#endif // _LOCALITY_H_

