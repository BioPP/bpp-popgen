// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _LOCALITY_H_
#define _LOCALITY_H_

// From std lib
#include <string>

#include <Bpp/Graphics/Point2D.h>

namespace bpp
{
/**
 * @brief The Locality class.
 *
 * This is a class derivated from the Point2D class.
 * It's a Point2D with a name.
 *
 * @author Sylvain Gaillard
 */
template<class T> class Locality :
  public bpp::Point2D<T>
{
protected:
  std::string name_;

public:
  // Constructors and destructor
  /**
   * @brief Build a new locality with name and coordinates.
   *
   * @param name The name of the locality.
   * @param x The longitude.
   * @param y The latitude.
   */
  Locality<T>(const std::string name, const T x = 0, const T y = 0) :
    bpp::Point2D<T>(x, y),
    name_(name) {}

  /**
   * @brief Build a new locality with name and coordinates.
   *
   * @param name The name of the locality.
   * @param coord The coordinates of the locality.
   */
  Locality<T>(const std::string name, const bpp::Point2D<T>& coord) :
    bpp::Point2D<T>(coord),
    name_(name) {}

  /**
   * @brief Destroy a locality.
   */
  virtual ~Locality<T>() {}

public:
  // Methodes
  /**
   * @brief Implements the Clonable interface.
   */
  Locality<T>* clone() const { return new Locality<T>(*this); }

  /**
   * @brief The == operator.
   *
   * returns true if both name and coordinates are identical between the two Locality objects.
   */
  virtual bool operator==(const Locality<T>& locality) const
  {
    return this->getX() == locality.getX() && this->getY() == locality.getY() && name_ == locality.name_;
  }

  /**
   * @brief The != operator.
   */
  virtual bool operator!=(const Locality<T>& locality) const
  {
    return !(locality == *this);
  }

  /**
   * @brief Set the name of the locality.
   */
  void setName(const std::string& name) { name_ = name; }

  /**
   * @brief Get the name of the locality.
   */
  const std::string& getName() const { return name_; }
};
} // end of namespace bpp;

#endif// _LOCALITY_H_
