// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef _ALLELEINFO_H_
#define _ALLELEINFO_H_

// From STL
#include <string>

#include <Bpp/Clonable.h>

namespace bpp
{
/**
 * @brief The AlleleInfo interface.
 *
 * An AlleleInfo is a data structure designed to store informations about
 * alleles in general like the size of the marker for example.
 *
 * @author Sylvain Gaillard
 */
class AlleleInfo :
  public virtual Clonable
{
public:
  // Destructor
  virtual ~AlleleInfo() {}

public:
  // Methodes
  /**
   * @brief Set the identifier of the allele.
   */
  virtual void setId(const std::string& allele_id) = 0;

  /**
   * @brief Get the identitier of the allele.
   */
  virtual const std::string& getId() const = 0;

  /**
   * @name The Clonable interface
   *
   * @{
   */
#ifdef NO_VIRTUAL_COV
  Clonable*
#else
  AlleleInfo*
#endif
  clone() const = 0;
  /** @} */
};
} // end of namespace bpp;

#endif // _ALLELEINFO_H_
