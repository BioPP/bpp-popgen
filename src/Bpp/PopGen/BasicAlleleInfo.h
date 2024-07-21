// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
class BasicAlleleInfo :
  public AlleleInfo
{
private:
  std::string id_;

public:
  // Constructors and destructor
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

public:
  // Methodes
  /**
   * @brief The assignation operator.
   */
  virtual BasicAlleleInfo& operator=(const BasicAlleleInfo& allele);

  /**
   * @brief The == operator.
   */
  virtual bool operator==(const BasicAlleleInfo& allele) const;

  /**
   * @brief The != operator.
   */
  virtual bool operator!=(const BasicAlleleInfo& allele) const;

  /**
   * @name The Clonable interface
   * @{
   */
#ifdef NO_VIRTUAL_COV
  Clonable*
#else
  BasicAlleleInfo*
#endif
  clone() const { return new BasicAlleleInfo(*this); }
  /** @} */

  /**
   * @name The AlleleInfo interface
   * @{
   */
  void setId(const std::string& allele_id);
  const std::string& getId() const;
  /** @} */
};
} // end of namespace bpp;

#endif // _BASICALLELEINFO_H_
