//
// File Date.h
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

#ifndef _DATE_H_
#define _DATE_H_

#include <Bpp/Exceptions.h>
#include <Bpp/Clonable.h>

namespace bpp
{
/**
 * @brief The Date class
 *
 * This is a little class to deal with dates.
 *
 * @author Sylvain Gaillard
 */
class Date : public Clonable
{
private:
  int day_;
  int month_;
  int year_;

public:
  // Constructors and destructor
  /**
   * @brief Build a new Date from three values.
   *
   * Build a new Date from three integers.
   * The default Date is set to 01-01-2000.
   *
   * @param day The day between 1 and 31.
   * @param month The month between 1 and 12.
   * @param year The year as a signed int.
   */
  Date(const int day = 1, const int month = 1, const int year = 2000);

  /**
   * @brief The Date copy constructor.
   */
  Date(const Date& date);

  /**
   * @brief Destroy the Date object.
   */
  ~Date();

public:
  // Methodes
  /**
   * @brief The Date copy operator.
   *
   * @return A ref toward the assigned Date.
   */
  Date& operator=(const Date& date);

  /**
   * @brief Set the Date.
   *
   * @param day The day as an integer between 1 and 31.
   * @param month The month as an integer between 1 and 12.
   * @param year The year as an integer.
   */
  void setDate(const int day, const int month, const int year);

  /**
   * @brief Set the year.
   *
   * @param year The year as an integer.
   */
  void setYear(const int year);

  /**
   * @brief Set the month.
   *
   * @param month The month as an integer between 1 and 12.
   */
  void setMonth(const int month);

  /**
   * @brief Set the day.
   *
   * @param day The day as an integer between 1 and 31.
   */
  void setDay(const int day);

  /**
   * @brief Get the Date as a string.
   *
   * @return The date as a string DDMMYYYY (i.e. January 1 2000 : 01012000).
   */
  std::string getDateStr() const;

  /**
   * @brief Get the Year as an int.
   */
  int getYear() const { return year_; }

  /**
   * @brief Get the month as an int.
   */
  int getMonth() const { return month_; }

  /**
   * @brief Get the day as an int.
   */
  int getDay() const { return day_; }

  /**
   * @brief The == operator.
   *
   * Test the numerical equality between to dates.
   */
  bool operator==(const Date& date) const;

  /**
   * @brief The < operator.
   *
   * Return true if the left Date is minor than the right Date.
   */
  bool operator<(const Date& date) const;

  /**
   * @brief The != operator.
   */
  bool operator!=(const Date& date) const { return !(*this == date); }

  /**
   * @brief The > operator.
   */
  bool operator>(const Date& date) const { return date < *this; }

  /**
   * @brief The <= operator.
   */
  bool operator<=(const Date& date) const { return !(date < *this); }

  /**
   * @brief The >= operator.
   */
  bool operator>=(const Date& date) const { return !(*this < date); }

  /**
   * @name The Clonable interface
   * @{
   */
#ifdef NO_VIRTUAL_COV
  Clonable*
#else
  Date*
#endif
  clone() const { return new Date(*this); }
};
} // end of namespace bpp;

#endif // _DATE_H_

