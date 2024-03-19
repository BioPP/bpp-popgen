// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
