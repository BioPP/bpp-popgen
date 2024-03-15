// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Text/TextTools.h>

// From Local
#include "Date.h"

using namespace bpp;
using namespace std;

// ** Class constructor: *******************************************************/

Date::Date(const int day, const int month, const int year) : day_(day),
  month_(month),
  year_(year)
{
  if (day < 1 || day > 31)
    throw (BadIntegerException("Date::Date: day must be in [1;31].", day));
  if (month < 1 || month > 12)
    throw (BadIntegerException("Date::Date: month must be in [1;12].", month));
}

Date::Date(const Date& date) : day_(date.getDay()),
  month_(date.getMonth()),
  year_(date.getYear()) {}

// ** Class destructor: ********************************************************/

Date::~Date() {}

// ** Other methodes: **********************************************************/

Date& Date::operator=(const Date& date)
{
  day_ = date.getDay();
  month_ = date.getMonth();
  year_ = date.getYear();
  return *this;
}

void Date::setDate(const int day, const int month, const int year)
{
  if (day >= 1 && day <= 31)
    day_ = day;
  else
    throw (BadIntegerException("Date::Date: day must be in [1;31].", day));
  if (month >= 1 && month <= 12)
    month_ = month;
  else
    throw (BadIntegerException("Date::Date: month must be in [1;12].", month));
  year_ = year;
}

void Date::setYear(const int year)
{
  year_ = year;
}

void Date::setMonth(const int month)
{
  if (month >= 1 && month <= 12)
    month_ = month;
  else
    throw (BadIntegerException("Date::Date: month must be in [1;12].", month));
}

void Date::setDay(const int day)
{
  if (day >= 1 && day <= 31)
    day_ = day;
  else
    throw (BadIntegerException("Date::Date: day must be in [1;31].", day));
}

std::string Date::getDateStr() const
{
  string date, uDay = "", uMonth = "";
  if (day_ < 10)
    uDay = "0";
  if (month_ < 10)
    uMonth = "0";
  date = uDay + TextTools::toString(day_) + uMonth + TextTools::toString(month_) + TextTools::toString(year_);
  return date;
}

bool Date::operator==(const Date& date) const
{
  if (day_ == date.getDay() && month_ == date.getMonth() && year_ == date.getYear())
    return true;
  else
    return false;
}

bool Date::operator<(const Date& date) const
{
  if (year_ < date.getYear() || (month_ < date.getMonth() && year_ == date.getYear()) || (day_ < date.getDay() && month_ == date.getMonth() && year_ == date.getYear()))
    return true;
  else
    return false;
}
