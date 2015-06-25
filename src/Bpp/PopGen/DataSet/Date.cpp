//
// File Date.cpp
// Author : Sylvain Gaillard <yragael2001@yahoo.fr>
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

#include <Bpp/Text/TextTools.h>

// From Local
#include "Date.h"

using namespace bpp;
using namespace std;

// ** Class constructor: *******************************************************/

Date::Date(const int day, const int month, const int year) throw (BadIntegerException) : day_(day),
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

void Date::setDate(const int day, const int month, const int year) throw (BadIntegerException)
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

void Date::setMonth(const int month) throw (BadIntegerException)
{
  if (month >= 1 && month <= 12)
    month_ = month;
  else
    throw (BadIntegerException("Date::Date: month must be in [1;12].", month));
}

void Date::setDay(const int day) throw (BadIntegerException)
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

