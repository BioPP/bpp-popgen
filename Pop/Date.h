/*
 * File Date.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday April 14 2004
 */

// Secured inclusion of header's file
#ifndef _DATE_H_
#define _DATE_h_

/**
 * @brief The Date class
 *
 * This is a littleclass to deal with the dates.
 */
class Date {
	public:
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
		Date(const int day=1, const int month=1, const int year=2000);

		/**
		 * @brief Destroy the Date object.
		 */
		~Date();

		/**
		 * @brief Set the Date.
		 */
		void setDate(const int day, const int month, const int year);

		/**
		 * @brief Get the Date as a string.
		 */
		string getDate() const;

		/**
		 * @brief Get the Year as an int.
		 */
		int getYear() const;

		/**
		 * @brief Get the month as an int.
		 */
		int getMonth() const;

		/**
		 * @brief Get the day as an int.
		 */
		int getDay() const;
	protected:
		int _day;
		int _month;
		int _year;
};

#endif // -DATE_H_
