/*
 * File Date.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday April 22 2004
 */

// Secured inclusion of header's file
#ifndef _DATE_H_
#define _DATE_h_

// From Utils
#include <Utils/Clonable.h>

/**
 * @brief The Date class
 *
 * This is a littleclass to deal with the dates.
 */
class Date : public Clonable {
	public: // Constructors and destructor
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
		 * @brief The Date copy constructor.
		 */
		Date(const Date & date);

		/**
		 * @brief Destroy the Date object.
		 */
		~Date();

	public: // Methodes

		/**
		 * @brief Implements the clonable interface.
		 */
		Clonable * clone() const;

		/**
		 * @brief The Date copy operator.
		 *
		 * @return A ref toward the assigned Date.
		 */
		Date & operator= (const Date & date);

		/**
		 * @brief Set the Date.
		 */
		void setDate(const int day, const int month, const int year);

		/**
		 * @brief Set the year.
		 */
		void setYear(const int year);

		/**
		 * @brief Set the month.
		 */
		void setMonth(const int month);

		/**
		 * @brief Set the day.
		 */
		void setDay(const int day);

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

		/**
		 * @brief The == operator.
		 *
		 * Test the numerical equality between to dates.
		 */
		bool operator== (const Date & date) const;

		/**
		 * @brief The < operator.
		 *
		 * Return true if the left Date is minor than the right Date.
		 */
		bool operator< (const Date & date) const;

		/**
		 * @brief The != operator.
		 */
		bool operator!= (const Date & date) { return !(*this == date); }

		/**
		 * @brief The > operator.
		 */
		bool operator> (const Date & date) { return date < *this; }

		/**
		 * @brief The <= operator.
		 */
		bool operator<= (const Date & date) { return !(date < *this); }

		/**
		 * @brief The >= operator.
		 */
		bool operator>= (const Date & date) { return !(*this < date); }

	protected:
		int _day;
		int _month;
		int _year;
};

#endif // -DATE_H_
