/*
 * File Individual.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday April 28 2004
 */

// Secured inclusion of header's file
#ifndef _INDIVIDUAL_H_
#define _INDIVIDUAL_H_

// From Utils
#include <Utils/Clonable.h>
#include <Utils/Exceptions.h>

// From PopLib
#include "Locality.h"
#include "Coord.h"
#include "Date.h"

/**
 * @brief The Individual class.
 *
 * This class is designed to store data on a single individual.
 */
class Individual : public Clonable {
	public: // Constructors and destructor :
		
		/**
		 * @brief Build a void new Individual.
		 */
		Individual();

		/**
		 * @brief Build a new Individual with parameters.
		 *
		 * @param date The date of the Individual as a Date object.
		 * @param coord The coordinates of the Individual as a Coord object.
		 * @param locality The locality of the Individual as a Locality object.
		 * @param sex The sex of the Individual as an unsigned short.
		 */
		Individual(const Date & date, const Coord<double> & coord, const Locality<double> & locality, const unsigned short sex);
		
		/**
		 * @brief The Individual copy constructor.
		 */
		Individual(const Individual &ind);

		/**
		 * @brief Destroy an Individual.
		 */
		virtual ~Individual();

	public: // Methodes
		
		/**
		 * @brief Implement the Clonable interface.
		 */
		Clonable * clone() const;

		/**
		 * @brief The Individual copy operator.
		 *
		 * @return A ref toward the assigned Individual.
		 */
		Individual & operator= (const Individual & ind);
		
		/**
		 * @brief Set the sex of the Individual.
		 *
		 * @param sex An unsigned short coding for the sex.
		 */
		void setSex(const unsigned short sex);

		/**
		 * @brief Get the sex of the Individual.
		 *
		 * @return The sex of the Individual as an unsigned short.
		 */
		unsigned short getSex() const;

		/**
		 * @brief Set the date of the Individual.
		 *
		 * @param date The date as a Date object.
		 */
		void setDate(const Date & date);

		/**
		 * @brief Get the date of the Individual.
		 *
		 * @return A pointer toward a Date object if the Individual has a date. Otherwise throw a NullPointerException.
		 */
		Date * getDate() const throw(NullPointerException);

		/**
		 * @brief Tell if this Individual has a date.
		 */
		bool hasDate() const;

		/**
		 * @brief Set the coodinates of the Individual.
		 *
		 * @param coord A Coord object.
		 */
		void setCoord(const Coord<double> & coord);

		/**
		 * @brief Set the coordinates of the Individual.
		 *
		 * @param x The X coordinate as a double.
		 * @param y The Y coordinate as a double.
		 */
		void setCoord(const double x, const double y);

		/**
		 * @brief Get the coordinates of the Induvidual.
		 *
		 * @return A pointer toward a Coord object if the Individual has coordinates. Otherwise throw a NullPointerException.
		 */
		Coord<double> * getCoord() const throw(NullPointerException);

		/**
		 * @brief Tell if this Individual has coordinates.
		 */
		bool hasCoord() const;

		/**
		 * @brief Set the X coordinate of the Individual.
		 *
		 * @param x The X coordinate as a double.
		 *
		 * Set the X coordinate if the Individual has coordinates. Otherwise throw a NullPointerException.
		 */
		void setX(const double x) throw(NullPointerException);

		/**
		 * @brief Set the Y coordinate of th Individual.
		 *
		 * @param The Y coordinate as a double.
		 *
		 * Set the Y coordinate if the Individual has coordinates. Otherwise throw a NullPointerException.
		 */
		void setY(const double y) throw(NullPointerException);

		/**
		 * @brief Get the X coordinate of the Individual.
		 *
		 * @return The X coordinate as a double if the Individual has coordinates. Otherwise throw a NullPointerException.
		 */
		double getX() const throw(NullPointerException);

		/**
		 * @brief Get the Y coordinate of the Individual.
		 *
		 * @return The Y coordinate as a double if the Individual has coordinates. Otherwise throw a NullPointerException.
		 */
		double getY() const throw(NullPointerException);

		/**
		 * @brief Set the locality of the Individual.
		 *
		 * @param locality A Locality object.
		 */
		void setLocality(const Locality<double> & locality);
		
		/**
		 * @brief Set the locality of the Individual.
		 *
		 * @param name The name of the locality as a string.
		 * @param coord The coordinates of the locality as a Coord object.
		 */
		void setLocality(const string name, const Coord<double> & coord);

		/**
		 * @brief Set the locality of the Individual.
		 *
		 * @param name The name of the locality as a string.
		 * @param x The X coordinate of the locality as a double.
		 * @param y The Y coordinate of the locality as a double.
		 */
		void setLocality(const string name, const double x, const double y);

		/**
		 * @brief Get the locality of the Individual.
		 *
		 * @return A pointer toward a Locality object if the Individual has a locality. Otherwise throw a NullPointerException.
		 */
		Locality<double> * getLocality() const throw(NullPointerException);

		/**
		 * @brief Tell if this Individual has a locality.
		 */
		bool hasLocality() const;
		
	protected:
		unsigned short _sex;
		Date * _date;
		Coord<double> * _coord;
		Locality<double> * _locality;
};
#endif // _INDIVIDUAL_H_
