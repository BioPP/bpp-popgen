/*
 * File Individual.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday April 27 2004
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
		 * @return The date of the Individual as a Date.
		 */
		Date getDate() const throw(NullPointerException);

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
		 * @return The coordinates of the Individual as a Coord object.
		 */
		Coord<double> getCoord() const throw(NullPointerException);

		/**
		 * @brief Set the X coordinate of the Individual.
		 */
		void setX(const double x) throw(NullPointerException);

		/**
		 * @brief Set the Y coordinate of th Individual.
		 */
		void setY(const double y) throw(NullPointerException);

		/**
		 * @brief Get the X coordinate of the Individual.
		 */
		double getX() const throw(NullPointerException);

		/**
		 * @brief Get the Y coordinate of the Individual.
		 */
		double getY() const throw(NullPointerException);

		/**
		 * @brief Set the locality of the Individual.
		 */
		void setLocality(const Locality<double> & locality);
		
		/**
		 * @brief Set the locality of the Individual.
		 */
		void setLocality(const string name, const Coord<double> & coord);

		/**
		 * @brief Set the locality of the Individual.
		 */
		void setLocality(const string name, const double x, const double y);

		/**
		 * @brief Get the locality of the Individual.
		 *
		 * @return The locality of the Individual as a Locality object.
		 */
		Locality<double> getLocality() const throw(NullPointerException);

	protected:
		unsigned short _sex;
		Date * _date;
		Coord<double> * _coord;
		Locality<double> * _locality;
};
#endif // _INDIVIDUAL_H_
