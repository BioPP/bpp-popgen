/*
 * File Individual.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday May 24 2004
 */

// Secured inclusion of header's file
#ifndef _INDIVIDUAL_H_
#define _INDIVIDUAL_H_

// From STL
#include <map>
#include <vector>

// From Utils
#include <Utils/Clonable.h>
#include <Utils/Exceptions.h>

// From SeqLib
#include <Seq/Sequence.h>
#include <Seq/VectorSequenceContainer.h>
#include <Seq/SequenceExceptions.h>

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
		 * @param id The id of the Individual as a string.
		 * @param date The date of the Individual as a Date object.
		 * @param coord The coordinates of the Individual as a Coord object.
		 * @param locality The locality of the Individual as a pointer to a Locality
		 * object.
		 * @param sex The sex of the Individual as an unsigned short.
		 */
		Individual(const string id,
		           const Date & date,
		           const Coord<double> & coord,
		           Locality<double> * locality,
		           const unsigned short sex);
		
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
		 * Make a copy of each atribute of the Individual.
		 */
		Individual & operator= (const Individual & ind);
		
		/**
		 * @brief Set the id of the Individual.
		 *
		 * @param id The id of the Individual as a string.
		 */
		void setId(const string id);

		/**
		 * @brief Get the id of the Individual.
		 *
		 * @return The id of the Individual as a string.
		 */
		string getId() const;

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
		 * @return A pointer toward a Date object if the Individual has a date.
		 * Otherwise throw a NullPointerException.
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
		 * @return A pointer toward a Coord object if the Individual has
		 * coordinates. Otherwise throw a NullPointerException.
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
		 * Set the X coordinate if the Individual has coordinates.
		 * Otherwise throw a NullPointerException.
		 */
		void setX(const double x) throw(NullPointerException);

		/**
		 * @brief Set the Y coordinate of th Individual.
		 *
		 * @param y The Y coordinate as a double.
		 *
		 * Set the Y coordinate if the Individual has coordinates.
		 * Otherwise throw a NullPointerException.
		 */
		void setY(const double y) throw(NullPointerException);

		/**
		 * @brief Get the X coordinate of the Individual.
		 *
		 * @return The X coordinate as a double if the Individual has coordinates.
		 * Otherwise throw a NullPointerException.
		 */
		double getX() const throw(NullPointerException);

		/**
		 * @brief Get the Y coordinate of the Individual.
		 *
		 * @return The Y coordinate as a double if the Individual has coordinates.
		 * Otherwise throw a NullPointerException.
		 */
		double getY() const throw(NullPointerException);

		/**
		 * @brief Set the locality of the Individual.
		 *
		 * @param locality A pointer to a Locality object.
		 */
		void setLocality(Locality<double> * locality);

		/**
		 * @brief Get the locality of the Individual.
		 *
		 * @return A pointer to the Locality of the Individual.
		 */
		Locality<double> * getLocality() const;

		/**
		 * @brief Tell if this Individual has a locality.
		 */
		bool hasLocality() const;

		/**
		 * @brief Add a sequence in a named sequence set.
		 *
		 * @param id The id of the sequence set.
		 * @param sequence The sequence to add.
		 */
		void addSequence(const string & id, const Sequence & sequence)
			throw (Exception);
		
		/**
		 * @brief Get a named sequence from a named sequence set.
		 *
		 * @param id The id of the sequence set.
		 * @param name The name of the sequence.
		 *
		 * @return A pointer to the sequence.
		 */
		const Sequence * getSequence(const string & id, const string & name)
			const throw(Exception);

		/**
		 * @brief Get an indexed sequence from a named sequence set.
		 *
		 * @param id The id of the sequence set.
		 * @param i The index of the sequence in the sequence set.
		 *
		 * @return A pointer tothe sequence.
		 */
		const Sequence * getSequence(const string & id, const unsigned int i)
			const throw(Exception);
		 
		/**
		 * @brief Get the sequence set ids.
		 *
		 * @return All the keys of the sequence sets in a vector.
		 */
		vector<string> getSequencesKeys() const;

		/**
		 * @brief Remove a named sequence from a named sequence set.
		 *
		 * @param id The id of the sequence set.
		 * @param name The name of the sequence.
		 *
		 * @return A pointer to a copy of the removed sequence.
		 */
		Sequence * removeSequence(const string & id, const string & name);

		/**
		 * @brief Delete a named sequence from a named sequence set.
		 *
		 * @param id The id of the sequence set.
		 * @param name The name of the sequence.
		 */
		void deleteSequence(const string & id, const string & name);

	protected:
		string _id;
		unsigned short _sex;
		Date * _date;
		Coord<double> * _coord;
		Locality<double> * _locality;
		map<string,VectorSequenceContainer *> _sequences;
};
#endif // _INDIVIDUAL_H_
