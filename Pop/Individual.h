/*
 * File Individual.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday August 03 2004
 *
 * Copyright (C) 2004 Sylvain Gaillard and the
 *                    PopGenLib Development Core Team
 *
 * This file is part of PopGenLib.
 *
 * PopGenLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopGenLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopGenLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
#include <Utils/TextTools.h>

// From SeqLib
#include <Seq/Sequence.h>
#include <Seq/OrderedSequenceContainer.h>
#include <Seq/MapSequenceContainer.h>
#include <Seq/SequenceExceptions.h>

// From PopGenLib
#include "Locality.h"
#include "Coord.h"
#include "Date.h"
#include "MultilocusGenotype.h"
#include "GeneralExceptions.h"

/**
 * @brief The Individual class.
 *
 * This class is designed to store data on a single individual.
 * This individual has only one sequence for each locus ... no information
 * about diploid sequence data.
 * See the no more in use MultiSeqIndividual documentation for an alternative.
 */
class Individual {
	public: // Constructors and destructor :
		
		/**
		 * @brief Build a void new Individual.
		 */
		Individual();

		/**
		 * @brief Build a new Individual with an identifier.
		 */
		Individual(const string & id);

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
		Individual(const string & id,
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
		const Date * getDate() const throw (NullPointerException);

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
		const Coord<double> * getCoord() const throw(NullPointerException);

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
		void setLocality(const Locality<double> * locality);

		/**
		 * @brief Get the locality of the Individual.
		 *
		 * @return A pointer to the Locality of the Individual.
		 */
		const Locality<double> * getLocality() const throw (NullPointerException);

		/**
		 * @brief Tell if this Individual has a locality.
		 */
		bool hasLocality() const;

		/**
		 * @brief Add a sequence to the Individual.
		 *
		 * Creates the sequence container when adding the first sequence.
		 * Otherwize add the sequence to the end of the sequence container.
		 *
		 * @param sequence_key the place where the sequence will be put.
		 * @param sequence The sequence to add.
		 * @throw AlphabetMismatchException if the sequence's alphabet doesn't match the container's alphabet.
		 * @throw BadIdentifierException if sequence's name is already in use.
		 * @throw BadIntegerException if sequence_position is already in use.
		 */
		void addSequence(unsigned int sequence_key, const Sequence & sequence)
			throw (Exception);
		
		/**
		 * @brief Get a sequence by its name.
		 *
		 * @param sequence_name The name of the sequence.
		 * @return A pointer to the sequence.
		 * @throw NullPointerException if there is no sequence container defined.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 */
		const Sequence * getSequenceByName(const string & sequence_name)
			const throw(Exception);

		/**
		 * @brief Get a sequence by its position.
		 *
		 * @param sequence_position The position of the sequence in the sequence set.
		 * @return A pointer to the sequence.
		 * @throw NullPointerException if there is no sequence container defined.
		 * @throw SequenceNotFoundException if sequence_position is not found (i.e. missing data or not used).
		 */
		const Sequence * getSequenceAtPosition(const unsigned int sequence_position)
			const throw(Exception);
		 
		/**
		 * @brief Delete a sequence.
		 *
		 * @param sequence_name The name of the sequence.
		 * @throw NullPointerException if there is no sequence container defined.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 */
		void deleteSequenceByName(const string & sequence_name) throw (Exception);

		/**
		 * @brief Delete a sequence.
		 *
		 * @param sequence_position The position of the sequence.
		 * @throw NullPointerException if there is no sequence container defined.
		 * @throw SequenceNotFoundException if sequence_postion is not found.
		 */
		void deleteSequenceAtPosition(unsigned int sequence_position) throw (Exception);
		
		/**
		 * @brief Tell if the Individual has some sequences.
		 *
		 * @return TRUE if the individual has at least one sequence.
		 * @return FALSE if the container is empty or undifined.
		 */
		bool hasSequences() const;

		/**
		 * @brief Tell if the Individual has a sequence at a given position.
		 */
		bool hasSequenceAtPosition(unsigned int position) const;

		/**
		 * @brief Return the alphabet of the sequences.
		 *
		 * @throw NullPointerException if there is no sequence container defined.
		 */
		const Alphabet * getSequenceAlphabet() const throw (NullPointerException);

		/**
		 * @brief Get the sequences' names.
		 *
		 * @return All the sequences' names of the individual in a vector of string.
		 * @throw NullPointerException if there is no sequence container defined.
		 */
		vector<string> getSequencesNames() const throw (NullPointerException);

		/**
		 * @brief Get the sequences' positions.
		 *
		 * @return All the positions where a sequence is found.
		 * @throw NullPointerException if there is no sequence container defined.
		 */
		vector<unsigned int> getSequencesPositions() const throw (NullPointerException);

		/**
		 * @brief Get the position of a sequence.
		 *
		 * @throw NullPointerException if there is no sequence container defined.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 */
		unsigned int getSequencePosition(const string & sequence_name) const throw (Exception);
		
		/**
		 * @brief Get the number of sequences.
		 */
		unsigned int getNumberOfSequences() const;

		/**
		 * @brief Set all the sequences with a map  sequence container.
		 */
		void setSequences(const MapSequenceContainer & msc);

		/**
		 * @brief Get a pointer to the sequence container.
		 *
		 * @throw NullPointerException if there is no sequence container defined.
		 */
		const OrderedSequenceContainer * getSequences() const throw (NullPointerException);
		
		/**
		 * @brief Set a genotype.
		 *
		 * @param genotype The MultilocusGenotype which will be copied.
		 */
		void setGenotype(const MultilocusGenotype & genotype);

		/**
		 * @brief Init the genotype.
		 *
		 * @throw Exception if the Individual already has a Genotype.
		 * @throw BadIntegerException if loci_number < 1.
		 */
		void initGenotype(unsigned int loci_number) throw (Exception);

		/**
		 * @brief Get the genotype.
		 */
		const MultilocusGenotype * getGenotype() const throw (NullPointerException);
		
		/**
		 * @brief Delete the genotype of the individual.
		 */
		void deleteGenotype();

		/**
		 * @brief Tell if the Individual has a MultilocusGenotype.
		 */
		bool hasGenotype() const;

		/**
		 * @brief Set a MonolocusGenotype.
		 *
		 * @throw NullPointerException if there is no genotype defined.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
		 */
		void setMonolocusGenotype(unsigned int locus_position, const MonolocusGenotype & monogen)
			throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype.
		 *
		 * @throw NullPointerException if there is no genotype defined.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
		 * @throw Exception if there is no key in allele_keys.
		 */
		void setMonolocusGenotypeByAlleleKey(unsigned int locus_position, const vector<unsigned int> allele_keys)
			throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype.
		 *
		 * @throw NullPointerException if there is no genotype defined.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
		 * @throw AlleleNotFoundException if at least one the id is not found in the LocusInfo.
		 */
		void setMonolocusGenotypeByAlleleId(unsigned int locus_position, const vector<string> allele_id, const LocusInfo & locus_info)
			throw (Exception);

		/**
		 * @brief Get a MonolocusGenotype.
		 *
		 * @throw NullPointerException if there is no genotype defined.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
		 */
		const MonolocusGenotype * getMonolocusGenotype(unsigned int locus_position) throw (Exception);

		/**
		 * @brief Count the number of non missing MonolocusGenotype.
		 *
		 * @throw NullPointerException if there is no genotype defined.
		 */
		unsigned int countNonMissingLoci() const throw (NullPointerException);

		/**
		 * @brief Count the number of homozygous MonolocusGenotype.
		 *
		 * @throw NullPointerException if there is no genotype defined.
		 */
		unsigned int countHomozygousLoci() const throw (NullPointerException);

		/**
		 * @brief Count the number of heterozygous MonolocusGenotype.
		 *
		 * @throw NullPointerException if there is no genotype defined.
		 */
		unsigned int countHeterozygousLoci() const throw (NullPointerException);
	
	protected:
		string _id;
		unsigned short _sex;
		Date * _date;
		Coord<double> * _coord;
		const Locality<double> * _locality;
		MapSequenceContainer * _sequences;
		MultilocusGenotype * _genotype;
};
#endif // _INDIVIDUAL_H_
