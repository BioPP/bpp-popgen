/*
 * File Individual.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday June 25 2004
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

// From PopLib
#include "Locality.h"
#include "Coord.h"
#include "Date.h"
#include "Genotype.h"
#include "GeneralExceptions.h"

/**
 * @brief The Individual class.
 *
 * This class is designed to store data on a single individual.
 * This individual has only one sequence for each locus ... no information
 * about diploid sequence data.
 * See the no more in use MultiSeqIndividual documentation for an alternative.
 */
class Individual : public Clonable {
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
		 * @throw BadIntegerException if sequence_index is already in use.
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
		 * @brief Get a sequence by its position (index).
		 *
		 * @param sequence_index The index of the sequence in the sequence set.
		 * @return A pointer to the sequence.
		 * @throw NullPointerException if there is no sequence container defined.
		 * @throw SequenceNotFoundException if sequence_index is not found (i.e. missing data or not used).
		 */
		const Sequence * getSequenceByIndex(const unsigned int sequence_index)
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
		 * @param sequence_index The index of the sequence.
		 * @throw NullPointerException if there is no sequence container defined.
		 * @throw SequenceNotFoundException if sequence_index is not found.
		 */
		void deleteSequenceByIndex(unsigned int sequence_index) throw (Exception);
		
		/**
		 * @brief Tell if the Individual has some sequences.
		 *
		 * @return TRUE if the individual has at least one sequence.
		 * @return FALSE if the container is empty or undifined.
		 */
		bool hasSequences() const;

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
		 * @brief Get the position (index) of a sequence.
		 *
		 * @throw NullPointerException if there is no sequence container defined.
		 * @throw SequenceNotFoundException if sequence_name is not found.
		 */
		unsigned int getSequencePosition(const string & sequence_name) const throw (Exception);
		
		/**
		 * @brief Get the number of sequences.
		 *
		 * @throw NullPointerException if there is no sequence container defined.
		 */
		unsigned int getNumberOfSequences() const throw (NullPointerException);

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
		 * @brief Add a genotype.
		 *
		 * @param genotype The Genotype to add.
		 */
		void addGenotype(const Genotype & genotype) throw (Exception);

		/**
		 * @brief Init the genotype.
		 *
		 * @throw NullPointerException if analyzed_loci is NULL.
		 */
		void initGenotype(const AnalyzedLoci * analyzed_loci) throw (Exception);

		/**
		 * @brief Get the genotype.
		 */
		const Genotype * getGenotype() const throw (NullPointerException);
		
		/**
		 * @brief Delete the genotype of the individual.
		 */
		void deleteGenotype();

		/**
		 * @brief Tell if the Individual has a Genotype.
		 */
		bool hasGenotype() const;

		/**
		 * @brief Set a MonolocusGenotype.
		 *
		 * @throw NullPointerException if there is no genotype defined.
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of loci.
		 */
		void setMonolocusGenotype(unsigned int locus_index, const MonolocusGenotype & monogen)
			throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype.
		 *
		 * @throw NullPointerException if there is no genotype defined.
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of loci.
		 * @throw Exception if the ploidy doesn't match.
		 */
		void setMonolocusGenotypeByAlleleKey(unsigned int locus_index, const vector<unsigned int> allele_keys)
			throw (Exception);

		/**
		 * @brief Set a MonolocusGenotype.
		 *
		 * @throw NullPointerException if there is no genotype defined.
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of loci.
		 * @throw Exception if the ploidy doesn' match.
		 */
		void setMonolocusGenotypeByAlleleId(unsigned int locus_index, const vector<unsigned int> allele_id)
			throw (Exception);

		/**
		 * @brief Get the ploidy of a locus.
		 *
		 * @throw NullPointerException if there is no genotype defined.
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of loci.
		 */
		unsigned int getPloidy(unsigned int locus_index) throw (Exception);

		/**
		 * @brief Get a MonolocusGenotype.
		 *
		 * @throw NullPointerException if there is no genotype defined.
		 * @throw IndexOutOfBoundsException if locus_index excedes the number of loci.
		 */
		const MonolocusGenotype * getMonolocusGenotype(unsigned int locus_index) throw (Exception);
	
	protected:
		string _id;
		unsigned short _sex;
		Date * _date;
		Coord<double> * _coord;
		const Locality<double> * _locality;
		MapSequenceContainer * _sequences;
		Genotype * _genotype;
};
#endif // _INDIVIDUAL_H_
