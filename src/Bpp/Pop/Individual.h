//
// File Individual.h
// Author : Sylvain Gaillard
// Last modification : Tuesday August 03 2004
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

#ifndef _INDIVIDUAL_H_
#define _INDIVIDUAL_H_

// From STL
#include <vector>
#include <memory>

#include <Bpp/Graphics/Point2D.h>
#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>

// From SeqLib
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/SequenceExceptions.h>
#include <Bpp/Seq/Container/OrderedSequenceContainer.h>
#include <Bpp/Seq/Container/MapSequenceContainer.h>

// From PopGenLib
#include "Locality.h"
#include "Date.h"
#include "MultilocusGenotype.h"
#include "GeneralExceptions.h"

namespace bpp
{

  /**
   * @brief The Individual class.
   *
   * This class is designed to store data on a single individual.
   * This individual has only one sequence for each locus ... no information
   * about diploid sequence data.
   * See the no more in use MultiSeqIndividual documentation for an alternative.
   *
   * @author Sylvain Gaillard
   */
  class Individual
  {
    protected:
      std::string id_;
      unsigned short sex_;
      std::auto_ptr<Date> date_;
      std::auto_ptr< Point2D<double> > coord_;
      const Locality<double>* locality_;
      std::auto_ptr<MapSequenceContainer> sequences_;
      std::auto_ptr<MultilocusGenotype> genotype_;

    public: // Constructors and destructor :

      /**
       * @brief Build a void new Individual.
       */
      Individual();

      /**
       * @brief Build a new Individual with an identifier.
       */
      Individual(const std::string& id);

      /**
       * @brief Build a new Individual with parameters.
       *
       * @param id The id of the Individual as a string.
       * @param date The date of the Individual as a Date object.
       * @param coord The coordinates of the Individual as a Point2D object.
       * @param locality The locality of the Individual as a pointer to a Locality
       * object.
       * @param sex The sex of the Individual as an unsigned short.
       */
      Individual(const std::string& id,
          const Date& date,
          const Point2D<double>& coord,
          Locality<double>* locality,
          const unsigned short sex);

      /**
       * @brief The Individual copy constructor.
       */
      Individual(const Individual& ind);

      /**
       * @brief Destroy an Individual.
       */
      virtual ~Individual();

    public: // Methods

      /**
       * @brief The Individual copy operator.
       *
       * @return A ref toward the assigned Individual.
       * Make a copy of each atribute of the Individual.
       */
      Individual& operator= (const Individual & ind);

      /**
       * @brief Set the id of the Individual.
       *
       * @param id The id of the Individual as a string.
       */
      void setId(const std::string& id);

      /**
       * @brief Get the id of the Individual.
       *
       * @return The id of the Individual as a string.
       */
      const std::string& getId() const { return id_; }

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
      unsigned short getSex() const { return sex_; }

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
      const Date& getDate() const throw (NullPointerException);

      /**
       * @brief Tell if this Individual has a date.
       */
      bool hasDate() const;

      /**
       * @brief Set the coodinates of the Individual.
       *
       * @param coord A Point2D object.
       */
      void setCoord(const Point2D<double> & coord);

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
       * @return A pointer toward a Point2D object if the Individual has
       * coordinates. Otherwise throw a NullPointerException.
       */
      const Point2D<double>& getCoord() const throw (NullPointerException);

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
      void setX(const double x) throw (NullPointerException);

      /**
       * @brief Set the Y coordinate of th Individual.
       *
       * @param y The Y coordinate as a double.
       *
       * Set the Y coordinate if the Individual has coordinates.
       * Otherwise throw a NullPointerException.
       */
      void setY(const double y) throw (NullPointerException);

      /**
       * @brief Get the X coordinate of the Individual.
       *
       * @return The X coordinate as a double if the Individual has coordinates.
       * Otherwise throw a NullPointerException.
       */
      double getX() const throw (NullPointerException);

      /**
       * @brief Get the Y coordinate of the Individual.
       *
       * @return The Y coordinate as a double if the Individual has coordinates.
       * Otherwise throw a NullPointerException.
       */
      double getY() const throw (NullPointerException);

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
      void addSequence(unsigned int sequence_key, const Sequence& sequence)
        throw (Exception);

      /**
       * @brief Get a sequence by its name.
       *
       * @param sequence_name The name of the sequence.
       * @return A reference to the sequence.
       * @throw NullPointerException if there is no sequence container defined.
       * @throw SequenceNotFoundException if sequence_name is not found.
       */
      const Sequence& getSequenceByName(const std::string& sequence_name)
        const throw(Exception);

      /**
       * @brief Get a sequence by its position.
       *
       * @param sequence_position The position of the sequence in the sequence set.
       * @return A reference to the sequence.
       * @throw NullPointerException if there is no sequence container defined.
       * @throw SequenceNotFoundException if sequence_position is not found (i.e. missing data or not used).
       */
      const Sequence& getSequenceAtPosition(const unsigned int sequence_position)
        const throw(Exception);

      /**
       * @brief Delete a sequence.
       *
       * @param sequence_name The name of the sequence.
       * @throw NullPointerException if there is no sequence container defined.
       * @throw SequenceNotFoundException if sequence_name is not found.
       */
      void deleteSequenceByName(const std::string & sequence_name) throw (Exception);

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
      std::vector<std::string> getSequencesNames() const throw (NullPointerException);

      /**
       * @brief Get the sequences' positions.
       *
       * @return All the positions where a sequence is found.
       * @throw NullPointerException if there is no sequence container defined.
       */
      std::vector<unsigned int> getSequencesPositions() const throw (NullPointerException);

      /**
       * @brief Get the position of a sequence.
       *
       * @throw NullPointerException if there is no sequence container defined.
       * @throw SequenceNotFoundException if sequence_name is not found.
       */
      unsigned int getSequencePosition(const std::string& sequence_name) const throw (Exception);

      /**
       * @brief Get the number of sequences.
       */
      unsigned int getNumberOfSequences() const;

      /**
       * @brief Set all the sequences with a MapSequenceContainer.
       */
      void setSequences(const MapSequenceContainer& msc);

      /**
       * @brief Get a reference to the sequence container.
       *
       * @throw NullPointerException if there is no sequence container defined.
       */
      const OrderedSequenceContainer& getSequences() const throw (NullPointerException);

      /**
       * @brief Set a genotype.
       *
       * @param genotype The MultilocusGenotype which will be copied.
       */
      void setGenotype(const MultilocusGenotype& genotype);

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
      const MultilocusGenotype& getGenotype() const throw (NullPointerException);

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
      void setMonolocusGenotype(unsigned int locus_position, const MonolocusGenotype& monogen)
        throw (Exception);

      /**
       * @brief Set a MonolocusGenotype.
       *
       * @throw NullPointerException if there is no genotype defined.
       * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
       * @throw Exception if there is no key in allele_keys.
       */
      void setMonolocusGenotypeByAlleleKey(unsigned int locus_position, const std::vector<unsigned int> allele_keys)
        throw (Exception);

      /**
       * @brief Set a MonolocusGenotype.
       *
       * @throw NullPointerException if there is no genotype defined.
       * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
       * @throw AlleleNotFoundException if at least one the id is not found in the LocusInfo.
       */
      void setMonolocusGenotypeByAlleleId(unsigned int locus_position, const std::vector<std::string> allele_id, const LocusInfo & locus_info)
        throw (Exception);

      /**
       * @brief Get a MonolocusGenotype.
       *
       * @throw NullPointerException if there is no genotype defined.
       * @throw IndexOutOfBoundsException if locus_position excedes the number of loci.
       */
      const MonolocusGenotype& getMonolocusGenotype(unsigned int locus_position) throw (Exception);

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
  };

} //end of namespace bpp;

#endif // _INDIVIDUAL_H_

