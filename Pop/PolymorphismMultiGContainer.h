/*
 * File PolymorphismMultiGContainer.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday July 20 2004
 */

// Secured inclusion of header's file
#ifndef _POLYMORPHYSMMULTIGCONTAINER_H_
#define _POLYMORPHYSMMULTIGCONTAINER_H_

// From STL
#include <vector>
#include <map>
using namespace std;

// From Utils
#include <Utils/Clonable.h>
#include <Utils/Exceptions.h>

// From poplib
#include "MultilocusGenotype.h"
#include "GeneralExceptions.h"

/**
 * @brief The PolymorphismMultiGContainer class
 *
 * This class is a container of MultilocusGenotype.
 */
class PolymorphismMultiGContainer : public Clonable {
	public: // Constructors and destructor
		/**
		 * @brief Build a new PolymorphismMultilocusGenotypeContainer.
		 */
		PolymorphismMultiGContainer();

		/**
		 * @brief The copy constructor.
		 */
		PolymorphismMultiGContainer(const PolymorphismMultiGContainer & pmgc);

		/**
		 * @brief Destroy a PolymorphismMultilocusGenotypeContainer.
		 */
		~PolymorphismMultiGContainer();

	public:
		/**
		 * @brief The assignation operator=.
		 */
		PolymorphismMultiGContainer & operator= (const PolymorphismMultiGContainer & pmgc);

		/**
		 * @name The clonable interface
		 * @{
		 */
		Clonable * clone() const;
		/** @} */

		/**
		 * @brief Add a MultilocusGenotype to the container.
		 */
		void addMultilocusGenotype(const MultilocusGenotype & mg, unsigned int group);

		/**
		 * @brief Get a MultilocusGenotype at a position.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		const MultilocusGenotype * getMultilocusGenotype(unsigned int position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Remove a MultilocusGenotype.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		MultilocusGenotype * removeMultilocusGenotype(unsigned int position) throw (IndexOutOfBoundsException);

		/**
		 * @brief Delete a MultilocusGenotype.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		void deleteMultilocusGenotype(unsigned int position) throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the Group of a MultilocusGenotype.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		unsigned int getGroup(unsigned int position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Tell if a group exists.
		 */
		bool groupExists(unsigned int group) const;
		
		/**
		 * @brief Get the number of MultilocusGenotype.
		 */
		unsigned int size() const;
		
		/**
		 * @brief Clear the container.
		 */
		void clear();

		/**
		 * @brief Count the different alleles at one locus.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		map<unsigned int, unsigned int> countAlleles(unsigned int locus_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Count the different alleles at one locus for one group.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw GroupNotFoundException if group is not found in the container.
		 */
		map<unsigned int, unsigned int> countAlleles(unsigned int locus_position, unsigned int group) const throw (Exception);

		/**
		 * @brief Count the number of non-missing data at a given locus.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		unsigned int countNonMissing(unsigned int locus_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Count the number of non-missing data at a given locus for one group.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw GroupNotFoundException if group is not found in the container.
		 */
		unsigned int countNonMissing(unsigned int locus_position, unsigned int group) const throw (Exception);

		/**
		 * @brief Count the number of bi-allelic MonolocusGenotype at a given locus.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		unsigned int countBiAllelic(unsigned int locus_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Count the number of bi-allelic MonolocusGenotype at a given locus for one group.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw GroupNotFoundException if group is not found in the container.
		 */
		unsigned int countBiAllelic(unsigned int locus_position, unsigned int group) const throw (Exception);
		
		/**
		 * @brief Count how many times each allele is found in an heterozygous MonolocusGenotype.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		map<unsigned int, unsigned int> countHeterozygous(unsigned int locus_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Count how many times each allele is found in an heterozygous MonolocusGenotype in one group.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw GroupNotFoundException if group is not found in the container.
		 */
		map<unsigned int, unsigned int> countHeterozygous(unsigned int locus_position, unsigned int group) const throw (Exception);

	protected:
		vector<MultilocusGenotype *> _multilocusGenotypes;
		vector<unsigned int> _groups;
};

#endif // _POLYMORPHYSMMULTIGCONTAINER_H_
