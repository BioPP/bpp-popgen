/*
 * File AnalyzedLoci.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 07 2004
 */

// Secured inclusion of header's file
#ifndef _ANALYZEDLOCI_H_
#define _ANALYZEDLOCI_H_

// From STL
#include <vector>
using namespace std;

//From Utils
#include <Utils/Exceptions.h>

// From local
#include "LocusInfo.h"
#include "GeneralExceptions.h"

/**
 * @brief The AnalyzedLoci class.
 *
 * This is a LocusInfo container.
 * Its instanciation requires a number of locus wich is fixed
 * and can't be modified.
 */
class AnalyzedLoci {
	public: // Constructors and Destructor
		/**
		 * @brief Build a void AnalyzedLoci with a specific number of loci.
		 */
		AnalyzedLoci(unsigned int number_of_loci);

		/**
		 * @brief Destroy the AnalyzedLoci.
		 */
		~AnalyzedLoci();

	public: // Other methodes
		/**
		 * @brief Add a LocusInfo.
		 */
		void setLocusInfo(unsigned int locus_index, const LocusInfo & locus)
			throw (IndexOutOfBoundsException);

		/**
		 * @brief Get a LocusInfo by name.
		 */
		const LocusInfo * getLocusInfoByName(const string & locus_name) const
			throw (BadIdentifierException);

		/**
		 * @brief Get a LocusInfo by index.
		 */
		const LocusInfo * getLocusInfoByIndex(unsigned int locus_index) const
			throw (Exception);

		/**
		 * @brief Add an AlleleInfo to a LocusInfo by LocusInfo name.
		 */
		void addAlleleInfoByLocusName(const string & locus_name,
				const AlleleInfo & allele)
			throw (Exception);

		/**
		 * @brief Add an AlleleInfo to a LocusInfo by index.
		 */
		void addAlleleInfoByLocusIndex(unsigned int locus_index,
				const AlleleInfo & allele)
			throw (Exception);

		/**
		 * @brief Get the number of loci.
		 */
		unsigned int getNumberOfLoci() const;
		
		/**
		 * @brief Get the number of alleles at each locus.
		 */
		vector<unsigned int> getNumberOfAlleles() const;

		/**
		 * @brief Get the ploidy of a locus by name.
		 */
		unsigned int getPloidyByLocusName(const string & locus_name) const
			throw (BadIdentifierException);

		/**
		 * @brief Get the ploidy of a locus by index.
		 */
		unsigned int getPloidyByLocusIndex(unsigned int locus_index) const
			throw (IndexOutOfBoundsException);
	protected:
		vector<LocusInfo *> _loci;
};
#endif // _ANALYZEDLOCI_H_
