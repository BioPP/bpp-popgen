/*
 * File AnalyzedLoci.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 01 2004
 */

// From STL
#include <vector>
using namespace std;

//From Utils
#include <Utils/Exceptions.h>

// From local
#include "LocusInfo.h"
#include "LocusInfoContainerExceptions.h"

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
		const LocusInfo * getLocusInfo(const string & locus_name)
			throw (LocusInfoNotFoundException);

		/**
		 * @brief Get a LocusInfo by index.
		 */
		const LocusInfo * getLocusInfo(unsigned int locus_index)
			throw (NullPointerException);

		/**
		 * @brief Add an AlleleInfo to a LocusInfo by LocusInfo name.
		 */
		void addAlleleInfo(const string & locus_name, const AlleleInfo & allele)
			throw (Exception);

		/**
		 * @brief Add an AlleleInfo to a LocusInfo by index.
		 */
		void addAlleleInfo(unsigned int locus_index, const AlleleInfo & allele)
			throw (Exception);
		
	protected:
		vector<LocusInfo *> _loci;
};
