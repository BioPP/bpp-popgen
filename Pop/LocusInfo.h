/*
 * File LocusInfo.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday June 07 2004
 */

// Secured inclusion of header's file
#ifndef _LOCUSINFO_H_
#define _LOCUSINFO_H_

//From STL
#include <map>
#include <vector>
using namespace std;

// From local Poplib
#include "AlleleInfo.h"
#include "GeneralExceptions.h"

// From Utils
#include <Utils/Exceptions.h>

/**
 * @brief The LocusInfo class.
 *
 * This is an AlleleInfo container with additionnal data like a name,
 * the ploidy and some comments.
 */
class LocusInfo {
	public: // Constantes
		static unsigned int HAPLODIPLOID;
		static unsigned int HAPLOID;
		static unsigned int DIPLOID;
		
	public: // Constructors and destructor
		/**
		 * @brief Build a new LocusInfo object.
		 *
		 * @param name The name of the locus.
		 * @param ploidy The ploidy of the locus.
		 */
		LocusInfo(const string &name, const unsigned int ploidy = DIPLOID);

		/**
		 * @brief Destroy the LocusInfo.
		 */
		virtual ~LocusInfo();
		
	public: // Methodes
		/**
		 * @brief Get the name of the locus.
		 */
		string getName() const;

		/**
		 * @brief Get the ploidy of the locus.
		 *
		 * @return The ploidy as an unsigned integer.
		 */
		unsigned int getPloidy() const;
		
		/**
		 * @brief Add an AlleleInfo to the LocusInfo.
		 */
		void addAlleleInfo(const AlleleInfo &allele) throw (Exception);

		/**
		 * @brief Retrieve an AlleleInfo object of the LocusInfo.
		 */
		AlleleInfo * getAlleleInfoById(unsigned int id) const throw (BadIdentifierException);

		/**
		 * @brief Retrieve an AlleleInfo object of the LocusInfo.
		 */
		AlleleInfo * getAlleleInfoByKey(unsigned int key) const throw (BadIntegerException);

		/**
		 * @brief Get the position of an AlleleInfo.
		 */
		unsigned int getAlleleInfoKey(unsigned int id) const
			throw (BadIdentifierException);
		
		/**
		 * @brief Get the number of alleles at this locus.
		 */
		unsigned int getNumberOfAlleles() const;

		/**
		 * @brief Get the identity number of each allele in the locus.
		 */
		vector<unsigned int> getAlleleInfosKeys() const throw (Exception);

		/**
		 * @brief Remove all alleles from the locus.
		 */
		void clear();

	protected:
		string _name;
		unsigned int _ploidy;
		map<unsigned int, AlleleInfo *> _alleles;
};

#endif // _LOCUSINFO_H_
