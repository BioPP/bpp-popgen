/*
 * File LocusInfo.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday June 18 2004
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
		 *
		 * @throw BadIdentifierException if the AlleleInfo's id already exists.
		 */
		void addAlleleInfo(const AlleleInfo &allele)
			throw (BadIdentifierException);

		/**
		 * @brief Retrieve an AlleleInfo object of the LocusInfo.
		 *
		 * @throw AlleleNotFoundException if the id is not found.
		 */
		AlleleInfo * getAlleleInfoById(unsigned int id) const
			throw (AlleleNotFoundException);

		/**
		 * @brief Retrieve an AlleleInfo object of the LocusInfo.
		 *
		 * @throw BadIntegerException if the key doesn't match any referenced key.
		 */
		AlleleInfo * getAlleleInfoByKey(unsigned int key) const
			throw (BadIntegerException);

		/**
		 * @brief Get the position of an AlleleInfo.
		 *
		 * @throw AlleleNotFoundException if the AlleleInfo's id is not found.
		 */
		unsigned int getAlleleInfoKey(unsigned int id) const
			throw (AlleleNotFoundException);
		
		/**
		 * @brief Get the number of alleles at this locus.
		 */
		unsigned int getNumberOfAlleles() const;

		/**
		 * @brief Get the identity number of each allele in the locus.
		 *
		 * @throw Exception if there is no AlleleInfo in the LocusInfo.
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
