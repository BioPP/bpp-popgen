/*
 * File MonolocusGenotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 22 2004
 */

// Secured inclusion of header's file
#ifndef _MONOLOCUSGENOTYPE_H_
#define _MONOLOCUSGENOTYPE_H_

//From Utils
#include <Utils/Clonable.h>

/**
 * @brief The MonolocusGenotype virtual class.
 *
 * A MonolocusGenotype containes the Alleles' keys defined in a Locus object.
 * This keys are returned as unsigned integers.
 * This class is an interface for all monolocus genotypes.
 */
class MonolocusGenotype : public Clonable {
	public: // Constructors and Destructor
		/**
		 * @brief Destroy a MonolocusGenotype.
		 */
		virtual ~MonolocusGenotype();
		
	public: // Methodes
		/**
		 * @brief Get the first allele's key.
		 */
		virtual unsigned int getAlleleKey() const = 0;

		/**
		 * @brief Get the ploidy of the locus.
		 */
		virtual unsigned int getPloidy() const = 0;

		/**
		 * @brief Get the maximum number of allele that con be stored.
		 */
		virtual unsigned int getSize() const  = 0;
};
#endif // _MONOLOCUSGENOTYPE_H_
