/*
 * File MonolocusGenotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Monday July 05 2004
 */

// Secured inclusion of header's file
#ifndef _MONOLOCUSGENOTYPE_H_
#define _MONOLOCUSGENOTYPE_H_

// From STL
#include <vector>
using namespace std;

// From Utils
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
		 * @brief Get the alleles' index.
		 *
		 * The alleles' index are the position of the AlleleInfo in a LocusInfo object.
		 * If no LocusInfo is used, the index are just numbers to identify the alleles.
		 *
		 * @return A vector of unsigned int.
		 *
		 * The size of the vector corresponds to the number of alleles at this locus.
		 */
		virtual vector<unsigned int> getAlleleIndex() const = 0;
};
#endif // _MONOLOCUSGENOTYPE_H_
