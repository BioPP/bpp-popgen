/*
 * File MonolocusGenotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday May 27 2004
 */

// From STL
#include <vector>
using namespace std;

// From Utils
#include <Utils/Exceptions.h>

/**
 * @brief The MonolocusGenotype class.
 *
 * A MonolocusGenotype containes the Alleles' keys defined in a Locus object.
 * This keys are unsigned integers.
 *
 * 0 (zero) is considered as missing data.
 */
class MonolocusGenotype {
	public: // Constructors and Destructor
		/**
		 * @brief Build a void MonolocusGenotype.
		 */
		MonolocusGenotype();

		/**
		 * @brief Destroy a MonolocusGenotype.
		 */
		~MonolocusGenotype();
		
	public: // Methodes
		/**
		 * @brief Add an Allele key to the MonolocusGenotype.
		 */
		void addKey(unsigned int key);

		/**
		 * @brief Get a special key.
		 */
		unsigned int getKey(unsigned int index) throw (IndexOutOfBoundsException);

		/**
		 * @brief Get all the valid Allele keys.
		 */
		vector<unsigned int> getKeys();
		
		/**
		 * @brief Get the number of alleles at this locus in an individual.
		 */
		unsigned int getNumberOfData();

		/**
		 * @brief Test the presence of missing data.
		 */
		bool hasMissingData();

	protected:
		vector<unsigned int> _allelekeys;
};
