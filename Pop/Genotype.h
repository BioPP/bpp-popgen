/*
 * File Genotype.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday May 28 2004
 */

// From STL
#include <vector>
using namespace std;

/**
 * @brief The Genotype class.
 *
 * This is a MonolocusGenotype containor.
 */
class Genotype {
	public: // Constructors and Destructor
		/**
		 * @brief Build a void Genotype.
		 */
		Genotype();

		/**
		 * @brief Build a Genotype linked to an AnalyzedLoci object.
		 */
		Gentype(AnalyzedLoci * analyzed_loci);

		/**
		 * @brief Destroy a Genotype.
		 */
		~Genotype();

	protected:
		AnalyzedLoci * _analyzedLoci;
		vector<MonolocusGenotype *> _loci;
};
