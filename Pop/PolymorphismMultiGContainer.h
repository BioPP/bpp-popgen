/*
 * File PolymorphismMultiGContainer.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday July 23 2004
 */

// Secured inclusion of header's file
#ifndef _POLYMORPHYSMMULTIGCONTAINER_H_
#define _POLYMORPHYSMMULTIGCONTAINER_H_

// From C library
#include <cmath>

// From STL
#include <vector>
#include <map>
#include <set>
#include <algorithm>
using namespace std;

// From Utils
#include <Utils/Clonable.h>
#include <Utils/Exceptions.h>
#include <Utils/MapTools.h>

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
		unsigned int getGroupId(unsigned int position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the groups' ids.
		 */
		set<unsigned int> getAllGroupsIds() const;

		/**
		 * @brief Tell if a group exists.
		 */
		bool groupExists(unsigned int group) const;

		/**
		 * @brief Get the number of groups.
		 */
		unsigned int getNumberOfGroups() const;

		/**
		 * @brief Get group size.
		 */
		unsigned int getGroupSize(unsigned int group) const;

		/**
		 * @brief Get the size of a group for a given locus.
		 */
		unsigned int getLocusGroupSize(unsigned int group, unsigned int locus_position) const;
		
		/**
		 * @brief Get the number of MultilocusGenotype.
		 */
		unsigned int size() const;
		
		/**
		 * @brief Clear the container.
		 */
		void clear();

		/**
		 * @brief Get all the alleles' id at one locus.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		vector<unsigned int> getAllelesIdsForAllGroups(unsigned int locus_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the alleles' id at one locus for a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		vector<unsigned int> getAllelesIdsForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Count the number of allele (gametes) at a locus for a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		unsigned int countGametesForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Count the different alleles at one locus.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		map<unsigned int, unsigned int> getAllelesMapForAllGroups(unsigned int locus_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Count the different alleles at one locus for one group.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw GroupNotFoundException if group is not found in the container.
		 */
		map<unsigned int, unsigned int> getAllelesMapForOneGroup(unsigned int locus_position, unsigned int group) const throw (Exception);

		/**
		 * @brief Get a map of allele count for a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		map<unsigned int, unsigned int> getAllelesMapForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the alleles frequencies at one locus for all groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		map<unsigned int, double> getAllelesFrqForAllGroups(unsigned int locus_position) const throw (Exception);

		/**
		 * @brief Get the alleles frequencies at one locus for one group.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw GroupNotFoundException if group is not found in the container.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		map<unsigned int, double> getAllelesFrqForOneGroup(unsigned int locus_position, unsigned int group) const throw (Exception);
		
		/**
		 * @brief Get the alleles frequencies at one locus for a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		map<unsigned int, double> getAllelesFrqForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (Exception);

		/**
		 * @brief Count the number of non-missing data at a given locus.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		unsigned int countNonMissingForAllGroups(unsigned int locus_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Count the number of non-missing data at a given locus for one group.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw GroupNotFoundException if group is not found in the container.
		 */
		unsigned int countNonMissingForOneGroup(unsigned int locus_position, unsigned int group) const throw (Exception);

		/**
		 * @brief Count the number of non-missing data at a given locus for a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		unsigned int countNonMissingForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Count the number of bi-allelic MonolocusGenotype at a given locus.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		unsigned int countBiAllelicForAllGroups(unsigned int locus_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Count the number of bi-allelic MonolocusGenotype at a given locus for one group.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw GroupNotFoundException if group is not found in the container.
		 */
		unsigned int countBiAllelicForOneGroup(unsigned int locus_position, unsigned int group) const throw (Exception);

		/**
		 * @brief Counr the number of bi-allelic MonolocusGenotype at a given locus for a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		unsigned int countBiAllelicForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Count how many times each allele is found in an heterozygous MonolocusGenotype.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		map<unsigned int, unsigned int> countHeterozygousForAllGroups(unsigned int locus_position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Count how many times each allele is found in an heterozygous MonolocusGenotype in one group.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw GroupNotFoundException if group is not found in the container.
		 */
		map<unsigned int, unsigned int> countHeterozygousForOneGroup(unsigned int locus_position, unsigned int group) const throw (Exception);

		/**
		 * @brief Count how many times each allele is found in an heterozygous MonolocusGenotype in a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		map<unsigned int, unsigned int> countHeterozygousForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the heterozygous frequencies for each allele at a locus.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		map<unsigned int, double> getHeterozygousFrqForAllGroups(unsigned int locus_position) const throw (Exception);

		/**
		 * @brief Get the heterozygous frequencies for each allele at a locus for one group.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw GroupNotFoundException if group is not found in the container.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		map<unsigned int, double> getHeterozygousFrqForOneGroup(unsigned int locus_position, unsigned int group) const throw (Exception);

		/**
		 * @brief Get the heterozygous frequencies for each allele at a locus in a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		map<unsigned int, double> getHeterozygousFrqForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (Exception);

		/**
		 * @brief Compute the observed heterozygosity for one locus.
		 *
		 * This is the mean value of the getHeterozygousFrqForGroups map.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		double getHobsForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (Exception);

		/**
		 * @brief Compute the expected heterozygosity for one locus.
		 *
		 * Nei 1977
		 * @f[
		 * H_{exp}=1-\sum_{i=1}^{n}x_i^2
		 * @f]
		 * where @f$x_i@f$ is the frequency of the i<sup>th</sup> allele and @f$n@f$ the number of alleles.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		double getHexpForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (Exception);

		/**
		 * @brief Compute the expected non biased heterozygosity for one locus.
		 *
		 * Nei 1978
		 * @f[
		 * H_{nb}=\frac{2n}{2n-1}\left(1-\sum_{i=1}^{n}x_i^2\right)=\frac{2n}{2n-1}H_{exp}
		 * @f]
		 * where @f$x_i@f$ is the frequency of the i<sup>th</sup> allele and @f$n@f$ the number of alleles.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		double getHnbForGroups(unsigned int locus_position, const set<unsigned int> & groups) const throw (Exception);

		/**
		 * @brief Compute the Nei distance between two groups at one locus.
		 *
		 * Nei 1972
		 * @f[
		 * \hat{D}_1=-\ln \left[\frac{\displaystyle\sum_{i=1}^{n}\left(x_i\times y_i\right)}
		 * {\sqrt{\displaystyle\sum_{i=1}^{n}x_i^2\times \displaystyle\sum_{i=1}^{n}y_i^2}}\right]
		 * @f]
		 * where @f$x_i@f$ and @f$y_i@f$ are respectively the i<sup>th</sup> allele's frequency of the first and second group
		 * and @f$n@f$ the total number of alleles of both groups.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		double getDnei72(unsigned int locus_position, unsigned int grp1, unsigned int grp2) const throw (Exception);

		/**
		 * @brief Compute the Nei unbiased distance between two groups at one locus.
		 *
		 * Nei 1978
		 * @f[
		 * \hat{D}=-\ln \left[\frac{\displaystyle\sum_{i=1}^{n}\left(x_i\times y_i\right)}
		 * {\sqrt{\frac{2n_XJ_X-1}{2n_X-1}\times\frac{2n_YJ_Y-1}{2n_YJ_Y}}}
		 * \right]
		 * @f]
		 * where @f$x_i@f$ and @f$y_i@f$ are respectively the i<sup>th</sup> allele's frequency of the first and second group,
		 * @f$n@f$ the total number of alleles of both groups, @f$n_X@f$ and @f$n_Y@f$ the number of alleles in the first and second group
		 * and
		 * @f[
		 * J_X=\sum_{i=1}^{n}x_i^2
		 * \qquad\textrm{and}\qquad
		 * J_Y=\sum_{i=1}^{n}y_i^2
		 * @f]
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		 double getDnei78(unsigned int locus_position, unsigned int grp1, unsigned int grp2) const throw (Exception);

	protected:
		vector<MultilocusGenotype *> _multilocusGenotypes;
		vector<unsigned int> _groups;
};

#endif // _POLYMORPHYSMMULTIGCONTAINER_H_
