/*
 * File MultilocusGenotypeStatistics.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday August 03 2004
 *
 * Copyright (C) 2004 Sylvain Gaillard and the
 *                    PopLib Development Core Team
 *
 * This file is part of PopLib.
 *
 * PopLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

// Secured inclusion of header's file
#ifndef _MULTILOCUSGENOTYPESTATISTICS_H_
#define _MULTILOCUSGENOTYPESTATISTICS_H_

// From C library
#include <cmath>

// From STL
#include <vector>
#include <map>
#include <set>
#include <algorithm>
using namespace std;

// From Utils
#include <Utils/Exceptions.h>
#include <Utils/MapTools.h>

// From poplib
#include "PolymorphismMultiGContainer.h"
#include "MultilocusGenotype.h"
#include "GeneralExceptions.h"

/**
 * @brief The MultilocusGenotypeStatistics class
 *
 * This class is a set of static method for PolymorphismMultiGContainer.
 */
class MultilocusGenotypeStatistics {
	public:
		struct VarComp {
			double a;
			double b;
			double c;
		};

		struct Fstats {
			double Fit;
			double Fst;
			double Fis;
		};

		/**
		 * @brief Get the alleles' id at one locus for a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		static vector<unsigned int> getAllelesIdsForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (IndexOutOfBoundsException);

		/**
		 * @brief Count the number of allele (gametes) at a locus for a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		static unsigned int countGametesForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (IndexOutOfBoundsException);

		/**
		 * @brief Get a map of allele count for a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		static map<unsigned int, unsigned int> getAllelesMapForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the alleles frequencies at one locus for a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		static map<unsigned int, double> getAllelesFrqForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);

		/**
		 * @brief Count the number of non-missing data at a given locus for a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		static unsigned int countNonMissingForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (IndexOutOfBoundsException);

		/**
		 * @brief Counr the number of bi-allelic MonolocusGenotype at a given locus for a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		static unsigned int countBiAllelicForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Count how many times each allele is found in an heterozygous MonolocusGenotype in a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 */
		static map<unsigned int, unsigned int> countHeterozygousForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the heterozygous frequencies for each allele at a locus in a set of groups.
		 *
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		static map<unsigned int, double> getHeterozygousFrqForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);

		/**
		 * @brief Compute the observed heterozygosity for one locus.
		 *
		 * This is the mean value of the getHeterozygousFrqForGroups map.
		 * @throw IndexOutOfBoundsException if locus_position excedes the number of loci of one MultilocusGenotype.
		 * @throw ZeroDivisionException if the number of considered alleles = 0.
		 */
		static double getHobsForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);

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
		static double getHexpForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);

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
		static double getHnbForGroups(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);

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
		static double getDnei72(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, unsigned int grp1, unsigned int grp2) throw (Exception);

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
		static double getDnei78(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, unsigned int grp1, unsigned int grp2) throw (Exception);

		/**
		 * @brief Compute the three F statistics of Weir and Cockerham.
		 */
		static map<unsigned int, Fstats>  getAllelesFstats(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);
		
		/**
		 * @brief Compute the Weir and Cockerham Fit on a set of groups.
		 */
		static map<unsigned int, double> getAllelesFit(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);

		/**
		 * @brief Compute the Weir and Cockerham Fst on a set of groups.
		 */
		static map<unsigned int, double> getAllelesFst(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);

		/**
		 * @brief Compute the Weir and Cockerham Fis on a set of groups.
		 */
		static map<unsigned int, double> getAllelesFis(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);
		
		/**
		 * @brief Get the variance components a, b and c (Weir and Cockerham, 1983).
		 */
		static map<unsigned int, VarComp> getVarianceComponents(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (ZeroDivisionException);
};

#endif // _MULTILOCUSGENOTYPESTATISTICS_H_
