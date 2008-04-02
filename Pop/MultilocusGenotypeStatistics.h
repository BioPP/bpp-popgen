//
// File MultilocusGenotypeStatistics.h
// Authors : Sylvain Gaillard
//           Khalid Belkhir
// Last modification : Wednesday August 04 2004
//

/*
   Copyright or © or Copr. CNRS, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for population genetics analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
   */

#ifndef _MULTILOCUSGENOTYPESTATISTICS_H_
#define _MULTILOCUSGENOTYPESTATISTICS_H_

// From STL
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

using namespace std;

// From Utils
#include <Utils/Exceptions.h>
#include <Utils/MapTools.h>

// From SeqLib
#include <Seq/DistanceMatrix.h>

// From popgenlib
#include "PolymorphismMultiGContainer.h"
#include "MultilocusGenotype.h"
#include "GeneralExceptions.h"

namespace bpp
{

  /**
   * @brief The MultilocusGenotypeStatistics class
   *
   * This class is a set of static method for PolymorphismMultiGContainer.
   */
  class MultilocusGenotypeStatistics
  {
    public:
      struct VarComp
      {
        double a;
        double b;
        double c;
      };

      struct Fstats
      {
        double Fit;
        double Fst;
        double Fis;
      };

      struct PermResults{
        double Statistic;
        double Percent_sup;
        double Percent_inf;
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
      static double getDnei72(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions, unsigned int grp1, unsigned int grp2) throw (Exception);

      /**
       * @brief Compute the Nei unbiased distance between two groups at a given number of loci.
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
      static double getDnei78(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions, unsigned int grp1, unsigned int grp2) throw (Exception);

      /**
       * @brief Compute the three F statistics of Weir and Cockerham for each allele of a given locus.
       */
      static map<unsigned int, Fstats>  getAllelesFstats(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);

      /**
       * @brief Compute the Weir and Cockerham Fit on a set of groups for each allele of a given locus.
       */
      static map<unsigned int, double> getAllelesFit(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);

      /**
       * @brief Compute the Weir and Cockerham @f$\theta@f$ on a set of groups for each allele of a given locus.
       */
      static map<unsigned int, double> getAllelesFst(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);

      /**
       * @brief Compute the Weir and Cockerham Fis on a set of groups for each allele of a given locus.
       */
      static map<unsigned int, double> getAllelesFis(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (Exception);

      /**
       * @brief Get the variance components a, b and c (Weir and Cockerham, 1983).
       */
      static map<unsigned int, VarComp> getVarianceComponents(const PolymorphismMultiGContainer & pmgc, unsigned int locus_position, const set<unsigned int> & groups) throw (ZeroDivisionException);

      /**
       * @brief Compute the Weir and Cockerham @f$\theta{wc}@f$ on a set of groups for a given set of loci.
       * The variance componenets for each allele are calculated and then combined over loci using Weir and Cockerham weighting.
       */
      static double getWCMultilocusFst(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions, const set<unsigned int> & groups) throw (Exception);

      /**
       * @brief Compute the Weir and Cockerham Fis on a set of groups for a given set of loci.
       * The variance componenets for each allele are calculated and then combined over loci using Weir and Cockerham weighting.
       */
      static double getWCMultilocusFis(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions, const set<unsigned int> & groups) throw (Exception);

      /**
       * @brief Compute the Weir and Cockerham @f$\theta_{wc}@f$ on a set of groups for a given set of loci and make a permutation test.
       * Multilocus @f$\theta@f$ is calculated as in getWCMultilocusFst on the original data set and on nb_perm data sets obtained after
       * a permutation of individuals between the different groups.
       * Return values are theta, % of values > theta and % of values < theta.
       */
      static PermResults getWCMultilocusFstAndPerm(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions,set<unsigned int> groups, int nb_perm) throw (Exception);

      /**
       * @brief Compute the Weir and Cockerham Fis on a set of groups for a given set of loci and make a permutation test.
       * Multilocus Fis is calculated as in getWCMultilocusFis on the original data set and on nb_perm data sets obtained after
       * a permutation of alleles between individual of each group.
       * Return values are Fis, % of values > Fis and % of values < Fis.
       */
      static PermResults getWCMultilocusFisAndPerm(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions,set<unsigned int> groups, int nb_perm) throw (Exception);


      /**
       * @brief Compute the @f$\theta_{RH}@f$ on a set of groups for a given set of loci.
       * The variance componenets for each allele are calculated and then combined over loci using RH weighting with alleles frequency.
       */
      static double getRHMultilocusFst(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions, const set<unsigned int> & groups) throw (Exception);

      /**
       * @brief Compute pairwise distances on a set of groups for a given set of loci.
       * distance is either Nei72, Nei78, Fst W&C or Fst Robertson & Hill, Nm,
       * D=-ln(1-Fst) of Reynolds et al. 1983, Rousset 1997 Fst/(1-Fst)
       */
      static DistanceMatrix * getDistanceMatrix(const PolymorphismMultiGContainer & pmgc, vector<unsigned int> locus_positions, const set<unsigned int> & groups, string distance_methode) throw (Exception);

  };

} //end of namespace bpp;

#endif // _MULTILOCUSGENOTYPESTATISTICS_H_

