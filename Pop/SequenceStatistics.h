/*
 * File SequenceStatistics.h
 * Author : Eric Bazin <bazin@univ-montp2.fr>
 *          Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Friday August 06 2004
*/
/*
Copyright or � or Copr. CNRS, (November 17, 2004)


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



//Secured inclusion of header's file
#ifndef _SEQUENCESTATISTICS_H_
#define _SEQUENCESTATISTICS_H_

//From the SeqLib library
#include <Seq/SiteIterator.h>
#include <Seq/SiteContainer.h>
#include <Seq/SymbolListTools.h>
#include <Seq/CodonAlphabet.h>
#include <Seq/GeneticCode.h>
#include <Seq/SiteContainerTools.h>

//From the PolyLib library
#include "PolymorphismSequenceContainer.h"

using namespace std;

class SequenceStatistics
{
	public:
		// Class destructor:
		~SequenceStatistics();

		/*******************************************************************************/
	public:

		/**
		 * @brief Compute the number of polymorphic site in an alignment
		 *
		 * Gaps are consider as mutations so if you want number of
		 * polymorphic site, you have to give a NonGapSiteIterator
		 *
		 * @param psc a PolymorphismSequenceContainer
		 * @param gapflag flag set by default to true if you don't want to
		 * take gap into account
		 */
		static unsigned int polymorphicSiteNumber(const PolymorphismSequenceContainer & psc, bool gapflag = true);



        /**
         *@brief Compute the number of parsimony informative sites in an alignment
         *@param a PolymorphicSequenceContainer, a boolean
         */
        static unsigned int parsimonyInformativeSiteNumber(const PolymorphismSequenceContainer & psc, bool gapflag =true);


		/**
		 * @brief Count the number of singleton nucleotides in an alignment.
		 *
		 * @param psc a PolymorphismSequenceContainer
		 * @param gapflag flag set by default to true if you don't want to
		 * take gap into account
		 */
		static unsigned int countSingleton(const PolymorphismSequenceContainer & psc, bool gapflag = true);

		/**
		 * @brief Count the total number of mutations in an alignment.
		 *
		 * This count is assumed to be under an infinite site model.
		 * @param psc a PolymorphismSequenceContainer
		 * @param gapflag flag set by default to true if you don't want to
		 * take gap into account
		 */
		static unsigned int totNumberMutations(const PolymorphismSequenceContainer & psc, bool gapflag = true);

        /**
		 * @brief Count the total number of mutations in external branchs.
                 * This is counted as the number of distinct singleton nucleotide  in the ingroup
                 * that are not shared with the outgroup
		 * @param  requires an ingroup and an outgroup
                 * @param gapflag flag set by default to true if you don't want to
		 * take gap into account
                 */
        static unsigned int totMutationsExternalBranchs(const PolymorphismSequenceContainer & ing,
                                                                 const PolymorphismSequenceContainer outg);


         /**
		 * @brief Compute the number of triplet in an alignment
		 *
		 * @param v a SiteContainer
		 * @param gapflag set by default to true if you don't want to take gap into account
		 */
        static unsigned int tripletNumber(const PolymorphismSequenceContainer & psc, bool gapflag = true);



        /**
         *@brief Compute the sum of per site heterozygosity in an alignment
         *@param a SiteContainer, a boolean
         */
        static double heterozygosity(const PolymorphismSequenceContainer & psc, bool gapflag=true);


        /**
         *@brief Compute the sum of per site squared heterozygosity in an alignment
         *@param a SiteContainer, a boolean
         */
        static double squaredHeterozygosity(const PolymorphismSequenceContainer & psc, bool gapflag=true);


        /**
		 * @brief Compute the mean GC content in an alignment
		 *
		 * @param v a PolymorphismSequenceContainer
		 */
		static double gcContent(const PolymorphismSequenceContainer & psc);


		/*
		 * @brief Return the number of GC alleles and the total number of alleles at polymorphic sites only
		 *
		 * G vs C and A vs T polymorphism are not taken into account
		 * @param psc a PolymorphismSequenceContainer
		 * @param stopflag a boolean set by default to true if you don't want
		 * to take stop codons into account
		 */
		static vector<unsigned int> gcPolymorphism(const PolymorphismSequenceContainer & psc, bool stopflag = true);



		/**
		 * @brief Compute diversity estimator Theta of Watterson (1975)
		 *
		 * @f[
		 * \hat{\theta}_S=\frac{S}{a_1}
		 * @f]
		 * where @f$S@f$ is the number of polymorphic sites and @f$a_1@f$ is describe in SequenceStatistics::_getUsefullValues().
		 * @param psc a PolymorphismSequenceContainer
		 * @param gapflag flag set by default to true if you don't want to
		 * take gap into account
		 */
		static double watterson75( const PolymorphismSequenceContainer & psc, bool gapflag = true );


		/**
		 * @brief Compute diversity estimator Theta of Tajima (1983)
		 *
		 * @f[
		 * \hat{\theta}_\pi=1-\sum_{i=1}^{S}\sum_{j=1}^{4}\frac{k_{j,i}\times\left(k_{j,i}-1\right)}
		 * {n_i\times\left(n_i-1\right)} \qquad \textrm{with }k_{j,i}>0
		 * @f]
		 * where @f$k_{j,i}@f$ is the count of the j<sup>th</sup> state at the i<sup>th</sup> site,
		 * @f$n_i@f$ the number of nucleotides and @f$S@f$ the number of polymorphic sites.
		 * @param psc a PolymorphismSequenceContainer
		 * @param gapflag flag set by default to true if you don't want to
		 * take gap into account
		 */
		static double tajima83( const PolymorphismSequenceContainer & psc, bool gapflag = true );


		/**
		 * @brief Return the number of haplotype in the sample.
		 * Depaulis and Veuille (1998)
		 *
		 * @param psc a PolymorphismSequenceContainer
		 * @param gapflag flag set by default to true if you don't want to
		 * take gap into account
		 */
		static unsigned int DVK ( const PolymorphismSequenceContainer & psc, bool gapflag = true );


		/**
		 * @brief Return the haplotype diversity of a sample.
		 * Depaulis and Veuille (1998)
		 *
		 * @param psc a PolymorphismSequenceContainer
		 * @param gapflag flag set by default to true if you don't want to
		 * take gaps into account
		 */
		static double DVH ( const PolymorphismSequenceContainer & psc, bool gapflag = true );


		/**
		 * @brief Return the number of transitions.
		 *
		 * @param psc a PolymorphismSequenceContainer
		 */
		static unsigned int getNumberOfTransitions( const PolymorphismSequenceContainer & psc );


		/**
		 * @brief Return the number of transversions.
		 *
		 * @param psc a PolymorphismSequenceContainer
		 */
		static unsigned int getNumberOfTransversions( const PolymorphismSequenceContainer & psc );


		/**
		 * @brief Return the ratio of transitions/transversions.
		 *
		 * @param psc a PolymorphismSequenceContainer
		 */
		 static double getTransitionsTransversionsRatio( const PolymorphismSequenceContainer & psc );




		/**
		 * @brief Compute the number of codon sites with stop codon
		 * @param v a SiteContainer
		 * @param gapfalg a boolean set by default to true if you don't want to take gaps into account
		 */
		static unsigned int stopCodonSiteNumber(const PolymorphismSequenceContainer & psc, bool gapflag = true);




		/**
		 * @brief Compute the number of polymorphic codon with only one mutated site
		 * @param v a SiteContainer
		 * @param stopflag a boolean set by default to true if you don't want
		 * to take stop codon into account
		 * @param gapflag a boolean set by default to true if you don't want
		 * to take gaps into account
		 */
		static unsigned int monoSitePolymorphicCodonNumber(const PolymorphismSequenceContainer & psc, bool stopflag = true, bool gapflag = true);



		/**
		 * @brief Compute the number of synonymous polymorphic codon sites
		 *
		 * Gaps are automatically excluded
		 * @param v a SiteContainer
		 * @param v a GeneticCode
		 * @param stopflag a boolean set by default to true if you don't want to take stop codons into account
		 */
		static unsigned int synonymousPolymorphicCodonNumber(const PolymorphismSequenceContainer & psc, const GeneticCode & gc,  bool stopflag = true);


		/**
		 * @brief Compute the Watterson(1975) estimator for synonymous positions
		 *
		 * Gaps are automatically excluded
		 * @param psc a PolymorphismSequenceContainer
		 * @param gc a GeneticCode
		 */
		static double watterson75Synonymous(const PolymorphismSequenceContainer & psc, const GeneticCode & gc);

		/**
		 * @brief Compute the Watterson(1975) estimator for non synonymous positions
		 *
		 * Gaps are automatically excluded
		 * @param psc a PolymorphismSequenceContainer
		 * @param gc a GeneticCode
		 */
		static double watterson75NonSynonymous(const PolymorphismSequenceContainer & psc, const GeneticCode & gc);

		/**
		  * @brief Compute the synonymous nucleotide diversity, pi
		  *
		  * Gaps are automatically excluded
		  * @param psc a PolymorphismSiteContainer
		  * @param gc a GeneticCode
		  * @param stopfalg a boolean set by default to true if you don't want
		  * to take gaps into account
		  * @param minchange a boolean set to false (see CodonSiteTools)
		  */
         static double piSynonymous(const PolymorphismSequenceContainer & psc, const GeneticCode & gc, bool stopflag = true, bool minchange=false);


		/**
		  * @brief Compute the non-synonymous nucleotide diversity, pi
		  *
		  * Gaps are automatically excluded
		  * @param v a SiteContainer
		  * @param gc a GeneticCode
		  * @param stopfalg a boolean set by default to true if you don't want
		  * to take gaps into account
		  * @param minchange a boolean set to false (see CodonSiteTools)
		  */
         static double piNonSynonymous(const PolymorphismSequenceContainer & psc, const GeneticCode & gc, bool stopflag = true, bool minchange=false);


		/**
		  * @brief compute the mean number of synonymous site in an alignment
		  *
		  * A site is x% synonymous if x% of possible mutations are synonymous
		  * The transition/transversion can be taken into account (use the variable ratio)
		  * Gaps are automatically excluded
		  * @param v a SiteContainer
		  * @param gc a GeneticCode
		  * @param ratio a double
		  * @param stopfalg a boolean set by default to true if you don't want
		  * to take stop codons into account
		  */
		static double meanSynonymousSitesNumber(const PolymorphismSequenceContainer & psc, const GeneticCode & gc, double ratio=1.0, bool stopflag=true);


		/**
		  * @brief compute the mean number of non-synonymous site in an alignment
		  *
		  * A site is x% synonymous if x% of possible mutations are synonymous
		  * The transition/transversion can be taken into account (use the variable ratio)
		  * Gaps are automatically excluded
		  * @param v a SiteContainer
		  * @param gc a GeneticCode
		  * @param ratio a double
		  * @param stopfalg a boolean set by default to true if you don't want
		  * to take stop codons into account
		  */
		static double meanNonSynonymousSitesNumber(const PolymorphismSequenceContainer & psc, const GeneticCode & gc, double ratio=1.0, bool stopflag=true);


		/**
		  * @brief compute the number of synonymous subsitutions in an alignment
		  *
		  * Gaps and unresolved sites are automatically excluded
		  * @param psc a PolymorphismSequenceContainer
		  * @param gc a GeneticCode
		  * @param freqmin a double, to exclude snp in frequency strictly lower than freqmin
		  */
		static unsigned int synonymousSubstitutionsNumber(const PolymorphismSequenceContainer & psc, const GeneticCode & gc, double freqmin = 0);


		/**
		  * @brief compute the number of  non synonymous subsitutions in an alignment
		  *
		  * Gaps and unresolved sites are automatically excluded
		  * @param sc a PolymorphismSequenceContainer
		  * @param gc a GeneticCode
		  * @param freqmin a double, to exclude snp in frequency strictly lower than freqmin
		  */
		static unsigned int nonSynonymousSubstitutionsNumber(const PolymorphismSequenceContainer & psc, const GeneticCode & gc, double freqmin = 0);


		/**
		  * @brief compute the number of fixed differences between two alignements
		  *
		  * @param pscin a PolymorphismSequenceContainer
		  * @param pscout a PolymorphismSequenceContainer
		  * @param psccons a PolymorphismSequenceContainer
		  * @param gc a GeneticCode
		  */
		static vector<unsigned int> fixedDifferences(const PolymorphismSequenceContainer & pscin, const PolymorphismSequenceContainer & pscout, PolymorphismSequenceContainer & psccons, const GeneticCode & gs);

		/**
		  * @brief return a vector containing Pa, Ps, Da, Ds
		  *
		  * Gaps and unresolved sites are automatically excluded
		  * @param ingroup a PolymorphismSequenceContainer
		  * @param outgroup a PolymorphismSequenceContainer
		  * @param gc a GeneticCode
		  * @param freqmin a double, to exclude snp in frequency strictly lower than freqmin
		  */
		static vector<unsigned int> MKtable(const PolymorphismSequenceContainer & ingroup, const PolymorphismSequenceContainer & outgroup , const GeneticCode & gc, double freqmin = 0);


		/**
		  * @brief return the neutrality index NI = (Pa/Ps)/(Da/Ds)
		  *
		  * Gaps and unresolved sites are automatically excluded
		  * @param ingroup a PolymorphismSequenceContainer
		  * @param outgroup a PolymorphismSequenceContainer
		  * @param gc a GeneticCode
		  * @param freqmin a double, to exclude snp in frequency strictly lower than freqmin
		  */
		static double neutralityIndex(const PolymorphismSequenceContainer & ingroup, const PolymorphismSequenceContainer & outgroup , const GeneticCode & gc, double freqmin = 0);



		/**
		 * @brief Return the Tajima's D test (Tajima 1989).
		 *
		 * Calculation using the number of polymorphic (segregating) sites.
		 * @f[
		 * D=\frac{\hat{\theta}_\pi-\hat{\theta}_S}{\sqrt{\textrm{V}\left(\hat{\theta}_\pi-\hat{\theta}_S\right)}}
		 * =\frac{\hat{\theta}_\pi-\hat{\theta}_S}{\sqrt{e_1S+e_2S(S-1)}}
		 * @f]
		 * @param psc a PolymorphismSequenceContainer
		 * @param gapflag flag set by default to true if you don't want to
		 * take gap into account
		 */
		static double tajimaDSS(const PolymorphismSequenceContainer & psc, bool gapflag = true);


		/**
		 * @brief Return the Tajima's D test (Tajima 1989).
		 *
		 * Calculation using the total number of mutation.
		 * @f[
		 * D=\frac{\hat{\theta}_\pi-\frac{\eta}{a_1}}{\sqrt{e_1\eta+e_2\eta(\eta-1)}}
		 * @f]
		 * @param psc a PolymorphismSequenceContainer
		 * @param gapflag flag set by default to true if you don't want to
		 * take gap into account
		 */
		static double tajimaDTNM(const PolymorphismSequenceContainer & psc, bool gapflag = true);


		/**
		 * @brief Return the Fu and Li D test (1993).
		 *
		 * @param ingroup a PolymorphismSequenceContainer
		 * @param outgroup a PolymorphismSequenceContainer
                 * This version uses the number of mutations
                 * If the outgroup contains more than one sequence the sites with more than one
                 * variant will not be considered for external branch mutations !
		 */
		static double fuliD(const PolymorphismSequenceContainer & ingroup, const PolymorphismSequenceContainer & outgroup);


		/**
		 * @brief Return the Fu and Li D<sup>*</sup> test (1993).
		 *
		 * @param group a PolymorphismSequenceContainer
		 */
		static double fuliDstar(const PolymorphismSequenceContainer & group);



		/**

		 * @brief Return the Fu and Li F test (1993).

		 *

		 * @param ingroup a PolymorphismSequenceContainer

		 * @param outgroup a PolymorphismSequenceContainer

                 * This version uses the number of mutations

		 * This version uses the number of mutations

                 * If the outgroup contains more than one sequence the sites with more than one

                 * variant will not be considered for externalbranch mutations !

                 */

		static double fuliF(const PolymorphismSequenceContainer & ingroup, const PolymorphismSequenceContainer & outgroup);



		/**

		 * @brief Return the Fu and Li F<sup>*</sup> test (1993).

		 *

		 * @param group a PolymorphismSequenceContainer

		 */

		static double fuliFstar(const PolymorphismSequenceContainer & group);





	    /**

	     * @brief generate a special PolymorphismSequenceContainer for linkage disequilbrium analysis

	     *

         * Only polymorphic sites with 2 alleles are kept

         * The value 1 is assigned to the more frequent allele, and 0 to the less frequent

	     * @param psc a PolymorphismSequenceContainer

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

		static PolymorphismSequenceContainer * generateLDContainer(const PolymorphismSequenceContainer & psc, bool keepsingleton=true, double freqmin=0) throw (Exception);

	    /**

	     * @brief give the vector of the pairwise distances between site positions corresponding to a LD SequencePolymorphismContainer

	     *

	     * @param psc a PolymorphismSequenceContainer

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

		static Vdouble pairwiseDistances1(const PolymorphismSequenceContainer & psc, bool keepsingleton=true, double freqmin=0);





	    /**

	     * @brief give the vector of all mean pairwise distance between two sites to a LD SequencePolymorphismContainer

	     *

	     * pairwise distances are computed for each sequence separately, excluding gaps. Then the mean is taken over all the sequences

	     * @param psc a PolymorphismSequenceContainer

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

		static Vdouble pairwiseDistances2(const PolymorphismSequenceContainer & psc, bool keepsingleton=true, double freqmin=0);





	    /**

	     * @brief give the vector of all mean pairwise D value between two sites (Lewontin & Kojima 1964)

	     *

	     * @param psc a PolymorphismSequenceContainer

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

		static Vdouble pairwiseD(const PolymorphismSequenceContainer & psc, bool keepsingleton=true, double freqmin=0);



	    /**

	     * @brief give the vector of all mean pairwise D' value between two sites (Lewontin 1964)

	     *

	     * @param psc a PolymorphismSequenceContainer

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

		static Vdouble pairwiseDprime(const PolymorphismSequenceContainer & psc, bool keepsingleton=true, double freqmin=0);



	    /**

	     * @brief give the vector of all mean pairwise R� value between two sites (Hill&Robertson 1968)

	     *

	     * @param psc a PolymorphismSequenceContainer

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

		static Vdouble pairwiseR2(const PolymorphismSequenceContainer & psc, bool keepsingleton=true, double freqmin=0);





	    /**

	     * @brief give mean D over all pairwise comparisons

	     *

	     * @param psc a PolymorphismSequenceContainer

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

		static double meanD(const PolymorphismSequenceContainer & psc, bool keepsingleton=true, double freqmin=0);



	    /**

	     * @brief give mean D' over all pairwise comparisons

	     *

	     * @param psc a PolymorphismSequenceContainer

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

		static double meanDprime(const PolymorphismSequenceContainer & psc, bool keepsingleton=true, double freqmin=0);



	    /**

	     * @brief give mean R� over all pairwise comparisons

	     *

	     * @param psc a PolymorphismSequenceContainer

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

		static double meanR2(const PolymorphismSequenceContainer & psc, bool keepsingleton=true, double freqmin=0);



	    /**

	     * @brief give mean pairwise distances between sites / method 1: differences between sequences are not taken into account

	     *

	     * @param psc a PolymorphismSequenceContainer

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

		static double meanDistance1(const PolymorphismSequenceContainer & psc, bool keepsingleton=true, double freqmin=0);



	    /**

	     * @brief give mean pairwise distances between sites / method 2: differences between sequences are taken into account

	     *

	     * @param psc a PolymorphismSequenceContainer

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

		static double meanDistance2(const PolymorphismSequenceContainer & psc, bool keepsingleton=true, double freqmin=0);



	    /**

	     * @brief give the slope of the regression |D| = 1+a*distance

	     *

             * The slope is given in |D| per kb

             *

	     * @param psc a PolymorphismSequenceContainer

             * @param distance1 a boolean (true to use distance1, false to use distance2, false by default)

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

        static double originRegressionD(const PolymorphismSequenceContainer & psc, bool distance1=false, bool keepsingleton=true, double freqmin=0);





	    /**

	     * @brief give the slope of the regression |D'| = 1+a*distance

	     *

             * The slope is given in |D'| per kb

             *

	     * @param psc a PolymorphismSequenceContainer

             * @param distance1 a boolean (true to use distance1, false to use distance2, false by default)

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

        static double originRegressionDprime(const PolymorphismSequenceContainer & psc, bool distance1=false, bool keepsingleton=true, double freqmin=0);



	    /**

	     * @brief give the slope of the regression R� = 1+a*distance

	     *

             * The slope is given in R� per kb

             *

	     * @param psc a PolymorphismSequenceContainer

             * @param distance1 a boolean (true to use distance1, false to use distance2, false by default)

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

        static double originRegressionR2(const PolymorphismSequenceContainer & psc, bool distance1=false, bool keepsingleton=true, double freqmin=0);



	    /**

	     * @brief give the slope and the origin of the regression |D| = a*distance+b

	     *

             * The slope is given in |D| per kb

             *

	     * @param psc a PolymorphismSequenceContainer

             * @param distance1 a boolean (true to use distance1, false to use distance2, false by default)

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

       static Vdouble linearRegressionD(const PolymorphismSequenceContainer & psc, bool distance1=false, bool keepsingleton=true, double freqmin=0);



	    /**

	     * @brief give the slope and the origin of the regression |D'| = a*distance+b

	     *

             * The slope is given in |D'| per kb

             *

	     * @param psc a PolymorphismSequenceContainer

             * @param distance1 a boolean (true to use distance1, false to use distance2, false by default)

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

       static Vdouble linearRegressionDprime(const PolymorphismSequenceContainer & psc, bool distance1=false, bool keepsingleton=true, double freqmin=0);



	    /**

	     * @brief give the slope and the origin of the regression R� = a*distance+b

	     *

             * The slope is given in R� per kb

             *

	     * @param psc a PolymorphismSequenceContainer

             * @param distance1 a boolean (true to use distance1, false to use distance2, false by default)

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

        static Vdouble linearRegressionR2(const PolymorphismSequenceContainer & psc, bool distance1=false, bool keepsingleton=true, double freqmin=0);



	    /**

	     * @brief give the slope of the regression R� = 1/(1+a*distance)

	     *

             * To fit the theoretical prediction R� = 1/(1+4Nr)

             * The slope is given in R� per kb

             *

	     * @param psc a PolymorphismSequenceContainer

             * @param distance1 a boolean (true to use distance1, false to use distance2, false by default)

	     * @param keepsingleton a boolean (true by default, false to exclude singleton)

	     * @param freqmin a float (to exlude site with the lowest allele frequency less than the threshold given by freqmin, 0 by default)

	     */

        static double inverseRegressionR2(const PolymorphismSequenceContainer & psc, bool distance1=false, bool keepsingleton=true, double freqmin=0);




	    /**

	     * @brief give estimate of C=4Nr using Hudson method (1987)
	     * @param psc a PolymorphismSequenceContainer
	     * @param precision
	     * @param cinf and csup, initial values
	     */
		static double hudson87(const PolymorphismSequenceContainer & psc, double precision = 0.000001, double cinf=0.001, double csup=10000);




	private:

		/**

		 * @brief Count the number of mutation for a site.

		 */

		static unsigned int _getMutationNumber(const Site & site);



		/**

		 * @brief Count the number of singleton for a site.

		 */

		static unsigned int _getSingletonNumber(const Site & site);

                /**

		 * @brief Count the number of singleton for a site.

                //khalid

                //will count singletons that are not in site_out (a site in outgroup)

                //site_in is a site from an ingroup

                */

                static unsigned int _getDerivedSingletonNumber(const Site & site_in,const Site & site_out );



		/**

		 * @brief Get usefull values for theta estimators.

		 *

		 * @return A map with 11 values. Keys are a1, a2, a1n, b1, b2, c1, c2, cn, dn, e1 and e2.

		 * The values are :

		 * @f[

		 * a_1=\sum_{i=1}^{n-1}\frac{1}{i} \qquad a_2=\sum_{i=1}^{n-1}\frac{1}{i^2}

		 * @f]

		 * @f[

		 * a_{1n}=\sum_{i=1}^{n}\frac{1}{i}

		 * @f]

		 * @f[

		 * b_1=\frac{n+1}{3(n-1)} \qquad b_2=\frac{2(n^2+n+3)}{9n(n-1)}

		 * @f]

		 * @f[

		 * c_1=b_1-\frac{1}{a_1} \qquad c_2=b_2-\frac{n+2}{a_1n}+\frac{a_2}{a_1^2}

		 * @f]

		 * @f[

		 * c_n=2\frac{na_1-2(n-1)}{(n-1)(n-2)}

		 * @f]

		 * @f[

		 * d_n=c_n+\frac{n-2}{(n-1)^2}+\frac{2}{n-1}\left(\frac{3}{2}-\frac{2a_{1n}-3}{n-2}-\frac{1}{n}\right)

		 * @f]

		 * @f[

		 * e_1=\frac{c_1}{a_1} \qquad e_2=\frac{c_2}{a_1^2+a_2}

		 * @f]

		 * where @f$n@f$ is the number of observed sequences.

		 */

		static map<string, double> _getUsefullValues(unsigned int n);


	    /**

	     * @brief give the left hand term of equation (4) in Hudson (1987)

         * This term is used in hudson87

	     * @param psc a PolymorphismSequenceContainer

	     */

	    /**

	     * @brief give the right hand term of equation (4) in Hudson (1987)

         * This term is used in hudson87

	     * @param psc a PolymorphismSequenceContainer

	     */

		static double _leftHandHudson(const PolymorphismSequenceContainer & psc);

		static double _rightHandHudson(double c, unsigned int n);


		/*******************************************************************************/

};

#endif // _SEQUENCESTATISTICS_H_

