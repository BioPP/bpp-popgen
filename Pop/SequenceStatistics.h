/*
 * File SequenceStatistics.h
 * Author : Eric Bazin <bazin@univ-montp2.fr>
 * Last modification : Wednesday July 28 2004
 */

//Secured inclusion of header's file
#ifndef _SEQUENCESTATISTICS_H_
#define _SEQUENCESTATISTICS_H_

//From the SeqLib library
#include <Seq/SiteIterator.h>
#include <Seq/SiteContainer.h>
#include <Seq/SymbolListTools.h>

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
		 * @brief Compute number of polymorphic site in an alignment
		 *
		 * Gaps are consider as mutations so if you want number of
		 * polymorphic site, you have to give a NonGapSiteIterator
		 *
		 * @param si a SiteIterator
		 */ 
		static unsigned int polymorphicSiteNumber( SiteIterator & si );

		/**
		 * @brief Compute number of polymorphic site in an alignment
		 *
		 * Gaps are consider as mutation variant
		 *
		 * @param v a SiteContainer
		 */
		static unsigned int polymorphicSiteNumber( const SiteContainer & v );

		/**
		 * @brief Count the number of singleton nucleotides in an alignment.
		 */
		static unsigned int countSingleton(SiteIterator & si);

		/**
		 * @brief Count the total number of mutatuons in an alignment.
		 *
		 * This count is assumed to be under an infinite site model.
		 */
		static unsigned int totNumberMutations(SiteIterator & si);

		/**
		 * @brief Compute diversity estimator Theta of Watterson (1975)
		 *
		 * @param v a SiteContainer
		 */
		static double watterson75( const SiteContainer & v );

		/**
		 * @brief Compute diversity estimator Theta of Tajima (1983)
		 *
		 * @param v a SiteContainer
		 */ 
		static double tajima83( const SiteContainer & v );

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
		 * take gap into account
		 */ 
		static double DVH ( const PolymorphismSequenceContainer & psc, bool gapflag = true );

		/**
		 * @brief Return the Tajima's D test (Tajima 1989).
		 *
		 * @f[
		 * D=\frac{\hat{\theta}_\pi-\hat{\theta}_S}{\sqrt{\textrm{V}\left(\hat{\theta}_\pi-\hat{\theta}_S\right)}}
		 * =\frac{\hat{\theta}_\pi-\hat{\theta}_S}{\sqrt{e_1S+e_2S(S-1)}}
		 * @f]
		 * @param psc a PolymorphismSequenceContainer
		 * @param gapflag flag set by default to true if you don't want to
		 * take gap into account
		 */
		static double tajimaD(const PolymorphismSequenceContainer & psc, bool gapflag = true);

	private:
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

		/*******************************************************************************/
};
#endif // _SEQUENCESTATISTICS_H_
