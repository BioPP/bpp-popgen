/*
 * File: PolymorphismSequenceContainerTools.h
 * Authors: Eric Bazin <bazin@univ-montp2.fr>
 *          Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
*/
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

#ifndef _POLYMORPHISMSEQUENCECONTAINERTOOL_H_
#define _POLYMORPHISMSEQUENCECONTAINERTOOL_H_

// from SeqLib
#include <Seq/Mase.h>
#include <Seq/VectorSiteContainer.h>
#include <Seq/MaseTools.h>
#include <Seq/SequenceContainerTools.h>
#include <Seq/SiteIterator.h>
#include <Seq/SiteTools.h>

// from Utils
#include <Utils/StringTokenizer.h>

// from STL
#include <map>

// From Local
#include "PolymorphismSequenceContainer.h"
#include "GeneralExceptions.h"

using namespace std;

class PolymorphismSequenceContainerTools
{
	public:
		// Class destructor:
		~PolymorphismSequenceContainerTools();

		/*******************************************************************************/
	public:

		/**
		 * @brief Read a Mase+ file and return a PolymorphismSequenceContainer. Toggle Sequence
		 * when selection tag begin with OUTGROUP (see Polymorphix)
		 *
		 * @param path Path to the Mase+ file
		 * @param alpha Sequence Alphabet
		 */
		static PolymorphismSequenceContainer * read(const string & path, const Alphabet * alpha) throw (Exception);

		/**
		 * @brief Extract ingroup sequences from a PolymorphismSequenceContainer and create a new one.
		 *
		 * @param psc a PolymorphismSequenceContainer reference
		 */
		static PolymorphismSequenceContainer * extractIngroup (const PolymorphismSequenceContainer & psc) throw (Exception);

		/**
		 * @brief Extract outgroup sequences from a PolymorphismSequenceContainer and create a new one.
		 *
		 * @param psc a PolymorphismSequenceContainer reference
		 */
		static PolymorphismSequenceContainer * extractOutgroup (const PolymorphismSequenceContainer & psc) throw (Exception);

		/**
		 * @brief Extract a special group from the PolymorphismSequenceContainer.
		 *
		 * @param psc a PolymorphismSequenceContainer reference.
		 * @param group_id the group identifier as an unsigned int.
		 *
		 * @throw GroupNotFoundException if group_id is not found.
		 */
		static PolymorphismSequenceContainer * extractGroup(const PolymorphismSequenceContainer & psc, unsigned int group_id) throw (Exception);

		/**
		 * @brief Retrieves sites without gaps from PolymorphismSequenceContainer.
		 *
		 * @param psc a PolymorphismSequenceContainer reference
		 */
		static PolymorphismSequenceContainer * getSitesWithoutGaps (const PolymorphismSequenceContainer & psc);

		/**
		 * @brief Return number of sites without gaps in a PolymorphismSequenceContainer.
		 *
		 * @param psc a PolymorphismSequenceContainer reference
		 * @param ingroup a boolean set to true if you want to take only ingroup sequences into account
		 */
		static unsigned int getNumberOfNonGapSites(const PolymorphismSequenceContainer & psc, bool ingroup) throw (Exception);

		/**
		 * @brief Return number of completely resolved sites in a PolymorphismSequenceContainer.
		 *
		 *
		 * @param psc a PolymorphismSequenceContainer reference
		 * @param ingroup a boolean set to true if you want to take only ingroup sequences into account
		 */
		static unsigned int getNumberOfCompleteSites(const PolymorphismSequenceContainer & psc, bool ingroup) throw (Exception);


		/**
		 * @brief Retrieves complete sites from a PolymorphismSequenceContainer.
		 *
		 * @param psc a PolymorphismSequenceContainer reference
		 */
		static PolymorphismSequenceContainer * getCompleteSites(const PolymorphismSequenceContainer & psc);


		 /**
          * @brief exclude flanking sites with gap but keep gap sites within the alignment
	      *
	      * @param psc a PolymorphismSequenceContainer reference
	      */
		static PolymorphismSequenceContainer * excludeFlankingGap(const PolymorphismSequenceContainer & psc) throw (Exception);

		/**
	     * @brief Retrieve coding sites defined in the mase file header
	     *
	     * @param psc a PolymorphismSequenceContainer reference
	     * Be carefull: to use before excluding gap
	     */
		static PolymorphismSequenceContainer * getCodingSites(const PolymorphismSequenceContainer & psc) throw (Exception);

	    /**
	     * @brief Retrieve non-coding sites defined in the mase file header
	     *
	     * @param psc a PolymorphismSequenceContainer reference
	     * Be carefull: to use before excluding gap
	     */
		static PolymorphismSequenceContainer * getNonCodingSites(const PolymorphismSequenceContainer & psc) throw (Exception);

	    /**
	     * @brief Retrieve sites at one codon position (1,2,3)
	     *
	     * @param psc a PolymorphismSequenceContainer reference  *
	     * Be carefull: to use before excluding gap
	     * The first position of the alignement is supposed to be 1 (to be checked before using this method)
	     */
	    static PolymorphismSequenceContainer * getOnePosition(const PolymorphismSequenceContainer & psc, unsigned int pos) throw (Exception);

	    /**
	     * @brief Retrieve intron sites
	     *
	     * Same as getNonCodgingSites but exclude 5' and 3' falnking regions if there are
         * Assumed that the first coding site correspond to the first position
	     * @param psc a PolymorphismSequenceContainer
	     */
        static PolymorphismSequenceContainer * getIntrons(const PolymorphismSequenceContainer & psc, unsigned int start) throw (Exception);

	    /**
	     * @brief Retrieve 5' sites
	     *
             * Assumed that the first coding site correspond to the first position
	     * @param psc a PolymorphismSequenceContainer
	     */
        static PolymorphismSequenceContainer * get5Prime(const PolymorphismSequenceContainer & psc) throw (Exception);

	    /**
	     * @brief Retrieve 3' sites
	     *
         * Assumed that the first coding site correspond to the first position
	     * @param psc a PolymorphismSequenceContainer
	     */
        static PolymorphismSequenceContainer * get3Prime(const PolymorphismSequenceContainer & psc, unsigned int start) throw (Exception);

	/**
	 * @brief Get a PolymorphismSequenceContainer corresponding to a site selection annotated in the mase comments
	 *
         * Assumed that the first coding site correspond to the first position
	 * @param psc a PolymorphismSequenceContainer.
	 * @param setName The name of the set to retrieve.
	 * @param phase a boolean set to true if you want to take the phase into account during the extraction. It removes the useless sites.
   	 */
	static PolymorphismSequenceContainer * getSelectedSites(const PolymorphismSequenceContainer & psc, const string &setName, bool phase) throw (Exception);

	/**
	 * @brief Get the species name of the ingroup
         *
	 * @param psc a PolymorphismSequenceContainer.
   	 */
	static string getIngroupSpeciesName(const PolymorphismSequenceContainer & psc) throw (Exception);

		/*******************************************************************************/
};
#endif // _POLYMORPHISMSEQUENCECONTAINERTOOL_H_
