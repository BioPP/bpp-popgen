/*
 * File: PolymorphismSequenceContainerTools.h
 * Authors: bazin <bazin@univ-montp2.fr>
 *          Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday July 28 2004
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
		 * Only complete sites are considered.
		 *
		 * @param psc a PolymorphismSequenceContainer reference
		 * @param ingroup a boolean set to true if you want to take only ingroup sequences into account
		 */
		static unsigned int getNumberOfNonGapSites(const PolymorphismSequenceContainer & psc, bool ingroup) throw (Exception);

		/*******************************************************************************/
};
#endif // _POLYMORPHISMSEQUENCECONTAINERTOOL_H_
