//
// File: PolymorphismSequenceContainerTools.h
// Authors: bazin <bazin@univ-montp2.fr>
//          Sylvain Gaillard <yragael2001@yahoo.fr>
// Last modification : Tuesday July 27 2004
//

#ifndef _POLYMORPHISMSEQUENCECONTAINERTOOL_H_
#define _POLYMORPHISMSEQUENCECONTAINERTOOL_H_

// from SeqLib
#include <Seq/Mase.h>
#include <Seq/VectorSiteContainer.h>
#include <Seq/MaseTools.h>
#include <Seq/SequenceContainerTools.h>
#include <Seq/SiteIterator.h>

// from STL
#include <map>

// From Local
#include "PolymorphismSequenceContainer.h"

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
		static PolymorphismSequenceContainer * read(const string & path,
				const Alphabet *alpha) throw (Exception);

		/**
		 * @brief Extract ingroup sequences from a PolymorphismSequenceContainer and create a new
		 * one.
		 *
		 * @param psc a PolymorphismSequenceContainer reference
		 */
		static PolymorphismSequenceContainer * extractIngroup (const PolymorphismSequenceContainer & psc ) throw (Exception);

		/**
		 * @brief Extract outgroup sequences from a PolymorphismSequenceContainer and create a new one.
		 *
		 * @param psc a PolymorphismSequenceContainer reference
		 */
		static PolymorphismSequenceContainer * extractOutgroup (const PolymorphismSequenceContainer & psc) throw (Exception);

		/**
		 * @brief Retrieves sites without gaps from PolymorphismSequenceContainer.
		 *
		 * @param psc a PolymorphismSequenceContainer reference
		 */
		static PolymorphismSequenceContainer * getSitesWithoutGaps (const PolymorphismSequenceContainer & psc ) throw (Exception);

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
