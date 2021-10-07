//
// File: PolymorphismSequenceContainerTools.h
// Authors: Eric Bazin
//          Sylvain Gaillard
// Created on: Thursday July 29 2004
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Text/StringTokenizer.h>

// from bpp-seq
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/Io/Mase.h>
#include <Bpp/Seq/Io/MaseTools.h>
#include <Bpp/Seq/Container/SiteContainerIterator.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>

// from STL
#include <string>

// From Local
#include "PolymorphismSequenceContainer.h"
#include "GeneralExceptions.h"

namespace bpp
{
/**
 * @brief Utilitary function to manipulate PolymorphismSequenceContainer
 *
 * @author Sylvain Gaillard
 */

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
   *
   * @throw Exception if the file is not in the specified format
   */
  static PolymorphismSequenceContainer* read(const std::string& path, const Alphabet* alpha);

  /**
   * @brief Extract ingroup sequences from a PolymorphismSequenceContainer and create a new one.
   *
   * @param psc a PolymorphismSequenceContainer reference
   *
   * @throw Exception if there is no ingroup sequence
   */
  static PolymorphismSequenceContainer* extractIngroup(const PolymorphismSequenceContainer& psc);

  /**
   * @brief Extract outgroup sequences from a PolymorphismSequenceContainer and create a new one.
   *
   * @param psc a PolymorphismSequenceContainer reference
   *
   * @throw Exception if there is no outgroup sequence
   */
  static PolymorphismSequenceContainer* extractOutgroup(const PolymorphismSequenceContainer& psc);

  /**
   * @brief Extract a special group from the PolymorphismSequenceContainer.
   *
   * @param psc a PolymorphismSequenceContainer reference.
   * @param group_id the group identifier as an size_t.
   *
   * @throw GroupNotFoundException if group_id is not found.
   */
  static PolymorphismSequenceContainer* extractGroup(const PolymorphismSequenceContainer& psc, size_t group_id);

  /**
   * @brief Extract selected sequences
   *
   * @param psc a PolymorphismSequenceContainer reference.
   * @param ss a sequence selection.
   *
   */
  static PolymorphismSequenceContainer* getSelectedSequences(const PolymorphismSequenceContainer& psc, const SequenceSelection& ss);

  /**
   * @brief Get a random set of sequences
   *
   * @param psc a PolymorphismSequenceContainer reference
   * @param n the number of sequence to get
   * @param replace a boolean flag true for sampling with replacement
   */
  static PolymorphismSequenceContainer* sample(const PolymorphismSequenceContainer& psc, size_t n, bool replace = true);

  /**
   * @brief Retrieves sites without gaps from PolymorphismSequenceContainer.
   *
   * @param psc a PolymorphismSequenceContainer reference
   */
  static PolymorphismSequenceContainer* getSitesWithoutGaps(const PolymorphismSequenceContainer& psc);

  /**
   * @brief Return number of sites without gaps in a PolymorphismSequenceContainer.
   *
   * @param psc a PolymorphismSequenceContainer reference
   * @param ingroup a boolean set to true if you want to take only ingroup sequences into account
   *
   * @throw Exception if there is no ingroup sequence
   */
  static size_t getNumberOfNonGapSites(const PolymorphismSequenceContainer& psc, bool ingroup);

  /**
   * @brief Return number of completely resolved sites in a PolymorphismSequenceContainer.
   *
   *
   * @param psc a PolymorphismSequenceContainer reference
   * @param ingroup a boolean set to true if you want to take only ingroup sequences into account
   *
   * @throw Exception if there is no ingroup sequence
   */
  static size_t getNumberOfCompleteSites(const PolymorphismSequenceContainer& psc, bool ingroup);


  /**
   * @brief Retrieves complete sites from a PolymorphismSequenceContainer.
   *
   * @param psc a PolymorphismSequenceContainer reference
   */
  static PolymorphismSequenceContainer* getCompleteSites(const PolymorphismSequenceContainer& psc);


  /**
   * @brief exclude flanking sites with gap but keep gap sites within the alignment
   *
   * @param psc a PolymorphismSequenceContainer reference
   */
  static PolymorphismSequenceContainer* excludeFlankingGap(const PolymorphismSequenceContainer& psc);

  /**
   * @brief Get a PolymorphismSequenceContainer corresponding to a site selection annotated in the mase comments
   *
   * Be carefull : in the new PolymorphismSequenceContainer the mase comments are lost
   * Information about cds positions and start codon is no more available
   *
   * @param psc a PolymorphismSequenceContainer.
   * @param setName The name of the set to retrieve.
   * @param phase a boolean set to true if you want to take the phase into account during the extraction. It removes the useless sites.
   */
  static PolymorphismSequenceContainer* getSelectedSites(const PolymorphismSequenceContainer& psc, const std::string& setName, bool phase);

  /**
   * @brief Retrieve non-coding sites defined in the mase file header
   *
   * Be carefull: to use before excluding gap
   *
   * @param psc a PolymorphismSequenceContainer reference
   * @param setName name of the CDS site selection
   */
  static PolymorphismSequenceContainer* getNonCodingSites(const PolymorphismSequenceContainer& psc, const std::string& setName);

  /**
   * @brief Retrieve sites at one codon position (1,2,3)
   *
   * Be carefull: to use before excluding gap
   * Be careful: if there is no phase information, the method catch an exception and set the phase to 1
   * This allows to use this method for PolymorphismSequenceContainer generated by getSelectedSequence
   *
   * @param psc a PolymorphismSequenceContainer reference
   * @param setName name of the CDS site selection
   * @param pos position index.
   */
  static PolymorphismSequenceContainer* getOnePosition(const PolymorphismSequenceContainer& psc, const std::string& setName, size_t pos);

  /**
   * @brief Retrieve intron sites
   *
   * Same as getNonCodgingSites but exclude 5' and 3' flanking regions if there are
   *
   * @param psc a PolymorphismSequenceContainer
   * @param setName name of the CDS site selection
   * @param gCode The genetic code to use
   */
  static PolymorphismSequenceContainer* getIntrons(const PolymorphismSequenceContainer& psc, const std::string& setName, const GeneticCode* gCode);

  /**
   * @brief Retrieve 5' sites
   *
   * @param psc a PolymorphismSequenceContainer
   * @param setName name of the CDS site selection
   */
  static PolymorphismSequenceContainer* get5Prime(const PolymorphismSequenceContainer& psc, const std::string& setName);

  /**
   * @brief Retrieve 3' sites
   *
   * @param psc a PolymorphismSequenceContainer
   * @param setName name of the CDS site selection
   * @param gCode The genetic code to use
   */
  static PolymorphismSequenceContainer* get3Prime(const PolymorphismSequenceContainer& psc, const std::string& setName, const GeneticCode* gCode);

  /**
   * @brief Get the species name of the ingroup
   *
   * @param psc a PolymorphismSequenceContainer.
   */
  static std::string getIngroupSpeciesName(const PolymorphismSequenceContainer& psc);

  /**
   * @brief Retrieve synonymous codon sites
   *
   * @param psc a PolymorphismSequenceContainer.
   * @param gCode The genetic code to use
   *
   * @return A new PolymorphismSequenceContainer with only synonymous sites.
   */
  static PolymorphismSequenceContainer* getSynonymousSites(const PolymorphismSequenceContainer& psc, const GeneticCode& gCode);

  /**
   * @brief Retrieve non-synonymous codon sites
   *
   * @param psc a PolymorphismSequenceContainer.
   * @param gCode The genetic code to use
   *
   * @return A new PolymorphismSequenceContainer with only non-synonymous sites.
   */
  static PolymorphismSequenceContainer* getNonSynonymousSites(const PolymorphismSequenceContainer& psc, const GeneticCode& gCode);
};
} // end of namespace bpp;

#endif// _POLYMORPHISMSEQUENCECONTAINERTOOL_H_
