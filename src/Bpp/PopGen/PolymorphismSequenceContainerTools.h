// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
using SimpleSiteContainerIterator = SimpleTemplateSiteContainerIterator<Site, Sequence, std::string>;
using CompleteSiteContainerIterator = CompleteTemplateSiteContainerIterator<Site, Sequence, std::string>;
using NoGapSiteContainerIterator = NoGapTemplateSiteContainerIterator<Site, Sequence, std::string>;

/**
 * @brief Utilitary function to manipulate PolymorphismSequenceContainer
 *
 * @author Sylvain Gaillard
 */
class PolymorphismSequenceContainerTools
{
public:
  // Class destructor:
  virtual ~PolymorphismSequenceContainerTools();

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
  static std::unique_ptr<PolymorphismSequenceContainer> read(
      const std::string& path,
      std::shared_ptr<const Alphabet> alpha);

  /**
   * @brief Extract ingroup sequences from a PolymorphismSequenceContainer and create a new one.
   *
   * @param psc a PolymorphismSequenceContainer reference
   *
   * @throw Exception if there is no ingroup sequence
   */
  static std::unique_ptr<PolymorphismSequenceContainer> extractIngroup(
      const PolymorphismSequenceContainer& psc);

  /**
   * @brief Extract outgroup sequences from a PolymorphismSequenceContainer and create a new one.
   *
   * @param psc a PolymorphismSequenceContainer reference
   *
   * @throw Exception if there is no outgroup sequence
   */
  static std::unique_ptr<PolymorphismSequenceContainer> extractOutgroup(
      const PolymorphismSequenceContainer& psc);

  /**
   * @brief Extract a special group from the PolymorphismSequenceContainer.
   *
   * @param psc a PolymorphismSequenceContainer reference.
   * @param group_id the group identifier as an size_t.
   *
   * @throw GroupNotFoundException if group_id is not found.
   */
  static std::unique_ptr<PolymorphismSequenceContainer> extractGroup(
      const PolymorphismSequenceContainer& psc,
      size_t groupId);

  /**
   * @brief Extract selected sequences
   *
   * @param psc a PolymorphismSequenceContainer reference.
   * @param ss a sequence selection.
   */
  static std::unique_ptr<PolymorphismSequenceContainer> getSelectedSequences(
      const PolymorphismSequenceContainer& psc,
      const SequenceSelection& ss);

  /**
   * @brief Get a random set of sequences
   *
   * @param psc a PolymorphismSequenceContainer reference
   * @param n the number of sequence to get
   * @param replace a boolean flag true for sampling with replacement
   */
  static std::unique_ptr<PolymorphismSequenceContainer> sample(
      const PolymorphismSequenceContainer& psc,
      size_t n,
      bool replace = true);

  /**
   * @brief Retrieves sites without gaps from PolymorphismSequenceContainer.
   *
   * @param psc a PolymorphismSequenceContainer reference
   */
  static std::unique_ptr<PolymorphismSequenceContainer> getSitesWithoutGaps(
      const PolymorphismSequenceContainer& psc);

  /**
   * @brief Return number of sites without gaps in a PolymorphismSequenceContainer.
   *
   * @param psc a PolymorphismSequenceContainer reference
   * @param ingroup a boolean set to true if you want to take only ingroup sequences into account
   *
   * @throw Exception if there is no ingroup sequence
   */
  static size_t getNumberOfNonGapSites(
      const PolymorphismSequenceContainer& psc,
      bool ingroup);

  /**
   * @brief Return number of completely resolved sites in a PolymorphismSequenceContainer.
   *
   *
   * @param psc a PolymorphismSequenceContainer reference
   * @param ingroup a boolean set to true if you want to take only ingroup sequences into account
   *
   * @throw Exception if there is no ingroup sequence
   */
  static size_t getNumberOfCompleteSites(
      const PolymorphismSequenceContainer& psc,
      bool ingroup);


  /**
   * @brief Retrieves complete sites from a PolymorphismSequenceContainer.
   *
   * @param psc a PolymorphismSequenceContainer reference
   */
  static std::unique_ptr<PolymorphismSequenceContainer> getCompleteSites(
      const PolymorphismSequenceContainer& psc);


  /**
   * @brief exclude flanking sites with gap but keep gap sites within the alignment
   *
   * @param psc a PolymorphismSequenceContainer reference
   */
  static std::unique_ptr<PolymorphismSequenceContainer> excludeFlankingGap(
      const PolymorphismSequenceContainer& psc);

  /**
   * @brief Get a PolymorphismSequenceContainer corresponding to a site selection annotated in the mase comments
   *
   * Warning: in the new PolymorphismSequenceContainer the mase comments are lost
   * Information about cds positions and start codon is no more available
   *
   * @param psc a PolymorphismSequenceContainer.
   * @param setName The name of the set to retrieve.
   * @param phase a boolean set to true if you want to take the phase into account during the extraction. It removes the useless sites.
   */
  static std::unique_ptr<PolymorphismSequenceContainer> getSelectedSites(
      const PolymorphismSequenceContainer& psc,
      const std::string& setName, bool phase);

  /**
   * @brief Retrieve non-coding sites defined in the mase file header
   *
   * Be carefull: to use before excluding gap
   *
   * @param psc a PolymorphismSequenceContainer reference
   * @param setName name of the CDS site selection
   */
  static std::unique_ptr<PolymorphismSequenceContainer> getNonCodingSites(
      const PolymorphismSequenceContainer& psc,
      const std::string& setName);

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
  static std::unique_ptr<PolymorphismSequenceContainer> getOnePosition(
      const PolymorphismSequenceContainer& psc,
      const std::string& setName,
      size_t pos);

  /**
   * @brief Retrieve intron sites
   *
   * Same as getNonCodgingSites but exclude 5' and 3' flanking regions if there are
   *
   * @param psc a PolymorphismSequenceContainer
   * @param setName name of the CDS site selection
   * @param gCode The genetic code to use
   */
  static std::unique_ptr<PolymorphismSequenceContainer> getIntrons(
      const PolymorphismSequenceContainer& psc,
      const std::string& setName,
      const GeneticCode& gCode);

  /**
   * @brief Retrieve 5' sites
   *
   * @param psc a PolymorphismSequenceContainer
   * @param setName name of the CDS site selection
   */
  static std::unique_ptr<PolymorphismSequenceContainer> get5Prime(
      const PolymorphismSequenceContainer& psc,
      const std::string& setName);

  /**
   * @brief Retrieve 3' sites
   *
   * @param psc a PolymorphismSequenceContainer
   * @param setName name of the CDS site selection
   * @param gCode The genetic code to use
   */
  static std::unique_ptr<PolymorphismSequenceContainer> get3Prime(
      const PolymorphismSequenceContainer& psc,
      const std::string& setName,
      const GeneticCode& gCode);

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
  static std::unique_ptr<PolymorphismSequenceContainer> getSynonymousSites(
      const PolymorphismSequenceContainer& psc,
      const GeneticCode& gCode);

  /**
   * @brief Retrieve non-synonymous codon sites
   *
   * @param psc a PolymorphismSequenceContainer.
   * @param gCode The genetic code to use
   *
   * @return A new PolymorphismSequenceContainer with only non-synonymous sites.
   */
  static std::unique_ptr<PolymorphismSequenceContainer> getNonSynonymousSites(
      const PolymorphismSequenceContainer& psc,
      const GeneticCode& gCode);
};
} // end of namespace bpp;

#endif // _POLYMORPHISMSEQUENCECONTAINERTOOL_H_
