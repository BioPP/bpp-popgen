// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "PolymorphismSequenceContainerTools.h"

#include <Bpp/Seq/CodonSiteTools.h>

using namespace bpp;
using namespace std;

PolymorphismSequenceContainerTools::~PolymorphismSequenceContainerTools() {}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::read(
    const std::string& path,
    shared_ptr<const Alphabet> alpha)
{
  Mase ms;
  string key;
  auto seqc = ms.readSequences(path, alpha);
  auto psc = make_unique<PolymorphismSequenceContainer>(*seqc);
  Comments maseFileHeader = seqc->getComments();
  auto groupMap = MaseTools::getAvailableSequenceSelections(maseFileHeader);
  for (auto& mi : groupMap)
  {
    key = mi.first;
    if (key.compare(0, 8, "OUTGROUP") == 0)
    {
      auto ss = MaseTools::getSequenceSet(maseFileHeader, key);
      for (size_t i = 0; i != ss.size(); ++i)
      {
        psc->setAsOutgroupMember(ss[i]);
      }
    }
  }
  return psc;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::extractIngroup(
    const PolymorphismSequenceContainer& psc)
{
  SequenceSelection ss;
  auto psci = make_unique<PolymorphismSequenceContainer>(psc);
  for (size_t i = 0; i < psc.getNumberOfSequences(); ++i)
  {
    if (!psc.isIngroupMember(i))
      ss.push_back(i);
  }
  if (ss.size() == psc.getNumberOfSequences())
  {
    throw Exception("PolymorphismSequenceContainerTools::extractIngroup: no Ingroup sequences found.");
  }
  for (size_t i = ss.size(); i > 0; --i)
  {
    psci->deleteSequence(ss[i - 1]);
  }
  return psci;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::extractOutgroup(
    const PolymorphismSequenceContainer& psc)
{
  SequenceSelection ss;
  auto psci = make_unique<PolymorphismSequenceContainer>(psc);
  for (size_t i = 0; i < psc.getNumberOfSequences(); ++i)
  {
    if (psc.isIngroupMember(i) )
      ss.push_back(i);
  }
  if (ss.size() == psc.getNumberOfSequences())
  {
    throw Exception("PolymorphismSequenceContainerTools::extractOutgroup: no Outgroup sequences found.");
  }
  for (size_t i = ss.size(); i > 0; i--)
  {
    psci->deleteSequence(ss[i - 1]);
  }
  return psci;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::extractGroup(
    const PolymorphismSequenceContainer& psc,
    size_t groupId)
{
  SequenceSelection ss;
  auto psci = make_unique<PolymorphismSequenceContainer>(psc);
  for (size_t i = 0; i < psc.getNumberOfSequences(); ++i)
  {
    if (psc.getGroupId(i) != groupId)
      ss.push_back(i);
  }
  if (ss.size() == psc.getNumberOfSequences())
  {
    throw GroupNotFoundException("PolymorphismSequenceContainerTools::extractGroup: group_id not found.", groupId);
  }
  for (size_t i = ss.size(); i > 0; i--)
  {
    psci->deleteSequence(ss[i - 1]);
  }
  return psci;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::getSelectedSequences(
  const PolymorphismSequenceContainer& psc,
  const SequenceSelection& ss)
{
  auto newpsc = make_unique<PolymorphismSequenceContainer>(psc.getAlphabet());
  for (size_t i = 0; i < ss.size(); ++i)
  {
    auto tmpSeq = make_unique<Sequence>(psc.sequence(ss[i]));
    newpsc->addSequenceWithFrequency(tmpSeq->getName(), tmpSeq, psc.getSequenceCount(i));
    if (psc.isIngroupMember(i))
      newpsc->setAsIngroupMember(i);
    else
    {
      newpsc->setAsOutgroupMember(i);
      newpsc->setGroupId(i, psc.getGroupId(i));
    }
  }
  newpsc->setComments(psc.getComments());
  return newpsc;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::sample(
    const PolymorphismSequenceContainer& psc,
    size_t n,
    bool replace)
{
  size_t nbSeq = psc.getNumberOfSequences();
  vector<size_t> v;
  for (size_t i = 0; i < nbSeq; ++i)
  {
    v.push_back(i);
  }
  vector<size_t> vv(n);
  RandomTools::getSample(v, vv, replace);
  auto newpsc = PolymorphismSequenceContainerTools::getSelectedSequences(psc, vv);
  return newpsc;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::getSitesWithoutGaps(
    const PolymorphismSequenceContainer& psc)
{
  auto seqNames = psc.getSequenceNames();
  auto noGapCont = make_unique<PolymorphismSequenceContainer>(psc.getNumberOfSequences(), psc.getAlphabet());
  noGapCont->setSequenceNames(seqNames, false);
  size_t nbSeq = psc.getNumberOfSequences();
  for (size_t i = 0; i < nbSeq; ++i)
  {
    noGapCont->setSequenceCount(i, psc.getSequenceCount(i));
    if (psc.isIngroupMember(i))
      noGapCont->setAsIngroupMember(i);
    else
    {
      noGapCont->setAsOutgroupMember(i);
      noGapCont->setGroupId(i, psc.getGroupId(i));
    }
  }
  NoGapSiteContainerIterator ngsi(psc);
  while (ngsi.hasMoreSites()) {
    auto tmpSite = make_unique<Site>(ngsi.nextSite());
    noGapCont->addSite(tmpSite);
  }
  return noGapCont;
}

/******************************************************************************/

size_t PolymorphismSequenceContainerTools::getNumberOfNonGapSites(
    const PolymorphismSequenceContainer& psc,
    bool ingroup)
{
  size_t count = psc.getNumberOfSites();
  unique_ptr<PolymorphismSequenceContainer> npsc = nullptr;
  unique_ptr<SimpleSiteContainerIterator> ssi = nullptr;
  if (ingroup)
  {
    npsc = extractIngroup(psc);
    ssi.reset(new SimpleSiteContainerIterator(*npsc));
  }
  else
    ssi.reset(new SimpleTemplateSiteContainerIterator<Site, Sequence, string>(psc));
  while (ssi->hasMoreSites())
    if (SiteTools::hasGap(ssi->nextSite()))
      count--;
  return count;
}

/******************************************************************************/

size_t PolymorphismSequenceContainerTools::getNumberOfCompleteSites(
    const PolymorphismSequenceContainer& psc,
    bool ingroup)
{
  size_t count = psc.getNumberOfSites();
  unique_ptr<PolymorphismSequenceContainer> npsc = nullptr;
  unique_ptr<SimpleSiteContainerIterator> ssi = nullptr;
  if (ingroup)
  {
    npsc = extractIngroup(psc);
    ssi.reset(new SimpleSiteContainerIterator(*npsc));
  }
  else
    ssi.reset(new SimpleSiteContainerIterator(psc));
  while (ssi->hasMoreSites())
    if (!SiteTools::isComplete(ssi->nextSite()))
      count--;
  return count;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::getCompleteSites(
    const PolymorphismSequenceContainer& psc)
{
  auto seqNames = psc.getSequenceNames();
  auto complete = make_unique<PolymorphismSequenceContainer>(psc.getNumberOfSequences(), psc.getAlphabet());
  complete->setSequenceNames(seqNames, false);
  size_t nbSeq = psc.getNumberOfSequences();
  for (size_t i = 0; i < nbSeq; ++i)
  {
    complete->setSequenceCount(i, psc.getSequenceCount(i));
    if (psc.isIngroupMember(i))
      complete->setAsIngroupMember(i);
    else
    {
      complete->setAsOutgroupMember(i);
      complete->setGroupId(i, psc.getGroupId(i));
    }
  }
  CompleteSiteContainerIterator csi(psc);
  while (csi.hasMoreSites()) {
    auto tmpSite = make_unique<Site>(csi.nextSite());
    complete->addSite(tmpSite);
  }
  return complete;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::excludeFlankingGap(
    const PolymorphismSequenceContainer& psc)
{
  auto psci = make_unique<PolymorphismSequenceContainer>(psc);
  while (SiteTools::hasGap(psci->site(0)))
    psci->deleteSite(0);
  size_t i = 0;
  size_t n = psci->getNumberOfSites();
  while (SiteTools::hasGap(psci->site(n - i - 1)))
  {
    psci->deleteSite(n - i - 1);
    i++;
  }
  return psci;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::getSelectedSites(
    const PolymorphismSequenceContainer& psc,
    const string& setName,
    bool phase)
{
  auto pscc = MaseTools::getSelectedSites(psc, setName);
  auto maseFileHeader = psc.getComments();
  if (phase)
  {
    for (size_t i = 1; i < MaseTools::getPhase(maseFileHeader, setName); ++i)
    {
      pscc->deleteSite(0);
    }
  }
  auto psci = make_unique<PolymorphismSequenceContainer>(*pscc);
  for (size_t i = 0; i < psc.getNumberOfSequences(); ++i)
  {
    if (psc.isIngroupMember(i))
      psci->setAsIngroupMember(i);
    else
    {
      psci->setAsOutgroupMember(i);
      psci->setGroupId(i, psc.getGroupId(i));
    }
  }
  psci->clearComments();
  return psci;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::getNonCodingSites(
    const PolymorphismSequenceContainer& psc,
    const string& setName)
{
  SiteSelection ss;
  auto maseFileHeader = psc.getComments();
  auto codss = MaseTools::getSiteSet(maseFileHeader, setName);
  for (size_t i = 0; i < psc.getNumberOfSites(); ++i)
  {
    if (find(codss.begin(), codss.end(), i) == codss.end())
      ss.push_back(i);
  }
  auto sc = SiteContainerTools::getSelectedSites(psc, ss);
  auto psci = make_unique<PolymorphismSequenceContainer>(*sc);
  for (size_t i = 0; i < psc.getNumberOfSequences(); ++i)
  {
    if (psc.isIngroupMember(i))
      psci->setAsIngroupMember(i);
    else
    {
      psci->setAsOutgroupMember(i);
      psci->setGroupId(i, psc.getGroupId(i));
    }
  }
  return psci;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::getOnePosition(
    const PolymorphismSequenceContainer& psc,
    const string& setName,
    size_t pos)
{
  auto maseFileHeader = psc.getComments();
  size_t start;
  try
  {
    start = MaseTools::getPhase(maseFileHeader, setName);
  }
  catch (Exception& e)
  {
    start = 1;
  }
  SiteSelection ss;
  size_t i;
  if (static_cast<int>(pos) - static_cast<int>(start) >= 0)
    i = pos - start;
  else
    i = pos - start + 3;
  while (i < psc.getNumberOfSites())
  {
    ss.push_back(i);
    i += 3;
  }
  auto sc = SiteContainerTools::getSelectedSites(psc, ss);
  auto newpsc = make_unique<PolymorphismSequenceContainer>(*sc);
  for (size_t j = 0; j < psc.getNumberOfSequences(); ++j)
  {
    if (psc.isIngroupMember(j))
      newpsc->setAsIngroupMember(j);
    else
    {
      newpsc->setAsOutgroupMember(j);
      newpsc->setGroupId(i, psc.getGroupId(j));
    }
  }
  return newpsc;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::getIntrons(
  const PolymorphismSequenceContainer& psc,
  const string& setName,
  const GeneticCode& gCode)
{
  auto maseFileHeader = psc.getComments();
  SiteSelection ss;
  auto codss = MaseTools::getSiteSet(maseFileHeader, setName);
  size_t start = MaseTools::getPhase(maseFileHeader, setName);
  size_t first = 0, last = psc.getNumberOfSites();
  // Check if the first codon is AUG
  if (start == 1 &&
      psc.site(codss[0]).getValue(0) == 0 &&
      psc.site(codss[1]).getValue(0) == 3 &&
      psc.site(codss[2]).getValue(0) == 2)
    first = codss[0];
  // Check if the last codon is a STOP one
  int c1 = psc.site(codss[codss.size() - 3]).getValue(0);
  int c2 = psc.site(codss[codss.size() - 2]).getValue(0);
  int c3 = psc.site(codss[codss.size() - 1]).getValue(0);
  if (gCode.isStop(gCode.codonAlphabet().getCodon(c1, c2, c3)))
    last = codss[codss.size() - 1];
  // Keep sites between AUG and STOP
  for (size_t i = first; i < last; i++)
  {
    if (find(codss.begin(), codss.end(), i) == codss.end())
    {
      ss.push_back(i);
    }
  }
  auto sc = SiteContainerTools::getSelectedSites(psc, ss);
  auto psci = make_unique<PolymorphismSequenceContainer>(*sc);
  for (size_t i = 0; i < psc.getNumberOfSequences(); ++i)
  {
    if (psc.isIngroupMember(i))
      psci->setAsIngroupMember(i);
    else
    {
      psci->setAsOutgroupMember(i);
      psci->setGroupId(i, psc.getGroupId(i));
    }
  }
  return psci;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::get5Prime(
    const PolymorphismSequenceContainer& psc,
    const string& setName)
{
  auto maseFileHeader = psc.getComments();
  SiteSelection ss;
  auto codss = MaseTools::getSiteSet(maseFileHeader, setName);
  size_t start = MaseTools::getPhase(maseFileHeader, setName);
  size_t last = 0;
  // Check if the first Codon is AUG
  if (start == 1 &&
      psc.site(codss[0]).getValue(0) == 0 &&
      psc.site(codss[1]).getValue(0) == 3 &&
      psc.site(codss[2]).getValue(0) == 2)
    last = codss[0];
  for (size_t i = 0; i < last; ++i)
  {
    if (find(codss.begin(), codss.end(), i) == codss.end())
    {
      ss.push_back(i);
    }
  }
  auto sc = SiteContainerTools::getSelectedSites(psc, ss);
  auto psci = make_unique<PolymorphismSequenceContainer>(*sc);
  for (size_t i = 0; i < psc.getNumberOfSequences(); ++i)
  {
    if (psc.isIngroupMember(i))
      psci->setAsIngroupMember(i);
    else
    {
      psci->setAsOutgroupMember(i);
      psci->setGroupId(i, psc.getGroupId(i));
    }
  }
  return psci;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::get3Prime(
  const PolymorphismSequenceContainer& psc,
  const string& setName,
  const GeneticCode& gCode)
{
  auto maseFileHeader = psc.getComments();
  SiteSelection ss;
  auto codss = MaseTools::getSiteSet(maseFileHeader, setName);
  size_t first = psc.getNumberOfSites() - 1;
  // Check if the last codon is a STOP one
  int c1 = psc.site(codss[codss.size() - 3]).getValue(0);
  int c2 = psc.site(codss[codss.size() - 2]).getValue(0);
  int c3 = psc.site(codss[codss.size() - 1]).getValue(0);
  if (gCode.isStop(gCode.codonAlphabet().getCodon(c1, c2, c3)))
    first = codss[codss.size() - 1];
  for (size_t i = first; i < psc.getNumberOfSites(); ++i)
  {
    if (find(codss.begin(), codss.end(), i) == codss.end())
    {
      ss.push_back(i);
    }
  }
  auto sc = SiteContainerTools::getSelectedSites(psc, ss);
  auto psci = make_unique<PolymorphismSequenceContainer>(*sc);
  for (size_t i = 0; i < psc.getNumberOfSequences(); ++i)
  {
    if (psc.isIngroupMember(i))
      psci->setAsIngroupMember(i);
    else
    {
      psci->setAsOutgroupMember(i);
      psci->setGroupId(i, psc.getGroupId(i));
    }
  }
  return psci;
}

/******************************************************************************/

string PolymorphismSequenceContainerTools::getIngroupSpeciesName(const PolymorphismSequenceContainer& psc)
{
  string key;
  string speciesName;
  auto maseFileHeader = psc.getComments();
  if (!maseFileHeader.size())
    return speciesName;
  auto groupMap = MaseTools::getAvailableSequenceSelections(maseFileHeader);
  for (auto& mi : groupMap)
  {
    key = mi.first;
    if (key.compare(0, 7, "INGROUP") == 0)
    {
      StringTokenizer sptk(key, "_");
      speciesName = sptk.getToken(1) + " " + sptk.getToken(2);
    }
  }
  return speciesName;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::getSynonymousSites(
    const PolymorphismSequenceContainer& psc,
    const GeneticCode& gCode)
{
  auto psco = make_unique<PolymorphismSequenceContainer>(psc.getSequenceNames(), psc.getAlphabet());
  for (size_t i = 0; i < psc.getNumberOfSites(); ++i)
  {
    const Site& site = psc.site(i);
    if (CodonSiteTools::isSynonymousPolymorphic(site, gCode))
    {
      auto tmpSite = make_unique<Site>(site);
      psco->addSite(tmpSite);
    }
  }
  for (size_t i = 0; i < psc.getNumberOfSequences(); ++i)
  {
    psco->setSequenceCount(i, psc.getSequenceCount(i));
    if (psc.isIngroupMember(i))
      psco->setAsIngroupMember(i);
    else
    {
      psco->setAsOutgroupMember(i);
      psco->setGroupId(i, psc.getGroupId(i));
    }
  }
  return psco;
}

/******************************************************************************/

unique_ptr<PolymorphismSequenceContainer> PolymorphismSequenceContainerTools::getNonSynonymousSites(
    const PolymorphismSequenceContainer& psc,
    const GeneticCode& gCode)
{
  auto psco = make_unique<PolymorphismSequenceContainer>(psc.getSequenceNames(), psc.getAlphabet());
  for (size_t i = 0; i < psc.getNumberOfSites(); ++i)
  {
    const Site& site = psc.site(i);
    if (!CodonSiteTools::isSynonymousPolymorphic(site, gCode))
    {
      auto tmpSite = make_unique<Site>(site);
      psco->addSite(tmpSite);
    }
  }
  for (size_t i = 0; i < psc.getNumberOfSequences(); ++i)
  {
    psco->setSequenceCount(i, psc.getSequenceCount(i));
    if (psc.isIngroupMember(i))
      psco->setAsIngroupMember(i);
    else
    {
      psco->setAsOutgroupMember(i);
      psco->setGroupId(i, psc.getGroupId(i));
    }
  }
  return psco;
}

/******************************************************************************/

