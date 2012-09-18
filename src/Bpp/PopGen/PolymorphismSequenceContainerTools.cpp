//
// File: PolymorphismSequenceContainerTools.cpp
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

#include "PolymorphismSequenceContainerTools.h"

using namespace bpp;
using namespace std;

PolymorphismSequenceContainerTools::~PolymorphismSequenceContainerTools() {}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::read(const std::string& path, const Alphabet* alpha) throw (Exception)
{
  Mase ms;
  string key;
  const OrderedSequenceContainer* seqc = 0;
  try
  {
    seqc = dynamic_cast<OrderedSequenceContainer*>(ms.readSequences(path, alpha ));
  }
  catch (Exception& e)
  {
    if (seqc != 0)
      delete seqc;
    throw e;
  }
  PolymorphismSequenceContainer* psc = new PolymorphismSequenceContainer(*seqc);
  Comments maseFileHeader = seqc->getGeneralComments();
  delete seqc;
  map<string, unsigned int> groupMap = MaseTools::getAvailableSequenceSelections(maseFileHeader);
  for (map<string, unsigned int>::iterator mi = groupMap.begin(); mi != groupMap.end(); mi++)
  {
    key = mi->first;
    if (key.compare(0, 8, "OUTGROUP") == 0)
    {
      SequenceSelection ss;
      try
      {
        ss = MaseTools::getSequenceSet(maseFileHeader, key);
      }
      catch (IOException& ioe)
      {
        delete psc;
        throw ioe;
      }
      for (unsigned int i = 0; i != ss.size(); i++)
      {
        try
        {
          psc->setAsOutgroupMember(ss[i]);
        }
        catch (SequenceNotFoundException& snfe)
        {
          delete psc;
          throw snfe;
        }
      }
    }
  }
  return psc;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::extractIngroup (const PolymorphismSequenceContainer& psc) throw (Exception)
{
  SequenceSelection ss;
  PolymorphismSequenceContainer* psci = dynamic_cast<PolymorphismSequenceContainer*>(psc.clone());
  for (unsigned int i = 0; i < psc.getNumberOfSequences(); i++)
  {
    if (!psc.isIngroupMember(i))
      ss.push_back(i);
  }
  if (ss.size() == psc.getNumberOfSequences())
  {
    delete psci;
    throw Exception("PolymorphismSequenceContainerTools::extractIngroup: no Ingroup sequences found.");
  }
  for (int i = ss.size() - 1; i >= 0; i--)
  {
    psci->deleteSequence(ss[i]);
  }
  return psci;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::extractOutgroup(const PolymorphismSequenceContainer& psc) throw (Exception)
{
  SequenceSelection ss;
  PolymorphismSequenceContainer* psci = dynamic_cast<PolymorphismSequenceContainer*>(psc.clone());
  for (unsigned int i = 0; i < psc.getNumberOfSequences(); i++)
  {
    if (psc.isIngroupMember(i) )
      ss.push_back(i);
  }
  if (ss.size() == psc.getNumberOfSequences())
  {
    delete psci;
    throw Exception("PolymorphismSequenceContainerTools::extractOutgroup: no Outgroup sequences found.");
  }
  for (unsigned int i = ss.size(); i > 0; i--)
  {
    psci->deleteSequence(ss[i - 1]);
  }
  return psci;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::extractGroup(const PolymorphismSequenceContainer& psc, unsigned int group_id) throw (Exception)
{
  SequenceSelection ss;
  PolymorphismSequenceContainer* psci = dynamic_cast<PolymorphismSequenceContainer*>(psc.clone());
  for (unsigned int i = 0; i < psc.getNumberOfSequences(); i++)
  {
    if (psc.getGroupId(i) != group_id)
      ss.push_back(i);
  }
  if (ss.size() == psc.getNumberOfSequences())
  {
    delete psci;
    throw GroupNotFoundException("PolymorphismSequenceContainerTools::extractGroup: group_id not found.", group_id);
  }
  for (unsigned int i = ss.size(); i > 0; i--)
  {
    psci->deleteSequence(ss[i - 1]);
  }
  return psci;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::getSelectedSequences(const PolymorphismSequenceContainer& psc, const SequenceSelection& ss)
{
  PolymorphismSequenceContainer* newpsc = new PolymorphismSequenceContainer(psc.getAlphabet());
  for (unsigned int i = 0; i < ss.size(); i++)
  {
    newpsc->addSequence(psc.getSequence(ss[i]), psc.getSequenceCount(i), false);
    if (psc.isIngroupMember(i))
      newpsc->setAsIngroupMember(i);
    else
    {
      newpsc->setAsOutgroupMember(i);
      newpsc->setGroupId(i, psc.getGroupId(i));
    }
  }
  newpsc->setGeneralComments(psc.getGeneralComments());
  return newpsc;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::sample(const PolymorphismSequenceContainer& psc, unsigned int n, bool replace)
{
  unsigned int nbSeq = psc.getNumberOfSequences();
  vector<unsigned int> v;
  for (unsigned int i = 0; i < nbSeq; ++i)
  {
    v.push_back(i);
  }
  vector<unsigned int> vv(n);
  RandomTools::getSample(v, vv, replace);
  PolymorphismSequenceContainer* newpsc = PolymorphismSequenceContainerTools::getSelectedSequences(psc, vv);
  return newpsc;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::getSitesWithoutGaps (const PolymorphismSequenceContainer& psc)
{
  vector<string> seqNames = psc.getSequencesNames();
  PolymorphismSequenceContainer* noGapCont = new PolymorphismSequenceContainer(psc.getNumberOfSequences(), psc.getAlphabet());
  noGapCont->setSequencesNames(seqNames, false);
  unsigned int nbSeq = psc.getNumberOfSequences();
  for (unsigned int i = 0; i < nbSeq; i++)
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
  NoGapSiteIterator ngsi(psc);
  while (ngsi.hasMoreSites())
    noGapCont->addSite(*ngsi.nextSite());
  return noGapCont;
}

/******************************************************************************/

unsigned int PolymorphismSequenceContainerTools::getNumberOfNonGapSites(const PolymorphismSequenceContainer& psc, bool ingroup) throw (Exception)
{
  unsigned int count = psc.getNumberOfSites();
  PolymorphismSequenceContainer* npsc = NULL;
  SimpleSiteIterator* ssi;
  if (ingroup)
  {
    try
    {
      npsc = extractIngroup(psc);
    }
    catch (Exception& e)
    {
      if (npsc != NULL)
        delete npsc;
      throw e;
    }
    ssi = new SimpleSiteIterator(*npsc);
  }
  else
    ssi = new SimpleSiteIterator(psc);
  while (ssi->hasMoreSites())
    if (SiteTools::hasGap(*ssi->nextSite()))
      count--;
  delete ssi;
  return count;
}

/******************************************************************************/

unsigned int PolymorphismSequenceContainerTools::getNumberOfCompleteSites(const PolymorphismSequenceContainer& psc, bool ingroup) throw (Exception)
{
  unsigned int count = psc.getNumberOfSites();
  PolymorphismSequenceContainer* npsc = NULL;
  SimpleSiteIterator* ssi;
  if (ingroup)
  {
    try
    {
      npsc = extractIngroup(psc);
    }
    catch (Exception& e)
    {
      if (npsc != NULL)
        delete npsc;
      throw e;
    }
    ssi = new SimpleSiteIterator(*npsc);
  }
  else
    ssi = new SimpleSiteIterator(psc);
  while (ssi->hasMoreSites())
    if (!SiteTools::isComplete(*ssi->nextSite()))
      count--;
  delete ssi;
  return count;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::getCompleteSites (const PolymorphismSequenceContainer& psc)
{
  vector<string> seqNames = psc.getSequencesNames();
  PolymorphismSequenceContainer* complete = new PolymorphismSequenceContainer(psc.getNumberOfSequences(), psc.getAlphabet());
  complete->setSequencesNames(seqNames, false);
  unsigned int nbSeq = psc.getNumberOfSequences();
  for (unsigned int i = 0; i < nbSeq; i++)
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
  CompleteSiteIterator csi(psc);
  while (csi.hasMoreSites())
    complete->addSite(*csi.nextSite());
  return complete;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::excludeFlankingGap(const PolymorphismSequenceContainer& psc)
{
  PolymorphismSequenceContainer* psci = dynamic_cast<PolymorphismSequenceContainer*>(psc.clone());
  while (SiteTools::hasGap(psci->getSite(0)))
    psci->deleteSite(0);
  unsigned int i = 0;
  unsigned int n = psci->getNumberOfSites();
  while (SiteTools::hasGap(psci->getSite(n - i - 1)))
  {
    psci->deleteSite(n - i - 1);
    i++;
  }
  return psci;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::getSelectedSites(const PolymorphismSequenceContainer& psc, const std::string& setName, bool phase)
{
  SiteContainer* pscc = MaseTools::getSelectedSites(psc, setName);
  Comments maseFileHeader = psc.getGeneralComments();
  if (phase)
  {
    for (unsigned int i = 1; i < MaseTools::getPhase(maseFileHeader, setName); i++)
    {
      pscc->deleteSite(0);
    }
  }
  PolymorphismSequenceContainer* psci = new PolymorphismSequenceContainer(*pscc);
  for (unsigned int i = 0; i < psc.getNumberOfSequences(); i++)
  {
    if (psc.isIngroupMember(i))
      psci->setAsIngroupMember(i);
    else
    {
      psci->setAsOutgroupMember(i);
      psci->setGroupId(i, psc.getGroupId(i));
    }
  }
  psci->deleteGeneralComments();
  delete pscc;
  return psci;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::getNonCodingSites(const PolymorphismSequenceContainer& psc, const std::string& setName)
{
  SiteSelection ss;
  Comments maseFileHeader = psc.getGeneralComments();
  SiteSelection codss = MaseTools::getSiteSet(maseFileHeader, setName);
  for (unsigned int i = 0; i < psc.getNumberOfSites(); i++)
  {
    if (find(codss.begin(), codss.end(), i) == codss.end())
      ss.push_back(i);
  }
  const SiteContainer* sc = SiteContainerTools::getSelectedSites(psc, ss);
  PolymorphismSequenceContainer* psci = new PolymorphismSequenceContainer(*sc);
  for (unsigned int i = 0; i < psc.getNumberOfSequences(); i++)
  {
    if (psc.isIngroupMember(i))
      psci->setAsIngroupMember(i);
    else
    {
      psci->setAsOutgroupMember(i);
      psci->setGroupId(i, psc.getGroupId(i));
    }
  }
  delete sc;
  return psci;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::getOnePosition(const PolymorphismSequenceContainer& psc, const std::string& setName, unsigned int pos)
{
  Comments maseFileHeader = psc.getGeneralComments();
  unsigned int start;
  try
  {
    start = MaseTools::getPhase(maseFileHeader, setName);
  }
  catch (Exception& e)
  {
    start = 1;
  }
  SiteSelection ss;
  unsigned int i;
  if ((int)pos - (int)start >= 0)
    i = pos - start;
  else
    i = pos - start + 3;
  while (i < psc.getNumberOfSites())
  {
    ss.push_back(i);
    i += 3;
  }
  const SiteContainer* sc = SiteContainerTools::getSelectedSites(psc, ss);
  PolymorphismSequenceContainer* newpsc = new PolymorphismSequenceContainer(*sc);
  for (unsigned int j = 0; j < psc.getNumberOfSequences(); j++)
  {
    if (psc.isIngroupMember(j))
      newpsc->setAsIngroupMember(j);
    else
    {
      newpsc->setAsOutgroupMember(j);
      newpsc->setGroupId(i, psc.getGroupId(j));
    }
  }
  delete sc;
  return newpsc;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::getIntrons(const PolymorphismSequenceContainer& psc, const std::string& setName, const CodonAlphabet* ca )
{
  Comments maseFileHeader = psc.getGeneralComments();
  SiteSelection ss;
  SiteSelection codss = MaseTools::getSiteSet(maseFileHeader, setName);
  unsigned int start;
  try
  {
    start = MaseTools::getPhase(maseFileHeader, setName);
  }
  catch (Exception& e)
  {
    throw e;
  }

  unsigned int first = 0, last = psc.getNumberOfSites();
  // Check if the first codon is AUG
  if (start == 1 &&
      psc.getSite(codss[0]).getValue(0) == 0 &&
      psc.getSite(codss[1]).getValue(0) == 3 &&
      psc.getSite(codss[2]).getValue(0) == 2)
    first = codss[0];
  // Check if the last codon is a STOP one
  int c1 = psc.getSite(codss[codss.size() - 3]).getValue(0);
  int c2 = psc.getSite(codss[codss.size() - 2]).getValue(0);
  int c3 = psc.getSite(codss[codss.size() - 1]).getValue(0);
  if (ca->isStop(ca->getCodon(c1, c2, c3)))
    last = codss[codss.size() - 1];
  // Keep sites between AUG and STOP
  for (unsigned int i = first; i < last; i++)
  {
    if (find(codss.begin(), codss.end(), i) == codss.end())
    {
      ss.push_back(i);
    }
  }
  const SiteContainer* sc = SiteContainerTools::getSelectedSites(psc, ss);
  PolymorphismSequenceContainer* psci = new PolymorphismSequenceContainer(*sc);
  for (unsigned int i = 0; i < psc.getNumberOfSequences(); i++)
  {
    if (psc.isIngroupMember(i))
      psci->setAsIngroupMember(i);
    else
    {
      psci->setAsOutgroupMember(i);
      psci->setGroupId(i, psc.getGroupId(i));
    }
  }
  delete sc;
  return psci;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::get5Prime(const PolymorphismSequenceContainer& psc, const std::string& setName)
{
  Comments maseFileHeader = psc.getGeneralComments();
  SiteSelection ss;
  SiteSelection codss = MaseTools::getSiteSet(maseFileHeader, setName);
  unsigned int start = MaseTools::getPhase(maseFileHeader, setName);
  unsigned int last = 0;
  // Check if the first Codon is AUG
  if (start == 1 &&
      psc.getSite(codss[0]).getValue(0) == 0 &&
      psc.getSite(codss[1]).getValue(0) == 3 &&
      psc.getSite(codss[2]).getValue(0) == 2)
    last = codss[0];
  for (unsigned int i = 0; i < last; i++)
  {
    if (find(codss.begin(), codss.end(), i) == codss.end())
    {
      ss.push_back(i);
    }
  }
  const SiteContainer* sc = SiteContainerTools::getSelectedSites(psc, ss);
  PolymorphismSequenceContainer* psci = new PolymorphismSequenceContainer(*sc);
  for (unsigned int i = 0; i < psc.getNumberOfSequences(); i++)
  {
    if (psc.isIngroupMember(i))
      psci->setAsIngroupMember(i);
    else
    {
      psci->setAsOutgroupMember(i);
      psci->setGroupId(i, psc.getGroupId(i));
    }
  }
  delete sc;
  return psci;
}

/******************************************************************************/

PolymorphismSequenceContainer* PolymorphismSequenceContainerTools::get3Prime(const PolymorphismSequenceContainer& psc, const std::string& setName, const CodonAlphabet* ca )
{
  Comments maseFileHeader = psc.getGeneralComments();
  SiteSelection ss;
  SiteSelection codss = MaseTools::getSiteSet(maseFileHeader, setName);
  unsigned int first = psc.getNumberOfSites() - 1;
  // Check if the last codon is a STOP one
  int c1 = psc.getSite(codss[codss.size() - 3]).getValue(0);
  int c2 = psc.getSite(codss[codss.size() - 2]).getValue(0);
  int c3 = psc.getSite(codss[codss.size() - 1]).getValue(0);
  if (ca->isStop(ca->getCodon(c1, c2, c3)))
    first = codss[codss.size() - 1];
  for (unsigned int i = first; i < psc.getNumberOfSites(); i++)
  {
    if (find(codss.begin(), codss.end(), i) == codss.end())
    {
      ss.push_back(i);
    }
  }
  const SiteContainer* sc = SiteContainerTools::getSelectedSites(psc, ss);
  PolymorphismSequenceContainer* psci = new PolymorphismSequenceContainer(*sc);
  for (unsigned int i = 0; i < psc.getNumberOfSequences(); i++)
  {
    if (psc.isIngroupMember(i))
      psci->setAsIngroupMember(i);
    else
    {
      psci->setAsOutgroupMember(i);
      psci->setGroupId(i, psc.getGroupId(i));
    }
  }
  delete sc;
  return psci;
}

/******************************************************************************/

string PolymorphismSequenceContainerTools::getIngroupSpeciesName(const PolymorphismSequenceContainer& psc)
{
  string key;
  string speciesName;
  Comments maseFileHeader = psc.getGeneralComments();
  if (!maseFileHeader.size())
    return speciesName;
  map<string, unsigned int> groupMap = MaseTools::getAvailableSequenceSelections(maseFileHeader);
  for (map<string, unsigned int>::iterator mi = groupMap.begin(); mi != groupMap.end(); mi++)
  {
    key = mi->first;
    if (key.compare(0, 7, "INGROUP") == 0)
    {
      StringTokenizer* sptk = new StringTokenizer(key, "_");
      speciesName = sptk->getToken(1) + " " + sptk->getToken(2);
    }
  }
  return speciesName;
}

/******************************************************************************/
