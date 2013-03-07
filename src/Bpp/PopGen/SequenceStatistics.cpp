//
// File SequenceStatistics.cpp
// Authors: Eric Bazin
//          Sylvain Gailard
//          Khalid Belkhir
//          Benoit Nabholz
// Created on: Wed Aug 04 2004
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

#include "SequenceStatistics.h" // class's header file
#include "PolymorphismSequenceContainerTools.h"
#include "PolymorphismSequenceContainer.h"

// From the STL:
#include <ctype.h>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

// From SeqLib:
#include <Bpp/Seq/Site.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/StringSequenceTools.h>
#include <Bpp/Seq/CodonSiteTools.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/StandardCodonAlphabet.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/VectorExceptions.h>

using namespace bpp;

// ******************************************************************************
// Basic statistics
// ******************************************************************************

size_t SequenceStatistics::polymorphicSiteNumber(const PolymorphismSequenceContainer& psc, bool gapflag, bool ignoreUnknown)
{
  size_t S = 0;
  const Site* site = 0;
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (!SiteTools::isConstant(*site, ignoreUnknown))
    {
      S++;
    }
  }
  delete si;
  return S;
}

size_t SequenceStatistics::parsimonyInformativeSiteNumber(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  size_t S = 0;
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (SiteTools::isParsimonyInformativeSite(*site))
    {
      S++;
    }
  }
  delete si;
  return S;
}

size_t SequenceStatistics::countSingleton(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  size_t nus = 0;
  const Site* site = 0;
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    nus += getSingletonNumber_(*site);
  }
  delete si;
  return nus;
}

size_t SequenceStatistics::tripletNumber(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  int S = 0;
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (SiteTools::isTriplet(*site))
    {
      S++;
    }
  }

  delete si;
  return S;
}

size_t SequenceStatistics::totNumberMutations(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  size_t tnm = 0;
  const Site* site = 0;
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    tnm += getMutationNumber_(*site);
  }
  delete si;
  return tnm;
}

size_t SequenceStatistics::totMutationsExternalBranchs(
  const PolymorphismSequenceContainer& ing,
  const PolymorphismSequenceContainer& outg) throw (Exception)
{
  if (ing.getNumberOfSites() != outg.getNumberOfSites())
    throw Exception("ing and outg must have the same size");
  size_t nmuts = 0;
  const Site* site_in = 0;
  const Site* site_out = 0;
  ConstSiteIterator* si = 0;
  ConstSiteIterator* so = 0;
  si = new SimpleSiteContainerIterator(ing);
  so = new SimpleSiteContainerIterator(outg);
  while (si->hasMoreSites())
  {
    site_in = si->nextSite();
    site_out = so->nextSite();
    // use fully resolved sites
    if (SiteTools::isComplete(*site_in) &&  SiteTools::isComplete(*site_out))
      nmuts += getDerivedSingletonNumber_(*site_in, *site_out);                                                                   // singletons that are not in outgroup
  }
  delete si;
  delete so;
  return nmuts;
}

double SequenceStatistics::heterozygosity(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  ConstSiteIterator* si = 0;
  const Site* site = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  double S = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    S += SiteTools::heterozygosity(*site);
  }
  delete si;
  return S;
}

double SequenceStatistics::squaredHeterozygosity(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  ConstSiteIterator* si = 0;
  const Site* site = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  double S = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    double h = SiteTools::heterozygosity(*site);
    S += h * h;
  }
  delete si;
  return S;
}

// ******************************************************************************
// GC statistics
// ******************************************************************************

double SequenceStatistics::gcContent(const PolymorphismSequenceContainer& psc)
{
  map<int, double> freqs;
  SequenceContainerTools::getFrequencies(psc, freqs);
  const Alphabet* alpha = psc.getAlphabet();
  return (freqs[alpha->charToInt("C")] + freqs[alpha->charToInt("G")]) / (freqs[alpha->charToInt("A")] + freqs[alpha->charToInt("C")] + freqs[alpha->charToInt("G")] + freqs[alpha->charToInt("T")]);
}

std::vector<size_t> SequenceStatistics::gcPolymorphism(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  size_t nbMut = 0;
  size_t nbGC = 0;
  const size_t nbSeq = psc.getNumberOfSequences();
  vector<size_t> vect(2);
  const Site* site = 0;
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new NoGapSiteContainerIterator(psc);
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (!SiteTools::isConstant(*site))
    {
      long double freqGC = SymbolListTools::getGCContent(*site);
      /*
       * Sylvain Gaillard 15/03/2010: realy unclear ...
       *   freqGC is always in [0,1] then why testing it ?
       *   why casting double into size_t ?
       *   is that method used by someone ?
       */
      if (freqGC > 0 && freqGC < 1)
      {
        nbMut += static_cast<size_t>(nbSeq);
        long double adGC = freqGC * nbSeq;
        nbGC += static_cast<size_t>(adGC);
      }
    }
  }
  vect[0] = nbMut;
  vect[1] = nbGC;
  delete si;
  return vect;
}

// ******************************************************************************
// Diversity statistics
// ******************************************************************************

double SequenceStatistics::watterson75(const PolymorphismSequenceContainer& psc, bool gapflag, bool ignoreUnknown)
{
  double ThetaW;
  size_t n = psc.getNumberOfSequences();
  size_t S = polymorphicSiteNumber(psc, gapflag, ignoreUnknown);
  map<string, double> values = getUsefullValues_(n);
  ThetaW = (double) S / values["a1"];
  return ThetaW;
}

double SequenceStatistics::tajima83(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  size_t alphabet_size = (psc.getAlphabet())->getSize();
  const Site* site = 0;
  ConstSiteIterator* si = 0;
  double value2 = 0.;
  if (gapflag)
    si = new CompleteSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (!SiteTools::isConstant(*site))
    {
      double value = 0.;
      map<int, size_t> count;
      SymbolListTools::getCounts(*site, count);
      map<int, size_t> tmp_k;
      size_t tmp_n = 0;
      for (map<int, size_t>::iterator it = count.begin(); it != count.end(); it++)
      {
        if (it->first >= 0 && it->first < static_cast<int>(alphabet_size))
        {
          tmp_k[it->first] = it->second * (it->second - 1);
          tmp_n += it->second;
        }
      }
      if (tmp_n == 0 || tmp_n == 1)
        continue;
      for (map<int, size_t>::iterator it = tmp_k.begin(); it != tmp_k.end(); it++)
      {
        value += static_cast<double>(it->second) / static_cast<double>(tmp_n * (tmp_n - 1));
      }
      value2 += 1. - value;
    }
  }
  delete si;
  return value2;
}

double SequenceStatistics::FayWu2000(const PolymorphismSequenceContainer& psc, const Sequence& ancestralSites)
{
  if (psc.getNumberOfSites() != ancestralSites.size())
    throw Exception("SequenceStatistics::FayWu2000: ancestralSites and psc don't have the same size!!!'" );

  const Sequence& tmps = psc.getSequence(0);

  size_t alphabet_size = (psc.getAlphabet())->getSize();
  double value = 0.;
  for (size_t i = 0; i < psc.getNumberOfSites(); i++)
  {
    const Site& site = psc.getSite(i);
    string ancB = ancestralSites.getChar(i);
    int ancV = ancestralSites.getValue(i);

    if (!SiteTools::isConstant(site) || tmps.getChar(i) != ancB)
    {
      if (ancV < 0)
        continue;

      map<int, size_t> count;
      SymbolListTools::getCounts(site, count);
      map<int, size_t> tmp_k;
      size_t tmp_n = 0;
      for (map<int, size_t>::iterator it = count.begin(); it != count.end(); it++)
      {
        if (it->first >= 0 && it->first < static_cast<int>(alphabet_size))
        {
          /* if derived allele */
          if (it->first != ancV)
          {
            tmp_k[it->first] = 2 * it->second * it->second;
          }
          tmp_n += it->second;
        }
      }
      if (tmp_n == 0 || tmp_n == 1)
        continue;
      for (map<int, size_t>::iterator it = tmp_k.begin(); it != tmp_k.end(); it++)
      {
        value += static_cast<double>(it->second) / static_cast<double>(tmp_n * (tmp_n - 1));
      }
    }
  }
  return value;
}

size_t SequenceStatistics::DVK(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  /*
   * Sylvain Gaillard 17/03/2010:
   * This implementation uses unneeded SequenceContainer recopy and works on
   * string. It needs to be improved.
   */
  PolymorphismSequenceContainer* sc = 0;
  if (gapflag)
    sc = PolymorphismSequenceContainerTools::getSitesWithoutGaps(psc);
  else
    sc = new PolymorphismSequenceContainer(psc);
  // int K = 0;
  vector<string> pscvector;
  pscvector.push_back(sc->toString(0));
  // K++;
  for (size_t i = 1; i < sc->getNumberOfSequences(); i++)
  {
    bool uniq = true;
    string query = sc->toString(i);
    for (vector<string>::iterator it = pscvector.begin(); it != pscvector.end(); it++)
    {
      if (query.compare(*it) == 0)
      {
        uniq = false;
        break;
      }
    }
    if (uniq)
    {
      // K++;
      pscvector.push_back(query);
    }
  }
  delete sc;
  // return K;
  return pscvector.size();
}

double SequenceStatistics::DVH(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  /*
   * Sylvain Gaillard 17/03/2010:
   * This implementation uses unneeded SequenceContainer recopy and works on
   * string. It needs to be improved.
   */
  PolymorphismSequenceContainer* sc = 0;
  if (gapflag)
    sc = PolymorphismSequenceContainerTools::getSitesWithoutGaps(psc);
  else
    sc = new PolymorphismSequenceContainer(psc);
  double H = 0.;
  size_t nbSeq;
  vector<string> pscvector;
  vector<size_t> effvector;
  pscvector.push_back(sc->toString(0));
  effvector.push_back(sc->getSequenceCount(0));
  nbSeq = sc->getSequenceCount(0);
  for (size_t i = 1; i < sc->getNumberOfSequences(); i++)
  {
    nbSeq += sc->getSequenceCount(i);
    bool uniq = true;
    string query = sc->toString(i);
    for (size_t j = 0; j < pscvector.size(); j++)
    {
      if (query.compare(pscvector[j]) == 0)
      {
        effvector[j] += sc->getSequenceCount(i);
        uniq = false;
        break;
      }
    }
    if (uniq)
    {
      pscvector.push_back(query);
      effvector.push_back(sc->getSequenceCount(i));
    }
  }
  for (size_t i = 0; i < effvector.size(); i++)
  {
    H -= (static_cast<double>(effvector[i]) / static_cast<double>(nbSeq)) * ( static_cast<double>(effvector[i]) / static_cast<double>(nbSeq));
  }
  H += 1.;
  delete sc;
  return H;
}

size_t SequenceStatistics::getNumberOfTransitions(const PolymorphismSequenceContainer& psc)
{
  size_t nbT = 0;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    // if (SiteTools::isConstant(*site) || SiteTools::isTriplet(*site)) continue;
    if (SiteTools::getNumberOfDistinctCharacters(*site) != 2)
      continue;
    vector<int> seq = site->getContent();
    int state1 = seq[0];
    int state2 = seq[0];
    for (size_t i = 1; i < seq.size(); i++)
    {
      if (state1 != seq[i])
      {
        state2 = seq[i];
        break;
      }
    }
    if (((state1 == 0 && state2 == 2) || (state1 == 2 && state2 == 0)) ||
        ((state1 == 1 && state2 == 3) || (state1 == 3 && state2 == 1)))
    {
      nbT++;
    }
  }
  delete si;
  return nbT;
}

size_t SequenceStatistics::getNumberOfTransversions(const PolymorphismSequenceContainer& psc)
{
  size_t nbTv = 0;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    // if (SiteTools::isConstant(*site) || SiteTools::isTriplet(*site)) continue;
    if (SiteTools::getNumberOfDistinctCharacters(*site) != 2)
      continue;
    vector<int> seq = site->getContent();
    int state1 = seq[0];
    int state2 = seq[0];
    for (size_t i = 1; i < seq.size(); i++)
    {
      if (state1 != seq[i])
      {
        state2 = seq[i];
        break;
      }
    }
    if (!(((state1 == 0 && state2 == 2) || (state1 == 2 && state2 == 0)) ||
          ((state1 == 1 && state2 == 3) || (state1 == 3 && state2 == 1))))
    {
      nbTv++;
    }
  }
  delete si;
  return nbTv;
}

double SequenceStatistics::getTransitionsTransversionsRatio(const PolymorphismSequenceContainer& psc) throw (Exception)
{
  // return (double) getNumberOfTransitions(psc)/getNumberOfTransversions(psc);
  size_t nbT = 0;
  size_t nbTv = 0;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  vector< int > state(2);
  while (si->hasMoreSites())
  {
    map<int, size_t> count;
    site = si->nextSite();
    SymbolListTools::getCounts(*site, count);
    if (count.size() != 2)
      continue;
    int i = 0;
    for (map<int, size_t>::iterator it = count.begin(); it != count.end(); it++)
    {
      state[i] = it->first;
      i++;
    }
    if (((state[0] == 0 && state[1] == 2) || (state[0] == 2 && state[1] == 0)) ||
        ((state[0] == 1 && state[1] == 3) || (state[0] == 3 && state[1] == 1)))
    {
      nbT++; // transitions
    }
    else
    {
      nbTv++; // transversion
    }
  }
  delete si;
  if (nbTv == 0)
    throw ZeroDivisionException("SequenceStatistics::getTransitionsTransversionsRatio.");
  return static_cast<double>(nbT) /  static_cast<double>(nbTv);
}

// ******************************************************************************
// Synonymous and non-synonymous polymorphism
// ******************************************************************************

size_t SequenceStatistics::stopCodonSiteNumber(const PolymorphismSequenceContainer& psc, bool gapflag)
{
  /*
   * Sylvain Gaillard 17/03/2010
   * What if the Alphabet is not a codon alphabet?
   */
  ConstSiteIterator* si = 0;
  if (gapflag)
    si = new NoGapSiteContainerIterator(psc);
  else
    si = new SimpleSiteContainerIterator(psc);
  size_t S = 0;
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (CodonSiteTools::hasStop(*site))
      S++;
  }
  delete si;
  return S;
}

size_t SequenceStatistics::monoSitePolymorphicCodonNumber(const PolymorphismSequenceContainer& psc, bool stopflag, bool gapflag)
{
  ConstSiteIterator* si = 0;
  if (stopflag)
    si = new CompleteSiteContainerIterator(psc);
  else
  {
    if (gapflag)
      si = new NoGapSiteContainerIterator(psc);
    else
      si = new SimpleSiteContainerIterator(psc);
  }
  size_t S = 0;
  const Site* site;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (CodonSiteTools::isMonoSitePolymorphic(*site))
      S++;
  }
  delete si;
  return S;
}

size_t SequenceStatistics::synonymousPolymorphicCodonNumber(const PolymorphismSequenceContainer& psc, const GeneticCode& gc)
{
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  size_t S = 0;
  const Site* site;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    if (CodonSiteTools::isSynonymousPolymorphic(*site, gc))
      S++;
  }
  delete si;
  return S;
}

double SequenceStatistics::watterson75Synonymous(const PolymorphismSequenceContainer& psc, const GeneticCode& gc)
{
  double ThetaW = 0.;
  size_t n = psc.getNumberOfSequences();
  size_t S = synonymousSubstitutionsNumber(psc, gc);
  map<string, double> values = getUsefullValues_(n);
  ThetaW = static_cast<double>(S) / values["a1"];
  return ThetaW;
}

double SequenceStatistics::watterson75NonSynonymous(const PolymorphismSequenceContainer& psc, const GeneticCode& gc)
{
  double ThetaW;
  size_t n = psc.getNumberOfSequences();
  size_t S = nonSynonymousSubstitutionsNumber(psc, gc);
  map<string, double> values = getUsefullValues_(n);
  ThetaW = static_cast<double>(S) / values["a1"];
  return ThetaW;
}

double SequenceStatistics::piSynonymous(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, bool minchange)
{
  double S = 0.;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    S += CodonSiteTools::piSynonymous(*site, gc, minchange);
  }
  delete si;
  return S;
}

double SequenceStatistics::piNonSynonymous(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, bool minchange)
{
  double S = 0.;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    S += CodonSiteTools::piNonSynonymous(*site, gc, minchange);
  }
  delete si;
  return S;
}

double SequenceStatistics::meanSynonymousSitesNumber(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, double ratio)
{
  double S = 0.;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    S += CodonSiteTools::meanNumberOfSynonymousPositions(*site, gc, ratio);
  }
  delete si;
  return S;
}

double SequenceStatistics::meanNonSynonymousSitesNumber(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, double ratio)
{
  double S = 0.;
  int n = 0;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    n = n + 3;
    S += CodonSiteTools::meanNumberOfSynonymousPositions(*site, gc, ratio);
  }
  delete si;
  return static_cast<double>(n - S);
}

size_t SequenceStatistics::synonymousSubstitutionsNumber(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, double freqmin)
{
  size_t St = 0, Sns = 0;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    St += CodonSiteTools::numberOfSubsitutions(*site, freqmin);
    Sns += CodonSiteTools::numberOfNonSynonymousSubstitutions(*site, gc, freqmin);
  }
  delete si;
  return St - Sns;
}

size_t SequenceStatistics::nonSynonymousSubstitutionsNumber(const PolymorphismSequenceContainer& psc, const GeneticCode& gc, double freqmin)
{
  size_t Sns = 0;
  ConstSiteIterator* si = new CompleteSiteContainerIterator(psc);
  const Site* site = 0;
  while (si->hasMoreSites())
  {
    site = si->nextSite();
    Sns += CodonSiteTools::numberOfNonSynonymousSubstitutions(*site, gc, freqmin);
  }
  delete si;
  return Sns;
}

vector<size_t> SequenceStatistics::fixedDifferences(const PolymorphismSequenceContainer& pscin, const PolymorphismSequenceContainer& pscout, PolymorphismSequenceContainer& psccons, const GeneticCode& gc)
{
  ConstSiteIterator* siIn = new CompleteSiteContainerIterator(pscin);
  ConstSiteIterator* siOut = new CompleteSiteContainerIterator(pscout);
  ConstSiteIterator* siCons = new CompleteSiteContainerIterator(psccons);
  const Site* siteIn = 0;
  const Site* siteOut = 0;
  const Site* siteCons = 0;
  size_t NfixS = 0;
  size_t NfixA = 0;
  while (siIn->hasMoreSites())
  {
    siteIn = siIn->nextSite();
    siteOut = siOut->nextSite();
    siteCons = siCons->nextSite();
    vector<size_t> v = CodonSiteTools::fixedDifferences(*siteIn, *siteOut, siteCons->getValue(0), siteCons->getValue(1), gc);
    NfixS += v[0];
    NfixA += v[1];
  }
  vector<size_t> v(2);
  v[0] = NfixS;
  v[1] = NfixA;
  delete siIn;
  delete siOut;
  delete siCons;
  return v;
}

vector<size_t> SequenceStatistics::MKtable(const PolymorphismSequenceContainer& ingroup, const PolymorphismSequenceContainer& outgroup, const GeneticCode& gc, double freqmin)
{
  PolymorphismSequenceContainer psctot(ingroup);
  for (size_t i = 0; i < outgroup.getNumberOfSequences(); i++)
  {
    psctot.addSequence(outgroup.getSequence(i));
    psctot.setAsOutgroupMember(i + ingroup.getNumberOfSequences());
  }
  const PolymorphismSequenceContainer* psccomplet = PolymorphismSequenceContainerTools::getCompleteSites(psctot);
  const PolymorphismSequenceContainer* pscin = PolymorphismSequenceContainerTools::extractIngroup(*psccomplet);
  const PolymorphismSequenceContainer* pscout = PolymorphismSequenceContainerTools::extractOutgroup(*psccomplet);
  const Sequence* consensusIn = SiteContainerTools::getConsensus(*pscin, "consensusIn");
  const Sequence* consensusOut = SiteContainerTools::getConsensus(*pscout, "consensusOut");
  PolymorphismSequenceContainer* consensus = new PolymorphismSequenceContainer(ingroup.getAlphabet());
  consensus->addSequence(*consensusIn);
  consensus->addSequence(*consensusOut);
  vector<size_t> u = SequenceStatistics::fixedDifferences(*pscin, *pscout, *consensus, gc);
  vector<size_t> v(4);
  v[0] = SequenceStatistics::nonSynonymousSubstitutionsNumber(*pscin, gc, freqmin);
  v[1] = SequenceStatistics::synonymousSubstitutionsNumber(*pscin, gc, freqmin);
  v[2] = u[1];
  v[3] = u[0];
  delete consensus;
  if (psccomplet)
  {
    delete psccomplet;
  }
  if (pscin)
  {
    delete pscin;
  }
  if (pscout)
  {
    delete pscout;
  }
  if (consensusIn)
  {
    delete consensusIn;
  }
  if (consensusOut)
  {
    delete consensusOut;
  }
  return v;
}

double SequenceStatistics::neutralityIndex(const PolymorphismSequenceContainer& ingroup, const PolymorphismSequenceContainer& outgroup, const GeneticCode& gc, double freqmin)
{
  vector<size_t> v = SequenceStatistics::MKtable(ingroup, outgroup, gc, freqmin);
  if (v[1] != 0 && v[2] != 0)
    return static_cast<double>(v[0] * v[3]) / static_cast<double>(v[1] * v[2]);
  else
    return -1;
}

// ******************************************************************************
// Statistical tests
// ******************************************************************************

double SequenceStatistics::tajimaDSS(const PolymorphismSequenceContainer& psc, bool gapflag) throw (ZeroDivisionException)
{
  double S = static_cast<double>(polymorphicSiteNumber(psc, gapflag));
  if (!S)
    throw ZeroDivisionException("S should not be null");
  double tajima = tajima83(psc, gapflag);
  double watterson = watterson75(psc, gapflag);
  size_t n = psc.getNumberOfSequences();
  map<string, double> values = getUsefullValues_(n);
  // if (S == 0)
  //  cout << "ARG S == 0" << endl;
  return (tajima - watterson) / sqrt((values["e1"] * S) + (values["e2"] * S * (S - 1)));
}

double SequenceStatistics::tajimaDTNM(const PolymorphismSequenceContainer& psc, bool gapflag) throw (ZeroDivisionException)
{
  double eta =  static_cast<double>(totNumberMutations(psc, gapflag));
  if (!eta)
    throw ZeroDivisionException("eta should not be null");
  double tajima = tajima83(psc, gapflag);
  size_t n = psc.getNumberOfSequences();
  map<string, double> values = getUsefullValues_(n);
  double eta_a1 = static_cast<double>(eta) / values["a1"];
  return (tajima - eta_a1) / sqrt((values["e1"] * eta) + (values["e2"] * eta * (eta - 1)));
}

double SequenceStatistics::fuliD(const PolymorphismSequenceContainer& ingroup, const PolymorphismSequenceContainer& outgroup, bool original) throw (ZeroDivisionException)
{
  size_t n = ingroup.getNumberOfSequences();
  map<string, double> values = getUsefullValues_(n);
  double vD = getVD_(n, values["a1"], values["a2"], values["cn"]);
  double uD = getUD_(values["a1"], vD);
  double eta = static_cast<double>(totNumberMutations(ingroup));
  if (eta == 0.)
    throw ZeroDivisionException("eta should not be null");
  double etae = 0.;
  if (original)
    etae = static_cast<double>(countSingleton(outgroup));
  else
    etae = static_cast<double>(totMutationsExternalBranchs(ingroup, outgroup));  // added by Khalid 13/07/2005
  return (eta - (values["a1"] * etae)) / sqrt((uD * eta) + (vD * eta * eta));
}

double SequenceStatistics::fuliDstar(const PolymorphismSequenceContainer& group) throw (ZeroDivisionException)
{
  size_t n = group.getNumberOfSequences();
  double nn = static_cast<double>(n);
  double _n = nn / (nn - 1.);
  map<string, double> values = getUsefullValues_(n);
  double vDs = getVDstar_(n, values["a1"], values["a2"], values["dn"]);
  double uDs = getUDstar_(n, values["a1"], vDs);
  double eta = static_cast<double>(totNumberMutations(group));
  if (eta == 0.)
    throw ZeroDivisionException("eta should not be null");
  double etas = static_cast<double>(countSingleton(group));

  // Fu & Li 1993
  return ((_n * eta) - (values["a1"] * etas)) / sqrt(uDs * eta + vDs * eta * eta);

  // Simonsen et al. 1995
  /*
     return ((eta / values["a1"]) - (etas * ((n - 1) / n))) / sqrt(uDs * eta + vDs * eta * eta);
   */
}

double SequenceStatistics::fuliF(const PolymorphismSequenceContainer& ingroup, const PolymorphismSequenceContainer& outgroup, bool original) throw (ZeroDivisionException)
{
  size_t n = ingroup.getNumberOfSequences();
  double nn = static_cast<double>(n);
  map<string, double> values = getUsefullValues_(n);
  double pi = tajima83(ingroup, true);
  double vF = (values["cn"] + values["b2"] - 2. / (nn - 1.)) / (pow(values["a1"], 2) + values["a2"]);
  double uF = ((1. + values["b1"] - (4. * ((nn + 1.) / ((nn - 1.) * (nn - 1.)))) * (values["a1n"] - (2. * nn) / (nn + 1.))) / values["a1"]) - vF;
  double eta = static_cast<double>(totNumberMutations(ingroup));
  if (eta == 0.)
    throw ZeroDivisionException("eta should not be null");
  double etae = 0.;
  if (original)
    etae = static_cast<double>(countSingleton(outgroup));
  else
    etae = static_cast<double>(totMutationsExternalBranchs(ingroup, outgroup));  // added by Khalid 13/07/2005
  return (pi - etae) / sqrt(uF * eta + vF * eta * eta);
}

double SequenceStatistics::fuliFstar(const PolymorphismSequenceContainer& group) throw (ZeroDivisionException)
{
  double n = static_cast<double>(group.getNumberOfSequences());
  map<string, double> values = getUsefullValues_(group.getNumberOfSequences());
  double pi = tajima83(group, true);

  // Fu & Li 1993
  //  double vFs = (values["dn"] + values["b2"] - (2. / (nn - 1.)) * (4. * values["a2"] - 6. + 8. / nn)) / (pow(values["a1"], 2) + values["a2"]);
  //  double uFs = (((nn / (nn - 1.)) + values["b1"] - (4. / (nn * (nn - 1.))) + 2. * ((nn + 1.) / (pow((nn - 1.), 2))) * (values["a1n"] - 2. * nn / (nn + 1.))) / values["a1"]) - vFs;

  // Simonsen et al. 1995
  double vFs = (((2 * n * n * n + 110 * n * n - 255 * n + 153) / (9 * n * n * (n - 1))) + ((2 * (n - 1) * values["a1"]) / (n * n)) - 8 * values["a2"] / n) / (pow(values["a1"], 2) + values["a2"]);
  double uFs = (((4 * n * n + 19 * n + 3 - 12 * (n + 1) * values["a1n"]) / (3 * n * (n - 1))) / values["a1"]) - vFs;
  double eta = static_cast<double>(totNumberMutations(group));
  if (eta == 0.)
    throw ZeroDivisionException("eta should not be null");
  double etas = static_cast<double>(countSingleton(group));
  // Fu & Li 1993
  // Simonsen et al. 1995
  return (pi - ((n - 1.) / n * etas)) / sqrt(uFs * eta + vFs * eta * eta);
}

double SequenceStatistics::FstHudson92(const PolymorphismSequenceContainer& psc, size_t id1, size_t id2)
{
  vector<double> vdiff;
  double piIntra1, piIntra2, meanPiIntra, piInter, Fst;

  PolymorphismSequenceContainer* Pop1 = PolymorphismSequenceContainerTools::extractGroup(psc, id1);
  PolymorphismSequenceContainer* Pop2 = PolymorphismSequenceContainerTools::extractGroup(psc, id2);

  piIntra1 = SequenceStatistics::tajima83(*Pop1, false);
  piIntra2 = SequenceStatistics::tajima83(*Pop2, false);

  meanPiIntra = (piIntra1 + piIntra2) / 2;

  double n = 0;
  for (size_t i = 0; i < Pop1->getNumberOfSequences(); i++)
  {
    const Sequence& s1 = Pop1->getSequence(i);
    for (size_t j = 0; j < Pop2->getNumberOfSequences(); j++)
    {
      n++;
      const Sequence& s2 = Pop2->getSequence(j);
      vdiff.push_back(SiteContainerTools::computeSimilarity(s1, s2, true, "no gap", true));
    }
  }
  piInter = (VectorTools::sum(vdiff) / n) * static_cast<double>(psc.getNumberOfSites());


  Fst = 1.0 - meanPiIntra / piInter;

  delete Pop1;
  delete Pop2;

  return Fst;
}

// ******************************************************************************
// Linkage disequilibrium statistics
// ******************************************************************************

/**********************/
/* Preliminary method */
/**********************/

PolymorphismSequenceContainer* SequenceStatistics::generateLDContainer(const PolymorphismSequenceContainer& psc, bool keepsingleton, double freqmin)
{
  SiteSelection ss;
  // Extract polymorphic site with only two alleles
  for (size_t i = 0; i < psc.getNumberOfSites(); i++)
  {
    if (keepsingleton)
    {
      if (SiteTools::isComplete(psc.getSite(i)) && !SiteTools::isConstant(psc.getSite(i)) && !SiteTools::isTriplet(psc.getSite(i)))
      {
        ss.push_back(i);
      }
    }
    else
    {
      if (SiteTools::isComplete(psc.getSite(i)) && !SiteTools::isConstant(psc.getSite(i)) && !SiteTools::isTriplet(psc.getSite(i)) && !SiteTools::hasSingleton(psc.getSite(i)))
      {
        ss.push_back(i);
      }
    }
  }

  const SiteContainer* sc = SiteContainerTools::getSelectedSites(psc, ss);
  Alphabet* alpha = new DNA(); // Sylvain Gaillard 17/03/2010: What if psc's Alphabet is not DNA
  PolymorphismSequenceContainer* ldpsc = new PolymorphismSequenceContainer(sc->getNumberOfSequences(), alpha);
  // Assign 1 to the more frequent and 0 to the less frequent alleles
  for (size_t i = 0; i < sc->getNumberOfSites(); i++)
  {
    const Site& site = sc->getSite(i);
    Site siteclone(site);
    bool deletesite = false;
    map<int, double> freqs;
    SymbolListTools::getFrequencies(siteclone, freqs);
    int first = 0;
    for (map<int, double>::iterator it = freqs.begin(); it != freqs.end(); it++)
    {
      if (it->second >= 0.5)
        first = it->first;
    }
    for (size_t j = 0; j < sc->getNumberOfSequences(); j++)
    {
      if (freqs[site.getValue(j)] >= 0.5 && site.getValue(j) == first)
      {
        if (freqs[site.getValue(j)] <= 1 - freqmin)
        {
          siteclone.setElement(j, 1);
          first = site.getValue(j);
        }
        else
          deletesite = true;
      }
      else
      {
        if (freqs[site.getValue(j)] >= freqmin)
          siteclone.setElement(j, 0);
        else
          deletesite = true;
      }
    }
    if (!deletesite)
      ldpsc->addSite(siteclone);
  }
  delete alpha;
  return ldpsc;
}

/*************************************/
/* Pairwise LD and distance measures */
/*************************************/

Vdouble SequenceStatistics::pairwiseDistances1(const PolymorphismSequenceContainer& psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  // get Positions with sites of interest
  SiteSelection ss;
  for (size_t i = 0; i < psc.getNumberOfSites(); i++)
  {
    if (keepsingleton)
    {
      if (SiteTools::isComplete(psc.getSite(i)) && !SiteTools::isConstant(psc.getSite(i)) && !SiteTools::isTriplet(psc.getSite(i)))
      {
        const Site& site = psc.getSite(i);
        bool deletesite = false;
        map<int, double> freqs;
        SymbolListTools::getFrequencies(site, freqs);
        for (int j = 0; j < static_cast<int>(site.getAlphabet()->getSize()); j++)
        {
          if (freqs[j] >= 1 - freqmin)
            deletesite = true;
        }
        if (!deletesite)
          ss.push_back(i);
      }
    }
    else
    {
      if (SiteTools::isComplete(psc.getSite(i)) && !SiteTools::isConstant(psc.getSite(i)) && !SiteTools::isTriplet(psc.getSite(i)) && !SiteTools::hasSingleton(psc.getSite(i)))
      {
        ss.push_back(i);
        const Site& site = psc.getSite(i);
        bool deletesite = false;
        map<int, double> freqs;
        SymbolListTools::getFrequencies(site, freqs);
        for (int j = 0; j < static_cast<int>(site.getAlphabet()->getSize()); j++)
        {
          if (freqs[j] >= 1 - freqmin)
            deletesite = true;
        }
        if (!deletesite)
          ss.push_back(i);
      }
    }
  }
  // compute pairwise distances
  if (ss.size() < 2)
    throw DimensionException("SequenceStatistics::pairwiseDistances1 : less than 2 sites are available", ss.size(), 2);
  Vdouble dist;
  for (size_t i = 0; i < ss.size() - 1; i++)
  {
    for (size_t j = i + 1; j < ss.size(); j++)
    {
      dist.push_back(static_cast<double>(ss[j] - ss[i]));
    }
  }
  return dist;
}

Vdouble SequenceStatistics::pairwiseDistances2(const PolymorphismSequenceContainer& psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  SiteSelection ss;
  for (size_t i = 0; i < psc.getNumberOfSites(); i++)
  {
    if (keepsingleton)
    {
      if (SiteTools::isComplete(psc.getSite(i)) && !SiteTools::isConstant(psc.getSite(i)) && !SiteTools::isTriplet(psc.getSite(i)))
      {
        const Site& site = psc.getSite(i);
        bool deletesite = false;
        map<int, double> freqs;
        SymbolListTools::getFrequencies(site, freqs);
        for (int j = 0; j < static_cast<int>(site.getAlphabet()->getSize()); j++)
        {
          if (freqs[j] >= 1 - freqmin)
            deletesite = true;
        }
        if (!deletesite)
          ss.push_back(i);
      }
    }
    else
    {
      if (SiteTools::isComplete(psc.getSite(i)) && !SiteTools::isConstant(psc.getSite(i)) && !SiteTools::isTriplet(psc.getSite(i)) && !SiteTools::hasSingleton(psc.getSite(i)))
      {
        ss.push_back(i);
        const Site& site = psc.getSite(i);
        bool deletesite = false;
        map<int, double> freqs;
        SymbolListTools::getFrequencies(site, freqs);
        for (int j = 0; j < static_cast<int>(site.getAlphabet()->getSize()); j++)
        {
          if (freqs[j] >= 1 - freqmin)
            deletesite = true;
        }
        if (!deletesite)
          ss.push_back(i);
      }
    }
  }
  size_t n = ss.size();
  if (n < 2)
    throw DimensionException("SequenceStatistics::pairwiseDistances1 : less than 2 sites are available", ss.size(), 2);
  Vdouble distance(n * (n - 1) / 2, 0);
  size_t nbsite = psc.getNumberOfSites();
  for (size_t k = 0; k < psc.getNumberOfSequences(); k++)
  {
    const Sequence& seq = psc.getSequence(k);
    SiteSelection gap, newss = ss;
    Vdouble dist;
    for (size_t i = 0; i < nbsite; i++)
    {
      if (seq.getValue(i) == -1)
        gap.push_back(i);
    }
    // Site positions are re-numbered to take gaps into account
    for (size_t i = 0; i < gap.size(); i++)
    {
      for (size_t j = 0; j < ss.size(); j++)
      {
        if (ss[j] > gap[i])
          newss[j]--;
      }
    }
    for (size_t i = 0; i < n - 1; i++)
    {
      for (size_t j = i + 1; j < n; j++)
      {
        dist.push_back(static_cast<double>(newss[j] - newss[i]));
      }
    }
    distance += dist;
  }
  distance = distance / static_cast<double>(psc.getNumberOfSequences());
  return distance;
}

Vdouble SequenceStatistics::pairwiseD(const PolymorphismSequenceContainer& psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  PolymorphismSequenceContainer* newpsc = SequenceStatistics::generateLDContainer(psc, keepsingleton, freqmin);
  Vdouble D;
  size_t nbsite = newpsc->getNumberOfSites();
  size_t nbseq = newpsc->getNumberOfSequences();
  if (nbsite < 2)
    throw DimensionException("SequenceStatistics::pairwiseD: less than two sites are available", nbsite, 2);
  if (nbseq < 2)
    throw DimensionException("SequenceStatistics::pairwiseD: less than two sequences are available", nbseq, 2);
  for (size_t i = 0; i < nbsite - 1; i++)
  {
    for (size_t j = i + 1; j < nbsite; j++)
    {
      double haplo = 0;
      const Site& site1 = newpsc->getSite(i);
      const Site& site2 = newpsc->getSite(j);
      map<int, double> freq1;
      map<int, double> freq2;
      SymbolListTools::getFrequencies(site1, freq1);
      SymbolListTools::getFrequencies(site2, freq2);
      for (size_t k = 0; k < nbseq; k++)
      {
        if (site1.getValue(k) + site2.getValue(k) == 2)
          haplo++;
      }
      haplo = haplo / static_cast<double>(nbseq);
      D.push_back(std::abs(haplo - freq1[1] * freq2[1]));
    }
  }
  return D;
}

Vdouble SequenceStatistics::pairwiseDprime(const PolymorphismSequenceContainer& psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  PolymorphismSequenceContainer* newpsc = SequenceStatistics::generateLDContainer(psc, keepsingleton, freqmin);
  Vdouble Dprime;
  size_t nbsite = newpsc->getNumberOfSites();
  size_t nbseq = newpsc->getNumberOfSequences();
  if (nbsite < 2)
    throw DimensionException("SequenceStatistics::pairwiseD: less than two sites are available", nbsite, 2);
  if (nbseq < 2)
    throw DimensionException("SequenceStatistics::pairwiseD: less than two sequences are available", nbseq, 2);
  for (size_t i = 0; i < nbsite - 1; i++)
  {
    for (size_t j = i + 1; j < nbsite; j++)
    {
      double haplo = 0;
      const Site& site1 = newpsc->getSite(i);
      const Site& site2 = newpsc->getSite(j);
      map<int, double> freq1;
      map<int, double> freq2;
      SymbolListTools::getFrequencies(site1, freq1);
      SymbolListTools::getFrequencies(site2, freq2);
      for (size_t k = 0; k < nbseq; k++)
      {
        if (site1.getValue(k) + site2.getValue(k) == 2)
          haplo++;
      }
      haplo = haplo / static_cast<double>(nbseq);
      double d, D = (haplo - freq1[1] * freq2[1]);
      if (D > 0)
      {
        if (freq1[1] * freq2[0] <= freq1[0] * freq2[1])
        {
          d = std::abs(D) / (freq1[1] * freq2[0]);
        }
        else
        {
          d = std::abs(D) / (freq1[0] * freq2[1]);
        }
      }
      else
      {
        if (freq1[1] * freq2[1] <= freq1[0] * freq2[0])
        {
          d = std::abs(D) / (freq1[1] * freq2[1]);
        }
        else
        {
          d = std::abs(D) / (freq1[0] * freq2[0]);
        }
      }
      Dprime.push_back(d);
    }
  }
  return Dprime;
}

Vdouble SequenceStatistics::pairwiseR2(const PolymorphismSequenceContainer& psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  PolymorphismSequenceContainer* newpsc = SequenceStatistics::generateLDContainer(psc, keepsingleton, freqmin);
  Vdouble R2;
  size_t nbsite = newpsc->getNumberOfSites();
  size_t nbseq = newpsc->getNumberOfSequences();
  if (nbsite < 2)
    throw DimensionException("SequenceStatistics::pairwiseD: less than two sites are available", nbsite, 2);
  if (nbseq < 2)
    throw DimensionException("SequenceStatistics::pairwiseD: less than two sequences are available", nbseq, 2);
  for (size_t i = 0; i < nbsite - 1; i++)
  {
    for (size_t j = i + 1; j < nbsite; j++)
    {
      double haplo = 0;
      const Site& site1 = newpsc->getSite(i);
      const Site& site2 = newpsc->getSite(j);
      map<int, double> freq1;
      map<int, double> freq2;
      SymbolListTools::getFrequencies(site1, freq1);
      SymbolListTools::getFrequencies(site2, freq2);
      for (size_t k = 0; k < nbseq; k++)
      {
        if (site1.getValue(k) + site2.getValue(k) == 2)
          haplo++;
      }
      haplo = haplo / static_cast<double>(nbseq);
      double r = ((haplo - freq1[1] * freq2[1]) * (haplo - freq1[1] * freq2[1])) / (freq1[0] * freq1[1] * freq2[0] * freq2[1]);
      R2.push_back(r);
    }
  }
  return R2;
}

/***********************************/
/* Global LD and distance measures */
/***********************************/

double SequenceStatistics::meanD(const PolymorphismSequenceContainer& psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  Vdouble D = pairwiseD(psc, keepsingleton, freqmin);
  return VectorTools::mean<double, double>(D);
}

double SequenceStatistics::meanDprime(const PolymorphismSequenceContainer& psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble Dprime = pairwiseDprime(psc, keepsingleton, freqmin);
    return VectorTools::mean<double, double>(Dprime);
  }
  catch (DimensionException& e)
  {
    throw e;
  }
}

double SequenceStatistics::meanR2(const PolymorphismSequenceContainer& psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble R2 = SequenceStatistics::pairwiseR2(psc, keepsingleton, freqmin);
    return VectorTools::mean<double, double>(R2);
  }
  catch (DimensionException& e)
  {
    throw e;
  }
}

double SequenceStatistics::meanDistance1(const PolymorphismSequenceContainer& psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble dist = pairwiseDistances1(psc, keepsingleton, freqmin);
    return VectorTools::mean<double, double>(dist);
  }
  catch (DimensionException& e)
  {
    throw e;
  }
}

double SequenceStatistics::meanDistance2(const PolymorphismSequenceContainer& psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble dist = pairwiseDistances2(psc, keepsingleton, freqmin);
    return VectorTools::mean<double, double>(dist);
  }
  catch (DimensionException& e)
  {
    throw e;
  }
}

/**********************/
/* Regression methods */
/**********************/

double SequenceStatistics::originRegressionD(const PolymorphismSequenceContainer& psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble D = pairwiseD(psc, keepsingleton, freqmin) - 1;
    Vdouble dist;
    if (distance1)
      dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else
      dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    return VectorTools::sum(D * dist) / VectorTools::sum(dist * dist);
  }
  catch (DimensionException& e)
  {
    throw e;
  }
}

double SequenceStatistics::originRegressionDprime(const PolymorphismSequenceContainer& psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble Dprime = pairwiseDprime(psc, keepsingleton, freqmin) - 1;
    Vdouble dist;
    if (distance1)
      dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else
      dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    return VectorTools::sum(Dprime * dist) / VectorTools::sum(dist * dist);
  }
  catch (DimensionException& e)
  {
    throw e;
  }
}

double SequenceStatistics::originRegressionR2(const PolymorphismSequenceContainer& psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble R2 = pairwiseR2(psc, keepsingleton, freqmin) - 1;
    Vdouble dist;
    if (distance1)
      dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else
      dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    return VectorTools::sum(R2 * dist) / VectorTools::sum(dist * dist);
  }
  catch (DimensionException& e)
  {
    throw e;
  }
}

Vdouble SequenceStatistics::linearRegressionD(const PolymorphismSequenceContainer& psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble D = pairwiseD(psc, keepsingleton, freqmin);
    Vdouble dist;
    Vdouble reg(2);
    if (distance1)
      dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else
      dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    reg[0] = VectorTools::cov<double, double>(dist, D) / VectorTools::var<double, double>(dist);
    reg[1] = VectorTools::mean<double, double>(D) - reg[0] * VectorTools::mean<double, double>(dist);
    return reg;
  }
  catch (DimensionException& e)
  {
    throw e;
  }
}

Vdouble SequenceStatistics::linearRegressionDprime(const PolymorphismSequenceContainer& psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble Dprime = pairwiseDprime(psc, keepsingleton, freqmin);
    Vdouble dist;
    Vdouble reg(2);
    if (distance1)
      dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else
      dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    reg[0] = VectorTools::cov<double, double>(dist, Dprime) / VectorTools::var<double, double>(dist);
    reg[1] = VectorTools::mean<double, double>(Dprime) - reg[0] * VectorTools::mean<double, double>(dist);
    return reg;
  }
  catch (DimensionException& e)
  {
    throw e;
  }
}

Vdouble SequenceStatistics::linearRegressionR2(const PolymorphismSequenceContainer& psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble R2 = pairwiseR2(psc, keepsingleton, freqmin);
    Vdouble dist;
    Vdouble reg(2);
    if (distance1)
      dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else
      dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    reg[0] = VectorTools::cov<double, double>(dist, R2) / VectorTools::var<double, double>(dist);
    reg[1] = VectorTools::mean<double, double>(R2) - reg[0] * VectorTools::mean<double, double>(dist);
    return reg;
  }
  catch (DimensionException& e)
  {
    throw e;
  }
}

double SequenceStatistics::inverseRegressionR2(const PolymorphismSequenceContainer& psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble R2 = pairwiseR2(psc, keepsingleton, freqmin);
    Vdouble unit(R2.size(), 1);
    Vdouble R2transformed = unit / R2 - 1;
    Vdouble dist;
    if (distance1)
      dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else
      dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    return VectorTools::sum(R2transformed * dist) / VectorTools::sum(dist * dist);
  }
  catch (DimensionException& e)
  {
    throw e;
  }
}

/**********************/
/*   Hudson method    */
/**********************/

double SequenceStatistics::hudson87(const PolymorphismSequenceContainer& psc, double precision, double cinf, double csup)
{
  double left = leftHandHudson_(psc);
  size_t n = psc.getNumberOfSequences();
  double dif = 1;
  double c1 = cinf;
  double c2 = csup;
  if (SequenceStatistics::polymorphicSiteNumber(psc) < 2)
    return -1;
  if (rightHandHudson_(c1, n) < left)
    return cinf;
  if (rightHandHudson_(c2, n) > left)
    return csup;
  while (dif > precision)
  {
    if (rightHandHudson_((c1 + c2) / 2, n) > left)
      c1 = (c1 + c2) / 2;
    else
      c2 = (c1 + c2) / 2;
    dif = std::abs(2 * (c1 - c2) / (c1 + c2));
  }
  return (c1 + c2) / 2;
}

/*****************/
/* Tests methods */
/*****************/

void SequenceStatistics::testUsefullValues(std::ostream& s, size_t n)
{
  map<string, double> v = getUsefullValues_(n);
  double vD = getVD_(n, v["a1"], v["a2"], v["cn"]);
  double uD = getUD_(v["a1"], vD);
  double vDs = getVDstar_(n, v["a1"], v["a2"], v["dn"]);
  double uDs = getUDstar_(n, v["a1"], vDs);

  s << n << "\t";
  s << v["a1"] << "\t";
  s << v["a2"] << "\t";
  s << v["a1n"] << "\t";
  s << v["b1"] << "\t";
  s << v["b2"] << "\t";
  s << v["c1"] << "\t";
  s << v["c2"] << "\t";
  s << v["cn"] << "\t";
  s << v["dn"] << "\t";
  s << v["e1"] << "\t";
  s << v["e2"] << "\t";
  s << uD << "\t";
  s << vD << "\t";
  s << uDs << "\t";
  s << vDs << endl;
}

// ******************************************************************************
// Private methods
// ******************************************************************************

size_t SequenceStatistics::getMutationNumber_(const Site& site)
{
  size_t tmp_count = 0;
  map<int, size_t> states_count;
  SymbolListTools::getCounts(site, states_count);

  for (map<int, size_t>::iterator it = states_count.begin(); it != states_count.end(); it++)
  {
    if (it->first >= 0)
      tmp_count++;
  }
  if (tmp_count > 0)
    tmp_count--;
  return tmp_count;
}

size_t SequenceStatistics::getSingletonNumber_(const Site& site)
{
  size_t nus = 0;
  map<int, size_t> states_count;
  SymbolListTools::getCounts(site, states_count);
  for (map<int, size_t>::iterator it = states_count.begin(); it != states_count.end(); it++)
  {
    if (it->second == 1)
      nus++;
  }
  return nus;
}

size_t SequenceStatistics::getDerivedSingletonNumber_(const Site& site_in, const Site& site_out)
{
  size_t nus = 0;
  map<int, size_t> states_count;
  map<int, size_t> outgroup_states_count;
  SymbolListTools::getCounts(site_in, states_count);
  SymbolListTools::getCounts(site_out, outgroup_states_count);
  // if there is more than one variant in the outgroup we will not be able to recover the ancestral state
  if (outgroup_states_count.size() == 1)
  {
    for (map<int, size_t>::iterator it = states_count.begin(); it != states_count.end(); it++)
    {
      if (it->second == 1)
      {
        if (outgroup_states_count.find(it->first) == outgroup_states_count.end())
          nus++;
      }
    }
  }
  return nus;
}

std::map<std::string, double> SequenceStatistics::getUsefullValues_(size_t n)
{
  double nn = static_cast<double>(n);
  map<string, double> values;
  values["a1"] = 0.;
  values["a2"] = 0.;
  values["a1n"] = 0.;
  values["b1"] = 0.;
  values["b2"] = 0.;
  values["c1"] = 0.;
  values["c2"] = 0.;
  values["cn"] = 0.;
  values["dn"] = 0.;
  values["e1"] = 0.;
  values["e2"] = 0.;
  if (n > 1)
  {
    for (double i = 1; i < nn; i++)
    {
      values["a1"] += 1. / i;
      values["a2"] += 1. / (i * i);
    }
    values["a1n"] = values["a1"] + (1. / nn);
    values["b1"] = (nn + 1.) / (3. * (nn - 1.));
    values["b2"] = 2. * ((nn * nn) + nn + 3.) / (9. * nn * (nn - 1.));
    values["c1"] = values["b1"] - (1. / values["a1"]);
    values["c2"] = values["b2"] - ((nn + 2.) / (values["a1"] * nn)) + (values["a2"] / (values["a1"] * values["a1"]));
    if (n == 2)
    {
      values["cn"] = 1.;
      values["dn"] = 2.;
    }
    else
    {
      values["cn"] = 2. * ((nn * values["a1"]) - (2. * (nn - 1.))) / ((nn - 1.) * (nn - 2.));
      values["dn"] =
        values["cn"]
        + ((nn - 2.) / ((nn - 1.) * (nn - 1.)))
        + (2. / (nn - 1.))
        * ((3. / 2.) - (((2. * values["a1n"]) - 3.) / (nn - 2.)) - (1. / nn));
    }
    values["e1"] = values["c1"] / values["a1"];
    values["e2"] = values["c2"] / ((values["a1"] * values["a1"]) + values["a2"]);
  }
  return values;
}

double SequenceStatistics::getVD_(size_t n, double a1, double a2, double cn)
{
  double nn = static_cast<double>(n);
  if (n < 3)
    return 0.;
  double vD = 1. + ((a1 * a1) / (a2 + (a1 * a1))) * (cn - ((nn + 1.) / (nn - 1.)));
  return vD;
}

double SequenceStatistics::getUD_(double a1, double vD)
{
  return a1 - 1. - vD;
}

double SequenceStatistics::getVDstar_(size_t n, double a1, double a2, double dn)
{
  double denom = (a1 * a1) + a2;
  if (n < 3 || denom == 0.)
    return 0.;
  double nn = static_cast<double>(n);
  double nnn = nn / (nn - 1.);
  // Fu & Li 1993
  double vDs = (
    (nnn * nnn * a2)
    + (a1 * a1 * dn)
    - (2. * (nn * a1 * (a1 + 1)) / ((nn - 1.) * (nn - 1.)))
    )
               /
               denom;
  // Simonsen et al. 1995
  /*
     double vDs = (
      (values["a2"] / pow(values["a1"], 2))
      - (2./nn) * (1. + 1./values["a1"] - values["a1"] + values["a1"]/nn)
      - 1./(nn*nn)
      )
      /
      (pow(values["a1"], 2) + values["a2"]);
   */
  return vDs;
}

double SequenceStatistics::getUDstar_(size_t n, double a1, double vDs)
{
  if (n < 3)
    return 0.;
  double nn = static_cast<double>(n);
  double nnn = nn / (nn - 1.);
  // Fu & Li 1993
  double uDs = (nnn * (a1 - nnn)) - vDs;
  // Simonsen et al. 1995
  /*
     double uDs = (((nn - 1.)/nn - 1./values["a1"]) / values["a1"]) - vDs;
   */
  return uDs;
}

double SequenceStatistics::leftHandHudson_(const PolymorphismSequenceContainer& psc)
{
  PolymorphismSequenceContainer* newpsc = PolymorphismSequenceContainerTools::getCompleteSites(psc);
  size_t nbseq = newpsc->getNumberOfSequences();
  double S1 = 0;
  double S2 = 0;
  for (size_t i = 0; i < nbseq - 1; i++)
  {
    for (size_t j = i + 1; j < nbseq; j++)
    {
      SequenceSelection ss(2);
      ss[0] = i;
      ss[1] = j;
      PolymorphismSequenceContainer* psc2 = PolymorphismSequenceContainerTools::getSelectedSequences(*newpsc, ss);
      S1 += SequenceStatistics::watterson75(*psc2, true);
      S2 += SequenceStatistics::watterson75(*psc2, true) * SequenceStatistics::watterson75(*psc2, true);
      delete psc2;
    }
  }
  double Sk = (2 * S2 - pow(2 * S1 / static_cast<double>(nbseq), 2.)) / pow(nbseq, 2.);
  double H = SequenceStatistics::heterozygosity(*newpsc);
  double H2 = SequenceStatistics::squaredHeterozygosity(*newpsc);
  delete newpsc;
  return static_cast<double>(Sk - H + H2) / pow(H * static_cast<double>(nbseq) / static_cast<double>(nbseq - 1), 2.);
}

double SequenceStatistics::rightHandHudson_(double c, size_t n)
{
  double nn = static_cast<double>(n);
  return 1. / (97. * pow(c, 2.) * pow(nn, 3.)) * ((nn - 1.) * (97. * (c * (4. + (c - 2. * nn) * nn) + (-2. * (7. + c) + 4. * nn + (c - 1.) * pow(nn, 2.)) * log((18. + c * (13. + c)) / 18.)) + sqrt(97.) * (110. + nn * (49. * nn - 52.) + c * (2. + nn * (15. * nn - 8.))) * log(-1. + (72. + 26. * c) / (36. + 13. * c - c * sqrt(97.)))));
}

