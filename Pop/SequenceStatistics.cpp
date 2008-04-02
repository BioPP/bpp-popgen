//
// File SequenceStatistics.cpp
// Author : Eric Bazin
//          Sylvain Gailard
//          Khalid Belkhir
// Last modification : Wednesday February 27 2008
//

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
#include <Seq/Site.h>
#include <Seq/SiteTools.h>
#include <Seq/StringSequenceTools.h>
#include <Seq/CodonSiteTools.h>
#include <Seq/DNA.h>
#include <Seq/StandardCodonAlphabet.h>
#include <Seq/StandardGeneticCode.h>

// from NumCalc
#include <NumCalc/VectorTools.h>
#include <NumCalc/VectorExceptions.h>

using namespace bpp;

SequenceStatistics::~SequenceStatistics() {}

//******************************************************************************************************************
//Basic statistics
//******************************************************************************************************************

unsigned int SequenceStatistics::polymorphicSiteNumber(const PolymorphismSequenceContainer & psc, bool gapflag)
{
  unsigned int S=0;
  const Site *site;
  SiteIterator * si = NULL;
  if (gapflag)
    si = new CompleteSiteIterator(psc);
  else
    si = new SimpleSiteIterator(psc);
  while ( si->hasMoreSites() ) {
    site=si->nextSite();
    if ( !SiteTools::isConstant(*site) ) {
      S++;
    }
  }
  delete si;
  return S;
}

unsigned int SequenceStatistics::parsimonyInformativeSiteNumber(const PolymorphismSequenceContainer & psc, bool gapflag)
{
  SiteIterator *si;
  if(gapflag) si = new CompleteSiteIterator(psc);
  else si = new SimpleSiteIterator(psc);
  unsigned int S=0;
  const Site *site;
  while ( si->hasMoreSites() ) {
    site=si->nextSite();
    if (SiteTools::isParsimonyInformativeSite(*site)) {
      S++;
    }
  }
  delete si;
  return S;
}

unsigned int SequenceStatistics::countSingleton(const PolymorphismSequenceContainer & psc, bool gapflag)
{
  unsigned int nus = 0;
  const Site * site;
  SiteIterator * si = NULL;
  if (gapflag)
    si = new CompleteSiteIterator(psc);
  else
    si = new SimpleSiteIterator(psc);
  while (si->hasMoreSites()) {
    site = si->nextSite();
    nus += _getSingletonNumber(* site);
  }
  delete si;
  return nus;
}

unsigned int SequenceStatistics::tripletNumber(const PolymorphismSequenceContainer & psc, bool gapflag)
{
  SiteIterator *si;
  if(gapflag) si = new CompleteSiteIterator(psc);
  else si = new SimpleSiteIterator(psc);
  int S=0;
  const Site *site;
  while ( si->hasMoreSites() ) {
    site=si->nextSite();
    if ( SiteTools::isTriplet(*site) ) {
      S++;
    }
  }

  delete si;
  return S;
}

unsigned int SequenceStatistics::totNumberMutations(const PolymorphismSequenceContainer & psc, bool gapflag)
{
  unsigned int tnm = 0;
  const Site * site;
  SiteIterator * si = NULL;
  if (gapflag)
    si = new CompleteSiteIterator(psc);
  else
    si = new SimpleSiteIterator(psc);
  while (si->hasMoreSites()) {
    site = si->nextSite();
    tnm += _getMutationNumber(* site);
  }
  delete si;
  return tnm;
}

unsigned int SequenceStatistics::totMutationsExternalBranchs(
    const PolymorphismSequenceContainer & ing,
    const PolymorphismSequenceContainer outg)
{
  unsigned int nmuts = 0;
  const Site * site_in;
  const Site * site_out;
  SiteIterator * si = NULL;
  SiteIterator * so = NULL;
  si = new SimpleSiteIterator(ing);
  so = new SimpleSiteIterator(outg);
  while (si->hasMoreSites()) {
    site_in = si->nextSite();
    site_out= so->nextSite();
    //use fully resolved sites
    if ( SiteTools::isComplete(*site_in) &&  SiteTools::isComplete(*site_out) )
      nmuts += _getDerivedSingletonNumber(* site_in, *site_out);//singletons that are not in outgroup
  }
  delete si;
  delete so;
  return nmuts;
}

double SequenceStatistics::heterozygosity(const PolymorphismSequenceContainer & psc, bool gapflag)
{
  SiteIterator *si;
  const Site * site;
  if(gapflag) si = new CompleteSiteIterator(psc);
  else si = new SimpleSiteIterator(psc);
  double S=0;
  while ( si->hasMoreSites() ) {
    site=si->nextSite();
    S+=SiteTools::heterozygosity(*site);
  }
  delete si;
  return S;
}

double SequenceStatistics::squaredHeterozygosity(const PolymorphismSequenceContainer & psc, bool gapflag)
{
  SiteIterator *si;
  const Site * site;
  if(gapflag) si = new CompleteSiteIterator(psc);
  else si = new SimpleSiteIterator(psc);
  double S=0;
  while ( si->hasMoreSites() ) {
    site=si->nextSite();
    double h = SiteTools::heterozygosity(*site);
    S+= h*h;
  }
  delete si;
  return S;
}

//******************************************************************************************************************
//GC statistics
//******************************************************************************************************************

double SequenceStatistics::gcContent(const PolymorphismSequenceContainer & psc)
{
  SiteContainer* sc = new VectorSiteContainer(psc);
  map<int, double> freqs = SequenceContainerTools::getFrequencies(*sc);
  delete sc;
  return (freqs[1] + freqs[2]) / (freqs[0] + freqs[1] + freqs[2] + freqs[3]);
}

vector<unsigned int> SequenceStatistics::gcPolymorphism(const PolymorphismSequenceContainer & psc, bool stopflag)
{
  unsigned int nbMut = 0;
  unsigned int nbGC = 0;
  const unsigned int nbSeq = psc.getNumberOfSequences();
  vector<unsigned int> vect(2);
  const Site * site;
  SiteIterator * si = NULL;
  if(stopflag) si = new CompleteSiteIterator(psc);
  else si = new NoGapSiteIterator(psc);
  while (si->hasMoreSites()) {
    site = si->nextSite();
    if(!SiteTools::isConstant(*site)) {
      double freqGC = StringSequenceTools::getGCcontent(site->toString(),0, nbSeq);
      if(freqGC > 0 && freqGC < 1) {
        nbMut += (unsigned int) nbSeq;
        double adGC = freqGC*nbSeq;
        nbGC += (unsigned int)adGC;
      }
    }
  }
  vect[0]=nbMut;
  vect[1]=nbGC;
  delete si;
  return vect;
}

//******************************************************************************************************************
//Diversity statistics
//******************************************************************************************************************

double SequenceStatistics::watterson75(const PolymorphismSequenceContainer & psc, bool gapflag)
{
  double ThetaW;
  unsigned int n = psc.getNumberOfSequences();
  unsigned int S = polymorphicSiteNumber(psc, gapflag);
  map<string, double> values = _getUsefullValues(n);
  ThetaW = (double) S / values["a1"];
  return ThetaW;
}

double SequenceStatistics::tajima83(const PolymorphismSequenceContainer & psc, bool gapflag)
{
  unsigned int alphabet_size = (psc.getAlphabet())->getSize();
  const Site * site;
  SiteIterator *si;
  double value2 = 0.;
  if (gapflag)
    si = new CompleteSiteIterator(psc);
  else
    si = new SimpleSiteIterator(psc);
  while (si->hasMoreSites()) {
    site = si->nextSite();
    if (! SiteTools::isConstant(* site)) {
      double value = 0.;
      map<int, unsigned int> count = SymbolListTools::getCounts(* site);
      map<int, unsigned int> tmp_k;
      unsigned int tmp_n = 0;
      for (map<int, unsigned int>::iterator it = count.begin() ; it != count.end() ; it++)
        if (it->first >= 0 && it->first < (int) alphabet_size) {
          tmp_k[it->first] = it->second * (it->second - 1);
          tmp_n += it->second;
        }
      if (tmp_n == 0 || tmp_n == 1) continue;
      for(map<int, unsigned int>::iterator it = tmp_k.begin() ; it != tmp_k.end() ; it++)
        value += (double) it->second / (tmp_n * (tmp_n - 1));
      value2 += 1. - value;
    }
  }
  delete si;
  return value2;
}

unsigned int SequenceStatistics::DVK( const PolymorphismSequenceContainer & psc, bool gapflag )
{
  PolymorphismSequenceContainer *sc;
  if (gapflag)
    sc = PolymorphismSequenceContainerTools::getSitesWithoutGaps(psc);
  else
    sc = new PolymorphismSequenceContainer(psc);
  //int K = 0;
  vector<string> pscvector;
  pscvector.push_back(sc->toString(0));
  //K++;
  for ( unsigned int i = 1; i < sc->getNumberOfSequences(); i++ ) {
    bool uniq = true;
    string query = sc->toString(i);
    for ( vector<string>::iterator it = pscvector.begin(); it != pscvector.end(); it++ ) {
      if ( query.compare(*it) == 0 ) {
        uniq = false;
        break;
      }
    }
    if (uniq) {
      //K++;
      pscvector.push_back(query);
    }
  }
  delete sc;
  //return K;
  return pscvector.size();
}

double SequenceStatistics::DVH( const PolymorphismSequenceContainer & psc, bool gapflag )
{
  PolymorphismSequenceContainer *sc;
  if (gapflag)
    sc = PolymorphismSequenceContainerTools::getSitesWithoutGaps(psc);
  else
    sc = new PolymorphismSequenceContainer(psc);
  double H = 0.0;
  unsigned int nbSeq;
  vector<string> pscvector;
  vector<int> effvector;
  pscvector.push_back(sc -> toString(0));
  effvector.push_back(sc -> getSequenceCount(0));
  nbSeq = sc -> getSequenceCount(0);
  for ( unsigned int i = 1; i < sc->getNumberOfSequences(); i++ ) {
    nbSeq += sc -> getSequenceCount(i);
    bool uniq = true;
    string query = sc -> toString(i);
    for (unsigned int j = 0; j < pscvector.size(); j++) {
      if ( query.compare(pscvector[j]) == 0) {
        effvector[j] += sc->getSequenceCount(i);
        uniq = false;
        break;
      }
    }
    if (uniq) {
      pscvector.push_back(query);
      effvector.push_back(sc -> getSequenceCount(i));
    }
  }
  for ( unsigned int i = 0; i < effvector.size(); i++ ) {
    H -= ( (double) effvector[i] / (double) nbSeq ) * ( (double) effvector[i] / (double) nbSeq );
  }
  H += 1.0;
  delete sc;
  return H;
}

unsigned int SequenceStatistics::getNumberOfTransitions( const PolymorphismSequenceContainer & psc )
{
  const Site *site;
  SiteIterator *si;
  unsigned int nbT = 0;
  si = new CompleteSiteIterator(psc);
  while (si->hasMoreSites()) {
    site = si->nextSite();
    //if(SiteTools::isConstant(*site) || SiteTools::isTriplet(*site)) continue;
    if (SiteTools::getNumberOfDistinctCharacters(*site) != 2) continue;
    vector<int> seq = site->getContent();
    int state1 = seq[0];
    int state2 = seq[0];
    for(unsigned int i = 1; i < seq.size(); i++) {
      if(state1 != seq[i]) {
        state2 = seq[i];
        break;
      }
    }
    if(((state1==0 && state2==2) || (state1==2 && state2==0)) ||
        ((state1==1 && state2==3) || (state1==3 && state2==1))) {
      nbT++;
    }
  }
  return nbT;
}

unsigned int SequenceStatistics::getNumberOfTransversions( const PolymorphismSequenceContainer & psc )
{
  const Site *site;
  SiteIterator *si;
  unsigned int nbTv = 0;
  si = new CompleteSiteIterator(psc);
  while (si->hasMoreSites()) {
    site = si->nextSite();
    //if(SiteTools::isConstant(*site) || SiteTools::isTriplet(*site)) continue;
    if (SiteTools::getNumberOfDistinctCharacters(*site) != 2) continue;
    vector<int> seq = site->getContent();
    int state1 = seq[0];
    int state2 = seq[0];
    for(unsigned int i = 1; i < seq.size(); i++) {
      if(state1 != seq[i]) {
        state2 = seq[i];
        break;
      }
    }
    if(!(((state1==0 && state2==2) || (state1==2 && state2==0)) ||
          ((state1==1 && state2==3) || (state1==3 && state2==1)))) {
      nbTv++;
    }

  }
  return nbTv;
}

double SequenceStatistics::getTransitionsTransversionsRatio( const PolymorphismSequenceContainer & psc ) throw(Exception)
{
  //return (double) getNumberOfTransitions(psc)/getNumberOfTransversions(psc);

  const Site *site;
  SiteIterator *si;
  unsigned int nbT = 0;
  unsigned int nbTv = 0;
  si = new CompleteSiteIterator(psc);
  map<int,unsigned int> count ;
  vector < int > state(2);

  while (si->hasMoreSites()) {
    site = si->nextSite();
    count = SymbolListTools::getCounts(*site);
    if (count.size() != 2) continue;
    int i = 0;
    for(map<int, unsigned int>::iterator it=count.begin(); it!=count.end(); it++)
    {
      state[i] = it->first;
      i++;
    }


    if(((state[0]==0 && state[1]==2) || (state[0]==2 && state[1]==0)) ||
        ((state[0]==1 && state[1]==3) || (state[0]==3 && state[1]==1))) {
      nbT++;//transitions
    }
    else
    {
      nbTv++;//transversion
    }

  }

  if (nbTv == 0)
    throw ZeroDivisionException("SequenceStatistics::getTransitionsTransversionsRatio.");

  return (double)nbT/nbTv;

}

//******************************************************************************************************************
//Synonymous and non-synonymous polymorphism
//******************************************************************************************************************

unsigned int SequenceStatistics::stopCodonSiteNumber(const PolymorphismSequenceContainer & psc, bool gapflag)
{
  SiteIterator *si = NULL;
  if(gapflag) si = new NoGapSiteIterator(psc);
  else si = new SimpleSiteIterator(psc);
  unsigned int S=0;
  const Site *site;
  while ( si->hasMoreSites() ) {
    site=si->nextSite();
    if (CodonSiteTools::hasStop(*site))S++;
  }
  delete si;
  return S;
}

unsigned int SequenceStatistics::monoSitePolymorphicCodonNumber(const PolymorphismSequenceContainer & psc, bool stopflag, bool gapflag)
{
  SiteIterator *si = NULL;
  if(stopflag) si = new CompleteSiteIterator(psc);
  else {
    if(gapflag) si = new NoGapSiteIterator(psc);
    else si = new SimpleSiteIterator(psc);
  }
  unsigned int S=0;
  const Site *site;
  while ( si->hasMoreSites() ) {
    site=si->nextSite();
    if (CodonSiteTools::isMonoSitePolymorphic(*site))S++;
  }
  delete si;
  return S;
}

unsigned int SequenceStatistics::synonymousPolymorphicCodonNumber(const PolymorphismSequenceContainer & psc, const GeneticCode & gc)
{
  SiteIterator* si = new CompleteSiteIterator(psc);
  unsigned int S=0;
  const Site *site;
  while ( si->hasMoreSites() ) {
    site=si->nextSite();
    if (CodonSiteTools::isSynonymousPolymorphic(*site,gc)) S++;
  }
  delete si;
  return S;
}

double SequenceStatistics::watterson75Synonymous(const PolymorphismSequenceContainer & psc, const GeneticCode & gc)
{
  double ThetaW;
  unsigned int n = psc.getNumberOfSequences();
  unsigned int S = synonymousSubstitutionsNumber(psc,gc);
  map<string, double> values = _getUsefullValues(n);
  ThetaW = (double) S / values["a1"];
  return ThetaW;
}

double SequenceStatistics::watterson75NonSynonymous(const PolymorphismSequenceContainer & psc, const GeneticCode & gc)
{
  double ThetaW;
  unsigned int n = psc.getNumberOfSequences();
  unsigned int S = nonSynonymousSubstitutionsNumber(psc,gc);
  map<string, double> values = _getUsefullValues(n);
  ThetaW = (double) S / values["a1"];
  return ThetaW;
}

double SequenceStatistics::piSynonymous(const PolymorphismSequenceContainer & psc, const GeneticCode & gc, bool minchange)
{
  SiteIterator* si = new CompleteSiteIterator(psc);
  double S=0.0;
  const Site *site;
  while(si->hasMoreSites()) {
    site=si->nextSite();
    S += CodonSiteTools::piSynonymous(*site,gc,minchange);
  }
  delete si;
  return S;
}

double SequenceStatistics::piNonSynonymous(const PolymorphismSequenceContainer & psc, const GeneticCode & gc, bool minchange)
{
  SiteIterator *si = new CompleteSiteIterator(psc);
  double S=0;
  const Site *site;
  while(si->hasMoreSites()) {
    site=si->nextSite();
    S += CodonSiteTools::piNonSynonymous(*site,gc,minchange);
  }
  delete si;
  return S;
}

double SequenceStatistics::meanSynonymousSitesNumber(const PolymorphismSequenceContainer & psc, const GeneticCode & gc, double ratio)
{
  SiteIterator *si =  new CompleteSiteIterator(psc);
  double S=0;
  const Site *site;
  while(si->hasMoreSites()) {
    site=si->nextSite();
    S += CodonSiteTools::meanNumberOfSynonymousPositions(*site,gc,ratio);
  }
  delete si;
  return S;
}

double SequenceStatistics::meanNonSynonymousSitesNumber(const PolymorphismSequenceContainer & psc, const GeneticCode & gc, double ratio)
{
  SiteIterator *si = new CompleteSiteIterator(psc);
  double S=0;
  int n=0;
  const Site *site;
  while(si->hasMoreSites()) {
    site=si->nextSite();
    n = n + 3;
    S += CodonSiteTools::meanNumberOfSynonymousPositions(*site,gc,ratio);
  }
  delete si;
  return ((double) n - S);
}

unsigned int SequenceStatistics::synonymousSubstitutionsNumber(const PolymorphismSequenceContainer & psc, const GeneticCode & gc, double freqmin)
{
  SiteIterator * si = new CompleteSiteIterator(psc);
  const Site * site;
  unsigned int St = 0, Sns = 0;
  while(si->hasMoreSites()) {
    site=si->nextSite();
    St += CodonSiteTools::numberOfSubsitutions(*site,freqmin);
    Sns += CodonSiteTools::numberOfNonSynonymousSubstitutions(*site,gc,freqmin);
  }
  delete si;
  return St - Sns;
}

unsigned int SequenceStatistics::nonSynonymousSubstitutionsNumber(const PolymorphismSequenceContainer & psc, const GeneticCode & gc, double freqmin)
{
  SiteIterator *si = new CompleteSiteIterator(psc);
  const Site * site;
  unsigned int Sns = 0;
  while(si->hasMoreSites()) {
    site=si->nextSite();
    Sns += CodonSiteTools::numberOfNonSynonymousSubstitutions(*site,gc,freqmin);
  }
  delete si;
  return Sns;
}

vector<unsigned int> SequenceStatistics::fixedDifferences(const PolymorphismSequenceContainer & pscin, const PolymorphismSequenceContainer & pscout, PolymorphismSequenceContainer & psccons, const GeneticCode & gc)
{
  SiteIterator *siIn = new CompleteSiteIterator(pscin);
  SiteIterator *siOut = new CompleteSiteIterator(pscout);
  SiteIterator *siCons = new CompleteSiteIterator(psccons);
  const Site *siteIn, *siteOut, *siteCons;
  unsigned int NfixS=0;
  unsigned int NfixA=0;
  while(siIn->hasMoreSites()){
    siteIn = siIn->nextSite();
    siteOut = siOut->nextSite();
    siteCons = siCons->nextSite();
    vector<unsigned int> v = CodonSiteTools::fixedDifferences(*siteIn,*siteOut,siteCons->getValue(0),siteCons->getValue(1),gc);
    NfixS += v[0];
    NfixA += v[1];
  }
  vector<unsigned int> v(2);
  v[0]=NfixS;
  v[1]=NfixA;
  delete siIn;
  delete siOut;
  delete siCons;
  return v;
}

vector<unsigned int> SequenceStatistics::MKtable(const PolymorphismSequenceContainer & ingroup, const PolymorphismSequenceContainer & outgroup , const GeneticCode & gc, double freqmin)
{
  PolymorphismSequenceContainer * psctot = new PolymorphismSequenceContainer(ingroup);
  for(unsigned int i = 0; i<outgroup.getNumberOfSequences();i++){
    psctot->addSequence(*outgroup.getSequence(i));
    psctot->setAsOutgroupMember(i+ingroup.getNumberOfSequences());
  }
  const PolymorphismSequenceContainer * psccomplet = PolymorphismSequenceContainerTools::getCompleteSites(*psctot);
  const PolymorphismSequenceContainer * pscin = PolymorphismSequenceContainerTools::extractIngroup(*psccomplet);
  const PolymorphismSequenceContainer * pscout = PolymorphismSequenceContainerTools::extractOutgroup(*psccomplet);
  const Sequence * consensusIn = SiteContainerTools::getConsensus(*pscin,"consensusIn");
  const Sequence * consensusOut = SiteContainerTools::getConsensus(*pscout,"consensusOut");
  PolymorphismSequenceContainer * consensus = new PolymorphismSequenceContainer(ingroup.getAlphabet());
  consensus->addSequence(*consensusIn);
  consensus->addSequence(*consensusOut);
  vector<unsigned int> u = SequenceStatistics::fixedDifferences(*pscin,*pscout,*consensus,gc);
  vector<unsigned int> v(4);
  v[0] = SequenceStatistics::nonSynonymousSubstitutionsNumber(*pscin,gc,freqmin);
  v[1] = SequenceStatistics::synonymousSubstitutionsNumber(*pscin,gc,freqmin);
  v[2] = u[1];
  v[3] = u[0];
  delete consensus;
  if(psccomplet) {delete psccomplet;}
  if(pscin) {delete pscin;}
  if(pscout) {delete pscout;}
  if(consensusIn) {delete consensusIn;}
  if(consensusOut) {delete consensusOut;}
  delete psctot;
  return v;
}

double SequenceStatistics::neutralityIndex(const PolymorphismSequenceContainer & ingroup, const PolymorphismSequenceContainer & outgroup , const GeneticCode & gc, double freqmin)
{
  vector<unsigned int> v = SequenceStatistics::MKtable(ingroup,outgroup,gc,freqmin);
  if(v[1]!=0 && v[2]!=0) return (double)(v[0]*v[3])/(v[1]*v[2]);
  else return -1;
}

//******************************************************************************************************************
//Statistical tests
//******************************************************************************************************************

double SequenceStatistics::tajimaDSS(const PolymorphismSequenceContainer & psc, bool gapflag)
{
  unsigned int S = polymorphicSiteNumber(psc, gapflag);
  double tajima = tajima83(psc, gapflag);
  double watterson = watterson75(psc, gapflag);
  unsigned int n = psc.getNumberOfSequences();
  map<string, double> values = _getUsefullValues(n);
  return (tajima - watterson) / sqrt((values["e1"] * S) + (values["e2"] * S * (S - 1)));
}

double SequenceStatistics::tajimaDTNM(const PolymorphismSequenceContainer & psc, bool gapflag)
{
  unsigned int eta = totNumberMutations(psc, gapflag);
  double tajima = tajima83(psc, gapflag);
  unsigned int n = psc.getNumberOfSequences();
  map<string, double> values = _getUsefullValues(n);
  double eta_a1 = (double) eta / values["a1"];
  return (tajima - eta_a1) / sqrt((values["e1"] * eta) + (values["e2"] * eta * (eta - 1)));
}

double SequenceStatistics::fuliD(const PolymorphismSequenceContainer & ingroup, const PolymorphismSequenceContainer & outgroup)
{
  unsigned int n = ingroup.getNumberOfSequences();
  double nn = (double) n;
  map<string, double> values = _getUsefullValues(n);
  double vD = 1. + (pow(values["a1"], 2) / (values["a2"] + pow(values["a1"], 2))) * (values["cn"] - ((nn + 1.) / (nn - 1.)));
  double uD = values["a1"] - 1. - vD;
  double eta = (double) totNumberMutations(ingroup);//using the number of mutations
  //double eta = (double)polymorphicSiteNumber(ingroup);
  double etae = (double) totMutationsExternalBranchs(ingroup,outgroup);
  return (eta - values["a1"] * etae) / sqrt((uD * eta) + (vD * eta * eta));
}

double SequenceStatistics::fuliDstar(const PolymorphismSequenceContainer & group)
{
  unsigned int n = group.getNumberOfSequences();
  double nn = (double) n;
  map<string, double> values = _getUsefullValues(n);
  // Fu & Li 1993
  double _n = nn / (nn - 1.);
  double vDs = (
      (_n * _n * values["a2"])
      + (values["a1"] * values["a1"] * values["dn"])
      - (2. * (nn * values["a1"] * (values["a1"] + 1.) / ((nn - 1.) * (nn - 1.))))
      )
    /
    (pow(values["a1"], 2) + values["a2"]);
  double uDs = _n * (values["a1"] - _n) - vDs;

  // Simonsen et al. 1995
  /*  double vDs = (
      (values["a2"] / pow(values["a1"], 2))
      - (2./nn) * (1. + 1./values["a1"] - values["a1"] + values["a1"]/nn)
      - 1./(nn*nn)
      )
      /
      (pow(values["a1"], 2) + values["a2"]);
      double uDs = (((nn - 1.)/nn - 1./values["a1"]) / values["a1"]) - vDs;
      */
  double eta = (double) totNumberMutations(group);
  double etas = (double) countSingleton(group);

  // Fu & Li 1993
  return ((_n * eta) - (values["a1"] * etas)) / sqrt(uDs * eta + vDs * eta * eta);

  // Simonsen et al. 1995
  //  return ((eta / values["a1"]) - (etas * ((n - 1) / n))) / sqrt(uDs * eta + vDs * eta * eta);
}

double SequenceStatistics::fuliF(const PolymorphismSequenceContainer & ingroup, const PolymorphismSequenceContainer & outgroup)
{
  unsigned int n = ingroup.getNumberOfSequences();
  double nn = (double) n;
  map<string, double> values = _getUsefullValues(n);
  double pi = tajima83(ingroup, true);
  double vF = (values["cn"] + values["b2"] - 2. / (nn - 1.)) / (pow(values["a1"], 2) + values["a2"]);
  double uF = ((1. + values["b1"] - (4. * ((nn + 1.) / ((nn - 1.) * (nn - 1.)))) * (values["a1n"] - (2. * nn) / (nn + 1.))) / values["a1"]) - vF;
  double eta = (double) totNumberMutations(ingroup);
  double etae = (double) totMutationsExternalBranchs(ingroup,outgroup);
  return (pi - etae) / sqrt(uF * eta + vF * eta * eta);
}

double SequenceStatistics::fuliFstar(const PolymorphismSequenceContainer & group) {
  unsigned int n = group.getNumberOfSequences();
  double nn = (double) n;
  map<string, double> values = _getUsefullValues(n);
  double pi = tajima83(group, true);

  // Fu & Li 1993
  //  double vFs = (values["dn"] + values["b2"] - (2. / (nn - 1.)) * (4. * values["a2"] - 6. + 8. / nn)) / (pow(values["a1"], 2) + values["a2"]);
  //  double uFs = (((nn / (nn - 1.)) + values["b1"] - (4. / (nn * (nn - 1.))) + 2. * ((nn + 1.) / (pow((nn - 1.), 2))) * (values["a1n"] - 2. * nn / (nn + 1.))) / values["a1"]) - vFs;

  // Simonsen et al. 1995
  double vFs = (((2*nn*nn*nn + 110*nn*nn - 255*nn + 153) / (9*nn*nn*(nn-1))) + ((2*(n-1)*values["a1"]) / (nn*nn)) - 8*values["a2"]/nn) / (pow(values["a1"], 2) + values["a2"]);
  double uFs = (((4*nn*nn + 19*nn + 3 - 12*(nn+1.)*values["a1n"]) / (3*nn*(n-1))) / values["a1"]) - vFs;
  double eta = (double) totNumberMutations(group);
  double etas = (double) countSingleton(group);
  // Fu & Li 1993
  // Simonsen et al. 1995
  return (pi - ((nn - 1.) / nn * etas)) / sqrt(uFs * eta + vFs * eta * eta);
}

//******************************************************************************************************************
//Linkage disequilibrium statistics
//******************************************************************************************************************

/**********************/
/* Preliminary method */
/**********************/

PolymorphismSequenceContainer * SequenceStatistics::generateLDContainer(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin)
{
  SiteSelection ss;
  // Extract polymorphic site with only two alleles
  for(unsigned int i=0; i<psc.getNumberOfSites(); i++){
    if(keepsingleton) {
      if(SiteTools::isComplete(*psc.getSite(i)) && !SiteTools::isConstant(*psc.getSite(i)) && !SiteTools::isTriplet(*psc.getSite(i))){
        ss.push_back(i);
      }
    }
    else{
      if(SiteTools::isComplete(*psc.getSite(i)) && !SiteTools::isConstant(*psc.getSite(i)) && !SiteTools::isTriplet(*psc.getSite(i)) && !SiteTools::hasSingleton(*psc.getSite(i))){
        ss.push_back(i);
      }
    }
  }

  const SiteContainer* sc = SiteContainerTools::getSelectedSites(psc,ss);
  Alphabet* alpha = new DNA();
  PolymorphismSequenceContainer *ldpsc = new PolymorphismSequenceContainer(sc->getNumberOfSequences(),alpha);
  // Assign 1 to the more frequent and 0 to the less frequent alleles
  for(unsigned int i=0; i<sc->getNumberOfSites(); i++){
    const Site* site = sc->getSite(i);
    Site* siteclone =  new Site(*site);
    bool deletesite = false;
    map<int, double> freqs = SymbolListTools::getFrequencies(*siteclone);
    int first = 0;
    for(map<int,double>::iterator it=freqs.begin(); it!=freqs.end(); it++){
      if(it->second>=0.5) first = it->first;
    }
    for(unsigned int j=0; j<sc->getNumberOfSequences(); j++){
      if(freqs[site->getValue(j)]>=0.5 && site->getValue(j)==first){
        if(freqs[site->getValue(j)]<=1-freqmin) {
          siteclone->setElement(j,1);
          first = site->getValue(j);
        }
        else deletesite = true;
      }
      else {
        if(freqs[site->getValue(j)]>=freqmin) siteclone->setElement(j,0);
        else deletesite = true;
      }
    }
    if(!deletesite)  ldpsc->addSite(*siteclone);
    delete siteclone;
  }
  delete alpha;
  return ldpsc;
}

/*************************************/
/* Pairwise LD and distance measures */
/*************************************/

Vdouble SequenceStatistics::pairwiseDistances1(const PolymorphismSequenceContainer & psc,bool keepsingleton, double freqmin) throw (DimensionException)
{
  //get Positions with sites of interest
  SiteSelection ss;
  for(unsigned int i=0; i<psc.getNumberOfSites(); i++){
    if(keepsingleton) {
      if(SiteTools::isComplete(*psc.getSite(i)) && !SiteTools::isConstant(*psc.getSite(i)) && !SiteTools::isTriplet(*psc.getSite(i))){
        const Site* site = psc.getSite(i);
        bool deletesite = false;
        map<int, double> freqs = SymbolListTools::getFrequencies(*site);
        for(unsigned int j=0; j<site->getAlphabet()->getSize(); j++){
          if(freqs[j]>=1-freqmin) deletesite = true;
        }
        if(!deletesite) ss.push_back(i);
      }
    }
    else{
      if(SiteTools::isComplete(*psc.getSite(i)) && !SiteTools::isConstant(*psc.getSite(i)) && !SiteTools::isTriplet(*psc.getSite(i)) && !SiteTools::hasSingleton(*psc.getSite(i))){
        ss.push_back(i);
        const Site* site = psc.getSite(i);
        bool deletesite = false;
        map<int, double> freqs = SymbolListTools::getFrequencies(*site);
        for(unsigned int j=0; j<site->getAlphabet()->getSize(); j++){
          if(freqs[j]>=1-freqmin) deletesite = true;
        }
        if(!deletesite) ss.push_back(i);
      }
    }
  }
  //compute pairwise distances
  if(ss.size()<2) throw DimensionException("SequenceStatistics::pairwiseDistances1 : less than 2 sites are available",ss.size(),2);
  Vdouble dist;
  for(unsigned int i=0; i<ss.size()-1; i++){
    for(unsigned int j=i+1; j<ss.size(); j++){
      dist.push_back(ss[j]-ss[i]);
    }
  }
  return dist;
}

Vdouble SequenceStatistics::pairwiseDistances2(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  SiteSelection ss;
  for(unsigned int i=0; i<psc.getNumberOfSites(); i++){
    if(keepsingleton) {
      if(SiteTools::isComplete(*psc.getSite(i)) && !SiteTools::isConstant(*psc.getSite(i)) && !SiteTools::isTriplet(*psc.getSite(i))){
        const Site* site = psc.getSite(i);
        bool deletesite = false;
        map<int, double> freqs = SymbolListTools::getFrequencies(*site);
        for(unsigned int j=0; j<site->getAlphabet()->getSize(); j++){
          if(freqs[j]>=1-freqmin) deletesite = true;
        }
        if(!deletesite) ss.push_back(i);
      }
    }
    else{
      if(SiteTools::isComplete(*psc.getSite(i)) && !SiteTools::isConstant(*psc.getSite(i)) && !SiteTools::isTriplet(*psc.getSite(i)) && !SiteTools::hasSingleton(*psc.getSite(i))){
        ss.push_back(i);
        const Site* site = psc.getSite(i);
        bool deletesite = false;
        map<int, double> freqs = SymbolListTools::getFrequencies(*site);
        for(unsigned int j=0; j<site->getAlphabet()->getSize(); j++){
          if(freqs[j]>=1-freqmin) deletesite = true;
        }
        if(!deletesite) ss.push_back(i);
      }
    }
  }
  unsigned int n = ss.size();
  if(n<2) throw DimensionException("SequenceStatistics::pairwiseDistances1 : less than 2 sites are available",ss.size(),2);
  Vdouble distance(n*(n-1)/2,0);
  unsigned int nbsite = psc.getNumberOfSites();
  for(unsigned int k=0; k<psc.getNumberOfSequences(); k++){
    const Sequence* seq = psc.getSequence(k);
    SiteSelection gap, newss = ss;
    Vdouble dist;
    for(unsigned int i=0; i<nbsite; i++){
      if(seq->getValue(i)==-1) gap.push_back(i);
    }
    //Site positions are re-numbered to take gaps into account
    for(unsigned int i=0; i<gap.size(); i++){
      for(unsigned int j=0; j<ss.size(); j++){
        if(ss[j]>gap[i]) newss[j]--;
      }
    }
    for(unsigned int i=0; i<n-1; i++){
      for(unsigned int j=i+1; j<n; j++){
        dist.push_back(newss[j]-newss[i]);
      }
    }
    distance += dist;
  }
  distance = distance/psc.getNumberOfSequences();
  return distance;
}

Vdouble SequenceStatistics::pairwiseD(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  PolymorphismSequenceContainer* newpsc = SequenceStatistics::generateLDContainer(psc, keepsingleton,  freqmin);
  Vdouble D;
  unsigned int nbsite = newpsc->getNumberOfSites();
  unsigned int nbseq = newpsc->getNumberOfSequences();
  if(nbsite<2) throw DimensionException("SequenceStatistics::pairwiseD: less than two sites are available",nbsite,2);
  if(nbseq<2) throw DimensionException("SequenceStatistics::pairwiseD: less than two sequences are available",nbseq,2);
  for(unsigned int i=0; i<nbsite-1; i++){
    for(unsigned int j=i+1; j<nbsite; j++){
      double haplo=0;
      const Site* site1 = newpsc->getSite(i);
      const Site* site2 = newpsc->getSite(j);
      map<int,double> freq1 = SymbolListTools::getFrequencies(*site1);
      map<int,double> freq2 = SymbolListTools::getFrequencies(*site2);
      for(unsigned int k=0; k<nbseq; k++){
        if(site1->getValue(k) + site2->getValue(k)==2) haplo++;
      }
      haplo = haplo/nbseq;
      D.push_back(std::abs(haplo-freq1[1]*freq2[1]));
    }
  }
  return D;
}

Vdouble SequenceStatistics::pairwiseDprime(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  PolymorphismSequenceContainer* newpsc = SequenceStatistics::generateLDContainer(psc, keepsingleton, freqmin);
  Vdouble Dprime;
  unsigned int nbsite = newpsc->getNumberOfSites();
  unsigned int nbseq = newpsc->getNumberOfSequences();
  if(nbsite<2) throw DimensionException("SequenceStatistics::pairwiseD: less than two sites are available",nbsite,2);
  if(nbseq<2) throw DimensionException("SequenceStatistics::pairwiseD: less than two sequences are available",nbseq,2);
  for(unsigned int i=0; i<nbsite-1; i++){
    for(unsigned int j=i+1; j<nbsite; j++){
      double haplo=0;
      const Site* site1 = newpsc->getSite(i);
      const Site* site2 = newpsc->getSite(j);
      map<int,double> freq1 = SymbolListTools::getFrequencies(*site1);
      map<int,double> freq2 = SymbolListTools::getFrequencies(*site2);
      for(unsigned int k=0; k<nbseq; k++){
        if(site1->getValue(k) + site2->getValue(k)==2) haplo++;
      }
      haplo = haplo/nbseq;
      double d, D = (haplo-freq1[1]*freq2[1]);
      if(D>0){
        if(freq1[1]*freq2[0]<=freq1[0]*freq2[1]){
          d=std::abs(D)/(freq1[1]*freq2[0]);
        }
        else{
          d=std::abs(D)/(freq1[0]*freq2[1]);
        }
      }
      else{
        if(freq1[1]*freq2[1]<=freq1[0]*freq2[0]){
          d=std::abs(D)/(freq1[1]*freq2[1]);
        }
        else{
          d=std::abs(D)/(freq1[0]*freq2[0]);
        }
      }
      Dprime.push_back(d);
    }
  }
  return Dprime;
}

Vdouble SequenceStatistics::pairwiseR2(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  PolymorphismSequenceContainer* newpsc = SequenceStatistics::generateLDContainer(psc, keepsingleton, freqmin);
  Vdouble R2;
  unsigned int nbsite = newpsc->getNumberOfSites();
  unsigned int nbseq = newpsc->getNumberOfSequences();
  if(nbsite<2) throw DimensionException("SequenceStatistics::pairwiseD: less than two sites are available",nbsite,2);
  if(nbseq<2) throw DimensionException("SequenceStatistics::pairwiseD: less than two sequences are available",nbseq,2);
  for(unsigned int i=0; i<nbsite-1; i++){
    for(unsigned int j=i+1; j<nbsite; j++){
      double haplo=0;
      const Site* site1 = newpsc->getSite(i);
      const Site* site2 = newpsc->getSite(j);
      map<int,double> freq1 = SymbolListTools::getFrequencies(*site1);
      map<int,double> freq2 = SymbolListTools::getFrequencies(*site2);
      for(unsigned int k=0; k<nbseq; k++){
        if(site1->getValue(k) + site2->getValue(k)==2) haplo++;
      }
      haplo = haplo/nbseq;
      double r = ((haplo-freq1[1]*freq2[1])*(haplo-freq1[1]*freq2[1]))/(freq1[0]*freq1[1]*freq2[0]*freq2[1]);
      R2.push_back(r);
    }
  }
  return R2;
}

/***********************************/
/* Global LD and distance measures */
/***********************************/

double SequenceStatistics::meanD(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  Vdouble D = pairwiseD(psc, keepsingleton, freqmin);
  return VectorTools::mean<double, double>(D);
}

double SequenceStatistics::meanDprime(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble Dprime = pairwiseDprime(psc, keepsingleton, freqmin);
    return VectorTools::mean<double, double>(Dprime);
  }
  catch (DimensionException & e) {throw e;}
}

double SequenceStatistics::meanR2(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble R2 = SequenceStatistics::pairwiseR2(psc, keepsingleton, freqmin);
    return VectorTools::mean<double, double>(R2);
  }
  catch (DimensionException & e) {throw e;}
}

double SequenceStatistics::meanDistance1(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble dist = pairwiseDistances1(psc, keepsingleton, freqmin);
    return VectorTools::mean<double, double>(dist);
  }
  catch (DimensionException & e) {throw e;}
}

double SequenceStatistics::meanDistance2(const PolymorphismSequenceContainer & psc, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble dist = pairwiseDistances2(psc, keepsingleton, freqmin);
    return VectorTools::mean<double, double>(dist);
  }
  catch (DimensionException & e) {throw e;}
}

/**********************/
/* Regression methods */
/**********************/

double SequenceStatistics::originRegressionD(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble D = pairwiseD(psc, keepsingleton, freqmin)-1;
    Vdouble dist;
    if(distance1) dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    return VectorTools::sum(D * dist) / VectorTools::sum(dist * dist);
  }
  catch (DimensionException & e) {throw e;}
}

double SequenceStatistics::originRegressionDprime(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble Dprime = pairwiseDprime(psc, keepsingleton, freqmin)-1;
    Vdouble dist;
    if(distance1) dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    return VectorTools::sum(Dprime * dist) / VectorTools::sum(dist * dist);
  }
  catch (DimensionException & e) {throw e;}
}

double SequenceStatistics::originRegressionR2(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble R2 = pairwiseR2(psc, keepsingleton, freqmin)-1;
    Vdouble dist;
    if(distance1) dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    return VectorTools::sum(R2 * dist) / VectorTools::sum(dist * dist);
  }
  catch (DimensionException & e) {throw e;}
}

Vdouble SequenceStatistics::linearRegressionD(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble D = pairwiseD(psc, keepsingleton, freqmin);
    Vdouble dist;
    Vdouble reg(2);
    if(distance1) dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    reg[0] = VectorTools::cov<double, double>(dist, D) / VectorTools::var<double, double>(dist);
    reg[1] = VectorTools::mean<double, double>(D) - reg[0] * VectorTools::mean<double, double>(dist);
    return reg;
  }
  catch (DimensionException & e) {throw e;}
}

Vdouble SequenceStatistics::linearRegressionDprime(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble Dprime = pairwiseDprime(psc, keepsingleton, freqmin);
    Vdouble dist;
    Vdouble reg(2);
    if(distance1) dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    reg[0] = VectorTools::cov<double, double>(dist, Dprime) / VectorTools::var<double, double>(dist);
    reg[1] = VectorTools::mean<double, double>(Dprime) - reg[0] * VectorTools::mean<double, double>(dist);
    return reg;
  }
  catch (DimensionException & e) {throw e;}
}

Vdouble SequenceStatistics::linearRegressionR2(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble R2 = pairwiseR2(psc, keepsingleton, freqmin);
    Vdouble dist;
    Vdouble reg(2);
    if(distance1) dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    reg[0] = VectorTools::cov<double, double>(dist, R2) / VectorTools::var<double, double>(dist);
    reg[1] = VectorTools::mean<double, double>(R2) - reg[0] * VectorTools::mean<double, double>(dist);
    return reg;
  }
  catch (DimensionException & e) {throw e;}
}

double SequenceStatistics::inverseRegressionR2(const PolymorphismSequenceContainer & psc, bool distance1, bool keepsingleton, double freqmin) throw (DimensionException)
{
  try
  {
    Vdouble R2 = pairwiseR2(psc, keepsingleton, freqmin);
    Vdouble unit(R2.size(),1);
    Vdouble R2transformed = unit/R2 -1;
    Vdouble dist;
    if(distance1) dist = pairwiseDistances1(psc, keepsingleton, freqmin) / 1000;
    else dist = pairwiseDistances2(psc, keepsingleton, freqmin) / 1000;
    return VectorTools::sum(R2transformed * dist) / VectorTools::sum(dist * dist);
  }
  catch (DimensionException & e) {throw e;}
}

/**********************/
/*   Hudson method    */
/**********************/

double SequenceStatistics::hudson87(const PolymorphismSequenceContainer & psc, double precision, double cinf, double csup){
  double left = _leftHandHudson(psc);
  unsigned int n = psc.getNumberOfSequences();
  double dif = 1;
  double c1=cinf;
  double c2=csup;
  if(SequenceStatistics::polymorphicSiteNumber(psc) < 2) return -1;
  if(_rightHandHudson(c1,n)<left) return cinf;
  if(_rightHandHudson(c2,n)>left) return csup;
  while(dif > precision){
    if(_rightHandHudson((c1+c2)/2,n)>left) c1=(c1+c2)/2;
    else c2=(c1+c2)/2;
    dif=std::abs(2*(c1-c2)/(c1+c2));
  }
  return (c1+c2)/2;
}

//******************************************************************************************************************
//Private methods
//******************************************************************************************************************

unsigned int SequenceStatistics::_getMutationNumber(const Site & site)
{
  unsigned int tmp_count = 0;

  map<int, unsigned int> states_count = SymbolListTools::getCounts(site);

  for (map<int, unsigned int>::iterator it = states_count.begin() ; it != states_count.end() ; it++)

    if (it->first >= 0)

      tmp_count++;

  if (tmp_count > 0)

    tmp_count--;

  return tmp_count;
}

unsigned int SequenceStatistics::_getSingletonNumber(const Site & site)
{
  unsigned int nus = 0;
  map<int, unsigned int> states_count = SymbolListTools::getCounts(site);
  for (map<int, unsigned int>::iterator it = states_count.begin() ; it != states_count.end() ; it++)
    if (it->second == 1)
      nus++;
  return nus;
}

unsigned int SequenceStatistics::_getDerivedSingletonNumber(const Site & site_in,const Site & site_out )
{
  unsigned int nus = 0;
  map<int, unsigned int> states_count = SymbolListTools::getCounts(site_in);
  map<int, unsigned int> outgroup_states_count = SymbolListTools::getCounts(site_out);
  //if there is more than one variant in the outgroup we will not be able to recover the ancestral state
  if (outgroup_states_count.size() == 1 )
  {
    for (map<int, unsigned int>::iterator it = states_count.begin() ; it != states_count.end() ; it++)
      if (it->second == 1)
      {       if ( outgroup_states_count.find(it->first) == outgroup_states_count.end() )
        nus++;
      }
  }
  return nus;
}

map<string, double> SequenceStatistics::_getUsefullValues(unsigned int n)
{
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
  if (n > 1) {
    for (unsigned int i = 1 ; i < n ; i++) {
      values["a1"] += 1. / i;
      values["a2"] += 1. / (i * i);
    }
    double nn = (double) n;
    values["a1n"] = values["a1"] + (1. / nn);
    values["b1"] = (nn + 1.) / (3. * (nn - 1.));
    values["b2"] = 2. * ((nn * nn) + nn + 3.) / (9. * nn * (nn - 1.));
    values["c1"] = values["b1"] - (1. / values["a1"]);
    values["c2"] = values["b2"] - ((nn + 2.) / (values["a1"] * nn)) + (values["a2"] / (values["a1"] * values["a1"]));
    values["cn"] = 2. * ((nn * values["a1"]) - (2. * (nn - 1.))) / ((nn - 1.) * (nn - 2.));
    values["dn"] = values["cn"] + ((nn - 2.) / ((nn - 1.) * (nn - 1.))) + ((2. / (nn - 1.)) * ((3. / 2.) - (((2. * values["a1n"]) - 3.) / (nn - 2.)) - (1. / nn)));
    values["e1"] = values["c1"] / values["a1"];
    values["e2"] = values["c2"] / ((values["a1"] * values["a1"]) + values["a2"]);
  }
  return values;
}

double SequenceStatistics::_leftHandHudson(const PolymorphismSequenceContainer & psc)
{
  PolymorphismSequenceContainer *newpsc = PolymorphismSequenceContainerTools::getCompleteSites(psc);
  unsigned int nbseq = newpsc->getNumberOfSequences();
  double S1 = 0;
  double S2 = 0;
  for(unsigned int i=0; i<nbseq-1; i++){
    for(unsigned int j=i+1; j<nbseq; j++){
      SequenceSelection ss(2);
      ss[0]=i;ss[1]=j;
      PolymorphismSequenceContainer *psc2 = PolymorphismSequenceContainerTools::getSelectedSequences(*newpsc,ss);
      S1+=SequenceStatistics::watterson75(*psc2,true);
      S2+=SequenceStatistics::watterson75(*psc2,true)*SequenceStatistics::watterson75(*psc2,true);
      delete psc2;
    }
  }
  double Sk = (2*S2-pow(2*S1/nbseq,2.))/pow(nbseq,2.);
  double H = SequenceStatistics::heterozygosity(*newpsc);
  double H2 = SequenceStatistics::squaredHeterozygosity(*newpsc);
  delete newpsc;
  return (double)(Sk-H+H2)/pow(H*nbseq/(nbseq-1),2.);
}

double SequenceStatistics::_rightHandHudson(double c, unsigned int n)
{
  return (double) (1/(97*pow(c,2.)*pow(n,3.))) * ((n-1)*(97*(c*(4+(c-2*n)*n)+(-2*(7+c)+4*n+(c-1)*pow(n,2.))*log((18+c*(13+c))/18))+sqrt(97.)*(110+n*(49*n-52)+c*(2+n*(15*n-8)))*log(-1+(72+26*c)/(36+13*c-c*sqrt(97.)))));
}

