/*
 * File: PolymorphismSequenceContainerTools.cpp
 * Authors: bazin <bazin@univ-montp2.fr>
 *          Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday July 28 2004
 */

// from PolyLib
#include "PolymorphismSequenceContainerTools.h"

using namespace std;

PolymorphismSequenceContainerTools::~PolymorphismSequenceContainerTools() {}

// Read a mase+ format
PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::read(const string & path, const Alphabet * alpha) throw (Exception) {
	Mase ms;
	string key;
	unsigned int n;
	const VectorSequenceContainer *seqc = NULL;
	try {
		seqc = ms.read( path, alpha );
	}
	catch (Exception & e) {
		if (seqc != NULL)
			delete seqc;
		throw e;
	}
	PolymorphismSequenceContainer * psc = new PolymorphismSequenceContainer(* seqc);
	Comments maseFileHeader = seqc->getGeneralComments();
	delete seqc;
	map<string, unsigned int> groupMap = MaseTools::getAvailableSequenceSelections(maseFileHeader);
	for(map<string, unsigned int>::iterator mi = groupMap.begin() ; mi != groupMap.end() ; mi++) {
		key = mi->first;
		n = mi->second;
		if (key.compare(0, 8, "OUTGROUP") == 0 ) {
			SequenceSelection ss;
			try {
				ss = MaseTools::getSequenceSet(maseFileHeader, key);
			}
			catch (IOException & ioe) {
				delete psc;
				throw ioe;
			}
			for (unsigned int i = 0 ; i != ss.size() ; i++) {
				try {
					psc->setAsOutgroupMember(ss[i]);
				}
				catch (SequenceNotFoundException & snfe) {
					delete psc;
					throw snfe;
				}
			}			
		}
	}
	return psc;
}
	 
PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::extractIngroup (const PolymorphismSequenceContainer & psc) throw (Exception) {
	SequenceSelection ss;	
	PolymorphismSequenceContainer * psci = dynamic_cast<PolymorphismSequenceContainer *>(psc.clone());
	for(unsigned int i = 0; i < psc.getNumberOfSequences(); i++)
		if (! psc.isIngroupMember(i))
			ss.push_back(i);
	if (ss.size() == psc.getNumberOfSequences()) {
		delete psci;
		throw Exception("PolymorphismSequenceContainerTools::extractIngroup: no Ingroup sequences found.");
	}
	for(unsigned int i = ss.size() - 1 ; i <= 0 ; i--)
		psci->deleteSequence(ss[i]);
	return psci;
}

PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::extractOutgroup(const PolymorphismSequenceContainer & psc) throw (Exception) {
	SequenceSelection ss;
	PolymorphismSequenceContainer * psci = dynamic_cast<PolymorphismSequenceContainer *>(psc.clone());
	for(unsigned int i = 0; i < psc.getNumberOfSequences(); i++)
		if (psc.isIngroupMember(i) )
			ss.push_back(i);
	if (ss.size() == psc.getNumberOfSequences()) {
		delete psci;
		throw Exception("PolymorphismSequenceContainerTools::extractOutgroup: no Outgroup sequences found.");
	}
	for(unsigned int i = ss.size() - 1; i <= 0; i--) {
		psci->deleteSequence(ss[i]);
	}
	return psci;
}

PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::extractGroup(const PolymorphismSequenceContainer & psc, unsigned int group_id) throw (Exception) {
	SequenceSelection ss;
	PolymorphismSequenceContainer * psci = dynamic_cast<PolymorphismSequenceContainer *>(psc.clone());
	for (unsigned int i = 0 ; i < psc.getNumberOfSequences() ; i++)
		if (psc.getGroupId(i) == group_id)
			ss.push_back(i);
	if (ss.size() == 0) {
		delete psci;
		throw GroupNotFoundException("PolymorphismSequenceContainerTools::extractGroup: group_id not found.", group_id);
	}
	for (unsigned int i = ss.size() - 1 ; i <= 0 ; i--)
		psci->deleteSequence(ss[i]);
	return psci;
}
	 
PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::getSitesWithoutGaps (const PolymorphismSequenceContainer & psc) {
	vector<string> seqNames = psc.getSequencesNames();
	PolymorphismSequenceContainer * noGapCont = new PolymorphismSequenceContainer(psc.getNumberOfSequences(), psc.getAlphabet());
	noGapCont->setSequencesNames(seqNames, false);
	unsigned int nbSeq = psc.getNumberOfSequences();
	for (unsigned int i = 0 ; i < nbSeq ; i++) {
		noGapCont->setSequenceCount(i, psc.getSequenceCount(i));
		if (! psc.isIngroupMember(i))
			noGapCont->setAsOutgroupMember(i);
	}
	NoGapSiteIterator ngsi(psc);
	while(ngsi.hasMoreSites())
		noGapCont->addSite(* ngsi.nextSite());
	return noGapCont;
}

unsigned int PolymorphismSequenceContainerTools::getNumberOfNonGapSites(const PolymorphismSequenceContainer & psc, bool ingroup) throw (Exception) {
	unsigned int count = 0;
	PolymorphismSequenceContainer * npsc = NULL;
	SimpleSiteIterator * ssi;
	if (ingroup) {
		try {
			npsc = extractIngroup(psc);
		}
		catch (Exception & e) {
			if (npsc != NULL)
				delete npsc;
			throw e;
		}
		ssi = new SimpleSiteIterator(* npsc);
	}
	else
		ssi = new SimpleSiteIterator(psc);
	while (ssi->hasMoreSites())
		if (SiteTools::hasGap(* ssi->nextSite()))
			count++;
	delete ssi;
	return count;	
}
