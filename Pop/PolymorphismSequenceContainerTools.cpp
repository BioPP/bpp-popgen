/*
 * File: PolymorphismSequenceContainerTools.cpp
 * Authors: Eric Bazin <bazin@univ-montp2.fr>
 *          Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 29 2004
 *
 * Copyright (C) 2004 Eric Bazin, Sylvain Gaillard and the
 *                    PopGenLib Development Core Team
 *
 * This file is part of PopGenLib.
 *
 * PopGenLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopGenLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopGenLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
	for(int i = ss.size() - 1; i >= 0 ; i--)
		{
		psci->deleteSequence(ss[i]);
		}
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
	for(int i = ss.size() - 1; i >= 0; i--) {
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
	for (int i = ss.size() - 1 ; i >= 0 ; i--)
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
	unsigned int count = psc.getNumberOfSites();
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
			count--;
	delete ssi;
	return count;	
}

unsigned int PolymorphismSequenceContainerTools::getNumberOfCompleteSites(const PolymorphismSequenceContainer & psc, bool ingroup) throw (Exception) {
	unsigned int count = psc.getNumberOfSites();
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
		if (!SiteTools::isComplete(* ssi->nextSite()))
			count--;
	delete ssi;
	return count;	
}
