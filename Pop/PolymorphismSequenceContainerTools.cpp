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


//*******************************************************************************************************************************


PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::getCompleteSites (const PolymorphismSequenceContainer & psc) {
	vector<string> seqNames = psc.getSequencesNames();
	PolymorphismSequenceContainer * complete = new PolymorphismSequenceContainer(psc.getNumberOfSequences(), psc.getAlphabet());
	complete->setSequencesNames(seqNames, false);
	unsigned int nbSeq = psc.getNumberOfSequences();
	for (unsigned int i = 0 ; i < nbSeq ; i++) {
		complete->setSequenceCount(i, psc.getSequenceCount(i));
		if (! psc.isIngroupMember(i))
			complete->setAsOutgroupMember(i);
	}
	CompleteSiteIterator csi(psc);
	while(csi.hasMoreSites())
		complete->addSite(* csi.nextSite());
	return complete;
}

PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::excludeFlankingGap(const PolymorphismSequenceContainer & psc) throw (Exception) {
	try {
        PolymorphismSequenceContainer *psci = new PolymorphismSequenceContainer(psc);
        while(SiteTools::hasGap(*psci->getSite(0))){
                psci->deleteSite(0);
        }
        int i=0;
        int n = psci->getNumberOfSites();
        while(SiteTools::hasGap(*psci->getSite(n-i-1))){
                psci->deleteSite(n-i-1);
                i++;
        }
	return psci;
	}
	catch(...) {}
}


// Be carefull: To use before excluding gap
PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::getCodingSites(const PolymorphismSequenceContainer & psc) throw (Exception) {
	try {
	Comments maseFileHeader = psc.getGeneralComments();
        SiteSelection ss = MaseTools::getCodingPositions(maseFileHeader);
        const SiteContainer *sc = SiteContainerTools::getSelectedSites(psc,ss);
        PolymorphismSequenceContainer *psci = new PolymorphismSequenceContainer(*sc);
        delete sc;
        return psci;
        }
	catch(...) {}

}


// Be carefull: To use before excluding gap

PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::getNonCodingSites(const PolymorphismSequenceContainer & psc) throw (Exception) {
try {
	Comments maseFileHeader = psc.getGeneralComments();
        SiteSelection ss;
        SiteSelection codss = MaseTools::getCodingPositions(maseFileHeader);
        for(unsigned int i=0; i<psc.getNumberOfSites(); i++) {
                if(find(codss.begin(),codss.end(),i)==codss.end()) {
                        ss.push_back(i);
                }
        }
        const SiteContainer *sc = SiteContainerTools::getSelectedSites(psc,ss);
        PolymorphismSequenceContainer *psci = new PolymorphismSequenceContainer(*sc);
        delete sc;
        return psci;
        }
catch(...) {}

}


// Be carefull : to use before excluding gap
// The first site of the alignment is supposed to correspond to the first position
PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::getOnePosition(const PolymorphismSequenceContainer & psc, unsigned int pos) throw (Exception) {
	try {
        SiteSelection ss;
        for(unsigned int i=0; i<psc.getNumberOfSites(); i++){
                if((i+1)%pos ==0) ss.push_back(i);
                }
        const SiteContainer *sc = SiteContainerTools::getSelectedSites(psc,ss);
        PolymorphismSequenceContainer *newpsc = new PolymorphismSequenceContainer(*sc);
        delete sc;
        return newpsc;
        }
	catch(...) {}
}


//Same as getNonCodgingSites but exclude 5' and 3' flanking regions if there are
//Assumed that the first coding site correspond to the first position
PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::getIntrons(const PolymorphismSequenceContainer & psc, unsigned int start) throw (Exception) {
        //Temporary: start (0,1,2) is used to indicated the begining of the first position of the firts codon
        try {
	Comments maseFileHeader = psc.getGeneralComments();
        SiteSelection ss;
        SiteSelection codss = MaseTools::getCodingPositions(maseFileHeader);
        unsigned int first=0, last=psc.getNumberOfSites();
        //Check the first Codon
	if(psc.getSite(codss[0])->getValue(0)==0 &&
	   psc.getSite(codss[1])->getValue(0)==3 &&
	   psc.getSite(codss[2])->getValue(0)==2) first = codss[0];
	//Check the last Codon
	if(psc.getSite(codss[codss.size()-3-start])->getValue(0)==3) {
		if(psc.getSite(codss[codss.size()-2-start])->getValue(0)==0) {
                        if(psc.getSite(codss[codss.size()-1-start])->getValue(0)==0 || psc.getSite(codss[codss.size()-1-start])->getValue(0)==2){
				last = codss[codss.size()-1-start];
			}
                }
		if(psc.getSite(codss[codss.size()-2-start])->getValue(0)==2 && psc.getSite(codss[codss.size()-1-start])->getValue(0)==0) {
                        last = codss[codss.size()-1-start];
                }
        }
        for(unsigned int i=first; i<last; i++) {
            if(find(codss.begin(),codss.end(),i)==codss.end()) {
               ss.push_back(i);
            }
        }
        const SiteContainer *sc = SiteContainerTools::getSelectedSites(psc,ss);
        PolymorphismSequenceContainer *psci = new PolymorphismSequenceContainer(*sc);
        delete sc;
        return psci;
        }
        catch(...) {}

}

PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::get5Prime(const PolymorphismSequenceContainer & psc) throw (Exception) {
        try {
	Comments maseFileHeader = psc.getGeneralComments();
        SiteSelection ss;
        SiteSelection codss = MaseTools::getCodingPositions(maseFileHeader);
        unsigned int last=0;
        //Check the first Codon
        if(psc.getSite(codss[0])->getValue(0)==0 &&
	   psc.getSite(codss[1])->getValue(0)==3 &&
	   psc.getSite(codss[2])->getValue(0)==2) last = codss[0];
        for(unsigned int i=0; i<last; i++) {
                if(find(codss.begin(),codss.end(),i)==codss.end()) {
                        ss.push_back(i);
                }
        }
        const SiteContainer *sc = SiteContainerTools::getSelectedSites(psc,ss);
        PolymorphismSequenceContainer *psci = new PolymorphismSequenceContainer(*sc);
        delete sc;
        return psci;
        }
        catch(...) {}

}

PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::get3Prime(const PolymorphismSequenceContainer & psc, unsigned int start) throw (Exception) {
        //Temporary: start (0,1,2) is used to indicated the begining of the first position of the firts codon
        try {
	Comments maseFileHeader = psc.getGeneralComments();
        SiteSelection ss;
        SiteSelection codss = MaseTools::getCodingPositions(maseFileHeader);
        unsigned int first = psc.getNumberOfSites()-1;
	//Check the last Codon
	if(psc.getSite(codss[codss.size()-3-start])->getValue(0)==3) {
		if(psc.getSite(codss[codss.size()-2-start])->getValue(0)==0) {
                        if(psc.getSite(codss[codss.size()-1-start])->getValue(0)==0 || psc.getSite(codss[codss.size()-1-start])->getValue(0)==2){
				first = codss[codss.size()-1-start];
			}
                }
		if(psc.getSite(codss[codss.size()-2-start])->getValue(0)==2 && psc.getSite(codss[codss.size()-1-start])->getValue(0)==0) {
                        first = codss[codss.size()-1-start];
                }
        }
        for(unsigned int i=first; i<psc.getNumberOfSites(); i++) {
                if(find(codss.begin(),codss.end(),i)==codss.end()) {
                        ss.push_back(i);
                }
        }
        const SiteContainer *sc = SiteContainerTools::getSelectedSites(psc,ss);
        PolymorphismSequenceContainer *psci = new PolymorphismSequenceContainer(*sc);
        delete sc;
        return psci;
        }
        catch(...) {}

}

// Be carefull: To use before excluding gap
PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::getSelectedSites(const PolymorphismSequenceContainer & psc, const string &setName, bool phase) throw (Exception) {
	try {
        SiteContainer *pscc = MaseTools::getSelectedSites(psc, setName);
	if (phase) {
		for (unsigned int i=1; i < psc.getPhase(setName); i++) {
			pscc -> deleteSite(i - 1);
		}		
	}		
        PolymorphismSequenceContainer *psci = new PolymorphismSequenceContainer(*pscc);
	for(unsigned int i = 0; i < psc.getNumberOfSequences(); i++) {
		if (! psc.isIngroupMember(i))
			psci -> setAsOutgroupMember(i);
		psci -> setGroupId(i, psc.getGroupId(i));
	}
        delete pscc;
        return psci;
        }
	catch(...) {}
}







