//
// File: PolymorphismSequenceContainerTools.cpp
// Authors: bazin <bazin@univ-montp2.fr>
//          Sylvain Gaillard <yragael2001@yahoo.fr>
// Last modification : Wednesday June 16 2004
//

// from PolyLib
#include "PolymorphismSequenceContainerTools.h"

using namespace std;

PolymorphismSequenceContainerTools::~PolymorphismSequenceContainerTools() {}

// Read a mase+ format
PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::read(const string & path, const Alphabet *
alpha) throw (Exception) {
	try {
	Mase ms;
	string clef;
	unsigned int n;
	const VectorSequenceContainer *seqc = ms.read( path, alpha );
	VectorSiteContainer *sitec = new VectorSiteContainer( *seqc );
	PolymorphismSequenceContainer *psc = new PolymorphismSequenceContainer( *sitec );
	Comments maseFileHeader = sitec -> getGeneralComments();
	map< string, unsigned int > groupMap = MaseTools::getAvailableSequenceSelections(maseFileHeader);
	for(map< string, unsigned int >::iterator mi = groupMap.begin(); mi != groupMap.end(); mi++) {
		clef = mi -> first;
		n = mi -> second;
		if ( clef.compare(0, 8, "OUTGROUP") == 0 ) {
			SequenceSelection ss = MaseTools::getSequenceSet(maseFileHeader, clef);
			for (unsigned int i = 0; i != ss.size(); i++) {
				psc -> setAsOutgroupMember(ss[i]);
			}			
		}
	}
	delete seqc, sitec;
	return(psc);
	} catch(...) {}
}
	 
PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::extractIngroup (const PolymorphismSequenceContainer & psc ) throw (Exception) {
	try {
	SequenceSelection ss;	
	PolymorphismSequenceContainer *psci = dynamic_cast<PolymorphismSequenceContainer *>(psc.clone());
	for(unsigned int i = 0; i < psc.getNumberOfSequences(); i++) {
		if ( psc.isIngroupMember(i) ) {
			ss.push_back(i);
		}
	}
	SequenceContainerTools::keepOnlySelectedSequences(
		* dynamic_cast<const OrderedSequenceContainer *>(
			dynamic_cast<const SiteContainer *>(psci)), ss );
	return( psci );
	} catch(...) {}
}
	 
PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::getSitesWithoutGaps (const PolymorphismSequenceContainer & psc ) throw (Exception) {
	try {
	vector<string> seqNames = psc.getSequencesNames();
	PolymorphismSequenceContainer * noGapCont = new PolymorphismSequenceContainer( psc.getNumberOfSequences(), psc.getAlphabet() );
	noGapCont -> setSequencesNames(seqNames, false);
	unsigned int nbSeq = psc.getNumberOfSequences();
	for (unsigned int i = 0; i < nbSeq; i++) {
		noGapCont -> setSequenceStrength( i, psc.getSequenceStrength(i) );
		if (! psc.isIngroupMember(i))
			noGapCont -> setAsOutgroupMember(i);
	}
	NoGapSiteIterator ngsi(psc);
	while(ngsi.hasMoreSites()) {
		noGapCont -> addSite(* ngsi.nextSite());
	}
	return noGapCont;
	}
	catch(...) {}
}
