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
VectorSiteContainer * PolymorphismSequenceContainerTools::read(const string & path, const Alphabet *
alpha) throw (Exception) {
	Mase ms;
	string clef;
	unsigned int n;
	unsigned int i;
	const VectorSequenceContainer *seqc = ms.read( path, alpha );
	VectorSiteContainer *sitec = new VectorSiteContainer( *seqc );
	PolymorphismSequenceContainer *psc = new PolymorphismSequenceContainer( *sitec );
	Comments maseFileHeader = sitec -> getGeneralComments();
	map< string, unsigned int > groupMap = MaseTools::getAvailableSequenceSelections(maseFileHeader);
	for(map< string, unsigned int >::iterator mi = groupMap.begin(); mi != groupMap.end(); mi++) {
		clef = mi -> first;
		n = mi -> second;
		if (! clef.compare(0, 8, "OUTGROUP")) {
			cout << "OUTGROUP" << '\t' << n << endl;
			SequenceSelection ss = MaseTools::getSequenceSet(maseFileHeader, clef);
			for (i = 0; i != ss.size(); i++)
				psc -> toggleIngroup(ss[i]);
		}
	}
	return(psc);
}
	 
PolymorphismSequenceContainer * PolymorphismSequenceContainerTools::extractIngroup (const
PolymorphismSequenceContainer & psc ) throw (Exception) {
	SequenceSelection s;
	PolymorphismSequenceContainer *psci = dynamic_cast<PolymorphismSequenceContainer *>(psc.clone());
	for(unsigned int i=0; i < psc.getNumberOfSequences() ; i++){
		if ( psc.isIngroup(i) ) {
			s.push_back(i);
		}
	}
	SequenceContainerTools::keepOnlySelectedSequences(
		dynamic_cast<const OrderedSequenceContainer &>(
			dynamic_cast<const SiteContainer &>(*psci))
		, s );
	return( psci );
}
