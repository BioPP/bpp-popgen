//
// File: PolymorphismSequenceContainer.h
// Authors: bazin <bazin@univ-montp2.fr>
//          Sylvain Gaillard <yragael2001@yahoo.fr>
// Last modification : Wednesday June 16 2004
//

// from PolyLib
#include "PolymorphismSequenceContainer.h"

//** Class constructor: *******************************************************/
PolymorphismSequenceContainer::PolymorphismSequenceContainer(const Alphabet *alpha) :	VectorSiteContainer(alpha) {}

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const OrderedSequenceContainer & sc) : VectorSiteContainer(sc) {
	for (unsigned int i = 0 ; i < sc.getNumberOfSequences() ; i++) {
		_ingroup.push_back(true);
		_strength.push_back(1);
	}
}

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const SiteContainer & sc) : VectorSiteContainer(sc) {
	for (unsigned int i = 0 ; i < sc.getNumberOfSequences() ; i++) {
		_ingroup.push_back(true);
		_strength.push_back(1);
	}
}

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const PolymorphismSequenceContainer & psc) : VectorSiteContainer(psc) {
	unsigned int nbSeq = psc.getNumberOfSequences();
	_strength.resize(nbSeq);
	_ingroup.resize(nbSeq);
	for(unsigned int i = 0; i < nbSeq; i++) {
		_strength[i] = getSequenceStrength(i);
		_ingroup[i] = isIngroupMember(i);
	}
}

PolymorphismSequenceContainer & PolymorphismSequenceContainer::operator= (const PolymorphismSequenceContainer & psc) {
	// Setting up alphabet
	_alphabet = psc.getAlphabet();
	// Setting up general comments
	setGeneralComments(psc.getGeneralComments());
	// Setting up the sequences names
	_names.resize(psc.getNumberOfSequences());
	setSequencesNames(psc.getSequencesNames(), true);
	// Setting up the sites
	for(unsigned int i = 0; i < psc.getNumberOfSites(); i++) {
		addSite(* psc.getSite(i), true);
	}
	// Setting up the sequences comments, numbers and ingroup state
	unsigned int nbSeq = psc.getNumberOfSequences();
	_comments.resize(nbSeq);
	_strength.resize(nbSeq);
	_ingroup.resize(nbSeq);
	for(unsigned int i = 0; i < nbSeq; i++) {
		_comments[i] = new Comments(psc.getComments(i));
		_strength[i] = getSequenceStrength(i);
		_ingroup[i] = isIngroupMember(i);
	}
	Sequence * s = NULL;
	_sequences = vector<Sequence *>(nbSeq, s);
	return * this;
}

//** Class destructor: *******************************************************/
PolymorphismSequenceContainer::~PolymorphismSequenceContainer() {
	clear();
}

Clonable * PolymorphismSequenceContainer::clone() const {
	return dynamic_cast<Clonable *>(dynamic_cast<OrderedSequenceContainer *>
			(dynamic_cast<SiteContainer *>
			 (new PolymorphismSequenceContainer(*this))));
}

//** Other methodes: *********************************************************/

Sequence * PolymorphismSequenceContainer::removeSequence(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::removeSequence: index out of bounds.", index, 0, getNumberOfSequences());
	_strength.erase(_strength.begin() + index);
	_ingroup.erase(_ingroup.begin() + index);
	return VectorSiteContainer::removeSequence(index);
}

Sequence * PolymorphismSequenceContainer::removeSequence(const string &name) throw (SequenceNotFoundException) {
	try {
		removeSequence(getSequencePosition(name));
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::removeSequence.", name);
	}
}

void PolymorphismSequenceContainer::deleteSequence(unsigned int index) throw (IndexOutOfBoundsException) {
	try {
		delete removeSequence(index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::deleteSequence.", index, 0, getNumberOfSequences());
	}
}

void PolymorphismSequenceContainer::deleteSequence(const string &name) throw (SequenceNotFoundException) {
	try {
		delete removeSequence(name);
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::deleteSequence.", name);
	}
}

void PolymorphismSequenceContainer::addSequence(const Sequence &sequence, unsigned int effectif, bool checkNames) throw (Exception) {
	try {
		VectorSiteContainer::addSequence(sequence, checkNames);
	}
	catch (Exception & e) {
		throw e;
	}
	_strength.push_back(effectif);
	_ingroup.push_back(true);
}

void PolymorphismSequenceContainer::clear() {
	VectorSiteContainer::clear();
	_strength.clear();
	_ingroup.clear();
}

bool PolymorphismSequenceContainer::isIngroupMember(unsigned int index) const throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::isIngroupMember: index out of bounds.", index, 0, getNumberOfSequences());
	return _ingroup[index];
}

bool PolymorphismSequenceContainer::isIngroupMember(const string &name) const throw (SequenceNotFoundException) {
	try {
		return _ingroup[getSequencePosition(name)];
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::isIngroupMember.", name);
	}
}

bool PolymorphismSequenceContainer::setAsIngroupMember(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setAsIngroupMember.", index, 0, getNumberOfSequences());
	_ingroup[index] = true;
}

bool PolymorphismSequenceContainer::setAsIngroupMember(const string &name) throw (SequenceNotFoundException) {
	try {
		unsigned int seqPos = getSequencePosition(name);
		_ingroup[seqPos] = true;
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::setAsIngroupMember.", name);
	}
}

bool PolymorphismSequenceContainer::setAsOutgroupMember(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setAsOutgroupMember.", index, 0, getNumberOfSequences());
	_ingroup[index] = false;
}

bool PolymorphismSequenceContainer::setAsOutgroupMember(const string &name) throw (SequenceNotFoundException) {
	try {
		unsigned int seqPos = getSequencePosition(name);
		_ingroup[seqPos] = false;
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::setAsOutgroupMember.", name);
	}
}

void PolymorphismSequenceContainer::setSequenceStrength(unsigned int index, unsigned int strength) throw (Exception) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setSequenceStrength.", index, 0, getNumberOfSequences());
	if (strength < 1)
		throw BadIntegerException("PolymorphismSequenceContainer::setSequenceStrength: strength can't be < 1.", strength);
	_strength[index] = strength;
}

void PolymorphismSequenceContainer::setSequenceStrength(const string &name, unsigned int strength) throw (Exception) {
	try {
		setSequenceStrength(getSequencePosition(name), strength);
	}
	catch (BadIntegerException & bie) {
		throw bie;
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::setSequenceStrength.", name);
	}
}

void PolymorphismSequenceContainer::incrementSequenceStrength(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::incrementSequenceStrength.", index, 0, getNumberOfSequences());
	_strength[index]++;
}

void PolymorphismSequenceContainer::incrementSequenceStrength(const string &name) throw (SequenceNotFoundException) {
	try {
		incrementSequenceStrength(getSequencePosition(name));
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::incrementSequenceStrength.", name);
	}
}

void PolymorphismSequenceContainer::decrementSequenceStrength(unsigned int index) throw (Exception) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::decrementSequenceStrength.", index, 0, getNumberOfSequences());
	if (_strength[index]-1 < 1)
		throw BadIntegerException("PolymorphismSequenceContainer::decrementSequenceStrength: strength can't be < 1.", _strength[index]-1);
	_strength[index]--;
}

void PolymorphismSequenceContainer::decrementSequenceStrength(const string &name) throw (Exception) {
	try {
		decrementSequenceStrength(getSequencePosition(name));
	}
	catch (BadIntegerException & bie) {
		throw bie;
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::decrementSequenceStrength.", name);
	}
}

unsigned int PolymorphismSequenceContainer::getSequenceStrength(unsigned int index) const throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::getSequenceStrength.", index, 0, getNumberOfSequences());
	return _strength[index];
}

unsigned int PolymorphismSequenceContainer::getSequenceStrength(const string &name) const throw (SequenceNotFoundException) {
	try {
		return getSequenceStrength(getSequencePosition(name));
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::getSequenceStrength.", name);
	}
}
