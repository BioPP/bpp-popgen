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
		_effectif.push_back(1);
	}
}

PolymorphismSequenceContainer::PolymorphismSequenceContainer(const SiteContainer & sc) : VectorSiteContainer(sc) {
	for (unsigned int i = 0 ; i < sc.getNumberOfSequences() ; i++) {
		_ingroup.push_back(true);
		_effectif.push_back(1);
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
	_effectif.resize(nbSeq);
	_ingroup.resize(nbSeq);
	for(unsigned int i = 0; i < nbSeq; i++) {
		_comments[i] = new Comments(psc.getComments(i));
		_effectif[i] = getEffectif(i);
		_ingroup[i] = isIngroup(i);
	}
	Sequence * s = NULL;
	_sequences = vector<Sequence *>(nbSeq, s);
}

//** Class destructor: *******************************************************/
PolymorphismSequenceContainer::~PolymorphismSequenceContainer() {
	clear();
}

Clonable* PolymorphismSequenceContainer::clone() const { new PolymorphismSequenceContainer(*this); }

//** Other methodes: *********************************************************/

Sequence * PolymorphismSequenceContainer::removeSequence(unsigned int index) throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::removeSequence: index out of bounds.", index, 0, getNumberOfSequences());
	_effectif.erase(_effectif.begin() + index);
	_ingroup.erase(_ingroup.begin() + index);
	return VectorSiteContainer::removeSequence(index);
}

Sequence * PolymorphismSequenceContainer::removeSequence(const string &name) throw (SequenceNotFoundException) {
	try {
		removeSequence(getSequencePosition(name));
	}
	catch (SequenceNotFoundException snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::removeSequence.", name);
	}
}

void PolymorphismSequenceContainer::deleteSequence(unsigned int index) throw (IndexOutOfBoundsException) {
	try {
		delete removeSequence(index);
	}
	catch (IndexOutOfBoundsException ioobe) {
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::deleteSequence.", index, 0, getNumberOfSequences());
	}
}

void PolymorphismSequenceContainer::deleteSequence(const string &name) throw (SequenceNotFoundException) {
	try {
		delete removeSequence(name);
	}
	catch (SequenceNotFoundException snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::deleteSequence.", name);
	}
}

void PolymorphismSequenceContainer::addSequence(const Sequence &sequence, unsigned int effectif, bool checkNames) throw (Exception) {
	try {
		VectorSiteContainer::addSequence(sequence, checkNames);
	}
	catch (Exception e) {
		throw e;
	}
	_effectif.push_back(effectif);
	_ingroup.push_back(true);
}

void PolymorphismSequenceContainer::clear() {
	VectorSiteContainer::clear();
	_effectif.clear();
	_ingroup.clear();
}

bool PolymorphismSequenceContainer::isIngroup(unsigned int index) const throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::isIngroup: index out of bounds.", index, 0, getNumberOfSequences());
	return _ingroup[index];
}

bool PolymorphismSequenceContainer::isIngroup(const string &name) const throw (SequenceNotFoundException) {
	try {
		return _ingroup[getSequencePosition(name)];
	}
	catch (SequenceNotFoundException snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::isIngroup.", name);
	}
}

bool PolymorphismSequenceContainer::toggleIngroup(unsigned int index) throw (IndexOutOfBoundsException) {
	try {
		if (isIngroup(index))
			_ingroup[index] = false;
		else
			_ingroup[index] = true;
	}
	catch (IndexOutOfBoundsException ioobe) {
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::toggleIngroup.", index, 0, getNumberOfSequences());
	}
}

bool PolymorphismSequenceContainer::toggleIngroup(const string &name) throw (SequenceNotFoundException) {
	try {
		unsigned int seqPos = getSequencePosition(name);
		if (isIngroup(seqPos))
			_ingroup[seqPos] = false;
		else
			_ingroup[seqPos] = true;
	}
	catch (SequenceNotFoundException snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::toggleIngroup.", name);
	}
}

void PolymorphismSequenceContainer::setEffectif(unsigned int index, unsigned int effectif) throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::setEffectif.", index, 0, getNumberOfSequences());
	_effectif[index] = effectif;
}

void PolymorphismSequenceContainer::setEffectif(const string &name, unsigned int effectif) throw (SequenceNotFoundException) {
	try {
		setEffectif(getSequencePosition(name), effectif);
	}
	catch (SequenceNotFoundException snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::setEffectif.", name);
	}
}

unsigned int PolymorphismSequenceContainer::getEffectif(unsigned int index) const throw (IndexOutOfBoundsException) {
	if (index >= getNumberOfSequences())
		throw IndexOutOfBoundsException("PolymorphismSequenceContainer::getEffectif.", index, 0, getNumberOfSequences());
	return _effectif[index];
}

unsigned int PolymorphismSequenceContainer::getEffectif(const string &name) const throw (SequenceNotFoundException) {
	try {
		return getEffectif(getSequencePosition(name));
	}
	catch (SequenceNotFoundException snfe) {
		throw SequenceNotFoundException("PolymorphismSequenceContainer::getEffectif.", name);
	}
}
