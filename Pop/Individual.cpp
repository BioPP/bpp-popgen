/*
 * File Individual.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday August 03 2004
 *
*/
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
#include "Individual.h"

//** Class constructor: *******************************************************/
Individual::Individual() {
	_id = "";
	_sex = 0;
	_date = NULL;
	_coord = NULL;
	_locality = NULL;
	_sequences = NULL;
	_genotype = NULL;
}

Individual::Individual(const string & id) {
	_id = id;
	_sex = 0;
	_date = NULL;
	_coord = NULL;
	_locality = NULL;
	_sequences = NULL;
	_genotype = NULL;
}

Individual::Individual(const string & id,
                       const Date & date,
                       const Coord<double> & coord,
                       Locality<double> * locality,
                       const unsigned short sex) {
	_id = id;
	_sex = sex;
	_date = new Date(date);
	_coord = new Coord<double>(coord);
	_locality = locality;
	_sequences = NULL;
	_genotype = NULL;
}

Individual::Individual(const Individual &ind) {
	setId(ind.getId());
	setSex(ind.getSex());
	_date = NULL;
	try {
		setDate(* ind.getDate());
	}
	catch (...) {}
	_coord = NULL;
	try {
		setCoord(* ind.getCoord());
	}
	catch (...) {}
	_locality = NULL;
	try {
		setLocality(ind.getLocality());
	}
	catch (...) {}
	_sequences = NULL;
	try {
		setSequences(* dynamic_cast<const MapSequenceContainer *>(ind.getSequences()));
	}
	catch (...) {}
	_genotype = NULL;
	if (ind.hasGenotype())
		_genotype = new MultilocusGenotype(* ind.getGenotype());
}

//** Class destructor: *******************************************************/
Individual::~Individual () {
	if (_date != NULL) {
		delete _date;
		_date = NULL;
	}
	if (_coord != NULL) {
		delete _coord;
		_coord = NULL;
	}
	if (_sequences != NULL) {
		delete _sequences;
		_sequences = NULL;
	}
	if (_genotype != NULL) {
		delete _genotype;
		_genotype = NULL;
	}
}

//** Other methodes: *********************************************************/
Individual & Individual::operator= (const Individual & ind) {
	setId(ind.getId());
	setSex(ind.getSex());
	try {
		setDate(* ind.getDate());
	}
	catch (NullPointerException) {
		_date = NULL;
	}
	try {
		setCoord(* ind.getCoord());
	}
	catch (NullPointerException) {
		_coord = NULL;
	}
	try {
		setLocality(ind.getLocality());
	}
	catch (NullPointerException) {
		_locality = NULL;
	}
	try {
		setSequences(* dynamic_cast<const MapSequenceContainer *>(ind.getSequences()));
	}
	catch (NullPointerException) {
		_sequences = NULL;
	}
	this->_genotype = ind.hasGenotype() ? new MultilocusGenotype(* ind.getGenotype()) : NULL;
	return * this;
}

// Id
void Individual::setId(const string id) {
	_id = id;
}

string Individual::getId() const {
	return _id;
}

// Sex
void Individual::setSex(const unsigned short sex) {
	_sex = sex;
}

unsigned short Individual::getSex() const {
	return _sex;
}

// Date
void Individual::setDate(const Date & date) {
	if (!hasDate()) {
		_date = new Date(date);
	}
	else if (* _date != date) {
		delete _date;
		_date = new Date(date);
	}
}

const Date * Individual::getDate() const throw (NullPointerException) {
	if (hasDate())
		return new Date(* _date);
	else
		throw(NullPointerException("Individual::getDate: no date associated to this individual."));
}

bool Individual::hasDate() const {
	return _date != NULL;
}

// Coord
void Individual::setCoord(const Coord<double> & coord) {
	if (!hasCoord()) {
		_coord = new Coord<double>(coord);
	}
	else if	(* _coord != coord) {
		delete _coord;
		_coord = new Coord<double>(coord);
	}
}

void Individual::setCoord(const double x, const double y) {
	if (!hasCoord()) {
		_coord = new Coord<double>(x, y);
	}
	else if (this->getX() != x || this->getY() != y) {
		delete _coord;
		_coord = new Coord<double>(x, y);
	}
}

const Coord<double> * Individual::getCoord() const throw(NullPointerException) {
	if (hasCoord())
		return new Coord<double>(* _coord);
	else
		throw(NullPointerException("Individual::getCoord: no coord associated to this individual."));
}

bool Individual::hasCoord() const {
	return _coord != NULL;
}

void Individual::setX(const double x) throw(NullPointerException) {
	if (hasCoord())
		_coord->setX(x);
	else
		throw(NullPointerException("Individual::setX: no coord associated to this individual."));
}

void Individual::setY(const double y) throw(NullPointerException) {
	if (hasCoord())
		_coord->setY(y);
	else
		throw(NullPointerException("Individual::setY: no coord associated to this individual."));
}

double Individual::getX() const throw(NullPointerException) {
	if (hasCoord())
		return _coord->getX();
	else
		throw(NullPointerException("Individual::getX: no coord associated to this individual."));
}

double Individual::getY() const throw(NullPointerException) {
	if (hasCoord())
		return _coord->getY();
	else
		throw(NullPointerException("Individual::getY: no coord associated to this individual."));
}

// Locality
void Individual::setLocality(const Locality<double> * locality) {
	_locality = locality;
}

const Locality<double> * Individual::getLocality() const  throw (NullPointerException) {
	if (hasLocality())
		return _locality;
	else
		throw(NullPointerException("Individual::getLocality: no locality associated to this individual."));
}

bool Individual::hasLocality() const {
	return _locality != NULL;
}

// Sequences
void Individual::addSequence(unsigned int sequence_key, const Sequence & sequence)
throw (Exception) {
	if (_sequences == NULL)
		_sequences = new MapSequenceContainer(sequence.getAlphabet());
	try {
		_sequences->addSequence(TextTools::toString(sequence_key), sequence);
	}
	catch (AlphabetMismatchException & ame)
	{
		throw(AlphabetMismatchException("Individual::addSequence: alphabets don't match.", ame.getAlphabets()[0], ame.getAlphabets()[1]));
	}
	catch (Exception & e)
	{
		if (string(e.what()).find("name") < string(e.what()).size())
			throw(BadIdentifierException("Individual::addSequence: sequence's name already in use.", sequence.getName()));
		// if (string(e.what()).find("key") < string(e.what()).size())
		else
			throw(BadIntegerException("Individual::addSequence: sequence_key already in use.", sequence_key));
	}
}

const Sequence * Individual::getSequenceByName(const string & sequence_name)
const throw (Exception) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::getSequenceByName: no sequence data.");
	try {
		return const_cast<const MapSequenceContainer *>(_sequences)->getSequence(sequence_name);
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Individual::getSequenceByName: sequence_name not found.", snfe.getSequenceId());
	}
}

const Sequence * Individual::getSequenceAtPosition(unsigned int sequence_position)
const throw (Exception) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::getSequenceAtPosition: no sequence data.");
	try {
		return const_cast<const MapSequenceContainer *>(_sequences)->getSequenceByKey(TextTools::toString(sequence_position));
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Individual::getSequenceAtPosition: sequence_position not found", snfe.getSequenceId());
	}
}

void Individual::deleteSequenceByName(const string & sequence_name) throw (Exception) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::deleteSequenceByName: no sequence data.");
	try {
		_sequences->deleteSequence(sequence_name);
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Individual::deleteSequenceByName: sequence_name not found.", snfe.getSequenceId());
	}
}

void Individual::deleteSequenceAtPosition(unsigned int sequence_position) throw (Exception) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::deleteSequenceAtPosition: no sequence data.");
	try {
		_sequences->deleteSequenceByKey(TextTools::toString(sequence_position));
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Individual::deleteSequenceAtPosition: sequence_position not found.", snfe.getSequenceId());
	}
}

vector<string> Individual::getSequencesNames() const throw (NullPointerException) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::getSequencesNames: no sequence data.");
	return _sequences->getSequencesNames();
}

vector<unsigned int> Individual::getSequencesPositions() const throw (NullPointerException) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::getSequencesPositions: no sequence data.");
	vector<unsigned int> seqpos;
	vector<string> seqkeys = _sequences->getKeys();
	for (unsigned int i = 0 ; i < seqkeys.size() ; i++)
		seqpos.push_back((unsigned int) TextTools::toInt(seqkeys[i]));
	return seqpos;
}

unsigned int Individual::getSequencePosition(const string & sequence_name) const throw (Exception) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::getSequencePosition: no sequence data.");
	try {
		return (unsigned int) TextTools::toInt(_sequences->getKey(getSequencePosition(sequence_name)));
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Individual::getSequencePosition: sequence_name not found.", snfe.getSequenceId());
	}
}
	
bool Individual::hasSequences() const {
	return !(getNumberOfSequences() == 0);
}

bool Individual::hasSequenceAtPosition(unsigned int position) const {
	if (hasSequences()) {
		vector<unsigned int> pos = getSequencesPositions();
		for (unsigned int i = 0 ; i < pos.size() ; i++)
			if (pos[i] == position)
				return true;
	}
	return false;
}

const Alphabet * Individual::getSequenceAlphabet() const throw (NullPointerException) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::getSequenceAlphabet: no sequence data.");
	return _sequences->getAlphabet();
}

unsigned int Individual::getNumberOfSequences() const {
	if (_sequences == NULL)
		return 0;
	return const_cast<const MapSequenceContainer *>(_sequences)->getNumberOfSequences();
}

void Individual::setSequences(const MapSequenceContainer & msc) {
	if (hasSequences()) {
		delete(_sequences);
		_sequences = NULL;
	}
	_sequences = new MapSequenceContainer(msc);
}

const OrderedSequenceContainer * Individual::getSequences() const throw (NullPointerException) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::getSequences: no sequence data.");
	return _sequences;
}

// MultilocusGenotype

void Individual::setGenotype(const MultilocusGenotype & genotype) {
	if (hasGenotype())
		delete _genotype;
	_genotype = new MultilocusGenotype(genotype);
}

void Individual::initGenotype(unsigned int loci_number) throw (Exception) {
	if (hasGenotype())
		throw Exception("Individual::initGenotype: individual already has a genotype.");
	try {
		_genotype = new MultilocusGenotype(loci_number);
	}
	catch (BadIntegerException & bie) {
		throw BadIntegerException("Individual::initGenotype: loci_number must be > 0.", bie.getBadInteger());
	}
}

const MultilocusGenotype * Individual::getGenotype() const throw (NullPointerException) {
	if (!hasGenotype())
		throw NullPointerException("Individual::getGenotype: individual has no genotype.");
	return _genotype;
}

void Individual::deleteGenotype() {
	if (hasGenotype()) delete _genotype;
}

bool Individual::hasGenotype() const {
	return _genotype != NULL;
}

void Individual::setMonolocusGenotype(unsigned int locus_position, const MonolocusGenotype & monogen) throw (Exception) {
	if (!hasGenotype())
		throw NullPointerException("Individual::setMonolocusGenotype: individual has no genotype.");
	try {
		_genotype->setMonolocusGenotype(locus_position, monogen);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Individual::setMonolocusGenotype: locus_position out of boubds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

void Individual::setMonolocusGenotypeByAlleleKey(unsigned int locus_position, const vector<unsigned int> allele_keys) throw (Exception) {
	if (!hasGenotype())
		throw NullPointerException("Individual::setMonolocusGenotypeByAlleleKey: individual has no genotype.");
	try {
		_genotype->setMonolocusGenotypeByAlleleKey(locus_position, allele_keys);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Individual::setMonolocusGenotypeByAlleleKey: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (Exception) {
		throw Exception("Individual::setMonolocusGenotypeByAlleleKey: no key in allele_keys.");
	}
}

void Individual::setMonolocusGenotypeByAlleleId(unsigned int locus_position, const vector<string> allele_id, const LocusInfo & locus_info) throw (Exception) {
	if (!hasGenotype())
		throw NullPointerException("Individual::setMonolocusGenotypeByAlleleId: individual has no genotype.");
	try {
		_genotype->setMonolocusGenotypeByAlleleId(locus_position, allele_id, locus_info);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Individual::setMonolocusGenotypeByAlleleId: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (AlleleNotFoundException & anfe) {
		throw AlleleNotFoundException("Individual::setMonolocusGenotypeByAlleleId: id not found.", anfe.getIdentifier());
	}
}

const MonolocusGenotype * Individual::getMonolocusGenotype(unsigned int locus_position) throw (Exception) {
	if (!hasGenotype())
		throw NullPointerException("Individual::getMonolocusGenotype: individual has no genotype.");
	try {
		return _genotype->getMonolocusGenotype(locus_position);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Individual::getMonolocusGenotype: locus_position out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

unsigned int Individual::countNonMissingLoci() const throw (NullPointerException) {
	if (!hasGenotype())
		throw NullPointerException("Individual::countNonMissingLoci: individual has no genotype.");
	return _genotype->countNonMissingLoci();
}

unsigned int Individual::countHomozygousLoci() const throw (NullPointerException) {
	if (!hasGenotype())
		throw NullPointerException("Individual::countHomozygousLoci: individual has no genotype.");
	return _genotype->countHomozygousLoci();
}

unsigned int Individual::countHeterozygousLoci() const throw (NullPointerException) {
	if (!hasGenotype())
		throw NullPointerException("Individual::countHeterozygousLoci: individual has no genotype.");
	return _genotype->countHeterozygousLoci();
}
