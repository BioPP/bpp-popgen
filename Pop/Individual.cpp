/*
 * File Individual.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday June 24 2004
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
	this->_genotype = ind.hasGenotype() ? new Genotype(* ind.getGenotype()) : NULL;
}

//** Class destructor: *******************************************************/
Individual::~Individual () {
	delete this->_date;
	delete this->_coord;
}

//** Other methodes: *********************************************************/
Clonable * Individual::clone() const {
	return new Individual(* this);
}

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
	this->_genotype = ind.hasGenotype() ? new Genotype(* ind.getGenotype()) : NULL;
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

const Sequence * Individual::getSequenceByIndex(unsigned int sequence_index)
const throw (Exception) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::getSequenceByIndex: no sequence data.");
	try {
		return const_cast<const MapSequenceContainer *>(_sequences)->getSequenceByKey(TextTools::toString(sequence_index));
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Individual::getSequenceByIndex: sequence_index not found", snfe.getSequenceId());
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

void Individual::deleteSequenceByIndex(unsigned int sequence_index) throw (Exception) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::deleteSequenceByIndex: no sequence data.");
	try {
		_sequences->deleteSequenceByKey(TextTools::toString(sequence_index));
	}
	catch (SequenceNotFoundException & snfe) {
		throw SequenceNotFoundException("Individual::deleteSequenceByIndex: sequence_index not found.", snfe.getSequenceId());
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
	return (_sequences != NULL && _sequences->getNumberOfSequences() != 0);
}

const Alphabet * Individual::getSequenceAlphabet() const throw (NullPointerException) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::getSequenceAlphabet: no sequence data.");
	return _sequences->getAlphabet();
}

unsigned int Individual::getNumberOfSequences() const throw (NullPointerException) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::getNumberOfSequences: no sequence data.");
	return const_cast<const MapSequenceContainer *>(_sequences)->getNumberOfSequences();
}

void Individual::setSequences(const MapSequenceContainer & msc) {
	if (_sequences != NULL)
		delete(_sequences);
	_sequences = new MapSequenceContainer(msc);
}

const OrderedSequenceContainer * Individual::getSequences() const throw (NullPointerException) {
	if (_sequences == NULL)
		throw NullPointerException("Individual::getSequences: no sequence data.");
	return _sequences;
}

// Genotype

void Individual::addGenotype(const Genotype & genotype) throw (Exception) {
	if (hasGenotype())
		throw Exception("Individual::addGenotype: individual already has a genotype.");
	_genotype = new Genotype(genotype);
}

void Individual::initGenotype(const AnalyzedLoci * analyzed_loci) throw (Exception) {
	if (hasGenotype())
		throw Exception("Individual::initGenotype: individual already has a genotype.");
	try {
		_genotype = new Genotype(analyzed_loci);
	}
	catch (NullPointerException) {
		throw NullPointerException("Individual::initGenotype: analyzed_loci is NULL.");
	}
}

const Genotype * Individual::getGenotype() const throw (NullPointerException) {
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

void Individual::setMonolocusGenotype(unsigned int locus_index, const MonolocusGenotype & monogen) throw (Exception) {
	if (!hasGenotype())
		throw NullPointerException("Individual::setMonolocusGenotype: individual has no genotype.");
	try {
		_genotype->setMonolocusGenotype(locus_index, monogen);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Individual::setMonolocusGenotype: locus_index out of boubds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

void Individual::setMonolocusGenotypeByAlleleKey(unsigned int locus_index, const vector<unsigned int> allele_keys) throw (Exception) {
	if (!hasGenotype())
		throw NullPointerException("Individual::setMonolocusGenotypeByAlleleKey: individual has no genotype.");
	try {
		_genotype->setMonolocusGenotypeByAlleleKey(locus_index, allele_keys);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Individual::setMonolocusGenotypeByAlleleKey: locus_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (Exception) {
		throw Exception("Individual::setMonolocusGenotypeByAlleleKey: allele_keys.size() doesn't match ploidy.");
	}
}

void Individual::setMonolocusGenotypeByAlleleId(unsigned int locus_index, const vector<unsigned int> allele_id) throw (Exception) {
	if (!hasGenotype())
		throw NullPointerException("Individual::setMonolocusGenotypeByAlleleId: individual has no genotype.");
	try {
		_genotype->setMonolocusGenotypeByAlleleId(locus_index, allele_id);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Individual::setMonolocusGenotypeByAlleleId: locus_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
	catch (Exception) {
		throw Exception("Individual::setMonolocusGenotypeByAlleleId: allele_keys.size() doesn't match ploidy.");
	}
}

unsigned int Individual::getPloidy(unsigned int locus_index) throw (Exception) {
	if (!hasGenotype())
		throw NullPointerException("Individual::getPloidy: individual has no genotype.");
	try {
		return _genotype->getPloidy(locus_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Individual::getPloidy: locus_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}

const MonolocusGenotype * Individual::getMonolocusGenotype(unsigned int locus_index) throw (Exception) {
	if (!hasGenotype())
		throw NullPointerException("Individual::getMonolocusGenotype: individual has no genotype.");
	try {
		return _genotype->getMonolocusGenotype(locus_index);
	}
	catch (IndexOutOfBoundsException & ioobe) {
		throw IndexOutOfBoundsException("Individual::getMonolocusGenotype: locus_index out of bounds.", ioobe.getBadInteger(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
	}
}
