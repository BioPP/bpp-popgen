/*
 * File Poplib.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday June 24 2004
 */

#include "Poplib.h"

Poplib::Poplib() {
	setMissingData();
	setDataSeparator();
}

Poplib::Poplib(const string & missing_data, const string & data_separator) throw (Exception) {
	try {
		setMissingData(missing_data);
		setDataSeparator(data_separator);
	}
	catch (Exception & e) {
		throw e;
	}
}

Poplib::~Poplib() {}

const string Poplib::WHITESPACE = "WHITESPACE";
const string Poplib::TAB = "TAB";

const string Poplib::getFormatName() {
	return "Poplib poplib format";
}

const string Poplib::getFormatDescription() {
	return "A very clear description coming soon ;-)";
}

void Poplib::setMissingData(const string & missing_data) throw (Exception) {
	if (missing_data.size() != 1 || isdigit(missing_data[0])
			|| TextTools::isWhiteSpaceCharacter(missing_data[0])
			|| missing_data == TextTools::toString(getDataSeparator()))
		throw Exception("Poplib::setMissingData: not expected value for missing_data.");
	
	_missing_data = missing_data.c_str()[0];
}

void Poplib::setDataSeparator(const string & data_separator) throw (Exception) {

	if (data_separator == WHITESPACE) _data_separator = ' ';
	else if (data_separator == TAB) _data_separator = '\t';
	else {
		if (isdigit(data_separator[0])
				|| data_separator == TextTools::toString(getMissingData())
			 )
			throw Exception("Poplib::setDataSeparator: not expected value for data_separator.");
		_data_separator = data_separator.c_str()[0];
	}
}

string Poplib::getMissingData() const {
	return TextTools::toString(_missing_data);
}

string Poplib::getDataSeparator() const {
	switch (_data_separator) {
		case (' '): return WHITESPACE;
		case ('\t'): return TAB;
		default: return TextTools::toString(_data_separator);
	}
}

char Poplib::getMissingDataChar() const {
	return _missing_data;
}

char Poplib::getDataSeparatorChar() const {
	return _data_separator;
}

void Poplib::read(istream & is, DataSet & data_set) throw (Exception) {
	;
}

void Poplib::read(const string & path, DataSet & data_set) throw (Exception) {
	AbstractIDataSet::read(path, data_set);
}
	
DataSet * Poplib::read(istream & is) throw (Exception) {
	return AbstractIDataSet::read(is);
}

DataSet * Poplib::read(const string & path) throw (Exception) {
	return AbstractIDataSet::read(path);
}

void Poplib::write(ostream & os, const DataSet & data_set) const throw (Exception) {
	unsigned int nbss = data_set.getNumberOfSequenceSets();
	unsigned int seqcpt = 1;
	// General section --------------------------------------
	os << "[General]" << endl;
	os << "MissingData = " << getMissingData() << endl;
	os << "DataSeparator = " << getDataSeparator() << endl;
	if (data_set.hasSequenceData())
		os << "SequenceType = " << data_set.getAlphabet()->getAlphabetType() << endl;
	// Localities section -----------------------------------
	if (data_set.hasLocality()) {
		os << endl << "[Localities]" << endl;
		for (unsigned int i = 0 ; i < data_set.getNumberOfLocalities() ; i++) {
			os << ">" << (data_set.getLocalityByIndex(i))->getName() << endl;
			os << "Coord = " << (data_set.getLocalityByIndex(i))->getX();
			os << " " << (data_set.getLocalityByIndex(i))->getY() << endl;
		}
	}

	// Sequences section ------------------------------------
	if (data_set.hasSequenceData()) {
		Fasta fasta(80);
		os << endl << "[Sequences]" << endl;
		for (unsigned int i = 0 ; i < data_set.getNumberOfGroups() ; i++)
			for (unsigned int j = 0 ; j < data_set.getNumberOfIndividualsInGroup(i) ; j++)
				fasta.write(os, * (data_set.getIndividualByIndexFromGroup(i,j))->getSequences());
	}

	// AllelicData section ----------------------------------
	if (data_set.hasAlleleicData()) {
		os << endl << "[Loci]" << endl;
	}

	// Individuals section ----------------------------------
	os << endl << "[Individuals]" << endl;
	for (unsigned int i = 0 ; i < data_set.getNumberOfGroups() ; i++)
		for (unsigned int j = 0 ; j < data_set.getNumberOfIndividualsInGroup(i) ; j++) {
			if (i>0 || j>0) os << endl;
			const Individual * tmp_ind = data_set.getIndividualByIndexFromGroup(i,j);
			os << ">" << tmp_ind->getId() << endl;
			os << "Group = " << TextTools::toString(i) << endl;
			if (tmp_ind->hasLocality()) os << "Locality = " << tmp_ind->getLocality()->getName() << endl;
			if (tmp_ind->hasCoord()) os << "Coord = " << tmp_ind->getX() << " " << tmp_ind->getY() << endl;
			if (tmp_ind->hasDate()) os << "Date = " << tmp_ind->getDate()->getDateStr() << endl;
			if (tmp_ind->hasSequences()) {
				os << "SequenceData = {" << endl;
				for (unsigned int k = 0 ; k < nbss ; k++) {
					try {
						tmp_ind->getSequenceByIndex(k);
						os << TextTools::toString(seqcpt++);
					}
					catch (SequenceNotFoundException) {
						os << getMissingDataChar();
					}
					if (k < nbss-1) os << getDataSeparatorChar();
					else os << endl;
				}
				os << "}" << endl;
			}
			if (tmp_ind->hasGenotype()) os << "AllelicData = { ... to be continued" << endl;
		}
}

void Poplib::write(const string & path, const DataSet & data_set, bool overwrite) const throw (Exception) {
	AbstractODataSet::write(path, data_set, overwrite);
}
