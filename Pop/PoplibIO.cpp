/*
 * File PoplibIO.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday July 07 2004
 */

#include "PoplibIO.h"

const string PoplibIO::WHITESPACE = "WHITESPACE";
const string PoplibIO::TAB = "TAB";
const string PoplibIO::COMA = "COMA";
const string PoplibIO::SEMICOLON = "SEMICOLON";

PoplibIO::PoplibIO() {
	setDataSeparator(WHITESPACE);
	setMissingDataSymbol("$");
}

PoplibIO::PoplibIO(const string & missing_data_symbol, const string & data_separator) throw (Exception) {
	try {
		setDataSeparator(data_separator);
		setMissingDataSymbol(missing_data_symbol);
	}
	catch (Exception & e) {
		throw e;
	}
}

PoplibIO::~PoplibIO() {}

const string PoplibIO::getFormatName() {
	return "PoplibIO ver 0.1";
}

const string PoplibIO::getFormatDescription() {
	return "A very clear description coming soon ;-)";
}

void PoplibIO::setMissingDataSymbol(const string & missing_data_symbol) throw (Exception) {
	if (missing_data_symbol.size() != 1 || isdigit(missing_data_symbol[0])
			|| TextTools::isWhiteSpaceCharacter(missing_data_symbol[0])
			|| missing_data_symbol[0] == _data_separator
			)
		throw Exception("PoplibIO::setMissingData: not expected value for missing_data_symbol.");
	
	_missing_data_symbol = missing_data_symbol[0];
}

void PoplibIO::setDataSeparator(const string & data_separator) throw (Exception) {

	if (data_separator == WHITESPACE) _data_separator = ' ';
	else if (data_separator == TAB) _data_separator = '\t';
	else if (data_separator == COMA) _data_separator = ',';
	else if (data_separator == SEMICOLON) _data_separator = ';';
	else {
		if (isdigit(data_separator[0])
				|| data_separator == getMissingDataSymbol()
			 )
			throw Exception("PoplibIO::setDataSeparator: not expected value for data_separator.");
		_data_separator = data_separator.c_str()[0];
	}
}

string PoplibIO::getMissingDataSymbol() const {
	return TextTools::toString(_missing_data_symbol);
}

string PoplibIO::getDataSeparator() const {
	switch (_data_separator) {
		case (' '): return WHITESPACE;
		case ('\t'): return TAB;
		case (','): return COMA;
		case (';'): return SEMICOLON;
		default: return TextTools::toString(_data_separator);
	}
}

char PoplibIO::getMissingDataChar() const {
	return _missing_data_symbol;
}

char PoplibIO::getDataSeparatorChar() const {
	return _data_separator;
}

void PoplibIO::read(istream & is, DataSet & data_set) throw (Exception) {
	if (!is)
		throw IOException("PoplibIO::read: fail to open stream.");
	string temp = "";
	unsigned int current_section = 0;
	// Main loop for all file lines
	while (!is.eof()) {
		getline(is, temp, '\n'); // Copy the line in a temporary string
		// Get the correct current section
		if (temp.find("[General]", 0) != string::npos)
			current_section = 1;
		if (temp.find("[Localities]", 0) != string::npos)
			current_section = 2;
		if (temp.find("[Sequences]", 0) != string::npos)
			current_section = 3;
		if (temp.find("[Loci]", 0) != string::npos)
			current_section = 4;
		if (temp.find("[Individuals]", 0) != string::npos)
			current_section = 5;
		// General section ------------------------------------
		if (current_section == 1) {
			cout << "General section" << endl;
			if (temp.find("MissingData", 0) != string::npos)
				setMissingDataSymbol(_getValues(temp, "=")[0]);
			if (temp.find("DataSeparator", 0) != string::npos)
				setDataSeparator(_getValues(temp, "=")[0]);
			if (temp.find("SequenceType", 0) != string::npos)
				// *** Must set the AnalyzedSequences alphabet ***
				cout << _getValues(temp, "=")[0] << endl;
		}
		// Localities section ---------------------------------
		if (current_section == 2) {
			cout << "Localities section" << endl;
			
		}
	}
}

void PoplibIO::read(const string & path, DataSet & data_set) throw (Exception) {
	AbstractIDataSet::read(path, data_set);
}
	
DataSet * PoplibIO::read(istream & is) throw (Exception) {
	return AbstractIDataSet::read(is);
}

DataSet * PoplibIO::read(const string & path) throw (Exception) {
	return AbstractIDataSet::read(path);
}

void PoplibIO::write(ostream & os, const DataSet & data_set) const throw (Exception) {
	unsigned int seqcpt = 1;
	// General section --------------------------------------
	os << "[General]" << endl;
	os << "MissingData = " << getMissingDataSymbol() << endl;
	os << "DataSeparator = " << getDataSeparator() << endl;
	if (data_set.hasSequenceData())
		os << "SequenceType = " << data_set.getAlphabet()->getAlphabetType() << endl;
	// Localities section -----------------------------------
	if (data_set.hasLocality()) {
		os << endl << "[Localities]" << endl;
		for (unsigned int i = 0 ; i < data_set.getNumberOfLocalities() ; i++) {
			os << ">" << (data_set.getLocalityAtPosition(i))->getName() << endl;
			os << "Coord = " << (data_set.getLocalityAtPosition(i))->getX();
			os << " " << (data_set.getLocalityAtPosition(i))->getY() << endl;
		}
	}

	// Sequences section ------------------------------------
	if (data_set.hasSequenceData()) {
		Fasta fasta(80);
		os << endl << "[Sequences]" << endl;
		for (unsigned int i = 0 ; i < data_set.getNumberOfGroups() ; i++)
			for (unsigned int j = 0 ; j < data_set.getNumberOfIndividualsInGroup(i) ; j++)
				fasta.write(os, * (data_set.getIndividualAtPositionFromGroup(i,j))->getSequences());
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
			const Individual * tmp_ind = data_set.getIndividualAtPositionFromGroup(i,j);
			os << ">" << tmp_ind->getId() << endl;
			os << "Group = " << TextTools::toString(i) << endl;
			if (tmp_ind->hasLocality()) os << "Locality = " << tmp_ind->getLocality()->getName() << endl;
			if (tmp_ind->hasCoord()) os << "Coord = " << tmp_ind->getX() << " " << tmp_ind->getY() << endl;
			if (tmp_ind->hasDate()) os << "Date = " << tmp_ind->getDate()->getDateStr() << endl;
			if (tmp_ind->hasSequences()) {
				unsigned int nbss = tmp_ind->getNumberOfSequences();
				os << "SequenceData = {" << endl;
				for (unsigned int k = 0 ; k < nbss ; k++) {
					try {
						tmp_ind->getSequenceAtPosition(k);
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

void PoplibIO::write(const string & path, const DataSet & data_set, bool overwrite) const throw (Exception) {
	AbstractODataSet::write(path, data_set, overwrite);
}

vector<string> PoplibIO::_getValues(string & param_line, const string & delim) {
	vector<string> values;
	int limit = param_line.find(delim, 0);
	if (limit >= 0)
		param_line = string(param_line.begin() + limit + delim.size(), param_line.end());
	param_line = TextTools::removeSurroundingWhiteSpaces(param_line);
	
	int bi = 0;
	int bs = param_line.find(getDataSeparatorChar(), bi);
	while (bs > 0) {
		values.push_back(string(param_line.begin() + bi, param_line.begin() + bs));
		bi = bs + 1;
		bs = param_line.find(getDataSeparatorChar(), bi);
	}
	values.push_back(string(param_line.begin() + bi, param_line.end()));
	return values;
}
