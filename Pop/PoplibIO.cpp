/*
 * File PoplibIO.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday July 21 2004
 */

#include "PoplibIO.h"

const string PoplibIO::WHITESPACE = string("WHITESPACE");
const string PoplibIO::TAB = string("TAB");
const string PoplibIO::COMA = string("COMA");
const string PoplibIO::SEMICOLON = string("SEMICOLON");

const string PoplibIO::DIPLOID = string("DIPLOID");
const string PoplibIO::HAPLOID = string("HAPLOID");
const string PoplibIO::HAPLODIPLOID = string("HAPLODIPLOID");

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
	vector<string> temp_v;
	stringstream tmp_ss;
	VectorSequenceContainer * tmp_vsc = NULL;
	Locality<double> tmp_locality("tmp");
	vector<LocusInfo> tmp_locinf;
	Individual tmp_indiv;
	bool section1 = true;
	bool section2 = true;
	bool section3 = true;
	bool section4 = true;
	bool section5 = true;
	unsigned int current_section = 0;
	unsigned int previous_section = 0;
	unsigned int linenum = 0;
	// Main loop for all file lines
	while (!is.eof()) {
		temp = FileTools::getNextLine(is);
		linenum++;
		// Get the correct current section
		if (temp.find("[General]", 0) != string::npos) {
			previous_section = current_section;
			current_section = 1;
			continue;
		}
		else if (temp.find("[Localities]", 0) != string::npos) {
			previous_section = current_section;
			current_section = 2;
			continue;
		}
		else if (temp.find("[Sequences]", 0) != string::npos) {
			previous_section = current_section;
			current_section = 3;
			continue;
		}
		else if (temp.find("[Loci]", 0) != string::npos) {
			previous_section = current_section;
			current_section = 4;
			continue;
		}
		else if (temp.find("[Individuals]", 0) != string::npos) {
			previous_section = current_section;
			current_section = 5;
			continue;
		}
		// General section ------------------------------------
		if (current_section == 1 && previous_section < 1) {
			temp_v.push_back(temp);
		}
		if (section1 && current_section != 1 && previous_section == 1) {
			section1 = false;
			_parseGeneral(temp_v, data_set);
			temp_v.clear();
			if (data_set.hasSequenceData() && tmp_vsc == NULL)
				tmp_vsc = new VectorSequenceContainer(data_set.getAlphabet());
		}

		// Localities section ---------------------------------
		if (current_section == 2 && previous_section < 2) {
			if (temp.find(">", 0) != string::npos) {
				_parseLocality(temp_v, data_set);
				temp_v.clear();
				temp_v.push_back(temp);
			}
			else
				temp_v.push_back(temp);
		}
		if (section2 && current_section !=2 && previous_section == 2) {
			section2 = false;
			_parseLocality(temp_v, data_set);
			temp_v.clear();
		}	

		// Sequences section ----------------------------------
		if (current_section == 3 && previous_section < 3) {
			if (temp.find(">", 0) != string::npos) {
				_parseSequence(temp_v, *tmp_vsc);
				temp_v.clear();
				temp_v.push_back(temp);
			}
			else
				temp_v.push_back(temp);
		}
		if (section3 && current_section !=3 && previous_section == 3) {
			section3 = false;
			_parseSequence(temp_v, *tmp_vsc);
			temp_v.clear();
		}

		// Loci section ---------------------------------------
		if (current_section == 4 && previous_section <4) {
			if (temp.find(">", 0) != string::npos) {
				_parseLoci(temp_v, tmp_locinf);
				temp_v.clear();
				temp_v.push_back(temp);
			}
			else
				temp_v.push_back(temp);
		}
		if (section4 && current_section != 4 && previous_section == 4) {
			section4 = false;
			_parseLoci(temp_v, tmp_locinf);
			temp_v.clear();
			AnalyzedLoci tmp_anloc(tmp_locinf.size());
			for (unsigned int i = 0 ; i < tmp_locinf.size() ; i++)
				tmp_anloc.setLocusInfo(i, tmp_locinf[i]);
			cout << "tmp_anloc.getNumberOfLoci() = " << tmp_anloc.getNumberOfLoci() << endl;
			for (unsigned int i = 0 ; i < tmp_anloc.getNumberOfLoci() ; i++)
				cout << "Locus " << TextTools::toString(i) << " : " << tmp_anloc.getLocusInfoAtPosition(i)->getName() << endl;
			data_set.setAnalyzedLoci(tmp_anloc);
			try {
				cout << "data_set.getNumberOfLoci() = " << data_set.getNumberOfLoci() << endl;
			}
			catch (Exception & e) {
				cerr << "### ERROR ###" << endl << e.what() << endl;
			}
			cout << "DataSet OK" << endl;
			cout << "_______________________________" << endl;
		}
		
		// Individuals section --------------------------------
		if (current_section == 5 && previous_section < 5) {
			if (temp.find(">", 0) != string::npos) {
				_parseIndividual(temp_v, data_set, *tmp_vsc);
				temp_v.clear();
				temp_v.push_back(temp);
			}
			else
				temp_v.push_back(temp);
		}
		if (section5 && current_section != 5 && previous_section == 5) {
			section5 = false;
			_parseIndividual(temp_v, data_set, *tmp_vsc);
			temp_v.clear();
		}
	}
	// Emptied the buffer if eof.
	if (section2 && current_section == 2)
		_parseLocality(temp_v, data_set);
	if (section3 && current_section == 3)
		_parseSequence(temp_v, *tmp_vsc);
	if (section5 && current_section == 5)
		_parseIndividual(temp_v, data_set, *tmp_vsc);
	temp_v.clear();
}

void PoplibIO::_parseGeneral(const vector<string> & in, DataSet & data_set) {
	stringstream is;
	for (unsigned int i = 0 ; i < in.size() ; i++)
		is << in[i] << endl;
	string temp;
	while (!is.eof() && in.size() != 0) {
		temp = FileTools::getNextLine(is);
		if (temp.find("MissingData", 0) != string::npos)
			setMissingDataSymbol(_getValues(temp, "=")[0]);
		if (temp.find("DataSeparator", 0) != string::npos)
			setDataSeparator(_getValues(temp, "=")[0]);
		if (temp.find("SequenceType", 0) != string::npos)
			data_set.setAlphabet(_getValues(temp, "=")[0]);
	}
}

void PoplibIO::_parseLocality(const vector<string> & in, DataSet & data_set) {
	stringstream is;
	for (unsigned int i = 0 ; i < in.size() ; i++)
		is << in[i] << endl;
	Locality<double> tmp_locality("");
	string temp;
	while (!is.eof() && in.size() != 0) {
		temp = FileTools::getNextLine(is);
//		cout << "_parseLocality: " << temp << endl;
		if (temp.find(">", 0) != string::npos) {
			tmp_locality.setName(TextTools::removeSurroundingWhiteSpaces(string(temp.begin() + 1, temp.end())));
		}
		if (temp.find("Coord", 0) != string::npos) {
			vector<string> v = _getValues(temp, "=");
			tmp_locality.setX(TextTools::toDouble(v[0]));
			tmp_locality.setY(TextTools::toDouble(v[1]));
		}
	}
	if (tmp_locality.getName() != "")
		data_set.addLocality(tmp_locality);
}

void PoplibIO::_parseSequence(const vector<string> & in, VectorSequenceContainer & vsc) {
	Fasta ifasta;
	stringstream is;
	for (unsigned int i = 0 ; i < in.size() ; i++)
		is << in[i] << endl;
	ifasta.read(is, vsc);
}

void PoplibIO::_parseLoci(const vector<string> & in, vector<LocusInfo> & locus_info) {
	stringstream is;
	for (unsigned int i = 0 ; i < in.size() ; i++)
		is << in[i] << endl;
	string locinf_name = "";
	unsigned int locinf_ploidy = LocusInfo::DIPLOID;
	string temp;
	while (!is.eof()) {
		temp = FileTools::getNextLine(is);
		if (temp.find(">", 0) != string::npos) {
			locinf_name = TextTools::removeSurroundingWhiteSpaces(string(temp.begin() + 1, temp.end()));
		}
		if (temp.find("Ploidy", 0) != string::npos) {
			vector<string> v = _getValues(temp, "=");
			string tmp_str_ploidy = TextTools::removeSurroundingWhiteSpaces(v[0]);
			tmp_str_ploidy = TextTools::toUpper(tmp_str_ploidy);
			if (tmp_str_ploidy == DIPLOID)
				locinf_ploidy = LocusInfo::DIPLOID;
			else if (tmp_str_ploidy == HAPLOID)
				locinf_ploidy = LocusInfo::HAPLOID;
			else if (tmp_str_ploidy == HAPLODIPLOID)
				locinf_ploidy = LocusInfo::HAPLODIPLOID;
		}
		if (temp.find("NbAlleles", 0) != string::npos) {
			// not used ...
		}
	}
	if (locinf_name != "")
		locus_info.push_back(LocusInfo(locinf_name, locinf_ploidy));
}

void PoplibIO::_parseIndividual(const vector<string> & in, DataSet & data_set, const VectorSequenceContainer & vsc) {
	Individual tmp_indiv;
	unsigned int tmp_group_pos = 0;
	string temp = "";
//	for (unsigned int i = 0 ; i < in.size() ; i++)
//		cout << "_parsIndiv: " << in[i] << endl;
	for (unsigned int i = 0 ; i < in.size() ; i++) {
		// Get Individual Id
		if (in[i].find(">", 0) != string::npos) {
			tmp_indiv.setId(TextTools::removeSurroundingWhiteSpaces(string(in[i].begin() + 1, in[i].end())));
		}
		// Get the Group
		if (in[i].find("Group", 0) != string::npos) {
			temp = in[i];
			tmp_group_pos = TextTools::toInt(_getValues(temp, "=")[0]);
			try {
				data_set.addEmptyGroup(tmp_group_pos);
			}
			catch (...) {}
/*			try {
			cout << "Group id : " << tmp_group_pos << endl;
			cout << "Group pos: " << data_set.getGroupPosition(tmp_group_pos) << endl;
			}
			catch (Exception & e) {
				cout << e.what() << endl;
			}
*/
		}
		// Find the locality
		if (in[i].find("Locality", 0) != string::npos) {
			temp = in[i];
			unsigned int sep_pos = temp.find("=", 0);
			string loc_name = TextTools::removeSurroundingWhiteSpaces(string(temp.begin()+sep_pos+1, temp.end()));
/*			cout << "### Locality: " << loc_name;
			try {
				cout << " - position: " << data_set.getLocalityPosition(loc_name) << endl;
			}
			catch (Exception & e) {
				cout << endl << e.what() << endl;
			}
*/
			try {
				tmp_indiv.setLocality(data_set.getLocalityByName(loc_name));
			}
			catch (Exception & e) {
//				cerr << "l.192:" << e.what() << endl;
			}
		}
		// Set the coord
		if (in[i].find("Coord", 0) != string::npos) {
			temp = in[i];
			tmp_indiv.setCoord(TextTools::toDouble(_getValues(temp, "=")[0]), TextTools::toDouble(_getValues(temp, "=")[1]));
		}
		// And the date
		if (in[i].find("Date", 0) != string::npos) {
			int d, m, y;
			temp = in[i];
			string tmp_date = _getValues(temp, "=")[0];
			d = TextTools::toInt(string(tmp_date.begin(), tmp_date.begin() + 2));
			m = TextTools::toInt(string(tmp_date.begin() + 2, tmp_date.begin() + 4));
			y = TextTools::toInt(string(tmp_date.begin() + 4, tmp_date.end()));
//			cout << "l.210:" << Date(d, m, y).getDateStr() << endl;
			tmp_indiv.setDate(Date(d, m, y));
		}
		// Now the sequences
		if (in[i].find("SequenceData", 0) != string::npos) {
			i++;
			temp = in[i];
			vector<string> seq_pos_str = _getValues(temp, "");
			for (unsigned int j = 0 ; j < seq_pos_str.size() ; j++) {
//				cout << "_parseIndiv l.308 " <<  seq_pos_str[j];
				try {
					if (seq_pos_str[j] != getMissingDataSymbol()) {
//						cout << endl;
						tmp_indiv.addSequence(j, *vsc.getSequence(TextTools::toInt(seq_pos_str[j])-1));
//						cout << tmp_indiv.getNumberOfSequences();
					}
//					else
//						cout << " is missing symbol !" << endl;
				}
				catch (Exception & e) {
//					cout << e.what() << endl;
				}
			}
		}
	}
//	cout << "nb grps : " << data_set.getNumberOfGroups() << endl;
	if (tmp_indiv.getId() != "") {
		try {
			data_set.addIndividualToGroup(data_set.getGroupPosition(tmp_group_pos), tmp_indiv);
		}
		catch (Exception & e) {
//			cout << "Et de une : " << e.what() << endl;
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
	if (data_set.hasSequenceData()) {
		string seq_type = data_set.getAlphabetType();
		os << "SequenceType = " << seq_type << endl;
	}
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
			os << "Group = " << TextTools::toString((data_set.getGroupAtPosition(i))->getGroupId()) << endl;
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
