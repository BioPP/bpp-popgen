//
// File PopgenlibIO.cpp
// Created by: Sylvain Gaillard
// Created on: Thursday July 29 2004
//

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

#include "PopgenlibIO.h"

using namespace bpp;
using namespace std;

const string PopgenlibIO::WHITESPACE = string("WHITESPACE");
const string PopgenlibIO::TAB = string("TAB");
const string PopgenlibIO::COMA = string("COMA");
const string PopgenlibIO::SEMICOLON = string("SEMICOLON");

const string PopgenlibIO::DIPLOID = string("DIPLOID");
const string PopgenlibIO::HAPLOID = string("HAPLOID");
const string PopgenlibIO::HAPLODIPLOID = string("HAPLODIPLOID");
const string PopgenlibIO::UNKNOWN = string("UNKNOWN");

PopgenlibIO::PopgenlibIO(): data_separator_(' '), missing_data_symbol_('$') {}

PopgenlibIO::PopgenlibIO(const std::string& missing_data_symbol,
    const std::string& data_separator)
throw (Exception): data_separator_(' '), missing_data_symbol_('$')
{
  try {
    setDataSeparator(data_separator);
    setMissingDataSymbol(missing_data_symbol);
  }
  catch (Exception& e) {
    throw e;
  }
}

PopgenlibIO::~PopgenlibIO() {}

void PopgenlibIO::setMissingDataSymbol(const std::string& missing_data_symbol) throw (Exception)
{
  if (missing_data_symbol.size() != 1 || isdigit(missing_data_symbol[0])
      || TextTools::isWhiteSpaceCharacter(missing_data_symbol[0])
      || missing_data_symbol[0] == data_separator_
     )
    throw Exception("PopgenlibIO::setMissingData: not expected value for missing_data_symbol.");

  missing_data_symbol_ = missing_data_symbol[0];
}

void PopgenlibIO::setDataSeparator(const std::string& data_separator) throw (Exception)
{
  if (data_separator == WHITESPACE) data_separator_ = ' ';
  else if (data_separator == TAB) data_separator_ = '\t';
  else if (data_separator == COMA) data_separator_ = ',';
  else if (data_separator == SEMICOLON) data_separator_ = ';';
  else {
    if (isdigit(data_separator[0])
        || data_separator == getMissingDataSymbol()
       )
      throw Exception("PopgenlibIO::setDataSeparator: not expected value for data_separator.");
    data_separator_ = data_separator.c_str()[0];
  }
}

std::string PopgenlibIO::getMissingDataSymbol() const
{
  return TextTools::toString(missing_data_symbol_);
}

std::string PopgenlibIO::getDataSeparator() const
{
  switch (data_separator_) {
    case (' '): return WHITESPACE;
    case ('\t'): return TAB;
    case (','): return COMA;
    case (';'): return SEMICOLON;
    default: return TextTools::toString(data_separator_);
  }
}

char PopgenlibIO::getMissingDataChar() const
{
  return missing_data_symbol_;
}

char PopgenlibIO::getDataSeparatorChar() const
{
  return data_separator_;
}

void PopgenlibIO::read(std::istream& is, DataSet& data_set) throw (Exception)
{
  if (!is)
    throw IOException("PopgenlibIO::read: fail to open stream.");
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
      parseGeneral_(temp_v, data_set);
      temp_v.clear();
      if (data_set.hasSequenceData() && tmp_vsc == NULL)
        tmp_vsc = new VectorSequenceContainer(data_set.getAlphabet());
    }

    // Localities section ---------------------------------
    if (current_section == 2 && previous_section < 2) {
      if (temp.find(">", 0) != string::npos) {
        parseLocality_(temp_v, data_set);
        temp_v.clear();
        temp_v.push_back(temp);
      }
      else
        temp_v.push_back(temp);
    }
    if (section2 && current_section !=2 && previous_section == 2) {
      section2 = false;
      parseLocality_(temp_v, data_set);
      temp_v.clear();
    }	

    // Sequences section ----------------------------------
    if (current_section == 3 && previous_section < 3) {
      if (temp.find(">", 0) != string::npos) {
        parseSequence_(temp_v, *tmp_vsc);
        temp_v.clear();
        temp_v.push_back(temp);
      }
      else
        temp_v.push_back(temp);
    }
    if (section3 && current_section !=3 && previous_section == 3) {
      section3 = false;
      parseSequence_(temp_v, *tmp_vsc);
      temp_v.clear();
    }

    // Loci section ---------------------------------------
    if (current_section == 4 && previous_section <4) {
      if (temp.find(">", 0) != string::npos) {
        parseLoci_(temp_v, tmp_locinf);
        temp_v.clear();
        temp_v.push_back(temp);
      }
      else
        temp_v.push_back(temp);
    }
    if (section4 && current_section != 4 && previous_section == 4) {
      section4 = false;
      parseLoci_(temp_v, tmp_locinf);
      temp_v.clear();
      AnalyzedLoci tmp_anloc(tmp_locinf.size());
      for (unsigned int i = 0 ; i < tmp_locinf.size() ; i++)
        tmp_anloc.setLocusInfo(i, tmp_locinf[i]);
      data_set.setAnalyzedLoci(tmp_anloc);
    }

    // Individuals section --------------------------------
    if (current_section == 5 && previous_section < 5) {
      if (temp.find(">", 0) != string::npos) {
        parseIndividual_(temp_v, data_set, * tmp_vsc);
        temp_v.clear();
        temp_v.push_back(temp);
      }
      else
        temp_v.push_back(temp);
    }
    if (section5 && current_section != 5 && previous_section == 5) {
      section5 = false;
      parseIndividual_(temp_v, data_set, * tmp_vsc);
      temp_v.clear();
    }
  }
  // Emptied the buffer if eof.
  if (section2 && current_section == 2)
    parseLocality_(temp_v, data_set);
  if (section3 && current_section == 3)
    parseSequence_(temp_v, * tmp_vsc);
  if (section5 && current_section == 5)
    parseIndividual_(temp_v, data_set, * tmp_vsc);
  temp_v.clear();
}

void PopgenlibIO::parseGeneral_(const std::vector<std::string>& in, DataSet& data_set)
{
  stringstream is;
  for (unsigned int i = 0 ; i < in.size() ; i++)
    is << in[i] << endl;
  string temp;
  while (!is.eof() && in.size() != 0) {
    temp = FileTools::getNextLine(is);
    if (temp.find("MissingData", 0) != string::npos)
      setMissingDataSymbol(getValues_(temp, "=")[0]);
    if (temp.find("DataSeparator", 0) != string::npos)
      setDataSeparator(getValues_(temp, "=")[0]);
    if (temp.find("SequenceType", 0) != string::npos)
      data_set.setAlphabet(getValues_(temp, "=")[0]);
  }
}

void PopgenlibIO::parseLocality_(const std::vector<std::string>& in, DataSet& data_set)
{
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
      vector<string> v = getValues_(temp, "=");
      tmp_locality.setX(TextTools::toDouble(v[0]));
      tmp_locality.setY(TextTools::toDouble(v[1]));
    }
  }
  if (tmp_locality.getName() != "")
    data_set.addLocality(tmp_locality);
}

void PopgenlibIO::parseSequence_(const std::vector<std::string>& in, VectorSequenceContainer& vsc) 
{
  Fasta ifasta;
  stringstream is;
  for (unsigned int i = 0 ; i < in.size() ; i++)
    is << in[i] << endl;
  ifasta.read(is, vsc);
}

void PopgenlibIO::parseLoci_(const std::vector<std::string>& in, std::vector<LocusInfo>& locus_info)
{
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
      vector<string> v = getValues_(temp, "=");
      string tmp_str_ploidy = TextTools::removeSurroundingWhiteSpaces(v[0]);
      tmp_str_ploidy = TextTools::toUpper(tmp_str_ploidy);
      //cout << "ploidy : " << tmp_str_ploidy << endl;
      if (tmp_str_ploidy == DIPLOID)
        locinf_ploidy = LocusInfo::DIPLOID;
      else if (tmp_str_ploidy == HAPLOID)
        locinf_ploidy = LocusInfo::HAPLOID;
      else if (tmp_str_ploidy == HAPLODIPLOID)
        locinf_ploidy = LocusInfo::HAPLODIPLOID;
      else if (tmp_str_ploidy == UNKNOWN)
        locinf_ploidy = LocusInfo::UNKNOWN;
    }
    if (temp.find("NbAlleles", 0) != string::npos) {
      // not used ...
    }
  }
  if (locinf_name != "")
    locus_info.push_back(LocusInfo(locinf_name, locinf_ploidy));
}

void PopgenlibIO::parseIndividual_(const std::vector<std::string>& in, DataSet& data_set, const VectorSequenceContainer& vsc)
{
  Individual tmp_indiv;
  unsigned int tmp_group_pos = 0;
  string temp = "";
  for (unsigned int i = 0 ; i < in.size() ; i++) {
    // Get Individual Id
    if (in[i].find(">", 0) != string::npos) {
      tmp_indiv.setId(TextTools::removeSurroundingWhiteSpaces(string(in[i].begin() + 1, in[i].end())));
    }
    // Get the Group
    if (in[i].find("Group", 0) != string::npos) {
      temp = in[i];
      tmp_group_pos = TextTools::toInt(getValues_(temp, "=")[0]);
      try {
        data_set.addEmptyGroup(tmp_group_pos);
      }
      catch (...) {}
    }
    // Find the locality
    if (in[i].find("Locality", 0) != string::npos) {
      temp = in[i];
      unsigned int sep_pos = temp.find("=", 0);
      string loc_name = TextTools::removeSurroundingWhiteSpaces(string(temp.begin()+sep_pos+1, temp.end()));
      try {
        tmp_indiv.setLocality(& data_set.getLocalityByName(loc_name));
      }
      catch (...) {}
    }
    // Set the coord
    if (in[i].find("Coord", 0) != string::npos) {
      temp = in[i];
      tmp_indiv.setCoord(TextTools::toDouble(getValues_(temp, "=")[0]), TextTools::toDouble(getValues_(temp, "=")[1]));
    }
    // And the date
    if (in[i].find("Date", 0) != string::npos) {
      int d, m, y;
      temp = in[i];
      string tmp_date = getValues_(temp, "=")[0];
      d = TextTools::toInt(string(tmp_date.begin(), tmp_date.begin() + 2));
      m = TextTools::toInt(string(tmp_date.begin() + 2, tmp_date.begin() + 4));
      y = TextTools::toInt(string(tmp_date.begin() + 4, tmp_date.end()));
      tmp_indiv.setDate(Date(d, m, y));
    }
    // Now the sequences
    if (in[i].find("SequenceData", 0) != string::npos) {
      i++;
      temp = in[i];
      vector<string> seq_pos_str = getValues_(temp, "");
      for (unsigned int j = 0 ; j < seq_pos_str.size() ; j++) {
        try {
          if (seq_pos_str[j] != getMissingDataSymbol())
            tmp_indiv.addSequence(j, vsc.getSequence(TextTools::toInt(seq_pos_str[j])-1));
        }
        catch (...) {}
      }
    }
    // Finally the loci
    if (in[i].find("AllelicData", 0) != string::npos)
    {
      string temp1 = in[++i];
      string temp2 = in[++i];
      vector<string> allele_pos_str1 = getValues_(temp1, "");
      vector<string> allele_pos_str2 = getValues_(temp2, "");
      try {
        tmp_indiv.initGenotype(data_set.getNumberOfLoci());
      }
      catch (...) {}
      if (allele_pos_str1.size() == allele_pos_str2.size()) {
        for (unsigned int j = 0 ; j < allele_pos_str1.size() ; j++) {
          const LocusInfo& locus_info = data_set.getLocusInfoAtPosition(j);
          allele_pos_str1[j] = TextTools::removeSurroundingWhiteSpaces(allele_pos_str1[j]);
          vector<string> tmp_alleles_id;
          if (allele_pos_str1[j] != getMissingDataSymbol()) {
            BasicAlleleInfo tmp_allele_info(allele_pos_str1[j]);
            try {
              data_set.addAlleleInfoByLocusPosition(j, tmp_allele_info);
            }
            catch (...) {}
            tmp_alleles_id.push_back(allele_pos_str1[j]);
          }
          allele_pos_str2[j] = TextTools::removeSurroundingWhiteSpaces(allele_pos_str2[j]);
          if (allele_pos_str2[j] != getMissingDataSymbol()) {
            BasicAlleleInfo tmp_allele_info(allele_pos_str2[j]);
            try {
              data_set.addAlleleInfoByLocusPosition(j, tmp_allele_info);
            }
            catch (...) {}
            tmp_alleles_id.push_back(allele_pos_str2[j]);
          }
          try {
            tmp_indiv.setMonolocusGenotypeByAlleleId(j, tmp_alleles_id, locus_info);
          }
          catch (...) {}
        }
      }
    }
  }
  if (tmp_indiv.getId() != "")
  {
    try {
      data_set.addIndividualToGroup(data_set.getGroupPosition(tmp_group_pos), tmp_indiv);
    }
    catch (...) {}
  }
}

void PopgenlibIO::read(const std::string& path, DataSet& data_set) throw (Exception)
{
  AbstractIDataSet::read(path, data_set);
}

DataSet* PopgenlibIO::read(std::istream& is) throw (Exception)
{
  return AbstractIDataSet::read(is);
}

DataSet* PopgenlibIO::read(const std::string& path) throw (Exception)
{
  return AbstractIDataSet::read(path);
}

void PopgenlibIO::write(std::ostream& os, const DataSet& data_set) const throw (Exception)
{
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
      os << ">" << (data_set.getLocalityAtPosition(i)).getName() << endl;
      os << "Coord = " << (data_set.getLocalityAtPosition(i)).getX();
      os << " " << (data_set.getLocalityAtPosition(i)).getY() << endl;
    }
  }

  // Sequences section ------------------------------------
  if (data_set.hasSequenceData()) {
    Fasta fasta(80);
    os << endl << "[Sequences]" << endl;
    for (unsigned int i = 0 ; i < data_set.getNumberOfGroups() ; i++)
      for (unsigned int j = 0 ; j < data_set.getNumberOfIndividualsInGroup(i) ; j++)
        fasta.write(os, data_set.getIndividualAtPositionFromGroup(i,j)->getSequences());
  }

  // AllelicData section ----------------------------------
  if (data_set.hasAlleleicData()) {
    os << endl << "[Loci]" << endl;
    for (unsigned int i = 0 ; i < data_set.getNumberOfLoci() ; i++) {
      const LocusInfo& tmp_locus_info = data_set.getLocusInfoAtPosition(i);
      os << ">" << tmp_locus_info.getName() << endl;
      os << "Ploidy = ";
      if (tmp_locus_info.getPloidy() == LocusInfo::HAPLOID)
        os << HAPLOID;
      else if (tmp_locus_info.getPloidy() == LocusInfo::DIPLOID)
        os << DIPLOID;
      else if (tmp_locus_info.getPloidy() == LocusInfo::HAPLODIPLOID)
        os << HAPLODIPLOID;
      else if (tmp_locus_info.getPloidy() == LocusInfo::UNKNOWN)
        os << UNKNOWN;
      os << endl;
      os << "NbAlleles = " << tmp_locus_info.getNumberOfAlleles() << endl;
    }
  }

  // Individuals section ----------------------------------
  os << endl << "[Individuals]" << endl;
  for (unsigned int i = 0 ; i < data_set.getNumberOfGroups() ; i++)
    for (unsigned int j = 0 ; j < data_set.getNumberOfIndividualsInGroup(i) ; j++) {
      if (i>0 || j>0) os << endl;
      const Individual * tmp_ind = data_set.getIndividualAtPositionFromGroup(i,j);
      os << ">" << tmp_ind->getId() << endl;
      os << "Group = " << TextTools::toString((data_set.getGroupAtPosition(i)).getGroupId()) << endl;
      if (tmp_ind->hasLocality()) os << "Locality = " << tmp_ind->getLocality()->getName() << endl;
      if (tmp_ind->hasCoord()) os << "Coord = " << tmp_ind->getX() << " " << tmp_ind->getY() << endl;
      if (tmp_ind->hasDate()) os << "Date = " << tmp_ind->getDate().getDateStr() << endl;
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
      if (tmp_ind->hasGenotype()) {
        const MultilocusGenotype& tmp_genotype = tmp_ind->getGenotype();
        vector<vector<string> > output(tmp_genotype.size());
        os << "AllelicData = {" << endl;
        for (unsigned int k = 0 ; k < tmp_genotype.size() ; k++) {
          output[k].resize(2);
          if (tmp_genotype.isMonolocusGenotypeMissing(k)) {
            output[k][0] = getMissingDataChar();
            output[k][1] = getMissingDataChar();
          }
          else {
            vector<unsigned int> tmp_all_ind = tmp_genotype.getMonolocusGenotype(k).getAlleleIndex();
            output[k][0] = data_set.getLocusInfoAtPosition(k).getAlleleInfoByKey(tmp_all_ind[0]).getId();
            if (tmp_all_ind.size() > 1)
              output[k][1] = data_set.getLocusInfoAtPosition(k).getAlleleInfoByKey(tmp_all_ind[1]).getId();
            else
              output[k][1] = getMissingDataChar();
          }
        }
        for (unsigned int k = 0 ; k < output.size() ; k++) {
          os << output[k][0];
          if (k < output.size() - 1)
            os << getDataSeparatorChar();
          else
            os << endl;
        }
        for (unsigned int k = 0 ; k < output.size() ; k++) {
          os << output[k][1];
          if (k < output.size() - 1)
            os << getDataSeparatorChar();
          else
            os << endl;
        }
        os << "}" << endl;
      }
    }
}

void PopgenlibIO::write(const std::string& path, const DataSet& data_set, bool overwrite) const throw (Exception) 
{
  AbstractODataSet::write(path, data_set, overwrite);
}

std::vector<std::string> PopgenlibIO::getValues_(std::string& param_line, const std::string& delim)
{
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

