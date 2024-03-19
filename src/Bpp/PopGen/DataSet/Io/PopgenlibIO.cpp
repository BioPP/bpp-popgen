// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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

PopgenlibIO::PopgenlibIO() : data_separator_(' '),
  missing_data_symbol_('$') {}

PopgenlibIO::PopgenlibIO(const std::string& missing_data_symbol,
    const std::string& data_separator) :
  data_separator_(' '),
  missing_data_symbol_('$')
{
  try
  {
    setDataSeparator(data_separator);
    setMissingDataSymbol(missing_data_symbol);
  }
  catch (Exception& e)
  {
    throw e;
  }
}

PopgenlibIO::~PopgenlibIO() {}

void PopgenlibIO::setMissingDataSymbol(const std::string& missing_data_symbol)
{
  if (missing_data_symbol.size() != 1 || isdigit(missing_data_symbol[0])
      || TextTools::isWhiteSpaceCharacter(missing_data_symbol[0])
      || missing_data_symbol[0] == data_separator_
      )
    throw Exception("PopgenlibIO::setMissingData: not expected value for missing_data_symbol.");

  missing_data_symbol_ = missing_data_symbol[0];
}

void PopgenlibIO::setDataSeparator(const std::string& data_separator)
{
  if (data_separator == WHITESPACE)
    data_separator_ = ' ';
  else if (data_separator == TAB)
    data_separator_ = '\t';
  else if (data_separator == COMA)
    data_separator_ = ',';
  else if (data_separator == SEMICOLON)
    data_separator_ = ';';
  else
  {
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
  switch (data_separator_)
  {
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

void PopgenlibIO::read(std::istream& is, DataSet& dataset)
{
  if (!is)
    throw IOException("PopgenlibIO::read: fail to open stream.");
  string temp = "";
  vector<string> temp_v;
  stringstream tmp_ss;
  VectorSequenceContainer* tmp_vsc = NULL;
  Locality<double> tmp_locality("tmp");
  vector<LocusInfo> tmp_locinf;
  Individual tmpIndiv;
  bool section1 = true;
  bool section2 = true;
  bool section3 = true;
  bool section4 = true;
  bool section5 = true;
  size_t current_section = 0;
  size_t previous_section = 0;
  // size_t linenum = 0;
  // Main loop for all file lines
  while (!is.eof())
  {
    temp = FileTools::getNextLine(is);
    // linenum++;
    // Get the correct current section
    if (temp.find("[General]", 0) != string::npos)
    {
      previous_section = current_section;
      current_section = 1;
      continue;
    }
    else if (temp.find("[Localities]", 0) != string::npos)
    {
      previous_section = current_section;
      current_section = 2;
      continue;
    }
    else if (temp.find("[Sequences]", 0) != string::npos)
    {
      previous_section = current_section;
      current_section = 3;
      continue;
    }
    else if (temp.find("[Loci]", 0) != string::npos)
    {
      previous_section = current_section;
      current_section = 4;
      continue;
    }
    else if (temp.find("[Individuals]", 0) != string::npos)
    {
      previous_section = current_section;
      current_section = 5;
      continue;
    }
    // General section ------------------------------------
    if (current_section == 1 && previous_section < 1)
    {
      temp_v.push_back(temp);
    }
    if (section1 && current_section != 1 && previous_section == 1)
    {
      section1 = false;
      parseGeneral_(temp_v, dataset);
      temp_v.clear();
      if (dataset.hasSequenceData() && tmp_vsc == NULL)
        tmp_vsc = new VectorSequenceContainer(dataset.getAlphabet());
    }

    // Localities section ---------------------------------
    if (current_section == 2 && previous_section < 2)
    {
      if (temp.find(">", 0) != string::npos)
      {
        parseLocality_(temp_v, dataset);
        temp_v.clear();
        temp_v.push_back(temp);
      }
      else
        temp_v.push_back(temp);
    }
    if (section2 && current_section != 2 && previous_section == 2)
    {
      section2 = false;
      parseLocality_(temp_v, dataset);
      temp_v.clear();
    }

    // Sequences section ----------------------------------
    if (current_section == 3 && previous_section < 3)
    {
      if (temp.find(">", 0) != string::npos)
      {
        parseSequence_(temp_v, *tmp_vsc);
        temp_v.clear();
        temp_v.push_back(temp);
      }
      else
        temp_v.push_back(temp);
    }
    if (section3 && current_section != 3 && previous_section == 3)
    {
      section3 = false;
      parseSequence_(temp_v, *tmp_vsc);
      temp_v.clear();
    }

    // Loci section ---------------------------------------
    if (current_section == 4 && previous_section < 4)
    {
      if (temp.find(">", 0) != string::npos)
      {
        parseLoci_(temp_v, tmp_locinf);
        temp_v.clear();
        temp_v.push_back(temp);
      }
      else
        temp_v.push_back(temp);
    }
    if (section4 && current_section != 4 && previous_section == 4)
    {
      section4 = false;
      parseLoci_(temp_v, tmp_locinf);
      temp_v.clear();
      AnalyzedLoci tmp_anloc(tmp_locinf.size());
      for (size_t i = 0; i < tmp_locinf.size(); i++)
      {
        tmp_anloc.setLocusInfo(i, tmp_locinf[i]);
      }
      dataset.setAnalyzedLoci(tmp_anloc);
    }

    // Individuals section --------------------------------
    if (current_section == 5 && previous_section < 5)
    {
      if (temp.find(">", 0) != string::npos)
      {
        parseIndividual_(temp_v, dataset, *tmp_vsc);
        temp_v.clear();
        temp_v.push_back(temp);
      }
      else
        temp_v.push_back(temp);
    }
    if (section5 && current_section != 5 && previous_section == 5)
    {
      section5 = false;
      parseIndividual_(temp_v, dataset, *tmp_vsc);
      temp_v.clear();
    }
  }
  // Emptied the buffer if eof.
  if (section2 && current_section == 2)
    parseLocality_(temp_v, dataset);
  if (section3 && current_section == 3)
    parseSequence_(temp_v, *tmp_vsc);
  if (section5 && current_section == 5)
    parseIndividual_(temp_v, dataset, *tmp_vsc);
  temp_v.clear();
}

void PopgenlibIO::parseGeneral_(const std::vector<std::string>& in, DataSet& dataset)
{
  stringstream is;
  for (size_t i = 0; i < in.size(); i++)
  {
    is << in[i] << endl;
  }
  string temp;
  while (!is.eof() && in.size() != 0)
  {
    temp = FileTools::getNextLine(is);
    if (temp.find("MissingData", 0) != string::npos)
      setMissingDataSymbol(getValues_(temp, "=")[0]);
    if (temp.find("DataSeparator", 0) != string::npos)
      setDataSeparator(getValues_(temp, "=")[0]);
    if (temp.find("SequenceType", 0) != string::npos)
      dataset.setAlphabet(getValues_(temp, "=")[0]);
  }
}

void PopgenlibIO::parseLocality_(const std::vector<std::string>& in, DataSet& dataset)
{
  stringstream is;
  for (size_t i = 0; i < in.size(); i++)
  {
    is << in[i] << endl;
  }
  Locality<double> tmp_locality("");
  string temp;
  while (!is.eof() && in.size() != 0)
  {
    temp = FileTools::getNextLine(is);
    //		cout << "_parseLocality: " << temp << endl;
    if (temp.find(">", 0) != string::npos)
    {
      tmp_locality.setName(TextTools::removeSurroundingWhiteSpaces(string(temp.begin() + 1, temp.end())));
    }
    if (temp.find("Coord", 0) != string::npos)
    {
      vector<string> v = getValues_(temp, "=");
      tmp_locality.setX(TextTools::toDouble(v[0]));
      tmp_locality.setY(TextTools::toDouble(v[1]));
    }
  }
  if (tmp_locality.getName() != "")
    dataset.addLocality(tmp_locality);
}

void PopgenlibIO::parseSequence_(const std::vector<std::string>& in, VectorSequenceContainer& vsc)
{
  Fasta ifasta;
  stringstream is;
  for (size_t i = 0; i < in.size(); i++)
  {
    is << in[i] << endl;
  }
  ifasta.readSequences(is, vsc);
}

void PopgenlibIO::parseLoci_(const std::vector<std::string>& in, std::vector<LocusInfo>& locus_info)
{
  stringstream is;
  for (size_t i = 0; i < in.size(); i++)
  {
    is << in[i] << endl;
  }
  string locinf_name = "";
  unsigned int locinf_ploidy = LocusInfo::DIPLOID;
  string temp;
  while (!is.eof())
  {
    temp = FileTools::getNextLine(is);
    if (temp.find(">", 0) != string::npos)
    {
      locinf_name = TextTools::removeSurroundingWhiteSpaces(string(temp.begin() + 1, temp.end()));
    }
    if (temp.find("Ploidy", 0) != string::npos)
    {
      vector<string> v = getValues_(temp, "=");
      string tmp_str_ploidy = TextTools::removeSurroundingWhiteSpaces(v[0]);
      tmp_str_ploidy = TextTools::toUpper(tmp_str_ploidy);
      // cout << "ploidy : " << tmp_str_ploidy << endl;
      if (tmp_str_ploidy == DIPLOID)
        locinf_ploidy = LocusInfo::DIPLOID;
      else if (tmp_str_ploidy == HAPLOID)
        locinf_ploidy = LocusInfo::HAPLOID;
      else if (tmp_str_ploidy == HAPLODIPLOID)
        locinf_ploidy = LocusInfo::HAPLODIPLOID;
      else if (tmp_str_ploidy == UNKNOWN)
        locinf_ploidy = LocusInfo::UNKNOWN;
    }
    if (temp.find("NbAlleles", 0) != string::npos)
    {
      // not used ...
    }
  }
  if (locinf_name != "")
    locus_info.push_back(LocusInfo(locinf_name, locinf_ploidy));
}

void PopgenlibIO::parseIndividual_(const std::vector<std::string>& in, DataSet& dataset, const VectorSequenceContainer& vsc)
{
  Individual tmpIndiv;
  size_t tmp_group_pos = 0;
  string temp = "";
  for (size_t i = 0; i < in.size(); i++)
  {
    // Get Individual Id
    if (in[i].find(">", 0) != string::npos)
    {
      tmpIndiv.setId(TextTools::removeSurroundingWhiteSpaces(string(in[i].begin() + 1, in[i].end())));
    }
    // Get the Group
    if (in[i].find("Group", 0) != string::npos)
    {
      temp = in[i];
      tmp_group_pos = TextTools::to<size_t>(getValues_(temp, "=")[0]);
      try
      {
        dataset.addEmptyGroup(tmp_group_pos);
      }
      catch (...)
      {}
    }
    // Find the locality
    if (in[i].find("Locality", 0) != string::npos)
    {
      temp = in[i];
      size_t sep_pos = temp.find("=", 0);
      string loc_name = TextTools::removeSurroundingWhiteSpaces(string(temp.begin() + static_cast<ptrdiff_t>(sep_pos + 1), temp.end()));
      try
      {
        tmpIndiv.setLocality(dataset.getLocalityByName(loc_name));
      }
      catch (...)
      {}
    }
    // Set the coord
    if (in[i].find("Coord", 0) != string::npos)
    {
      temp = in[i];
      tmpIndiv.setCoord(TextTools::toDouble(getValues_(temp, "=")[0]), TextTools::toDouble(getValues_(temp, "=")[1]));
    }
    // And the date
    if (in[i].find("Date", 0) != string::npos)
    {
      int d, m, y;
      temp = in[i];
      string tmp_date = getValues_(temp, "=")[0];
      d = TextTools::toInt(string(tmp_date.begin(), tmp_date.begin() + 2));
      m = TextTools::toInt(string(tmp_date.begin() + 2, tmp_date.begin() + 4));
      y = TextTools::toInt(string(tmp_date.begin() + 4, tmp_date.end()));
      tmpIndiv.setDate(Date(d, m, y));
    }
    // Now the sequences
    if (in[i].find("SequenceData", 0) != string::npos)
    {
      i++;
      temp = in[i];
      vector<string> seq_pos_str = getValues_(temp, "");
      for (size_t j = 0; j < seq_pos_str.size(); ++j)
      {
        try
        {
          if (seq_pos_str[j] != getMissingDataSymbol())
          {
            auto tmpSeq = unique_ptr<Sequence>(vsc.sequence(TextTools::to<size_t>(seq_pos_str[j]) - 1).clone());
            tmpIndiv.addSequence(j, tmpSeq);
          }
        }
        catch (...)
        {}
      }
    }
    // Finally the loci
    if (in[i].find("AllelicData", 0) != string::npos)
    {
      string temp1 = in[++i];
      string temp2 = in[++i];
      vector<string> allele_pos_str1 = getValues_(temp1, "");
      vector<string> allele_pos_str2 = getValues_(temp2, "");
      try
      {
        tmpIndiv.initGenotype(dataset.getNumberOfLoci());
      }
      catch (...)
      {}
      if (allele_pos_str1.size() == allele_pos_str2.size())
      {
        for (size_t j = 0; j < allele_pos_str1.size(); j++)
        {
          const LocusInfo& locus_info = dataset.getLocusInfoAtPosition(j);
          allele_pos_str1[j] = TextTools::removeSurroundingWhiteSpaces(allele_pos_str1[j]);
          vector<string> tmp_alleles_id;
          if (allele_pos_str1[j] != getMissingDataSymbol())
          {
            BasicAlleleInfo tmp_allele_info(allele_pos_str1[j]);
            try
            {
              dataset.addAlleleInfoByLocusPosition(j, tmp_allele_info);
            }
            catch (...)
            {}
            tmp_alleles_id.push_back(allele_pos_str1[j]);
          }
          allele_pos_str2[j] = TextTools::removeSurroundingWhiteSpaces(allele_pos_str2[j]);
          if (allele_pos_str2[j] != getMissingDataSymbol())
          {
            BasicAlleleInfo tmp_allele_info(allele_pos_str2[j]);
            try
            {
              dataset.addAlleleInfoByLocusPosition(j, tmp_allele_info);
            }
            catch (...)
            {}
            tmp_alleles_id.push_back(allele_pos_str2[j]);
          }
          try
          {
            tmpIndiv.setMonolocusGenotypeByAlleleId(j, tmp_alleles_id, locus_info);
          }
          catch (...)
          {}
        }
      }
    }
  }
  if (tmpIndiv.getId() != "")
  {
    try
    {
      dataset.addIndividualToGroup(dataset.getGroupPosition(tmp_group_pos), tmpIndiv);
    }
    catch (...)
    {}
  }
}

void PopgenlibIO::read(const std::string& path, DataSet& dataset)
{
  AbstractIDataSet::read(path, dataset);
}

DataSet* PopgenlibIO::read(std::istream& is)
{
  return AbstractIDataSet::read(is);
}

DataSet* PopgenlibIO::read(const std::string& path)
{
  return AbstractIDataSet::read(path);
}

void PopgenlibIO::write(std::ostream& os, const DataSet& dataset) const
{
  size_t seqcpt = 1;
  // General section --------------------------------------
  os << "[General]" << endl;
  os << "MissingData = " << getMissingDataSymbol() << endl;
  os << "DataSeparator = " << getDataSeparator() << endl;
  if (dataset.hasSequenceData())
  {
    string seq_type = dataset.getAlphabetType();
    os << "SequenceType = " << seq_type << endl;
  }
  // Localities section -----------------------------------
  if (dataset.hasLocality())
  {
    os << endl << "[Localities]" << endl;
    for (size_t i = 0; i < dataset.getNumberOfLocalities(); i++)
    {
      os << ">" << (dataset.localityAtPosition(i)).getName() << endl;
      os << "Coord = " << (dataset.localityAtPosition(i)).getX();
      os << " " << (dataset.localityAtPosition(i)).getY() << endl;
    }
  }

  // Sequences section ------------------------------------
  if (dataset.hasSequenceData())
  {
    Fasta fasta(80);
    os << endl << "[Sequences]" << endl;
    for (size_t i = 0; i < dataset.getNumberOfGroups(); i++)
    {
      for (size_t j = 0; j < dataset.getNumberOfIndividualsInGroup(i); j++)
      {
        fasta.writeSequences(os, dataset.getIndividualAtPositionFromGroup(i, j).sequences());
      }
    }
  }

  // AllelicData section ----------------------------------
  if (dataset.hasAlleleicData())
  {
    os << endl << "[Loci]" << endl;
    for (size_t i = 0; i < dataset.getNumberOfLoci(); i++)
    {
      const LocusInfo& tmp_locus_info = dataset.getLocusInfoAtPosition(i);
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
  for (size_t i = 0; i < dataset.getNumberOfGroups(); i++)
  {
    for (size_t j = 0; j < dataset.getNumberOfIndividualsInGroup(i); j++)
    {
      if (i > 0 || j > 0)
        os << endl;
      const auto& tmpInd = dataset.getIndividualAtPositionFromGroup(i, j);
      os << ">" << tmpInd.getId() << endl;
      os << "Group = " << TextTools::toString((dataset.getGroupAtPosition(i)).getGroupId()) << endl;
      if (tmpInd.hasLocality())
        os << "Locality = " << tmpInd.locality().getName() << endl;
      if (tmpInd.hasCoord())
        os << "Coord = " << tmpInd.getX() << " " << tmpInd.getY() << endl;
      if (tmpInd.hasDate())
        os << "Date = " << tmpInd.date().getDateStr() << endl;
      if (tmpInd.hasSequences())
      {
        size_t nbss = tmpInd.getNumberOfSequences();
        os << "SequenceData = {" << endl;
        for (size_t k = 0; k < nbss; k++)
        {
          try
          {
            tmpInd.sequenceAtPosition(k);
            os << TextTools::toString(seqcpt++);
          }
          catch (SequenceNotFoundException&)
          {
            os << getMissingDataChar();
          }
          if (k < nbss - 1)
            os << getDataSeparatorChar();
          else
            os << endl;
        }
        os << "}" << endl;
      }
      if (tmpInd.hasGenotype())
      {
        const MultilocusGenotype& tmp_genotype = tmpInd.getGenotype();
        vector<vector<string>> output(tmp_genotype.size());
        os << "AllelicData = {" << endl;
        for (size_t k = 0; k < tmp_genotype.size(); k++)
        {
          output[k].resize(2);
          if (tmp_genotype.isMonolocusGenotypeMissing(k))
          {
            output[k][0] = getMissingDataChar();
            output[k][1] = getMissingDataChar();
          }
          else
          {
            vector<size_t> tmp_all_ind = tmp_genotype.monolocusGenotype(k).getAlleleIndex();
            output[k][0] = dataset.getLocusInfoAtPosition(k).getAlleleInfoByKey(tmp_all_ind[0]).getId();
            if (tmp_all_ind.size() > 1)
              output[k][1] = dataset.getLocusInfoAtPosition(k).getAlleleInfoByKey(tmp_all_ind[1]).getId();
            else
              output[k][1] = getMissingDataChar();
          }
        }
        for (size_t k = 0; k < output.size(); k++)
        {
          os << output[k][0];
          if (k < output.size() - 1)
            os << getDataSeparatorChar();
          else
            os << endl;
        }
        for (size_t k = 0; k < output.size(); k++)
        {
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
}

void PopgenlibIO::write(const std::string& path, const DataSet& dataset, bool overwrite) const
{
  AbstractODataSet::write(path, dataset, overwrite);
}

std::vector<std::string> PopgenlibIO::getValues_(std::string& param_line, const std::string& delim)
{
  vector<string> values;
  size_t limit = param_line.find(delim, 0);
  if (limit != string::npos)
    param_line = string(param_line.begin() + static_cast<ptrdiff_t>(limit + delim.size()), param_line.end());
  param_line = TextTools::removeSurroundingWhiteSpaces(param_line);

  size_t bi = 0;
  size_t bs = param_line.find(getDataSeparatorChar(), bi);
  while (bs > 0)
  {
    values.push_back(string(param_line.begin() + static_cast<ptrdiff_t>(bi), param_line.begin() + static_cast<ptrdiff_t>(bs)));
    bi = bs + 1;
    bs = param_line.find(getDataSeparatorChar(), bi);
  }
  values.push_back(string(param_line.begin() + static_cast<ptrdiff_t>(bi), param_line.end()));
  return values;
}
