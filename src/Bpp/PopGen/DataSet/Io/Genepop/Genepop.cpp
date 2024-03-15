// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "Genepop.h"

using namespace bpp;
using namespace std;

Genepop::Genepop() {}

Genepop::~Genepop() {}

void Genepop::read(istream& is, DataSet& data_set)
{
  if (!is)
    throw IOException("Genepop::read: fail to open stream.");
  // Skip first line
  FileTools::getNextLine(is);
  ios::pos_type entry_point = is.tellg();
  bool eof_ok = false;
  bool loc_def_ok = false;
  bool loc_nbr_ok = false;
  size_t grp_nbr = 0;
  vector<LocusInfo> tmp_loc;
  vector<set<string> > al_ids;
  map<string, size_t> ind_id_count;
  map<string, size_t> ind_id_index;

  string temp("");
  // First read : file structure
  while (!eof_ok)
  {
    if (is.peek() == EOF && !eof_ok)
    {
      // If eof rewind to entry_point
      is.seekg(entry_point);
      eof_ok = true;
    }
    else
    {
      // Count everything
      temp = FileTools::getNextLine(is);
      string cp_temp = TextTools::removeSurroundingWhiteSpaces(temp);
      cp_temp = TextTools::toUpper(cp_temp);
      if (cp_temp == string("POP"))
      {
        loc_def_ok = true;
        grp_nbr++;
        data_set.addEmptyGroup(grp_nbr);
      }
      if (!loc_def_ok)
      {
        StringTokenizer st(temp, string(", "), true);
        while (st.hasMoreToken())
          tmp_loc.push_back(LocusInfo(TextTools::removeSurroundingWhiteSpaces(st.nextToken())));
      }
      if (loc_def_ok && !loc_nbr_ok)
      {
        al_ids.resize(tmp_loc.size());
        loc_nbr_ok = true;
      }
      if (loc_def_ok)
      {
        string alleles;
        StringTokenizer st(temp, string(","));
        if (st.numberOfRemainingTokens() == 2)
        {
          ind_id_count[TextTools::removeSurroundingWhiteSpaces(st.nextToken())]++;
          alleles = st.nextToken();
        }
        StringTokenizer st2(alleles);
        if ((size_t)st2.numberOfRemainingTokens() == tmp_loc.size())
        {
          size_t i = 0;
          while (st2.hasMoreToken())
          {
            string ids = TextTools::removeSurroundingWhiteSpaces(st2.nextToken());
            string tmp_id = string(ids.begin(), ids.begin() + (ids.size() / 2));
            if (tmp_id != string("00") && tmp_id != string("000"))
              al_ids[i].insert(tmp_id);
            tmp_id = string(ids.begin() + (ids.size() / 2), ids.end());
            if (tmp_id != string("00") && tmp_id != string("000"))
              al_ids[i].insert(tmp_id);
            i++;
          }
        }
      }
    }
  }

  // Set AnalyzedLoci
  data_set.initAnalyzedLoci(tmp_loc.size());
  for (size_t i = 0; i < tmp_loc.size(); i++)
  {
    data_set.setLocusInfo(i, tmp_loc[i]);
    for (set<string>::iterator it = al_ids[i].begin(); it != al_ids[i].end(); it++)
    {
      data_set.addAlleleInfoByLocusPosition(i, BasicAlleleInfo(*it));
    }
  }

  // Second read : file data
  grp_nbr = 0;
  size_t grp_pos = 0;
  loc_def_ok = false;
  while (!is.eof())
  {
    temp = FileTools::getNextLine(is);
    string cp_temp = TextTools::removeSurroundingWhiteSpaces(temp);
    cp_temp = TextTools::toUpper(cp_temp);
    if (cp_temp == string("POP"))
    {
      grp_nbr++;
      loc_def_ok = true;
      grp_pos = data_set.getGroupPosition(grp_nbr);
    }
    else
    {
      if (loc_def_ok)
      {
        string alleles;
        StringTokenizer st(temp, string(","));
        size_t ind_pos = 0;
        if (st.numberOfRemainingTokens() == 2)
        {
          string ind_id = TextTools::removeSurroundingWhiteSpaces(st.nextToken());
          if (ind_id_count[ind_id] > 1)
            ind_id = ind_id + string("_") + TextTools::toString(++ind_id_index[ind_id]);
          data_set.addEmptyIndividualToGroup(grp_pos, ind_id);
          ind_pos = data_set.getIndividualPositionInGroup(grp_pos, ind_id);
          data_set.initIndividualGenotypeInGroup(grp_pos, ind_pos);
          alleles = st.nextToken();
        }
        StringTokenizer st2(alleles);
        if ((size_t)st2.numberOfRemainingTokens() == tmp_loc.size())
        {
          size_t i = 0;
          while (st2.hasMoreToken())
          {
            string ids = TextTools::removeSurroundingWhiteSpaces(st2.nextToken());
            vector<string> tmp_ids;
            tmp_ids.push_back(string(ids.begin(), ids.begin() + (ids.size() / 2)));
            tmp_ids.push_back(string(ids.begin() + (ids.size() / 2), ids.end()));
            if (tmp_ids[0] != string("00") && tmp_ids[0] != string("000")
                && tmp_ids[1] != string("00") && tmp_ids[1] != string("000"))
            {
              data_set.setIndividualMonolocusGenotypeByAlleleIdInGroup(grp_pos, ind_pos, i, tmp_ids);
            }
            i++;
            tmp_ids.clear();
          }
        }
      }
    }
  }
}

void Genepop::read(const string& path, DataSet& data_set)
{
  AbstractIDataSet::read(path, data_set);
}

DataSet* Genepop::read(istream& is)
{
  return AbstractIDataSet::read(is);
}

DataSet* Genepop::read(const string& path)
{
  return AbstractIDataSet::read(path);
}
