//
// File PolymorphismMultiGContainerTools.cpp
// Author : Sylvain Gailard
//          Khalid Belkhir
// Last modification : june 15 2006
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

#include "PolymorphismMultiGContainerTools.h"
#include <algorithm>

using namespace std;
using namespace bpp;

/******************************************************************************/

PolymorphismMultiGContainer PolymorphismMultiGContainerTools::permutMultiG(const PolymorphismMultiGContainer & pmgc)
{
  PolymorphismMultiGContainer permuted_pmgc(pmgc);
  vector<unsigned int> groups;
  for (unsigned int i = 0 ; i < permuted_pmgc.size() ; i++)
    groups.push_back(permuted_pmgc.getGroupId(i));
  // use std::random_shuffle instead of RandomTools::getSampl
  //groups = RandomTools::getSample(groups, groups.size());
  std::random_shuffle(groups.begin(), groups.end());
  for (unsigned int i = 0 ; i < permuted_pmgc.size() ; i++)
    permuted_pmgc.setGroupId(i, groups[i]);
  return permuted_pmgc;
}

/******************************************************************************/

PolymorphismMultiGContainer PolymorphismMultiGContainerTools::permutMonoG(const PolymorphismMultiGContainer & pmgc, const std::set<unsigned int> & groups)
{
  PolymorphismMultiGContainer permuted_pmgc;
  unsigned int loc_num = pmgc.getNumberOfLoci();
  vector<vector<const MonolocusGenotype *> > mono_gens;
  mono_gens.resize(loc_num);
  // Get all the MonolocusGenotypes to permut
  for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
    if (groups.find(pmgc.getGroupId(i)) != groups.end()) {
      for (unsigned int j = 0 ; j < loc_num ; j++)
        mono_gens[j].push_back(pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j));
    }
  }
  // Permut the MonolocusGenotypes
  for (unsigned int i = 0 ; i < loc_num ; i++)
    //mono_gens[i] = RandomTools::getSample(mono_gens[i], mono_gens[i].size());
    std::random_shuffle(mono_gens[i].begin(), mono_gens[i].end());
  // Build the new PolymorphismMultiGContainer
  unsigned int k = 0;
  for (unsigned int i = 0 ; i < pmgc.size() ; i++) {
    if (groups.find(pmgc.getGroupId(i)) != groups.end()) {
      MultilocusGenotype tmp_mg(loc_num);
      for (unsigned int j = 0 ; j < loc_num ; j++) {
        if (mono_gens[j][k] != NULL)
          tmp_mg.setMonolocusGenotype(j, * (mono_gens[j][k]));
      }
      permuted_pmgc.addMultilocusGenotype(tmp_mg, pmgc.getGroupId(i));
      k++;
    }
    else {
      permuted_pmgc.addMultilocusGenotype(* (pmgc.getMultilocusGenotype(i)), pmgc.getGroupId(i));
    }
  }

  //update groups names
  set<unsigned int> grp_ids = pmgc.getAllGroupsIds();
  for (set<unsigned int>::iterator it = grp_ids.begin(); it != grp_ids.end(); it++)
  {
    unsigned int id = *it;
    string name = pmgc.getGroupName(id);
    permuted_pmgc.setGroupName(id, name);
  }

  return permuted_pmgc;
}

/******************************************************************************/

PolymorphismMultiGContainer PolymorphismMultiGContainerTools::permutIntraGroupMonoG(const PolymorphismMultiGContainer & pmgc, const std::set<unsigned int> & groups)
{
  PolymorphismMultiGContainer permuted_pmgc;
  unsigned int loc_num = pmgc.getNumberOfLoci();
  vector<vector<const MonolocusGenotype *> > mono_gens;
  mono_gens.resize(loc_num);

  for(set<unsigned int>::const_iterator g = groups.begin(); g != groups.end(); g++)//for each group
  {
    unsigned int nb_ind_in_group = 0;
    // Get all the MonolocusGenotypes of group g to permut
    for (unsigned int i = 0 ; i < pmgc.size() ; i++)
    {
      int indiv_grp = pmgc.getGroupId(i);
      if (groups.find(indiv_grp) != groups.end())
      {
        if(indiv_grp == (int)(*g))
        {
          nb_ind_in_group ++;

          for(unsigned int j = 0; j < loc_num ; j++)
            mono_gens[j].push_back(pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j));
        }
      }
      else //insert as is
      {
        permuted_pmgc.addMultilocusGenotype(* (pmgc.getMultilocusGenotype(i)),indiv_grp);
      }
    }//for i

    // Permut the MonolocusGenotypes
    if(nb_ind_in_group > 0)
    {
      for(unsigned int j = 0 ; j < loc_num ; j++)
        //mono_gens[j] = RandomTools::getSample(mono_gens[j], mono_gens[j].size());
        std::random_shuffle(mono_gens[j].begin(), mono_gens[j].end());

      // Build the new multilocus genotypes
      MultilocusGenotype tmp_mg(loc_num);
      for(unsigned int k = 0; k < nb_ind_in_group; k++)
      {
        for (unsigned int j = 0 ; j < loc_num ; j++)
        {
          if (mono_gens[j][k] != NULL) tmp_mg.setMonolocusGenotype(j, * (mono_gens[j][k]));
        }//for j

        permuted_pmgc.addMultilocusGenotype(tmp_mg, (*g));
      }//for k
    }//if nb_ind_in_group

  }//for g

  //update groups names
  set<unsigned int> grp_ids = pmgc.getAllGroupsIds();
  for (set<unsigned int>::iterator it = grp_ids.begin(); it != grp_ids.end(); it++)
  {
    unsigned int id = *it;
    string name = pmgc.getGroupName(id);
    permuted_pmgc.setGroupName(id, name);
  }

  return permuted_pmgc;
}

/******************************************************************************/

PolymorphismMultiGContainer PolymorphismMultiGContainerTools::permutAlleles(const PolymorphismMultiGContainer & pmgc, const std::set<unsigned int> & groups)
{
  PolymorphismMultiGContainer permuted_pmgc;
  unsigned int loc_num = pmgc.getNumberOfLoci();
  vector<vector<unsigned int> > alleles;
  alleles.resize(loc_num);
  // Get all the alleles to permut
  for (unsigned int i = 0 ; i < pmgc.size() ; i++)
  {
    if (groups.find(pmgc.getGroupId(i)) != groups.end())
    {
      for (unsigned int j = 0 ; j < loc_num ; j++)
        if (pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j) != NULL)
          for (unsigned int k = 0 ; k < pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j)->getAlleleIndex().size() ; k++)
            alleles[j].push_back(pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j)->getAlleleIndex()[k]);
    }
  }
  // Permut the alleles
  for (unsigned int i = 0 ; i < loc_num ; i++)
    //alleles[i] = RandomTools::getSample(alleles[i], alleles[i].size());
    std::random_shuffle(alleles[i].begin(), alleles[i].end());
  // Build the new PolymorphismMultiGContainer
  vector<unsigned int> k(loc_num,0);
  for (unsigned int i = 0 ; i < pmgc.size() ; i++)
  {
    if (groups.find(pmgc.getGroupId(i)) != groups.end())
    {
      MultilocusGenotype tmp_mg(loc_num);
      for (unsigned int j = 0 ; j < loc_num ; j++)
      {
        if (pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j) != NULL)
        {
          if (pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j)->getAlleleIndex().size() == 1)
            tmp_mg.setMonolocusGenotype(j, MonoAlleleMonolocusGenotype(alleles[j][k[j]++]));
          if (pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j)->getAlleleIndex().size() == 2)
            tmp_mg.setMonolocusGenotype(j, BiAlleleMonolocusGenotype(alleles[j][k[j]++], alleles[j][k[j]++]));
        }
      }
      permuted_pmgc.addMultilocusGenotype(tmp_mg, pmgc.getGroupId(i));
    }
    else
    {
      permuted_pmgc.addMultilocusGenotype(* (pmgc.getMultilocusGenotype(i)), pmgc.getGroupId(i));
    }
  }

  //update groups names
  set<unsigned int> grp_ids = pmgc.getAllGroupsIds();
  for (set<unsigned int>::iterator it = grp_ids.begin(); it != grp_ids.end(); it++)
  {
    unsigned int id = *it;
    string name = pmgc.getGroupName(id);
    permuted_pmgc.setGroupName(id, name);
  }

  return permuted_pmgc;
}

/******************************************************************************/

PolymorphismMultiGContainer PolymorphismMultiGContainerTools::permutIntraGroupAlleles(const PolymorphismMultiGContainer & pmgc, const std::set<unsigned int> & groups)
{
  PolymorphismMultiGContainer permuted_pmgc;
  unsigned int loc_num = pmgc.getNumberOfLoci();
  vector<vector<unsigned int> > alleles;
  alleles.resize(loc_num);

  for (set<unsigned int>::const_iterator g = groups.begin(); g != groups.end(); g++)//for each group
  {
    int nb_ind_in_group = 0;

    vector< vector<unsigned int> > nb_alleles_for_inds;
    nb_alleles_for_inds.resize(loc_num);
    // Get all the alleles to permut
    for (unsigned int i = 0; i < pmgc.size(); i++)
    {
      int indiv_grp = pmgc.getGroupId(i);
      if (groups.find(indiv_grp) != groups.end() )
      {
        if (indiv_grp == (int)(*g) )
        {
          nb_ind_in_group++;
          for (unsigned int j = 0 ; j < loc_num ; j++)
          {
            if (pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j) != NULL) //? données manquantes
            {
              unsigned int nb_alls = pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j)->getAlleleIndex().size();
              nb_alleles_for_inds[j].push_back(nb_alls);
              for (unsigned int k = 0 ; k < nb_alls ; k++)
                alleles[j].push_back(pmgc.getMultilocusGenotype(i)->getMonolocusGenotype(j)->getAlleleIndex()[k]);
            }//if
          }//for j
        }
      }
      else //inserer tel quel
      {
        permuted_pmgc.addMultilocusGenotype(* (pmgc.getMultilocusGenotype(i)), indiv_grp);
      }
    }//for i

    // Permut the alleles
    if (nb_ind_in_group > 0)
    {
      for (unsigned int i = 0 ; i < loc_num ; i++)
        //alleles[i] = RandomTools::getSample(alleles[i], alleles[i].size());
        std::random_shuffle(alleles[i].begin(), alleles[i].end());

      // Build the new PolymorphismMultiGContainer
      vector<unsigned int> k(loc_num,0);

      for (int ind=0; ind < nb_ind_in_group; ind++)
      {
        MultilocusGenotype tmp_mg(loc_num);
        for (unsigned int j = 0 ; j < loc_num ; j++)
        {

          if (nb_alleles_for_inds[j][ind] == 1)
            tmp_mg.setMonolocusGenotype(j, MonoAlleleMonolocusGenotype(alleles[j][k[j]++]));
          if (nb_alleles_for_inds[j][ind] == 2)
            tmp_mg.setMonolocusGenotype(j, BiAlleleMonolocusGenotype(alleles[j][k[j]++], alleles[j][k[j]++]));

        }//for j

        permuted_pmgc.addMultilocusGenotype(tmp_mg, (*g));

      }//for ind

    }//if nb_ind_in_group

  }//for g


  //update groups names
  set<unsigned int> grp_ids = pmgc.getAllGroupsIds();
  for (set<unsigned int>::iterator it = grp_ids.begin(); it != grp_ids.end(); it++)
  {
    unsigned int id = *it;
    string name = pmgc.getGroupName(id);
    permuted_pmgc.setGroupName(id, name);
  }

  return permuted_pmgc;
}

/******************************************************************************/

PolymorphismMultiGContainer PolymorphismMultiGContainerTools::extractGroups(const PolymorphismMultiGContainer & pmgc, const std::set<unsigned int> & groups)
{
  PolymorphismMultiGContainer sub_pmgc;
  for (set<unsigned int>::const_iterator g = groups.begin(); g != groups.end() ; g++)//for each group
  {

    // Get all the MonolocusGenotypes of group g to extract
    for(unsigned int i = 0 ; i < pmgc.size() ; i++)
    {
      int indiv_grp = pmgc.getGroupId(i);
      if (groups.find(indiv_grp) != groups.end() )
      {
        if(indiv_grp == (int)(*g))
        {
          sub_pmgc.addMultilocusGenotype(* (pmgc.getMultilocusGenotype(i)),indiv_grp);
        }
      }

    }//for i  

  }//for g

  //update groups names
  set<unsigned int> grp_ids = sub_pmgc.getAllGroupsIds();
  for (set<unsigned int>::iterator it = grp_ids.begin(); it != grp_ids.end(); it++)
  {
    unsigned int id = *it;
    string name = pmgc.getGroupName(id);
    sub_pmgc.setGroupName(id, name);
  }

  return sub_pmgc;
}

/******************************************************************************/

