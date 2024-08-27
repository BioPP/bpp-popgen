/*
 * File MultilocusGenotypeStatistics.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday August 04 2004
 *
 */

// Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)
// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Utils/MapTools.h>

#include "MultilocusGenotypeStatistics.h"
#include "PolymorphismMultiGContainerTools.h"

using namespace bpp;

// From STL

#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

vector<size_t> MultilocusGenotypeStatistics::getAllelesIdsForGroups(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  map<size_t, size_t> tmp_alleles;
  try
  {
    tmp_alleles = getAllelesMapForGroups(pmgc, locusPosition, groups);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getAllelesIdsForGroups: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  return MapTools::getKeys(tmp_alleles);
}

size_t MultilocusGenotypeStatistics::countGametesForGroups(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  map<size_t, size_t> allele_count;
  size_t nb_tot_allele = 0;
  try
  {
    allele_count = getAllelesMapForGroups(pmgc, locusPosition, groups);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::countGametesForGroups: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  vector<size_t> counter = MapTools::getValues(allele_count);
  for (size_t i = 0; i < counter.size(); i++)
  {
    nb_tot_allele += counter[i];
  }
  return nb_tot_allele;
}

map<size_t, size_t> MultilocusGenotypeStatistics::getAllelesMapForGroups(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  map<size_t, size_t> alleles_count;
  for (size_t i = 0; i < pmgc.size(); i++)
  {
    try
    {
      if (!pmgc.multilocusGenotype(i).isMonolocusGenotypeMissing(locusPosition) &&  (groups.find(pmgc.getGroupId(i)) != groups.end()) )
      {
        // if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locusPosition) &&  (find(groups.begin(), groups.end(), pmgc.getGroupId(i)) != groups.end()) ) {
        vector<size_t> tmp_alleles = pmgc.multilocusGenotype(i).monolocusGenotype(locusPosition).getAlleleIndex();
        for (size_t j = 0; j < tmp_alleles.size(); j++)
        {
          alleles_count[tmp_alleles[j]]++;
        }
      }
    }
    catch (IndexOutOfBoundsException& ioobe)
    {
      throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getAllelesMapForGroups: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
    }
  }
  return alleles_count;
}

map<size_t, double> MultilocusGenotypeStatistics::getAllelesFrqForGroups(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  map<size_t, double> alleles_frq;
  size_t nb_tot_allele = 0;
  map<size_t, size_t> tmp_alleles;
  try
  {
    tmp_alleles = getAllelesMapForGroups(pmgc, locusPosition, groups);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getAllelesFrqForGroups: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  vector<size_t> counter = MapTools::getValues(tmp_alleles);
  for (size_t i = 0; i < counter.size(); i++)
  {
    nb_tot_allele += counter[i];
  }
  if (nb_tot_allele == 0)
    throw ZeroDivisionException("MultilocusGenotypeStatistics::getAllelesFrqForGroups.");
  for (map<size_t, size_t>::iterator it = tmp_alleles.begin(); it != tmp_alleles.end(); it++)
  {
    alleles_frq[it->first] = static_cast<double>(it->second) / static_cast<double>(nb_tot_allele);
  }
  return alleles_frq;
}

size_t MultilocusGenotypeStatistics::countNonMissingForGroups(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  size_t counter = 0;
  for (size_t i = 0; i < pmgc.size(); i++)
  {
    try
    {
      // if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locusPosition) &&  (find(groups.begin(), groups.end(), pmgc.getGroupId(i)) != groups.end()) )
      if (!pmgc.multilocusGenotype(i).isMonolocusGenotypeMissing(locusPosition) &&  (groups.find(pmgc.getGroupId(i) ) != groups.end()) )
        counter++;
    }
    catch (IndexOutOfBoundsException& ioobe)
    {
      throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::countNonMissing: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
    }
  }
  return counter;
}

size_t MultilocusGenotypeStatistics::countBiAllelicForGroups(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  size_t counter = 0;
  for (size_t i = 0; i < pmgc.size(); i++)
  {
    try
    {
      // if (! pmgc.getMultilocusGenotype(i)->isMonolocusGenotypeMissing(locusPosition) && (find(groups.begin(), groups.end(), pmgc.getGroupId(i)) != groups.end()) )
      if (!pmgc.multilocusGenotype(i).isMonolocusGenotypeMissing(locusPosition) && (groups.find(pmgc.getGroupId(i)) != groups.end()) )
        if ((pmgc.multilocusGenotype(i).monolocusGenotype(locusPosition).getAlleleIndex()).size() == 2)
          counter++;
    }
    catch (IndexOutOfBoundsException& ioobe)
    {
      throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::countBiAllelic: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
    }
  }
  return counter;
}

map<size_t, size_t> MultilocusGenotypeStatistics::countHeterozygousForGroups(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  map<size_t, size_t> counter;
  for (size_t i = 0; i < pmgc.size(); i++)
  {
    try
    {
      if (!pmgc.multilocusGenotype(i).isMonolocusGenotypeMissing(locusPosition) && (groups.find(pmgc.getGroupId(i)) != groups.end() ))
      {
        const auto& tmpMg = pmgc.multilocusGenotype(i).monolocusGenotype(locusPosition);
        if ((tmpMg.getAlleleIndex()).size() == 2)
        {
          if (!dynamic_cast<const BiAlleleMonolocusGenotype&>(tmpMg).isHomozygous())
          {
            vector<size_t> tmpAlleles = tmpMg.getAlleleIndex();
            for (size_t j = 0; j < tmpAlleles.size(); j++)
            {
              counter[tmpAlleles[j]]++;
            }
          }
        }
      }
    }
    catch (IndexOutOfBoundsException& ioobe)
    {
      throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::countHeterozygous: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
    }
  }
  return counter;
}

map<size_t, double> MultilocusGenotypeStatistics::getHeterozygousFrqForGroups(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  map<size_t, double> freq;
  size_t counter = 0;
  for (size_t i = 0; i < pmgc.size(); i++)
  {
    try
    {
      if (!pmgc.multilocusGenotype(i).isMonolocusGenotypeMissing(locusPosition) && (groups.find(pmgc.getGroupId(i)) != groups.end()) )
      {
        const auto& tmpMg = pmgc.multilocusGenotype(i).monolocusGenotype(locusPosition);
        if ((tmpMg.getAlleleIndex()).size() == 2)
        {
          counter++;
          if (!dynamic_cast<const BiAlleleMonolocusGenotype&>(tmpMg).isHomozygous())
          {
            vector<size_t> tmpAlleles = tmpMg.getAlleleIndex();
            for (size_t j = 0; j < tmpAlleles.size(); j++)
            {
              freq[tmpAlleles[j]]++;
            }
          }
        }
      }
    }
    catch (IndexOutOfBoundsException& ioobe)
    {
      throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getHeterozygousFrqForGroups: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
    }
  }
  if (counter == 0)
    throw ZeroDivisionException("MultilocusGenotypeStatistics::getHeterozygousFrqForGroups.");
  for (map<size_t, double>::iterator i = freq.begin(); i != freq.end(); i++)
  {
    i->second = (double) i->second / (double) counter;
  }
  return freq;
}

double MultilocusGenotypeStatistics::getHobsForGroups(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  map<size_t, double> heterozygous_frq;
  double frq = 0.;
  try
  {
    heterozygous_frq = getHeterozygousFrqForGroups(pmgc, locusPosition, groups);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getHobsForGroups: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (ZeroDivisionException& zde)
  {
    throw ZeroDivisionException("MultilocusGenotypeStatistics::getHobsForGroups.");
  }
  for (map<size_t, double>::iterator it = heterozygous_frq.begin(); it != heterozygous_frq.end(); it++)
  {
    frq += it->second;
  }
  return frq / static_cast<double>(heterozygous_frq.size());
}

double MultilocusGenotypeStatistics::getHexpForGroups(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  map<size_t, double> allele_frq;
  double frqsqr = 0.;
  try
  {
    allele_frq = getAllelesFrqForGroups(pmgc, locusPosition, groups);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getHexpForGroups: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (ZeroDivisionException& zde)
  {
    throw ZeroDivisionException("MultilocusGenotypeStatistics::getHexpForGroups.");
  }
  for (map<size_t, double>::iterator it = allele_frq.begin(); it != allele_frq.end(); it++)
  {
    frqsqr += it->second * it->second;
  }
  return 1 - frqsqr;
}

double MultilocusGenotypeStatistics::getHnbForGroups(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  size_t nb_alleles;
  double Hexp;
  try
  {
    nb_alleles = countGametesForGroups(pmgc, locusPosition, groups);
    Hexp = getHexpForGroups(pmgc, locusPosition, groups);
  }
  catch (IndexOutOfBoundsException& ioobe)
  {
    throw IndexOutOfBoundsException("MultilocusGenotypeStatistics::getHnbForGroups: locusPosition out of bounds.", ioobe.getBadIndex(), ioobe.getBounds()[0], ioobe.getBounds()[1]);
  }
  catch (ZeroDivisionException& zde)
  {
    throw ZeroDivisionException("MultilocusGenotypeStatistics::getHnbForGroups.");
  }
  return 2 * static_cast<double>(nb_alleles) * Hexp  / static_cast<double>((2 * nb_alleles) - 1);
}

double MultilocusGenotypeStatistics::getDnei72(const PolymorphismMultiGContainer& pmgc, vector<size_t> locusPositions, size_t grp1, size_t grp2)
{
  map<size_t, double> allele_frq1, allele_frq2;
  vector<size_t> allele_ids;
  set<size_t> group1_id;
  set<size_t> group2_id;
  set<size_t> groups_id;
  double Jx = 0.;
  double Jy = 0.;
  double Jxy = 0.;
  group1_id.insert(grp1);
  group2_id.insert(grp2);
  groups_id.insert(grp1);
  groups_id.insert(grp2);
  for (size_t i = 0; i < locusPositions.size(); i++)
  {
    allele_ids.clear();
    allele_frq1.clear();
    allele_frq2.clear();
    try
    {
      allele_ids = getAllelesIdsForGroups(pmgc, locusPositions[i], groups_id);
      allele_frq1 = getAllelesFrqForGroups(pmgc, locusPositions[i], group1_id);
      allele_frq2 = getAllelesFrqForGroups(pmgc, locusPositions[i], group2_id);
    }
    catch (Exception& e)
    {
      throw e;
    }
    for (size_t j = 0; j < allele_ids.size(); j++)
    {
      map<size_t, double>::iterator it1 = allele_frq1.find(allele_ids[j]);
      map<size_t, double>::iterator it2 = allele_frq2.find(allele_ids[j]);
      double tmp_frq1 = (it1 != allele_frq1.end()) ? it1->second : 0.;
      double tmp_frq2 = (it2 != allele_frq2.end()) ? it2->second : 0.;
      Jx += tmp_frq1 * tmp_frq1;
      Jy += tmp_frq2 * tmp_frq2;
      Jxy += tmp_frq1 * tmp_frq2;
    }
  }
  if (Jx * Jy == 0.)
    throw ZeroDivisionException("MultilocusGenotypeStatistics::getDnei72.");
  return -log(Jxy / sqrt(Jx * Jy));
}

double MultilocusGenotypeStatistics::getDnei78(const PolymorphismMultiGContainer& pmgc, vector<size_t> locusPositions, size_t grp1, size_t grp2)
{
  map<size_t, double> allele_frq1, allele_frq2;
  vector<size_t> allele_ids;
  set<size_t> group1_id;
  set<size_t> group2_id;
  set<size_t> groups_id;
  double Jx = 0.;
  double Jy = 0.;
  double Jxy = 0.;
  size_t nx = 0, ny = 0;
  group1_id.insert(grp1);
  group2_id.insert(grp2);
  groups_id.insert(grp1);
  groups_id.insert(grp2);
  for (size_t i = 0; i < locusPositions.size(); i++)
  {
    allele_ids.clear();
    allele_frq1.clear();
    allele_frq2.clear();
    try
    {
      allele_ids = getAllelesIdsForGroups(pmgc, locusPositions[i], groups_id);
      allele_frq1 = getAllelesFrqForGroups(pmgc, locusPositions[i], group1_id);
      allele_frq2 = getAllelesFrqForGroups(pmgc, locusPositions[i], group2_id);
      nx = countBiAllelicForGroups(pmgc, locusPositions[i], group1_id);
      ny = countBiAllelicForGroups(pmgc, locusPositions[i], group2_id);
    }
    catch (Exception& e)
    {
      throw e;
    }
    double tmp_Jx = 0.;
    double tmp_Jy = 0.;
    for (size_t j = 0; j < allele_ids.size(); j++)
    {
      map<size_t, double>::iterator it1 = allele_frq1.find(allele_ids[j]);
      map<size_t, double>::iterator it2 = allele_frq2.find(allele_ids[j]);
      double tmp_frq1 = (it1 != allele_frq1.end()) ? it1->second : 0.;
      double tmp_frq2 = (it2 != allele_frq2.end()) ? it2->second : 0.;
      tmp_Jx += tmp_frq1 * tmp_frq1;
      tmp_Jy += tmp_frq2 * tmp_frq2;
      Jxy += tmp_frq1 * tmp_frq2;
    }
    Jx += ((2. * (double) nx * tmp_Jx) - 1.) / ((2. * (double) nx) - 1.);
    Jy += ((2. * (double) ny * tmp_Jy) - 1.) / ((2. * (double) ny) - 1.);
  }
  double denom = Jx * Jy;
  if (denom == 0.)
    throw ZeroDivisionException("MultilocusGenotypeStatistics::getDnei78.");
  return -log(Jxy / sqrt(denom));
}

map<size_t, MultilocusGenotypeStatistics::Fstats> MultilocusGenotypeStatistics::getAllelesFstats(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  map<size_t, MultilocusGenotypeStatistics::VarComp> vc = getVarianceComponents(pmgc, locusPosition, groups);
  map<size_t, MultilocusGenotypeStatistics::Fstats> f_stats;
  for (map<size_t, MultilocusGenotypeStatistics::VarComp>::iterator it = vc.begin(); it != vc.end(); it++)
  {
    double abc = it->second.a + it->second.b + it->second.c;
    double bc  = it->second.b + it->second.c;

    if (abc == 0)
    {
      f_stats[it->first].Fit = NAN;
      f_stats[it->first].Fst = NAN;
    }
    {
      f_stats[it->first].Fit = 1. - it->second.c / abc;
      f_stats[it->first].Fst = it->second.a / abc;
    }
    if (bc == 0)
      f_stats[it->first].Fis = NAN;
    else
      f_stats[it->first].Fis = 1. - it->second.c / bc;
  }
  return f_stats;
}

map<size_t, double> MultilocusGenotypeStatistics::getAllelesFit(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  map<size_t, MultilocusGenotypeStatistics::VarComp> values = getVarianceComponents(pmgc, locusPosition, groups);
  map<size_t, double> Fit;
  for (map<size_t, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin(); it != values.end(); it++)
  {
    Fit[it->first] = it->second.a + it->second.b + it->second.c;
    if (Fit[it->first] == 0.)
      throw ZeroDivisionException("MultilocusGenotypeStatistics::getAllelesFit.");
    Fit[it->first] = 1. - it->second.c / Fit[it->first];
  }
  return Fit;
}

map<size_t, double> MultilocusGenotypeStatistics::getAllelesFst(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  if (groups.size() <= 1)
    throw BadIntegerException("MultilocusGenotypeStatistics::getAllelesFst: groups must be >= 2.", static_cast<int>(groups.size()));
  map<size_t, MultilocusGenotypeStatistics::VarComp> values = getVarianceComponents(pmgc, locusPosition, groups);
  map<size_t, double> Fst;
  for (map<size_t, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin(); it != values.end(); it++)
  {
    Fst[it->first] = it->second.a + it->second.b + it->second.c;
    if (Fst[it->first] == 0.)
      throw ZeroDivisionException("MultilocusGenotypeStatistics::getAllelesFst.");
    Fst[it->first] = it->second.a / Fst[it->first];
  }
  return Fst;
}

map<size_t, double> MultilocusGenotypeStatistics::getAllelesFis(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  map<size_t, MultilocusGenotypeStatistics::VarComp> values = getVarianceComponents(pmgc, locusPosition, groups);
  map<size_t, double> Fis;
  for (map<size_t, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin(); it != values.end(); it++)
  {
    Fis[it->first] = it->second.b + it->second.c;
    if (Fis[it->first] == 0.)
      throw ZeroDivisionException("MultilocusGenotypeStatistics::getAllelesFis.");
    Fis[it->first] = 1. - it->second.c / Fis[it->first];
  }
  return Fis;
}

map<size_t, MultilocusGenotypeStatistics::VarComp> MultilocusGenotypeStatistics::getVarianceComponents(const PolymorphismMultiGContainer& pmgc, size_t locusPosition, const set<size_t>& groups)
{
  map<size_t, MultilocusGenotypeStatistics::VarComp> values;
  // Base values computation
  double nbar = 0.;
  double nc = 0.;
  vector<size_t> ids = getAllelesIdsForGroups(pmgc, locusPosition, groups);
  map<size_t, double> pbar;
  map<size_t, double> s2;
  map<size_t, double> hbar;
  for (size_t i = 0; i < ids.size(); i++)
  {
    pbar[ids[i]] = 0.;
    s2[ids[i]] = 0.;
    hbar[ids[i]] = 0.;
  }
  double r = static_cast<double>(groups.size());
  for (set<size_t>::iterator set_it = groups.begin(); set_it != groups.end(); set_it++)
  {
    size_t i  = (*set_it);
    double ni = static_cast<double>(pmgc.getLocusGroupSize(i, locusPosition));
    set<size_t> group_id;
    group_id.insert( i );
    map<size_t, double> pi = getAllelesFrqForGroups(pmgc, locusPosition, group_id);
    map<size_t, double> hi = getHeterozygousFrqForGroups(pmgc, locusPosition, group_id);
    nbar += ni;
    if (r > 1)
      nc += ni * ni;

    for (map<size_t, double>::iterator it = pi.begin(); it != pi.end(); it++)
    {
      pbar[it->first] += ni * it->second;
    }
    for (map<size_t, double>::iterator it = hi.begin(); it != hi.end(); it++)
    {
      hbar[it->first] += ni * it->second;
    }

    group_id.clear();
  }
  nbar = nbar / r;
  if (nbar <= 1)
    throw ZeroDivisionException("MultilocusGenotypeStatistics::getVarianceComponents.");
  if (r > 1)
    nc = (r * nbar) - (nc / (r * nbar)) / (r - 1.);
  for (map<size_t, double>::iterator it = pbar.begin(); it != pbar.end(); it++)
  {
    it->second = it->second / (r * nbar);
  }
  for (map<size_t, double>::iterator it = hbar.begin(); it != hbar.end(); it++)
  {
    it->second = it->second / ( r * nbar);
  }

  for (set<size_t>::iterator set_it = groups.begin(); set_it != groups.end(); set_it++)
  {
    size_t i  = (*set_it);
    double ni = static_cast<double>(pmgc.getLocusGroupSize( i, locusPosition));
    set<size_t> group_id;
    group_id.insert( i );
    map<size_t, double> pi = getAllelesFrqForGroups(pmgc, locusPosition, group_id);
    for (size_t j = 0; j < ids.size(); j++)
    {
      pi[ids[j]];
    }
    for (map<size_t, double>::iterator it = pi.begin(); it != pi.end(); it++)
    {
      s2[it->first] += ni * (it->second - pbar[it->first]) * (it->second - pbar[it->first]);
    }
    group_id.clear();
  }
  for (map<size_t, double>::iterator it = s2.begin(); it != s2.end(); it++)
  {
    it->second = it->second / ((r - 1.) * nbar);
  }

  // a, b, c computation
  for (size_t i = 0; i < ids.size(); i++)
  {
    values[ids[i]];
  }

  for (map<size_t, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin(); it != values.end(); it++)
  {
    it->second.a = (nbar / nc) * (s2[it->first] - ((1. / (nbar - 1.)) * ((pbar[it->first] * (1. - pbar[it->first])) - (s2[it->first] * ((double) r - 1.) / r) - ((1. / 4.) * hbar[it->first]))));
    it->second.b = (nbar / (nbar - 1.)) * ((pbar[it->first] * (1. - pbar[it->first])) - (s2[it->first] * ((double) r - 1.) / (double) r) - ((((2. * nbar) - 1.) / (4. * nbar)) * hbar[it->first]));
    it->second.c = hbar[it->first] / 2.;
  }
  return values;
}

double MultilocusGenotypeStatistics::getWCMultilocusFst(const PolymorphismMultiGContainer& pmgc, vector<size_t> locusPositions, const set<size_t>& groups)
{
  double A, B, C;
  A = B = C = 0.0;
  for (size_t i = 0; i < locusPositions.size(); i++)
  {
    // count total number of individuals without missing data
    size_t ni = 0;
    for (set<size_t>::iterator setIt = groups.begin(); setIt != groups.end(); setIt++)
    {
      ni += pmgc.getLocusGroupSize( (*setIt), i);
    }

    // reduce computation for polymorphic loci for that groups
    vector<size_t> ids = getAllelesIdsForGroups(pmgc, i, groups);
    if (ids.size() >= 2 && ni >= 1)
    {
      map<size_t, MultilocusGenotypeStatistics::VarComp> values = getVarianceComponents(pmgc, locusPositions[i], groups);
      for (map<size_t, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin(); it != values.end(); it++)
      {
        A += it->second.a;
        B += it->second.b;
        C += it->second.c;
      }
    }
  }
  if ((A + B + C) == 0)
    throw ZeroDivisionException("MultilocusGenotypeStatistics::getWCMultilocusFst.");
  return A / (A + B + C);
}

double MultilocusGenotypeStatistics::getWCMultilocusFis(const PolymorphismMultiGContainer& pmgc, vector<size_t> locusPositions, const set<size_t>& groups)
{
  double B, C;
  B = C = 0.0;
  for (size_t i = 0; i < locusPositions.size(); i++)
  {
    // count total number of individuals without missing data
    size_t ni = 0;
    for (set<size_t>::iterator setIt = groups.begin(); setIt != groups.end(); setIt++)
    {
      ni += pmgc.getLocusGroupSize( (*setIt), i);
    }

    // reduce computation for polymorphic loci for that groups
    vector<size_t> ids = getAllelesIdsForGroups(pmgc, i, groups);
    if (ids.size() >= 2 && ni >= 1)
    {
      map<size_t, MultilocusGenotypeStatistics::VarComp> values = getVarianceComponents(pmgc, locusPositions[i], groups);
      for (map<size_t, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin(); it != values.end(); it++)
      {
        B += it->second.b;
        C += it->second.c;
      }
    }
  }
  if ((B + C) == 0)
    throw ZeroDivisionException("MultilocusGenotypeStatistics::getWCMultilocusFis.");
  return 1.0 - C / (B + C);
}

MultilocusGenotypeStatistics::PermResults MultilocusGenotypeStatistics::getWCMultilocusFstAndPerm(
    const PolymorphismMultiGContainer& pmgc,
    vector<size_t> locusPositions,
    set<size_t> groups,
    unsigned int nbPerm)
{
  // extract a PolymorphismMultiGContainer with only those groups
  auto subPmgc = PolymorphismMultiGContainerTools::extractGroups(pmgc, groups);
  double nbSup = 0.0;
  double nbInf = 0.0;
  PermResults results;
  results.statistic = getWCMultilocusFst(*subPmgc, locusPositions, groups);
  if (nbPerm > 0)
  {
    for (size_t i = 0; i < nbPerm; ++i)
    {
      auto permutedPmgc = PolymorphismMultiGContainerTools::permuteMultiG(*subPmgc);
      double FstPerm = getWCMultilocusFst(*permutedPmgc, locusPositions, groups);
      // cout << Fst_perm << endl;
      if (FstPerm > results.statistic)
        nbSup++;
      if (FstPerm < results.statistic)
        nbInf++;
    }

    nbSup /= static_cast<double>(nbPerm);
    nbInf /= static_cast<double>(nbPerm);
  }

  results.percentSup = nbSup;
  results.percentInf = nbInf;
  return results;
}

MultilocusGenotypeStatistics::PermResults MultilocusGenotypeStatistics::getWCMultilocusFisAndPerm(
    const PolymorphismMultiGContainer& pmgc,
    vector<size_t> locusPositions,
    set<size_t> groups,
    unsigned int nbPerm)
{
  // extract a PolymorphismMultiGContainer with only those groups
  auto subPmgc =  PolymorphismMultiGContainerTools::extractGroups(pmgc, groups);
  double nbSup = 0.0;
  double nbInf = 0.0;
  PermResults results;
  results.statistic =  getWCMultilocusFis(*subPmgc, locusPositions, groups);
  if (nbPerm > 0)
  {
    for (unsigned int i = 0; i < nbPerm; ++i)
    {
      auto permutedPmgc = PolymorphismMultiGContainerTools::permuteIntraGroupAlleles(*subPmgc, groups);
      double FisPerm = getWCMultilocusFis(*permutedPmgc, locusPositions, groups);

      if (FisPerm > results.statistic)
        nbSup++;
      if (FisPerm < results.statistic)
        nbInf++;
    }

    nbSup /= static_cast<double>(nbPerm);
    nbInf /= static_cast<double>(nbPerm);
  }

  results.percentSup = nbSup;
  results.percentInf = nbInf;
  return results;
}

double MultilocusGenotypeStatistics::getRHMultilocusFst(
    const PolymorphismMultiGContainer& pmgc,
    vector<size_t> locusPositions,
    const set<size_t>& groups)
{
  double Au, Bu, Cu;
  double RH = 0.0;
  int nb_alleles = 0;
  int total_alleles = 0;

  for (size_t i = 0; i < locusPositions.size(); i++)
  {
    // reduce computation for polymorphic loci for that groups
    vector<size_t> ids = getAllelesIdsForGroups(pmgc, i, groups);
    if (ids.size() >= 2)
    {
      nb_alleles = 0;
      // mean allelic frequencies
      map<size_t, double> P = MultilocusGenotypeStatistics::getAllelesFrqForGroups (pmgc, locusPositions[i], groups);
      // variance components from W&C
      map<size_t, MultilocusGenotypeStatistics::VarComp> values = getVarianceComponents(pmgc, locusPositions[i], groups);
      for (map<size_t, MultilocusGenotypeStatistics::VarComp>::iterator it = values.begin(); it != values.end(); it++)
      {
        Au = it->second.a;
        Bu = it->second.b;
        Cu = it->second.c;
        if ((Au + Bu + Cu) != 0)
        {
          double Pu = P[it->first]; // it->first is the allele number
          RH += (1 - Pu) * Au / (Au + Bu + Cu);
          nb_alleles++;
        }
      }
      total_alleles += (nb_alleles - 1);
    }
  }
  if (total_alleles == 0)
    throw ZeroDivisionException("MultilocusGenotypeStatistics::getRHMultilocusFst.");
  return RH / double(total_alleles);
}

std::unique_ptr<DistanceMatrix> MultilocusGenotypeStatistics::getDistanceMatrix(const PolymorphismMultiGContainer& pmgc, vector<size_t> locusPositions, const set<size_t>& groups, string distance_methode)
{
  vector<string> names = pmgc.getAllGroupsNames();
  vector<size_t> grp_ids_vect;
  for (set<size_t>::iterator i = groups.begin(); i != groups.end(); i++)
  {
    grp_ids_vect.push_back(*i);
  }

  unique_ptr<DistanceMatrix> _dist(new DistanceMatrix(names));
  for (size_t i = 0; i < groups.size(); i++)
  {
    (*_dist)(i, i) = 0;
  }

  set<size_t> pairwise_grp;

  for (size_t j = 0; j < groups.size () - 1; j++)
  {
    for (size_t k = j + 1; k < groups.size (); k++)
    {
      double distance = 0;
      if (distance_methode ==  "nei72")
        distance = MultilocusGenotypeStatistics::getDnei72( pmgc, locusPositions, grp_ids_vect[j], grp_ids_vect[k] );
      else if  (distance_methode == "nei78")
        distance = MultilocusGenotypeStatistics::getDnei78( pmgc, locusPositions, grp_ids_vect[j], grp_ids_vect[k] );
      else if (distance_methode == "WC") // Fst multilocus selon W&C
      {
        pairwise_grp.insert(grp_ids_vect[j] );
        pairwise_grp.insert(grp_ids_vect[k] );
        distance = MultilocusGenotypeStatistics::getWCMultilocusFst( pmgc, locusPositions, pairwise_grp);
        pairwise_grp.clear();
      }
      else if (distance_methode == "RH") // Fst multilocus selon ponderation Robertson & Hill
      {
        pairwise_grp.insert(grp_ids_vect[j] );
        pairwise_grp.insert(grp_ids_vect[k] );
        distance = MultilocusGenotypeStatistics::getRHMultilocusFst( pmgc, locusPositions, pairwise_grp);
        pairwise_grp.clear();
      }
      else if (distance_methode == "Nm") // Nm déduit des Fst multilocus selon W&C modèle en îles Fst = 1/(1+4Nm)
      {
        pairwise_grp.insert(grp_ids_vect[j] );
        pairwise_grp.insert(grp_ids_vect[k] );
        distance = MultilocusGenotypeStatistics::getWCMultilocusFst( pmgc, locusPositions, pairwise_grp);
        if (distance != 0)
          distance = 0.25 * (1 - distance) / distance;
        else
          distance = NAN;
        pairwise_grp.clear();
      }
      else if (distance_methode == "D") // D=-ln(1-Fst) of Reynolds, Weir and Cockerham, 1983
      {
        pairwise_grp.insert(grp_ids_vect[j] );
        pairwise_grp.insert(grp_ids_vect[k] );
        distance = MultilocusGenotypeStatistics::getWCMultilocusFst( pmgc, locusPositions, pairwise_grp);
        if (distance != 1)
          distance =  -log(1 - distance);
        else
          distance = NAN;
        pairwise_grp.clear();
      }
      else if (distance_methode == "Rousset") // Calcul de Fst/(1-Fst). Rousset F. 1997
      {
        pairwise_grp.insert(grp_ids_vect[j] );
        pairwise_grp.insert(grp_ids_vect[k] );
        distance = MultilocusGenotypeStatistics::getWCMultilocusFst( pmgc, locusPositions, pairwise_grp);
        if (distance != 1)
          distance = distance / (1 - distance);
        else
          distance = NAN;
        pairwise_grp.clear();
      }

      (*_dist)(k, j) =  distance;
      (*_dist)(j, k) =  distance;
    } // for k
  } // for j

  return _dist;
}
