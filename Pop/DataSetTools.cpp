/*
 * File DataSetTools.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday August 04 2004
 *
 * Copyright (C) 2004 Sylvain Gaillard and the
 *                    PopGenLib Development Core Team
 *
 * This file is part of PopGenLib.
 *
 * PopGenLib is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * PopGenLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PopGenLib; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "DataSetTools.h"

DataSet * DataSetTools::buildDataSet(const OrderedSequenceContainer & osc) throw (Exception) {
	DataSet * d_s = new DataSet();
	d_s->addEmptyGroup(0);
	for (unsigned int i = 0 ; i < osc.getNumberOfSequences() ; i++) {
		d_s->addEmptyIndividualToGroup(0, string("Individual_") + TextTools::toString(i + 1));
		try {
			d_s->addIndividualSequenceInGroup(0, i, 0, * osc.getSequence(i));
		}
		catch (Exception & e) {
			throw e;
		}
	}
	return d_s;
}

DataSet * DataSetTools::buildDataSet(const PolymorphismSequenceContainer & psc) throw (Exception) {
	DataSet * d_s = new DataSet();
	set<unsigned int> grp_ids = psc.getAllGroupsIds();
	for (set<unsigned int>::iterator it = grp_ids.begin() ; it != grp_ids.end() ; it++)
		d_s->addEmptyGroup(* it);
	unsigned int ind_count = 0;
	for (unsigned int i = 0 ; i < psc.getNumberOfSequences() ; i++) {
		for (unsigned int j = 0 ; j < psc.getSequenceCount(i) ; j++) {
			d_s->addEmptyIndividualToGroup(psc.getGroupId(i), string("Individual_") + TextTools::toString(ind_count++));
			try {
				d_s->addIndividualSequenceInGroup(psc.getGroupId(i), i, 0, * psc.getSequence(i));
			}
			catch (Exception & e) {
				throw e;
			}
		}
	}
	return d_s;
}
