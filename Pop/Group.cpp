/*
 * File Group.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Wednesday June 01 2004
 */

#include "Group.h"

//** Class constructor: *******************************************************/
Group::Group() {}

//** Class destructor: *******************************************************/
Group::~Group () {}

//** Assignation opperators: *************************************************/
void Group::addIndividual(const Individual & ind) {
	_individuals.push_back(new Individual(ind));
}

void Group::clear() {
	_individuals.clear();
}
//** Consultation opperators: ************************************************/
int Group::getNumberOfIndividuals() {
	return _individuals.size();
}
