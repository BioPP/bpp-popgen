/*
 * File MonolocusGenotype.cpp
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday May 27 2004
 */

#include "MonolocusGenotype.h"

//** Class constructor: *******************************************************/
MonolocusGenotype::MonolocusGenotype() {};

//** Class destructor: ********************************************************/
MonolocusGenotype::~MonolocusGenotype() {};

//** Other methodes: **********************************************************/
void MonolocusGenotype::addKey(const unsigned int key) {
	_allelekeys.push_back(key);
}

unsigned int MonolocusGenotype::getKey(unsigned int index)
throw (IndexOutOfBoundsException) {
	if (index >= 0 && index < _allelekeys.size())
		return _allelekeys[index];
	else
		throw IndexOutOfBoundsException("MonolocusGenotype::getKey: index out of bounds",
				index, 0, _allelekeys.size()-1);
}

vector<unsigned int> MonolocusGenotype::getKeys() {
	vector<unsigned int> keys;
	for (unsigned int i = 0 ; i < _allelekeys.size() ; i++)
		if (_allelekeys[i] != 0) keys.push_back(_allelekeys[i]);
	return keys;
}

unsigned int MonolocusGenotype::getNumberOfData() {
	return _allelekeys.size();
}

bool MonolocusGenotype::hasMissingData() {
	for (unsigned int i = 0 ; i < _allelekeys.size() ; i++)
		if (_allelekeys[i] == 0) return true;
	return false;
}
