/*
 * File PolymorphismMultiGContainer.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Tuesday July 20 2004
 */

// Secured inclusion of header's file
#ifndef _POLYMORPHYSMMULTIGCONTAINER_H_
#define _POLYMORPHYSMMULTIGCONTAINER_H_

// From STL
#include <vector>
using namespace std;

// From Utils
#include <Utils/Clonable.h>
#include <Utils/Exceptions.h>

// From poplib
#include "MultilocusGenotype.h"

/**
 * @brief The PolymorphismMultiGContainer class
 *
 * This class is a container of MultilocusGenotype.
 */
class PolymorphismMultiGContainer : public Clonable {
	public: // Constructors and destructor
		/**
		 * @brief Build a new PolymorphismMultilocusGenotypeContainer.
		 */
		PolymorphismMultiGContainer();

		/**
		 * @brief The copy constructor.
		 */
		PolymorphismMultiGContainer(const PolymorphismMultiGContainer & pmgc);

		/**
		 * @brief Destroy a PolymorphismMultilocusGenotypeContainer.
		 */
		~PolymorphismMultiGContainer();

	public:
		/**
		 * @brief The assignation operator=.
		 */
		PolymorphismMultiGContainer & operator= (const PolymorphismMultiGContainer & pmgc);

		/**
		 * @name The clonable interface
		 * @{
		 */
		Clonable * clone() const;
		/** @} */

		/**
		 * @brief Add a MultilocusGenotype to the container.
		 */
		void addMultilocusGenotype(const MultilocusGenotype & mg, unsigned int group);

		/**
		 * @brief Get a MultilocusGenotype at a position.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		const MultilocusGenotype * getMultilocusGenotype(unsigned int position) const throw (IndexOutOfBoundsException);

		/**
		 * @brief Remove a MultilocusGenotype.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		MultilocusGenotype * removeMultilocusGenotype(unsigned int position) throw (IndexOutOfBoundsException);

		/**
		 * @brief Delete a MultilocusGenotype.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		void deleteMultilocusGenotype(unsigned int position) throw (IndexOutOfBoundsException);

		/**
		 * @brief Get the Group of a MultilocusGenotype.
		 *
		 * @throw IndexOutOfBoundsException if position excedes the size of the container.
		 */
		unsigned int getGroup(unsigned int position) const throw (IndexOutOfBoundsException);
		
		/**
		 * @brief Get the number of MultilocusGenotype.
		 */
		unsigned int size() const;
		
		/**
		 * @brief Clear the container.
		 */
		void clear();

	protected:
		vector<MultilocusGenotype *> _multilocusGenotypes;
		vector<unsigned int> _groups;
};

#endif // _POLYMORPHYSMMULTIGCONTAINER_H_
