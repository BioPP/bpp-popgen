/*
 * File AnalyzedSequences.h
 * Author : Sylvain Gaillard <yragael2001@yahoo.fr>
 * Last modification : Thursday July 08 2004
 */

// Secured inclusion of header's file
#ifndef _ANALYZEDSEQUENCES_H_
#define _ANALYZEDSEQUENCES_H_

// From Seq
#include <Seq/Alphabet.h>
#include <Seq/DNA.h>
#include <Seq/RNA.h>
#include <Seq/ProteicAlphabet.h>

/**
 * @brief The AnalyzedSequences class.
 *
 * This is a class to store info about the sequences.
 */
class AnalyzedSequences {
	public: // Constructor and destructor
		AnalyzedSequences();
		~AnalyzedSequences();

	public:
		/**
		 * @brief Set the alphabet used for the sequences.
		 */
		void setAlphabet(const Alphabet * alpha);

		/**
		 * @brief Set the alphabet used for the sequences by alphabet type.
		 */
		void setAlphabet(const string & alpha_type) throw (Exception);

		/**
		 * @brief Get the alphabet.
		 */
		const Alphabet * getAlphabet() const;

		/**
		 * @brief Get the alphabet type as a string.
		 */
		string getAlphabetType() const;

	protected:
		const Alphabet * _alphabet;
};

#endif // _ANALYZEDSEQUENCES_H_
