/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and 
Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
(c) Copyright 2020 by Joseph Guscott, Daniel Barker, Miguel Pinheiro,
Fernando Guntoro, Maximilian Strobl, Chang Sik Kim, Martyn Winn and Chris Wood.

All rights reserved.
 
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


/* ========== ReadFiles.h - interface for ReadFiles.cpp ========== */

#ifndef CREADFILES_H_ 
#define CREADFILES_H_

#include <string>
#include <string.h>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <algorithm>

using namespace std;

class CReadFiles {

public:
	static const int FORMAT_PHYLIP_ = 0;
	static const int FORMAT_FASTA_ = 1;
	static const int FORMAT_NEXUS_ = 2;
	static const int FORMAT_CLUSTAL_ = 3;


public:
	/// this is because some compilers has a problem with to_string method
	template < typename T > std::string to_string_( const T& n )
	{
		std::ostringstream stm ;
		stm << n ;
		return stm.str() ;
	}

public:
	CReadFiles();
	virtual ~CReadFiles();
	static int exit_error(int ierr, std::string sz_error);
	static bool is_file_exist(std::string file_name);

	void save_file(std::string sz_file_name_temp);
	int read_file(std::string file_name_out, int n_file_type);
	unsigned int get_length_sequences() {
		if ((int) lst_sequences.size() > 0) return lst_sequences[0].length();
		return 0;
	}
	int get_number_seqs() { return (int) lst_names_seq.size(); }
	int get_max_length_seq_name(){ return n_max_length_name_seq; }
	int get_length_seq_name(int n_seq) { return (int) lst_names_seq[n_seq].length(); }
	char get_char_seq_name(int n_seq, int n_pos_char) { return (char) lst_names_seq[n_seq].at(n_pos_char); }
	char get_char_sequences(int n_seq, int n_pos_char) { return (char) lst_sequences[n_seq].at(n_pos_char); }

/// data structure
private:
	std::vector< std::string > lst_sequences;
	std::vector< std::string > lst_names_seq;
	std::string sz_extension;			/// extension of the file
	std::string sz_file_name;			/// file name possible all path
	std::string sz_only_file_name;			/// only yhe file name
	std::string sz_accept_chars;			// chars to pass on filter
	int n_max_length_name_seq;

	int clean_data();

private:
	/// several read file methods
	int read_clustal(int filetype);
	int read_phylip();
	int read_fasta();
	int read_nexus();

	bool is_only_one_sequence_in_array(std::vector< std::string > lst_strings);
	std::string get_string_from_list(std::vector< std::string > lst_strings, bool b_last_one);

	std::string trim(std::string const& str);
	std::string trim(std::string const& str, char c_char);
	std::vector<std::string> split(const std::string &s, char delim);
	std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);

	/// methods used in phylip files
	bool is_phylip_interleaved();
	bool is_phylip_line_sequential(int &n_total_lines);
	std::string clean_phylip_dna_sequence(std::string sz_sequence);

	/// used for phylip files, is the max length of the namess
	int n_nmlngth_phylip_names;
	std::string sz_phylip_accept_chars;

	//// trim chars form the strings
	std::vector<char> vect_trim_char;
	bool b_debug;

	/// members for nexus format
	void get_dimensions_nexus_format(std::string sz_line, int &n_seqs, int &n_length_seq);

};

#endif /* CREADFILES_H_ */ 
