/* LVB

(c) Copyright 2003-2012 by Daniel Barker
(c) Copyright 2013, 2014 by Daniel Barker and Maximilian Strobl
(c) Copyright 2014 by Daniel Barker, Miguel Pinheiro, and Maximilian Strobl
(c) Copyright 2015 by Daniel Barker, Miguel Pinheiro, Maximilian Strobl,
and Chris Wood.
(c) Copyright 2019 by Daniel Barker, Miguel Pinheiro, Joseph Guscott,
Fernando Guntoro, Maximilian Strobl and Chris Wood.
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

/* ========== CReadFiles.cpp - data matrix inputs ========== */

#include "MSAInput.h"

CReadFiles::CReadFiles() {

	//// trim chars form the strings
	vect_trim_char.push_back('\n');
	vect_trim_char.push_back('\r');
	vect_trim_char.push_back('\n');
//	vect_trim_char.push_back('\t');
	vect_trim_char.push_back(' ');

	/// these lines are getted from phylip code version 3.695 at 18-Jun-2014
	/// src/phylip.h:#define nmlngth         10   /* number of characters in species name    */
	/// src/phylip.h:#define MAXNCH          20   /* must be greater than or equal to nmlngth */
	n_nmlngth_phylip_names = 10;
	sz_phylip_accept_chars = "ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz*?-";
	sz_accept_chars = "ACGTUYRWSKMBDHVNX?O-";
	n_max_length_name_seq = 0;
	b_debug = false;		/// print some messages
}

CReadFiles::~CReadFiles() {
	// TODO Auto-generated destructor stub
}

bool CReadFiles::is_file_exist(std::string file_name)
{
	if (FILE *file = fopen(file_name.c_str(), "r")) {
		fclose(file);
		return true;
	} else { return false; }
}

int CReadFiles::exit_error(int ierr, std::string sz_error)
{
	if (ierr == 1) cout << "FATAL ERROR: " << sz_error << endl;
	else cout << "Warning: " << sz_error << endl;

	/* Exit if necessary */
	#ifndef NP_Implementation
	if (ierr==1 || ierr==-1) { return 1; }
	#else
	if (ierr==1 || ierr==-1) { exit(1); }
	#endif
	return 0;
}

std::string CReadFiles::trim(std::string const& str)
{
	string sz_return = str;
	for (std::vector<char >::iterator it = vect_trim_char.begin(); it != vect_trim_char.end(); ++it){
		sz_return = trim(sz_return, *it);
	}
    return sz_return;
}

std::string CReadFiles::trim(std::string const& str, char c_char)
{
    std::size_t first = str.find_first_not_of(c_char);
    if(first == std::string::npos) return "";
    std::size_t last  = str.find_last_not_of(c_char);
    return str.substr(first, last-first+1);
}

std::vector<std::string> &CReadFiles::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s); std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> CReadFiles::split(const std::string &s, char delim) {
    std::vector<std::string> elems, vect_out;
    split(s, delim, elems);
    for (std::vector<std::string>::iterator it = elems.begin() ; it != elems.end(); ++it) vect_out.push_back(trim(*it));
    return elems;
}

int CReadFiles::read_clustal(int filetype)
{
	std::string sz_line = "";
	std::vector< std::string > lst_strings;
	int header = 0, length = 0, n_count_sequence = 0;
	bool b_is_name_seq = false;

	ifstream filein;
	filein.open(sz_file_name.c_str());
	while (std::getline(filein, sz_line)) {

		length = (int) sz_line.length();
		if (header == 0 && length > 1 && filetype == 1) {
		  if ((int) sz_line.find("CLUSTAL") != -1){
			  header = 1;
			  continue;
		  }
		}
		else if (header == 0 && length > 1 && filetype == 2) {
		  if ((int) sz_line.find("//") != -1){
			  header = 1;
			  continue;
		  }
		}
		else if(length > 0 and header == 1) { /// header is here now

			lst_strings = split(sz_line, ' ');
			if ((int) lst_strings.size() == 0) continue; /// nothing to do

			if (is_only_one_sequence_in_array(lst_strings) || lst_strings[0].empty()){
				/// escape for the **** **************.*********.**.*************** ***
				if (! b_is_name_seq){
					b_is_name_seq = true;
					continue;
				}
				/// escape for the **** **************.*********.**.********** ********
				if (n_count_sequence == (int) lst_names_seq.size()){
					n_count_sequence = 0;
					continue;
				}
			}
			if (! b_is_name_seq){
				if ((int) lst_strings.size() > 1){
					lst_names_seq.push_back(get_string_from_list(lst_strings, false));
					lst_sequences.push_back(get_string_from_list(lst_strings, true));
				}
			}
			else{
				#ifndef NP_Implementation
				if (n_count_sequence >= (int) lst_sequences.size()) return exit_error(1 , "Some problem reading the file.");
				#else
				if (n_count_sequence >= (int) lst_sequences.size()) exit_error(1 , "Some problem reading the file.");
				#endif
				lst_sequences[n_count_sequence] += get_string_from_list(lst_strings, true);
				n_count_sequence ++;
			}
		}
		else{
			n_count_sequence = 0;
			if ((int) lst_names_seq.size() > 0) b_is_name_seq = true;
		}
	}
	filein.close();
	#ifndef NP_Implementation
	return 0;
	#endif
}

bool CReadFiles::is_only_one_sequence_in_array(std::vector< std::string > lst_strings){

	int n_number_seqs = 0;
	for (int i = 0; i < (int) lst_strings.size(); ++i){
		if ((int) lst_strings[i].length() > 0) n_number_seqs += 1;
	}
	if (n_number_seqs == 1) return true;
	return false;
}

std::string CReadFiles::get_string_from_list(std::vector< std::string > lst_strings, bool b_last_one){

	int n_last_one_with_sequence = -1;
	for (int i = (int) lst_strings.size() - 1; i > -1; --i){
		if ((int) lst_strings[i].length() > 0){
			if (b_last_one){
				if (b_debug) cout << "last seq: " << lst_strings[i] << endl;
				return lst_strings[i];
			}
			n_last_one_with_sequence = i;
			break;
		}
	}
	if (b_last_one) return "";

	//// test the first one
	std::string sz_return = "";
	for (int i = 0; i < n_last_one_with_sequence; ++i){
		if ((int) lst_strings[i].length() > 0){
			if ((int) sz_return.length() > 0) sz_return += " ";
			sz_return += lst_strings[i];
		}
	}
	if (b_debug) cout << "first seq: " << sz_return << endl;
	return sz_return;
}

/// these lines are getted from phylip code version 3.695 at 18-Jun-2014
/// src/phylip.h:#define nmlngth         10   /* number of characters in species name    */
/// src/phylip.h:#define MAXNCH          20   /* must be greater than or equal to nmlngth */
/// need to remove numbers from the sequences...
int CReadFiles::read_phylip()
{
/// interleaved example, examples from the phylip web site
//	  5    42
//	Turkey    AAGCTNGGGC ATTTCAGGGT
//	Salmo gairAAGCCTTGGC AGTGCAGGGT
//	H. SapiensACCGGTTGGC CGTTCAGGGT
//	Chimp     AAACCCTTGC CGTTACGCTT
//	Gorilla   AAACCCTTGC CGGTACGCTT
//
//	GAGCCCGGGC AATACAGGGT AT
//	GAGCCGTGGC CGGGCACGGT AT
//	ACAGGTTGGC CGTTCAGGGT AA
//	AAACCGAGGC CGGGACACTC AT
//	AAACCATTGC CGGTACGCTT AA
	bool b_is_interlevead = is_phylip_interleaved();

/// interleaved sequential example, examples from the phylip web site
//	5    42
//	Turkey    AAGCTNGGGC ATTTCAGGGT
//	GAGCCCGGGC AATACAGGGT AT
//	Salmo gairAAGCCTTGGC AGTGCAGGGT
//	GAGCCGTGGC CGGGCACGGT AT
//	H. SapiensACCGGTTGGC CGTTCAGGGT
//	ACAGGTTGGC CGTTCAGGGT AA
//	Chimp     AAACCCTTGC CGTTACGCTT
//	AAACCGAGGC CGGGACACTC AT
//	Gorilla   AAACCCTTGC CGGTACGCTT
//	AAACCATTGC CGGTACGCTT AA
	int n_total_lines;
	bool b_is_line_sequencial = is_phylip_line_sequential(n_total_lines);

	//// get the number of sequences and length of sequences...
	int n_seqs = 0, n_length_seq = 0, n_count_line = 0, length;
	std::string sz_line = "", error;
	bool b_is_start_read = false;

	ifstream filein;
	filein.open(sz_file_name.c_str());
	while (std::getline(filein, sz_line)) {

		sz_line = trim(sz_line);
		if (sz_line.empty()) continue;

		/// try to get the number of sequences and its length
		if (n_seqs == 0 && n_length_seq == 0){
			int scan_value = std::sscanf(sz_line.c_str(), "%d%d", &n_seqs, &n_length_seq);
			if (scan_value != 2){
				n_seqs = 0;
				n_length_seq = 0;
			}
			continue;
		}

		length = (int) sz_line.length();
		if(b_is_interlevead || (!b_is_line_sequencial && !b_is_interlevead)){
			if (!b_is_start_read){
				if (sz_phylip_accept_chars.find(sz_line[length - 1]) != string::npos) b_is_start_read = true;
				else continue;
			}

			if ((int)lst_names_seq.size() < n_seqs){		///
				/// the clustal put other headers in phylip format sometimes...
				if (length < n_nmlngth_phylip_names){
					error = "Failed to open alignment file: " + sz_file_name + "\nThis line has less than " + to_string_(n_nmlngth_phylip_names) + " chars...";
					#ifndef NP_Implementation
					return exit_error(1 , error);
					#else
					exit_error(1, error);
					#endif
					
				}
				lst_names_seq.push_back(trim(sz_line.substr(0, n_nmlngth_phylip_names)));
				lst_sequences.push_back(sz_line.substr(n_nmlngth_phylip_names, (int) sz_line.length() -  n_nmlngth_phylip_names));
			}
			else{
				if (sz_line.empty()) continue;
				lst_sequences[n_count_line] += sz_line;
				n_count_line += 1;
				if (n_count_line == n_seqs) n_count_line = 0;
			}
		}
		else if (b_is_line_sequencial){
			if (!b_is_start_read){
				if (sz_phylip_accept_chars.find(sz_line[length - 1]) != string::npos) b_is_start_read = true;
				else continue;
			}


			if (n_count_line == 0){		///
				/// the clustal put other headers in phylip format sometimes...
				if (length < n_nmlngth_phylip_names){
					error = "Failed to open alignment file: " + sz_file_name + "\nThis line has less than " + to_string_(n_nmlngth_phylip_names) + " chars...";
					#ifndef NP_Implementation
					return exit_error(1 , error);
					#else
					exit_error(1, error);
					#endif
				}
				lst_names_seq.push_back(trim(sz_line.substr(0, n_nmlngth_phylip_names)));
				lst_sequences.push_back(sz_line.substr(n_nmlngth_phylip_names, (int) sz_line.length() -  n_nmlngth_phylip_names));
				if ((int) lst_names_seq.size() > n_seqs){
					error = "There are more sequences than the predict ones. Please, check the file.";
					#ifndef NP_Implementation
					return exit_error(1 , error);
					#else
					exit_error(1, error);
					#endif
				}
				n_count_line += 1;
			}
			else{
				lst_sequences[(int) lst_sequences.size() - 1] += sz_line;
				n_count_line += 1;
				if ((n_total_lines / n_seqs) == n_count_line) n_count_line = 0;
			}
		}
		/// try to find the first line with sequence data...
	}


	if (n_seqs == 0){
		error = "Some problem reading the file. Please, check the file format.\nMore information in http://evolution.genetics.washington.edu/phylip/doc/sequence.html";
		filein.close();
		#ifndef NP_Implementation
				return exit_error(1 , error);
			#else
				exit_error(1, error);
			#endif
	}
	if ((int) lst_names_seq.size() != (int) lst_sequences.size()){
		error = "The number of sequences read are different from the number of name samples.";
		#ifndef NP_Implementation
				return exit_error(1 , error);
			#else
				exit_error(1, error);
			#endif
	}

	/// clean the numbers and spaces...
	if ((int) lst_names_seq.size() != n_seqs){
		error = "The file has a different number of sequences.\nRead: " + to_string_((int) lst_names_seq.size()) + "\nIn the header: " + to_string_(n_seqs);
		#ifndef NP_Implementation
				return exit_error(1 , error);
			#else
				exit_error(1, error);
			#endif
	}
	for(int i = 0; i < (int) lst_sequences.size(); i++){
		lst_sequences[i] = clean_phylip_dna_sequence(lst_sequences[i]);
		if ((int) lst_sequences[i].size() != n_length_seq){
			cout << lst_sequences[i] << endl;
			error = "This sequence " + lst_names_seq[i] + " has a different length " + to_string_((int) lst_sequences[i].size()) + " from the one read in the header: " + to_string_(n_length_seq);
			#ifndef NP_Implementation
				return exit_error(1 , error);
			#else
				exit_error(1, error);
			#endif
		}
	}
	filein.close();
	#ifndef NP_Implementation
	return 0;
	#endif
}


std::string CReadFiles::clean_phylip_dna_sequence(std::string sz_sequence){
	char chars[] = "0123456789 ";
	for (unsigned int i = 0; i < strlen(chars); ++i){
		sz_sequence.erase (std::remove(sz_sequence.begin(), sz_sequence.end(), chars[i]), sz_sequence.end());
	}
	return sz_sequence;
}


bool CReadFiles::is_phylip_interleaved(){

	std::string sz_line;
	int n_seqs = 0, n_length_seq = 0;
	bool b_is_start_read = false, b_empty_line = false;
	int n_total_lines = 0;
	ifstream filein;
	filein.open(sz_file_name.c_str());
	while (std::getline(filein, sz_line)) {

		sz_line = trim(sz_line);
		if (sz_line.empty() && !b_is_start_read) continue;

		if (n_seqs == 0 && n_length_seq == 0){
			int scan_value = std::sscanf( sz_line.c_str(), "%d%d", &n_seqs, &n_length_seq);
			if (scan_value != 2){
				n_seqs = 0;
				n_length_seq = 0;
			}
			continue;
		}
		if (!b_is_start_read){
			if (sz_phylip_accept_chars.find(sz_line[sz_line.length() - 1]) != string::npos){
				b_is_start_read = true;
				n_total_lines += 1;
				continue;
			}
			else continue;
		}
		if (b_is_start_read && sz_line.empty()){
			b_empty_line = true;
			continue;
		}
		else if (b_is_start_read && !b_empty_line){
			n_total_lines += 1;
		}
		else if (b_is_start_read && b_empty_line){
			return true;
		}
	}
	return false;
}


bool CReadFiles::is_phylip_line_sequential(int &n_total_lines){

	std::string sz_line;
	int n_seqs = 0, n_length_seq = 0;
	bool b_is_start_read = false;
	n_total_lines = 0;
	ifstream filein;
	filein.open(sz_file_name.c_str());
	while (std::getline(filein, sz_line)) {

		sz_line = trim(sz_line);
		if (sz_line.empty() && !b_is_start_read) continue;

		if (n_seqs == 0 && n_length_seq == 0){
			int scan_value = std::sscanf( sz_line.c_str(), "%d%d", &n_seqs, &n_length_seq);
			if (scan_value != 2){
				n_seqs = 0;
				n_length_seq = 0;
			}
			continue;
		}
		if (!b_is_start_read){
			if (sz_phylip_accept_chars.find(sz_line[sz_line.length() - 1]) != string::npos){
				b_is_start_read = true;
				n_total_lines += 1;
				continue;
			}
			else continue;
		}
		if (b_is_start_read && !sz_line.empty()) n_total_lines += 1;
		if (sz_line.empty() && n_total_lines == n_seqs) return false;
	}

	if (n_seqs > 0 && n_length_seq > 0 && (n_total_lines % n_seqs) == 0 && n_total_lines > n_seqs) return true;
	return false;
}



int CReadFiles::read_fasta()
{
	#ifndef NP_Implementation
	std::string sz_line = "";
	#else
	std::string sz_line;
	#endif
	int n_count_sequence = -1;

	ifstream filein;
	filein.open(sz_file_name.c_str());
	// line by line
	while (std::getline(filein, sz_line)) {

		sz_line = trim(sz_line);
		if (b_debug) cout << sz_line << endl;
		#ifndef NP_Implementation
		if ((int) sz_line.length() == 0) continue;
		#else
		if ((int) sz_line.length() == 0) { continue; };
		#endif
		if (sz_line[0] == '>'){
			lst_names_seq.push_back(sz_line.substr(1, sz_line.size() - 1));
			n_count_sequence += 1;
		}
		else{
			#ifdef NP_Implementation
			if (n_count_sequence == -1) { continue; }
			#endif
			if (n_count_sequence == (int) lst_sequences.size()) lst_sequences.push_back(sz_line);
			else lst_sequences[n_count_sequence] += sz_line;
		}
	}
	filein.close();
	#ifndef NP_Implementation
	return 0;
	#endif
}

/// get the dimensions of nexus files...
// dimensions ntax=10 nchar=705;
void CReadFiles::get_dimensions_nexus_format(std::string sz_line, int &n_seqs, int &n_length_seq){

	n_seqs = 0;
	n_length_seq = 0;
	int scan_value = std::sscanf( sz_line.c_str(), "dimensions ntax=%d nchar=%d;", &n_seqs, &n_length_seq);
	if (scan_value == 0){
		std::sscanf( sz_line.c_str(), "dimensions nchar=%d ntax=%d;", &n_seqs, &n_length_seq);
	}
}

/// read nexus format
int CReadFiles::read_nexus()
{

	//// get the number of sequences and length of sequences...
	int n_seqs = 0, n_length_seq = 0;
	std::string sz_line = "", error;
	std::vector< std::string > lst_strings;
	bool b_is_start_read = false;

	std::string sz_dimensions = "dimensions";
	std::string sz_matix = "matrix";
	ifstream filein;
	filein.open(sz_file_name.c_str());
	bool b_is_name_seq = false;
	int n_count_sequence = 0;

	while (std::getline(filein, sz_line)) {

		sz_line = trim(sz_line);
		if (sz_line.empty()) continue;

		if (sz_line.find(sz_dimensions) != string::npos){
			get_dimensions_nexus_format(sz_line, n_seqs, n_length_seq);
			if (n_seqs == 0 || n_length_seq == 0){
				error = "Wans't possible to get the dimensions of the matrix.\nLine: " + sz_line;
				#ifndef NP_Implementation
				return exit_error(1 , error);
				#else
				exit_error(1, error);
				#endif
			}
		}
		if (sz_line.find(sz_matix) != string::npos){
			b_is_start_read = true;
			continue;
		}

		if (b_is_start_read){
			if (sz_line.find(";") != string::npos){
				b_is_start_read = false;	// end of reading
				continue;
			}

			lst_strings = split(sz_line, ' ');
			if ((int) lst_strings.size() == 0) continue; /// nothing to do

			if (! b_is_name_seq){
				if ((int) lst_strings.size() > 1){
					lst_names_seq.push_back(get_string_from_list(lst_strings, false));
					lst_sequences.push_back(get_string_from_list(lst_strings, true));
				}
				if ((int) lst_names_seq.size() == n_seqs) b_is_name_seq = true;
			}
			else{
				#ifndef NP_Implementation
				if (n_count_sequence >= (int) lst_sequences.size()) return exit_error(1 , "Some problem reading the file.");
				#else
				if (n_count_sequence >= (int) lst_sequences.size()) exit_error(1 , "Some problem reading the file.");
				#endif
				lst_sequences[n_count_sequence] += get_string_from_list(lst_strings, true);
				//cout << lst_sequences[n_count_sequence].length() << endl;
				n_count_sequence ++;
				if (n_count_sequence == n_seqs) n_count_sequence = 0;
			}
		}
	}


	if (n_seqs == 0){
		error = "Some problem reading the file. Please, check the file format.\nMore information in http://evolution.genetics.washington.edu/phylip/doc/sequence.html";
		filein.close();
		#ifndef NP_Implementation
		return exit_error(1 , error);
		#else
		exit_error(1, error);
		#endif
	}
	if ((int) lst_names_seq.size() != (int) lst_sequences.size()){
		error = "The number of sequences read are different from the number of name samples.";
		#ifndef NP_Implementation
		return exit_error(1 , error);
		#else
		exit_error(1, error);
		#endif
	}

	/// clean the numbers and spaces...
	if ((int) lst_names_seq.size() != n_seqs){
		error = "The file has a different number of sequences.\nRead: " + to_string_((int) lst_names_seq.size()) + "\nIn the header: " + to_string_(n_seqs);
		#ifndef NP_Implementation
		return exit_error(1 , error);
		#else
		exit_error(1, error);
		#endif
	}
	for(int i = 0; i < (int) lst_sequences.size(); i++){
		if ((int) lst_sequences[i].size() != n_length_seq){
			error = "This sequence " + lst_names_seq[i] + " has a different length " + to_string_((int) lst_sequences[i].size()) + " from the one read in the header: " + to_string_(n_length_seq);
			#ifndef NP_Implementation
		return exit_error(1 , error);
		#else
		exit_error(1, error);
		#endif
		}
	}
	filein.close();
	#ifndef NP_Implementation
	return 0;
	#endif
}


int CReadFiles::clean_data(){
	std::string error;

	/// test if the file exist
	if (!is_file_exist(sz_file_name)){
		error = "File doesn't exist: " + sz_file_name;
		#ifndef NP_Implementation
		return exit_error(1 , error);
		#else
		exit_error(1, error);
		#endif
	}

	ifstream filein;
	filein.open(sz_file_name.c_str());
	if (!filein.good()){
		error = "Failed to open alignment file: " + sz_file_name;
		#ifndef NP_Implementation
		return exit_error(1 , error);
		#else
		exit_error(1, error);
		#endif
	}
	filein.close();

	lst_sequences.clear();
	lst_names_seq.clear();
	n_max_length_name_seq = 0;

	/// get the extension from the file name
#ifdef WINDOWS_KEY_WORD
	if (sz_file_name.find_last_of("\\") != string::npos){
		sz_only_file_name = sz_file_name.substr(sz_file_name.find_last_of("\\") + 1, sz_file_name.length() - sz_file_name.find_last_of("\\"));
	}
#else
	// linux and mac
	if (sz_file_name.find_last_of("/") != string::npos){
		sz_only_file_name = sz_file_name.substr(sz_file_name.find_last_of("/") + 1, sz_file_name.length() - sz_file_name.find_last_of("/"));
	}
#endif
	else sz_only_file_name = sz_file_name;


	if (sz_only_file_name.find_last_of(".") != string::npos)
		sz_extension = sz_only_file_name.substr(sz_only_file_name.find_last_of(".") + 1, sz_only_file_name.length() - sz_only_file_name.find_last_of("."));
	else sz_extension = "";
	#ifndef NP_Implementation
	return 0;
	#endif
}


int CReadFiles::read_file(std::string sz_file_name_temp, int n_file_type){

	std::string error;
	#ifndef NP_Implementation
	int filetype = 0, n_error_code = EXIT_SUCCESS;
	#else
	int filetype = 0;
	#endif
	sz_file_name = sz_file_name_temp;

	//// setup initial data...
	#ifndef NP_Implementation
	n_error_code = clean_data();
	if (n_error_code != EXIT_SUCCESS) return n_error_code;

	if (n_file_type == CReadFiles::FORMAT_CLUSTAL_) {           // Clustal format
		filetype=1;
		n_error_code = read_clustal(filetype);
		if (n_error_code != EXIT_SUCCESS) return n_error_code;
	}
	else if (n_file_type == CReadFiles::FORMAT_MSF_) {          // MSF format
		filetype=2;
		n_error_code = read_clustal(filetype);
		if (n_error_code != EXIT_SUCCESS) return n_error_code;
	}
	else if (n_file_type == CReadFiles::FORMAT_PHYLIP_) {          // Phylip format
		n_error_code = read_phylip();
		if (n_error_code != EXIT_SUCCESS) return n_error_code;
	}
	else if (n_file_type == CReadFiles::FORMAT_FASTA_) {          // Fasta format
		n_error_code = read_fasta();
		if (n_error_code != EXIT_SUCCESS) return n_error_code;
	}

	else if (n_file_type == CReadFiles::FORMAT_NEXUS_) {          // nexus format
		n_error_code = read_nexus();
		if (n_error_code != EXIT_SUCCESS) return n_error_code;
	}
	else {
		error = "Unrecognized file..." + sz_file_name;
		return exit_error(1 , error);
	}
	#else
	clean_data();

	if (n_file_type == CReadFiles::FORMAT_CLUSTAL) {               // Clustal format
		filetype=1;
		read_clustal(filetype);
	}
	else if (n_file_type == CReadFiles::FORMAT_PHYLIP) {          // Phylip format
		read_phylip();
	}
	else if (n_file_type == CReadFiles::FORMAT_FASTA) {          // Fasta format
		read_fasta();
	}
	else if (n_file_type == CReadFiles::FORMAT_NEXUS) {          // nexus format
		read_nexus();
	}
	else {
		error = "Unrecognized file..." + sz_file_name;
		exit_error(1 , error);
	}
	#endif


	if (b_debug) cout << "testing sequences..." << endl;
	/// test if there is any sequence in the arrays, need to have more than one
	if ((int) lst_sequences.size() < 2){
		if ((int) lst_sequences.size() == 0){
			error = "Zero sequences were read from the file: " + sz_file_name;
			#ifndef NP_Implementation
			return exit_error(1 , error);
			#else
			exit_error(1, error);
			#endif
		}
		error = "Only one sequence was read from the file: " + sz_file_name;
		#ifndef NP_Implementation
			return exit_error(1 , error);
			#else
			exit_error(1, error);
			#endif
	}

	/// test the length of the two arrays
	if (lst_sequences.size() != lst_names_seq.size()){
		error = "Something wrong with the file.\nThe number of sequences names are different from the number of sequences in the file: " + sz_file_name;
		#ifndef NP_Implementation
			return exit_error(1 , error);
			#else
			exit_error(1, error);
			#endif
	}

	/// get the max length of size name
	n_max_length_name_seq = 0;
	for (unsigned int i = 1; i < lst_names_seq.size(); i++){
		if ((int) lst_names_seq[i].length() > n_max_length_name_seq) n_max_length_name_seq = (int) lst_names_seq[i].length();
	}

	/// test the size of the sequences, must be equal for all of them
	unsigned int n_max_size = lst_sequences[0].length();
	if (b_debug) cout << lst_names_seq[0] << ": " << n_max_size << endl;
	for (unsigned int i = 1; i < lst_sequences.size(); i++){
		if (b_debug) cout << lst_names_seq[i] << ": " << lst_sequences[i].length() << endl;
		if (lst_sequences[i].length() != n_max_size){
			error = "Something wrong with the file.\nThe sequence lengths are different in the file: " + sz_file_name;
			#ifndef NP_Implementation
			return exit_error(1 , error);
			#else
			exit_error(1, error);
			#endif
		}
	}

	/// check the chars available to read
	unsigned int n_max_seq = get_length_sequences();
	for (unsigned int i = 0; i < lst_sequences.size(); i++){
		std::transform(lst_sequences[i].begin(), lst_sequences[i].end(), lst_sequences[i].begin(), ::toupper);	// make ti upper case
		for(unsigned int x = 0; x < n_max_seq; x++){
			if (sz_accept_chars.find(get_char_sequences(i, x)) == string::npos){
				std::string a = "This char is not allowed (";
				a += get_char_sequences(i, x);
				a += ")\nThe char is in this line: ";
				a += lst_names_seq[i].c_str();
				a += ": ";
				a += lst_sequences[i].c_str();
				a += "\n";
				#ifndef NP_Implementation
				return exit_error(1 , a.c_str());
				#else
				exit_error(1, a.c_str());
				#endif
			}
		}

	}

	if (b_debug) cout << "Read: " << lst_names_seq.size() << " sequences with size " << lst_sequences.size() << endl;
	#ifndef NP_Implementation
	return EXIT_SUCCESS;
	#endif
}

void CReadFiles::save_file(std::string file_name_out){

	//// save sequences....
	cout << "save sequences..." << endl;
	std::string sz_out;
	int n_line_size = 60;
	FILE *fileout = fopen(file_name_out.c_str(), "w");
	for (unsigned int i = 0; i < lst_sequences.size(); i++){
		sz_out = ">" + lst_names_seq[i] + "\n";
		fprintf(fileout, "%s", sz_out.c_str());
		/// print the sequences
		for (unsigned int x = 0; x < lst_sequences[i].size(); x += n_line_size){
			sz_out = lst_sequences[i].substr(x, n_line_size) + "\n";
			fprintf(fileout, "%s", sz_out.c_str());
		}
	}
	fclose(fileout);
}