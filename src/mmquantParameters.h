/*
mmquant, a multi-mapping quantification tool.
Copyright (C) 2016-2017 Matthias Zytnicki

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MMQUANTPARAMETERS_H
#define MMQUANTPARAMETERS_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <queue>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <thread> 
#include <mutex> 
#include <atomic> 
#include <functional> 
#include <utility>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <locale>
#include <libgen.h>
#include <zlib.h>

#ifdef MMSTANDALONE
#   define MMOUT  std::cout
#   define MMERR  std::cerr
#   define MMEXIT exit(EXIT_FAILURE)
#else
#   define MMOUT  Rcpp::Rcout
#   define MMERR  Rcpp::Rcerr
#   define MMEXIT Rcpp::stop("Halting now.")
#endif

static const char VERSION[] = "1.2";

enum class Strandedness {U, F, R, FF, FR, RF};
enum class ReadsFormat {unknown, sam, bam};

static inline void lower (std::string &s) {
	transform(s.begin(), s.end(), s.begin(), ::tolower);
}

inline bool strandF (bool strand, bool isFirst, bool isPaired) {
	if (isPaired) {
		MMERR << "Error! Strandedness 'F' should be used for single-end reads.\nExiting." << std::endl;
		MMEXIT;
	}
	return strand;
}
inline bool strandR (bool strand, bool isFirst, bool isPaired) {
	if (isPaired) {
		MMERR << "Error! Strandedness 'R' should be used for single-end reads.\nExiting." << std::endl;
		MMEXIT;
	}
	return ! strand;
}
inline bool strandFF (bool strand, bool isFirst, bool isPaired) {
	if (! isPaired) {
		MMERR << "Error! Strandedness 'FF' should be used for paired-end reads.\nExiting." << std::endl;
	    MMEXIT;
	}
	return strand;
}
inline bool strandFR (bool strand, bool isFirst, bool isPaired) {
	if (! isPaired) {
		MMERR << "Error! Strandedness 'FR' should be used for paired-end reads.\nExiting." << std::endl;
	    MMEXIT;
	}
	return (strand == isFirst);
}
inline bool strandRF (bool strand, bool isFirst, bool isPaired) {
	if (! isPaired) {
		MMERR << "Error! Strandedness 'RF' should be used for paired-end reads.\nExiting." << std::endl;
	    MMEXIT;
	}
	return (strand != isFirst);
}
inline bool strandU (bool strand, bool isFirst, bool isPaired) {
	return true;
}


class Gene;
class Read;
struct MmquantParameters;

bool geneInclusion (MmquantParameters &parameters, Gene &g, Read &r);
bool geneOverlapPc (MmquantParameters &parameters, Gene &g, Read &r);
bool geneOverlap (MmquantParameters &parameters, Gene &g, Read &r);

typedef unsigned long Position;
typedef bool(*StrandednessFunction)(bool, bool, bool);
typedef bool(*GeneOverlapFunction)(MmquantParameters &, Gene &, Read &);

static const Position UNKNOWN = std::numeric_limits<Position>::max();


//' @export MmquantParameters
struct MmquantParameters {
	std::vector <std::string> args;

	std::vector <bool>                 sortednesses;
	std::vector <Strandedness>         strandednesses;
	std::vector <StrandednessFunction> strandednessFunctions;
	std::vector <ReadsFormat>          formats;

	GeneOverlapFunction       geneOverlapFunction;
	std::string               gtfFileName;
	std::string               outputFileName;
	std::vector <std::string> readsFileNames;
	std::vector <std::string> names;
	//std::filebuf              fileBuffer;
	std::ostream             *outputFile;

	unsigned int nInputs;
	unsigned int countThreshold      {     0   };
	float        overlap             {    -1.0 };
	float        mergeThreshold      {     0.0 };
	unsigned int nThreads            {     1   };
	unsigned int nOverlapDifference  {    30   };
	float        pcOverlapDifference {     2.0 };
	Position     binSize             { 16384   };
	bool         featureCountStyle   { false };
	bool         allSorted           {  true };
	bool         printGeneName       { false };
	bool         progress            { false };
	bool         quiet               { false };
	bool         printStructure      { false };
	
	inline void printUsage () {
#ifdef MMSTANDALONE
		MMERR << "Usage: mmquant [options]\n";
		MMERR <<   "\tCompulsory options:\n";
		MMERR <<     "\t\t-a file: annotation file in GTF format\n";
		MMERR <<     "\t\t-r file1 [file2 ...]: reads in BAM/SAM format\n";
		MMERR << "\tMain options:\n";
		MMERR <<     "\t\t-o output: output file (default: stdout)\n";
		MMERR <<     "\t\t-n name1 name2...: short name for each of the reads files\n";
		MMERR <<     "\t\t-s strand: string (U, F, R, FR, RF, FF, defaut: U) (use several strand types if the library strategies differ)\n";
		MMERR <<     "\t\t-e sorted: string (Y if reads are position-sorted, N otherwise, defaut: Y) (use several times if reads are not consistently (un)sorted)\n";
		MMERR <<     "\t\t-f format (SAM or BAM): format of the read files (default: guess from file extension)\n";
		MMERR <<     "\t\t-l integer: overlap type (<0: read is included, <1: % overlap, otherwise: # nt, default: " << overlap << ")\n";
		MMERR << "\tAmbiguous reads options:\n";
		MMERR <<     "\t\t-c integer: count threshold (default: " << countThreshold << ")\n";
		MMERR <<     "\t\t-m float: merge threshold (default: " << mergeThreshold << ")\n";
		MMERR <<     "\t\t-d integer: number of overlapping bp between the best matches and the other matches (default: " << nOverlapDifference << ")\n";
		MMERR <<     "\t\t-D float: ratio of overlapping bp between the best matches and the other matches (default: " << pcOverlapDifference << ")\n";
		MMERR << "\tOutput options:\n";
		MMERR <<     "\t\t-g: print gene name instead of gene ID in the output file\n";
		MMERR <<     "\t\t-0: use featureCounts output style\n";
		MMERR <<     "\t\t-p: print progress\n";
		MMERR <<     "\t\t-t integer: # threads (default: " << nThreads << ")\n";
		MMERR <<     "\t\t-v: version" << std::endl;
		MMERR <<     "\t\t-h: this help" << std::endl;
#endif
	}
	
	void printState() {
		MMOUT << "GTF file: " << gtfFileName << "\n";
		MMOUT << "Read(s) file:";
		for (std::string &f: readsFileNames) {
			MMOUT << " " << f;
		}
		MMOUT << "\n";
		MMOUT << "Output file: " << outputFileName << " (" << outputFile << ")\n";
		MMOUT << "Overlap function: " << geneOverlapFunction << "\n";
	}

	void setGtfFileName(std::string &s) {
		gtfFileName = s;
	}

	void addReadsFileName(std::string &s) {
		readsFileNames.push_back(s);
	}

	void addName(std::string &s) {
		names.push_back(s);
	}

	void setOutputFileName(std::string &s) {
		outputFileName = s;
	}

	void setOverlap(float f) {
		overlap = f;
		if      (overlap < 0.0) geneOverlapFunction = geneInclusion;
		else if (overlap < 1.0) geneOverlapFunction = geneOverlapPc;
		else                    geneOverlapFunction = geneOverlap;
	}

	int addStrand(std::string &s) {
	    if      (s == "U")  { strandednesses.push_back(Strandedness::U);  strandednessFunctions.push_back(strandU); }
	    else if (s == "F")  { strandednesses.push_back(Strandedness::F);  strandednessFunctions.push_back(strandF); }
	    else if (s == "R")  { strandednesses.push_back(Strandedness::R);  strandednessFunctions.push_back(strandR); }
	    else if (s == "FR") { strandednesses.push_back(Strandedness::FR); strandednessFunctions.push_back(strandFR); }
	    else if (s == "FF") { strandednesses.push_back(Strandedness::FF); strandednessFunctions.push_back(strandFF); }
	    else if (s == "RF") { strandednesses.push_back(Strandedness::RF); strandednessFunctions.push_back(strandRF); }
	    else {
	        MMERR << "Do not understand strandedness " << s << "\n" << "Exiting." << std::endl;
	        printUsage();
	        return EXIT_FAILURE;
	    }
	    return 0;
	}

	int addSort(std::string &s) {
			if      (s == "Y")  { sortednesses.push_back(true); }
			else if (s == "N")  { sortednesses.push_back(false); }
			else {
				MMERR << "Do not understand sortedness " << s << "\n" << "Exiting." << std::endl;
				printUsage();
				return EXIT_FAILURE;
			}
			return 0;
	}

	void setCountThreshold(unsigned int u) {
		countThreshold = u;
	}

	void setMergeThreshold(float f) {
		mergeThreshold = f;
	}

	void setPrintGeneName(bool b) {
		printGeneName = b;
	}

	void setFeatureCountStyle(bool b) {
		featureCountStyle = b;
	}
	
	void setQuiet(bool b) {
		quiet = b;
	}

	void setProgress(bool b) {
		progress = b;
	}

	void setNThreads(int n) {
		nThreads = n;
	}

	int addFormat(std::string &s) {
			lower(s);
			if      (s == "sam")  { formats.push_back(ReadsFormat::sam); }
			else if (s == "bam")  { formats.push_back(ReadsFormat::bam); }
			else {
				MMERR << "Do not understand reads format " << s << "\n" << "Exiting." << std::endl;
				printUsage();
				return EXIT_FAILURE;
			}
			return 0;
	}

	void setNOverlapDifference(int n) {
		nOverlapDifference = n;
	}

	void setPcOverlapDifference(float f) {
		pcOverlapDifference = f;
	}
	
	int parse(std::vector < std::string > &a) {
		int code;
		size_t nArgs = a.size();
		args = a;
		if (nArgs == 1) {
			printUsage();
			return EXIT_FAILURE;
		}
		for (size_t i = 1; i < nArgs; i++) {
			std::string &s = args[i];
			if (! s.empty()) {
				if (s == "-a") {
					setGtfFileName(args[++i]);
				}
				else if (s == "-r") {
					for (++i; i < nArgs; ++i) {
						s = args[i];
						if (s[0] == '-') {--i; break;}
						else addReadsFileName(s);
					}
				}
				else if (s == "-n") {
					for (++i; i < nArgs; ++i) {
						s = args[i];
						if (s[0] == '-') {--i; break;}
						else addName(s);
					}
				}
				else if (s == "-o") {
					setOutputFileName(args[++i]);
				}
				else if (s == "-l") {
					setOverlap(stof(args[++i]));
				}
				else if (s == "-s") {
					for (++i; i < nArgs; ++i) {
						s = args[i];
						if (s.empty())        {--i; break;}
						else if (s[0] == '-') {--i; break;}
						else if ((code = addStrand(s)) != 0) {
							return EXIT_FAILURE;
						}
					}
				}
				else if (s == "-e") {
					for (++i; i < nArgs; ++i) {
						s = args[i];
						if (s.empty())        {--i; break;}
						else if (s[0] == '-') {--i; break;}
						else if ((code = addSort(s)) != 0) {
							return EXIT_FAILURE;
						}
					}
				}
				else if (s == "-c") {
					setCountThreshold(stoul(args[++i]));
				}
				else if (s == "-m") {
					setMergeThreshold(stof(args[++i]));
				}
				else if (s == "-g") {
					setPrintGeneName(true);
				}
				else if (s == "-O") {
					setFeatureCountStyle(true);
				}
				else if (s == "-p") {
					setProgress(true);
				}
				else if (s == "-t") {
					setNThreads(stoi(args[++i]));
				}
				else if (s == "-f") {
					for (++i; i < nArgs; ++i) {
						s = args[i];
						lower(s);
						if (s.empty())        {--i; break;}
						else if (s[0] == '-') {--i; break;}
						else if ((code = addFormat(s)) != 0) {
							return EXIT_FAILURE;
						}
					}
				}
				else if (s == "-d") {
					setNOverlapDifference(stoi(args[++i]));
				}
				else if (s == "-D") {
					setPcOverlapDifference(stof(args[++i]));
				}
				else if (s == "-v") {
					MMERR << "mmquant version " << VERSION << std::endl;
					return EXIT_FAILURE;
				}
				else if (s == "-h") {
					printUsage();
					return EXIT_FAILURE;
				}
				else {
					MMERR << "Error: wrong parameter '" << s << "'.\nExiting." << std::endl;
					printUsage();
					return EXIT_FAILURE;
				}
			}
		}
		return 0;
	}

	int check () {
		if (gtfFileName.empty()) {
			MMERR << "Missing input GTF file.\nExiting." << std::endl;
			printUsage();
			return EXIT_FAILURE;
		}
		if (readsFileNames.empty()) {
			MMERR << "Missing input BAM file.\nExiting." << std::endl;
			printUsage();
			return EXIT_FAILURE;
		}
		nInputs = readsFileNames.size();
		if (names.empty()) {
			for (std::string &fileName: readsFileNames) {
				std::string n = fileName;
				size_t p = n.find_last_of("/");
				if (p != std::string::npos) n = n.substr(p+1); 
				p = n.find_last_of(".");
				if (p != std::string::npos) n = n.substr(0, p); 
				names.push_back(n);
			}
		}
		else if (names.size() != nInputs) {
			MMERR << "Number of names is not equal to number of file names.\nExiting." << std::endl;
			printUsage();
			return EXIT_FAILURE;
		}
		if (strandednesses.size() == 0) {
			strandednesses        = std::vector <Strandedness>         (nInputs, Strandedness::U);
			strandednessFunctions = std::vector <StrandednessFunction> (nInputs, strandU);
		}
		else if (strandednesses.size() == 1) {
			if (nInputs != 1) {
				strandednesses        = std::vector <Strandedness>         (nInputs, strandednesses.front());
				strandednessFunctions = std::vector <StrandednessFunction> (nInputs, strandednessFunctions.front());
			}
		}
		else if (strandednesses.size() != nInputs) {
			MMERR << "Number of strandedness is not equal to number of file names.\nExiting." << std::endl;
			printUsage();
			return EXIT_FAILURE;
		}
		if (sortednesses.size() == 0) {
			sortednesses = std::vector <bool> (nInputs, true);
		}
		else if (sortednesses.size() == 1) {
			if (nInputs != 1) {
				sortednesses = std::vector <bool> (nInputs, sortednesses.front());
			}
		}
		else if (sortednesses.size() != nInputs) {
			MMERR << "Number of sortedness is not equal to number of file names.\nExiting." << std::endl;
			printUsage();
			return EXIT_FAILURE;
		}
		allSorted = all_of(sortednesses.begin(), sortednesses.end(), [](bool v) {return v;});
		if (formats.size() == 0) {
			formats = std::vector <ReadsFormat> (nInputs, ReadsFormat::unknown);
		}
		else if (formats.size() == 1) {
			if (nInputs != 1) {
				formats = std::vector <ReadsFormat> (nInputs, formats.front());
			}
		}
		else if (formats.size() != nInputs) {
			MMERR << "Number of reads formats is not equal to number of file names.\nExiting." << std::endl;
			printUsage();
			return EXIT_FAILURE;
		}
		if (outputFileName.empty()) {
			outputFile = new std::ostream(MMOUT.rdbuf());
		}
		else {
		    std::filebuf *fileBuffer = new std::filebuf;
			fileBuffer->open(outputFileName.c_str(), std::ios::out);
			outputFile = new std::ostream(fileBuffer);
			//fileBuffer.open(outputFileName.c_str(), std::ios::out);
			//outputFile = new std::ostream(&fileBuffer);
		}
		return EXIT_SUCCESS;
	}

	std::ostream &getOS () {
		return *outputFile;
	}
};

#endif