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

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <queue>
#include <unordered_map>
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

static const char BAM_CIGAR_LOOKUP[] = "MIDNSHP=X";

enum class Strandedness {U, F, R, FF, FR, RF};
enum class ReadsFormat {unknown, sam, bam};

static inline void lower (std::string &s) {
	transform(s.begin(), s.end(), s.begin(), ::tolower);
}

namespace std {
	template <>
		struct hash<std::vector<unsigned int>> {
			size_t operator() (const std::vector<unsigned int> &vc) const {
				size_t seed = 0;
				for(auto& i : vc) seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
				return seed;
			}
		};
}

class comma_numpunct : public std::numpunct<char> {
	protected:
		virtual char   do_thousands_sep () const { return ','; }
		virtual std::string do_grouping ()      const { return "\03"; }
};


static void split(std::string &inString, char delim, std::vector <std::string> &outStrings) {
	std::stringstream ss(inString);
	std::string item;
	while (getline(ss, item, delim)) {
		outStrings.push_back(item);
	}
}

static void join(std::vector <std::string> &inStrings, std::string &outString, const char *delim) {
	typename std::vector <std::string>::iterator it = inStrings.begin();
	std::stringstream ss;
	ss << *it;
	for (++it; it != inStrings.end(); it++) ss << delim << *it;
	outString += ss.str();
}

static inline void ltrim(std::string &s) {
	s.erase(s.begin(), find_if(s.begin(), s.end(), not1(std::ptr_fun<int, int>(isspace))));
}

static inline void rtrim(std::string &s) {
	s.erase(find_if(s.rbegin(), s.rend(), not1(std::ptr_fun<int, int>(isspace))).base(), s.end());
}

static inline void trim(std::string &s) {
	ltrim(s); rtrim(s);
}

static inline void printStats(unsigned int n, const char *s, unsigned int nHits) {
	unsigned int size = log10(static_cast<double>(nHits))+1;
	size += static_cast<unsigned int>(size / 3.0);
	MMERR << "\t" << s << std::setw(size) << n << " (" << std::setw(5) << std::fixed << std::setprecision(1) << (static_cast<float>(n)/nHits*100) << "%)\n";
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
struct Parameters;

bool geneInclusion (Parameters &parameters, Gene &g, Read &r);
bool geneOverlapPc (Parameters &parameters, Gene &g, Read &r);
bool geneOverlap (Parameters &parameters, Gene &g, Read &r);

typedef unsigned long Position;
typedef bool(*StrandednessFunction)(bool, bool, bool);
typedef bool(*GeneOverlapFunction)(Parameters &, Gene &, Read &);

static const Position UNKNOWN = std::numeric_limits<Position>::max();


/*
struct Globals {
	static std::vector < std::string > chromosomes;
};
std::vector < std::string > allChromosomes {};
*/
std::vector < std::string > allChromosomes {};
std::vector < std::pair < std::string, std::vector < unsigned int > > > outputTable {};

struct Parameters {
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
	std::filebuf              fileBuffer;
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
	std::mutex   printMutex;

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
			fileBuffer.open(outputFileName.c_str(), std::ios::out);
			outputFile = new std::ostream(&fileBuffer);
		}
		return EXIT_SUCCESS;
	}

	std::ostream &getOS () {
		return *outputFile;
	}
};


class GtfLineParser {
	protected:
		std::string type, chromosome, geneId, transcriptId, geneName;
		Position start, end;
		bool strand;

	public:
		GtfLineParser (std::string line) {
			std::vector <std::string> splittedLine;
			split(line, '\t', splittedLine);
			assert(splittedLine.size() == 9);
			type       = splittedLine[2];
			chromosome = splittedLine[0];
			strand     = (splittedLine[6] == "+");
			start      = stoul(splittedLine[3]);
			end        = stoul(splittedLine[4]);
			std::string remaining = splittedLine[8];
			trim(remaining);
			while (! remaining.empty()) {
				size_t splitPosSpace = remaining.find(" "), splitPosEq = remaining.find("="), splitPos, endVal, endTag;
				if      (splitPosEq    == std::string::npos) splitPos = splitPosSpace;
				else if (splitPosSpace == std::string::npos) splitPos = splitPosEq;
				else                                         splitPos = std::min<size_t>(splitPosSpace, splitPosEq);
				std::string tag = remaining.substr(0, splitPos), value;
				rtrim(tag);
				remaining = remaining.substr(splitPos+1);
				ltrim(remaining);
				if (remaining[0] == '"') {
					remaining = remaining.substr(1);
					endVal    = remaining.find("\"");
					value     = remaining.substr(0, endVal);
					remaining = remaining.substr(endVal+1);
				}
				else {
					endVal = remaining.find(";");
					if (endVal == std::string::npos) endVal = remaining.size();
					value  = remaining.substr(0, endVal);
					rtrim(value);
				}
				if      (tag == "gene_id")       geneId       = value;
				else if (tag == "gene_name")     geneName     = value;
				else if (tag == "transcript_id") transcriptId = value;
				endTag = remaining.find(";");
				if (endTag == std::string::npos) {
					remaining.clear();
				}
				else {
					remaining = remaining.substr(endTag+1);
					ltrim(remaining);
				}
			}
		}
		std::string &getType ()               { return type; }
		std::string &getChromosome ()         { return chromosome; }
		bool         getStrand ()       const { return strand; }
		Position     getStart ()        const { return start; }
		Position     getEnd ()          const { return end; }
		std::string  getGeneId ()       const { return geneId; }
		std::string  getTranscriptId () const { return transcriptId; }
		std::string  getGeneName ()     const { return geneName; }
};

class XamRecord {
	protected:
		std::string name, chromosome;
		unsigned int flags, nHits;
		Position start;
		std::vector <std::pair <char, int>> cigar;
		size_t size;
		bool over;

	public:
		XamRecord (): start(UNKNOWN), over(false) { }
		void setChromosome (const std::string &c) { chromosome = c; }
		void setStart (Position s) { start = s; }
		void setName (const std::string &n) {
			name = n;
			if (isPaired()) {
				size_t pos = name.rfind('_');
				if ((pos != std::string::npos) && (pos < name.size()-1) && ((name[pos+1] == '1') || (name[pos+1] == '2'))) {
					name.resize(pos);
				}
			}
		}
		void setSize  (size_t s)                                         { size = s; }
		void setFlags (unsigned int f)                                   { flags = f; }
		void setCigar (const std::vector < std::pair < char, int > > &g) { cigar = g; }
		void setNHits (unsigned int n)                                   { nHits = n; }
		void setOver  ()                                                 { over = true; }
		std::string                         &getChromosome()             { return chromosome; }
		std::string                         &getName ()                  { return name; }
		std::vector <std::pair <char, int>> &getCigar ()                 { return cigar; }
		Position                             getStart()          const   { return start; }
		size_t                               getSize ()          const   { return size; }
		bool                                 isMapped ()         const   { return ((flags & 0x4) == 0);}
		bool                                 isProperlyMapped () const   { return ((flags & 0x2) == 0x2);}
		bool                                 getStrand ()        const   { return ((flags & 0x10) == 0);}
		bool                                 isPaired ()         const   { return ((flags & 0x1) == 0x1);}
		bool                                 isFirst ()          const   { return ((flags & 0x40) == 0x40);}
		bool                                 isOver ()           const   { return over; }
		unsigned int                         getNHits ()         const   { return nHits; }
};

class Reader {
	protected:
		std::ifstream    file;
		XamRecord   record;
		Parameters &parameters;

	public:
		Reader (Parameters &p, std::string &fileName): file(fileName.c_str()), parameters(p) {
			if (! file.good()) {
				MMERR << "Error, file '" << fileName << "' does not exists!" << std::endl;
				MMEXIT;
			}
		}
		virtual ~Reader () {}
		virtual XamRecord &getRecord() {
			getNextRecord();
			return record;
		}
		virtual void getNextRecord() = 0;
};

class SamReader: public Reader {
	public:
		SamReader (Parameters &p, std::string &fileName): Reader(p, fileName) {
			std::lock_guard<std::mutex> lock(parameters.printMutex);
			if (! parameters.quiet) MMERR << "Reading SAM file " << fileName << std::endl;
		}
		virtual void getNextRecord () {
			std::string line;
			std::vector <std::string> splittedLine;
			std::vector <std::pair <char, int>> cigar;
			do {
				if (! getline(file, line)) {
					record.setOver();
					return;
				}
			}
			while ((line.empty()) || (line[0] == '@') || (line[0] == '#'));
			split(line, '\t', splittedLine);
			assert(splittedLine.size() >= 12);
			record.setFlags(stoi(splittedLine[1]));
			record.setChromosome(splittedLine[2]);
			record.setStart(stoul(splittedLine[3]));
			record.setName(splittedLine[0]);
			record.setSize(splittedLine[9].size());
			int value = 0;
			for (char c: splittedLine[5]) {
				if ((c >= '0') && (c <= '9')) {
					value *= 10;
					value += (c - '0');
				}
				else {
					cigar.push_back(std::make_pair(c, value));
					value = 0;
				}
			}
			record.setCigar(cigar);
			std::string key;
			for (unsigned int i = 11; i < splittedLine.size(); i++) {
				std::string &part = splittedLine[i];
				size_t pos   = part.find(':');
				key = part.substr(0, pos);
				if (key == "NH") {
					record.setNHits(stoul(part.substr(part.find(':', pos+1)+1)));
					return;
				}
			}
			record.setNHits(1);
		}
};

class BamReader: public Reader {
	protected:
		std::vector <std::string> chromosomes;
		gzFile file;
	public:
		BamReader (Parameters &p, std::string &fileName): Reader(p, fileName) {
			std::lock_guard<std::mutex> lock(parameters.printMutex);
			if (! parameters.quiet) MMERR << "Reading BAM file " << fileName << std::endl;
			char buffer[1000000];
			int32_t tlen, nChrs, size;
			file = gzopen(fileName.c_str(), "rb");
			if (! file) {
				if (! parameters.quiet) MMERR << "Cannot open file '" << fileName << "'." << std::endl;
				return;
			}
			gzread(file, buffer, 4);
			buffer[4] = 0;
			gzread(file, reinterpret_cast<char*>(&tlen), 4);
			gzread(file, buffer, tlen); // text
			gzread(file, reinterpret_cast<char*>(&nChrs), 4);
			for (int i = 0; i < nChrs; i++) {
				gzread(file, reinterpret_cast<char*>(&size), 4);
				gzread(file, buffer, size);
				chromosomes.push_back(buffer);
				gzread(file, reinterpret_cast<char*>(&buffer), 4);
			}
			chromosomes.push_back("*");
		}
		virtual void getNextRecord () {
			if (gzeof(file)) {
				record.setOver();
				gzclose(file);
				return;
			}
			char     buffer[10000], v_c;
			int8_t   v_8;
			uint8_t  v_u8;
			int16_t  v_16;
			uint16_t v_u16;
			int32_t  v_32, size, chrId, pos, lReadName, lSeq, flags, nCigar;
			uint32_t v_u32, binMqNl, flagNc;
			float    v_f;
			std::vector <std::pair <char, int>> cigar;
			gzread(file, reinterpret_cast<char*>(&size), 4);
			if (gzeof(file)) {
				record.setOver();
				gzclose(file);
				return;
			}
			gzread(file, reinterpret_cast<char*>(&chrId), 4);
			if (chrId == -1) record.setChromosome(chromosomes.back());
			else record.setChromosome(chromosomes[chrId]);
			gzread(file, reinterpret_cast<char*>(&pos), 4);
			record.setStart(++pos);
			gzread(file, reinterpret_cast<char*>(&binMqNl), 4);
			lReadName = binMqNl & 0xff;
			gzread(file, reinterpret_cast<char*>(&flagNc), 4);
			flags  = flagNc >> 16;
			nCigar = flagNc & 0xffff;
			record.setFlags(flags);
			gzread(file, reinterpret_cast<char*>(&lSeq), 4);
			record.setSize(lSeq);
			gzread(file, buffer, 4);
			gzread(file, buffer, 4);
			gzread(file, buffer, 4);
			gzread(file, buffer, lReadName);
			record.setName(buffer);
			cigar.reserve(nCigar);
			for (int i = 0; i < nCigar; i++) {
				gzread(file, reinterpret_cast<char*>(&v_u32), 4);
				uint32_t s = v_u32 >> 4;
				char op = BAM_CIGAR_LOOKUP[v_u32 & ((1 << 4) - 1)];
				cigar.push_back(std::make_pair(op, s));
			}
			record.setCigar(cigar);
			gzread(file, buffer, (lSeq+1)/2);
			gzread(file, buffer, lSeq);
			record.setNHits(1);
			std::string key(2, 0);
			char c;
			for (int32_t i = 33+lReadName+4*nCigar+(lSeq+1)/2+lSeq; i < size; ) {
				for (unsigned int j = 0; j < 2; j++) {
					gzread(file, &c, 1);
					key[j] = c;
				}
				gzread(file, &c, 1);
				i += 3;
				int8_t n = 1;
				switch(c) {
					case 'H':
						gzread(file, reinterpret_cast<char*>(&n), 1);
						c = 'C';
						i++;
						break;
					case 'B':
						int8_t s = 0, m = 1;
						gzread(file, &c, 1);
						n = 0;
						for (unsigned int j = 0; j < 4; j++) {
							gzread(file, reinterpret_cast<char*>(&s), 1);
							n += s * m;
							m *= 16;
						}
						i += 5;
						break;
				}
				for (int j = 0; j < n; j++) {
					switch(c) {
						case 'A':
							gzread(file, reinterpret_cast<char*>(&v_c), 1);
							i++;
							break;
						case 'c':
							gzread(file, reinterpret_cast<char*>(&v_8), 1);
							i++;
							break;
						case 'C':
							gzread(file, reinterpret_cast<char*>(&v_u8), 1);
							i++;
							break;
						case 's':
							gzread(file, reinterpret_cast<char*>(&v_16), 2);
							i += 2;
							break;
						case 'S':
							gzread(file, reinterpret_cast<char*>(&v_u16), 2);
							i += 2;
							break;
						case 'i':
							gzread(file, reinterpret_cast<char*>(&v_32), 4);
							i += 4;
							break;
						case 'I':
							gzread(file, reinterpret_cast<char*>(&v_u32), 4);
							i += 4;
							break;
						case 'f':
							gzread(file, reinterpret_cast<char*>(&v_f), 4);
							i += 4;
							break;
						case 'Z':
							while ((v_c = gzgetc(file)) != 0) {
								//ss.put(v_c);
								i++;
							}
							i++;
							break;
						default:
							MMERR << "Problem with tag type '" << c << "'" << std::endl;
							return;
					}
				}
				if (key == "NH") record.setNHits(v_u8);
			}
		}
};

class Interval {
	protected:
		Position start, end;
	public:
		Interval (): start(UNKNOWN), end(UNKNOWN) {}
		Interval (Position s, Position e): start(s), end(e) { }
			Interval (GtfLineParser &line): Interval(line.getStart(), line.getEnd()) {}
			Position getStart () const { return start; }
			Position getEnd   () const { return end; }
			Position getSize  () const { return end - start + 1; }
			size_t overlaps (Interval &i) {
				Position s = std::max<Position>(start, i.start), e = std::min<Position>(end, i.end);
				if (s >= e) return 0;
				return e - s;
			}
			bool includes (const Interval &i) const {
				return ((i.start >= start) && (i.end <= end));
			}
			bool isBefore (const Interval &i) const {
				return (end < i.start);
			}
			bool isAfter (const Interval &i) const {
				return (i.isBefore(*this));
			}
			void merge (const Interval &i) {
				start = std::min<Position>(start, i.start);
				end   = std::max<Position>(end, i.end);
			}
			bool operator< (const Interval &i) { return (this->isBefore(i)); }
};

class Transcript: public Interval {
	protected:
		std::string            name;
		std::vector <Interval> exons;
		std::vector <Interval> introns;
	public:
		Transcript () {}
		Transcript (std::string n, Position s, Position e = UNKNOWN): Interval(s, e), name(n) {}
		Transcript (Interval e): Interval(e), name("unnamed_transcript") {}
		Transcript (GtfLineParser &line): Transcript(line.getTranscriptId(), line.getStart(), line.getEnd()) {}
		std::string &getName () { return name; }
		std::vector <Interval> &getExons () { return exons; }
		void addExon (Interval &e) {
			exons.push_back(e);
			start = std::min<Position>(start, e.getStart());
			end   = std::max<Position>(end,   e.getEnd());
		}
		void addExon (Position s, Position e) {
			exons.push_back(Interval(s, e));
		}
		void checkStructure () {
			sort(exons.begin(), exons.end());
			if (exons.empty()) exons.push_back(*this);
			for (size_t i = 1; i < exons.size(); i++) {
				introns.push_back(Interval(exons[i-1].getEnd()+1, exons[i].getStart()-1));
			}
		}
		size_t overlaps (Transcript &t, bool isIntrons = false) {
			if (Interval::overlaps(t) == 0) return 0;
			std::vector <Interval> &theseElement = ((isIntrons)? introns: exons), &thoseElement = ((isIntrons)? t.introns: t.exons);
			size_t o = 0;
			for (Interval &i1: theseElement) {
				for (Interval &i2: thoseElement) {
					o += i1.overlaps(i2);
				}
			}
			return o;
		}
		bool includes (Transcript &t) {
			if (! Interval::includes(t)) return false;
			for (Interval &e2: t.exons) {
				if (! any_of(exons.begin(), exons.end(), [&e2](Interval &e1){return e1.includes(e2);})) return false;
			}
			return true;
		}
};

class Read: public Transcript {
	protected:
		std::string chromosome;
		bool strand, paired, first;
		unsigned int nHits;
		size_t size;
	public:
		Read () {}
		Read (XamRecord &record, StrandednessFunction strandednessFunction): Transcript(record.getName(), record.getStart()), chromosome(record.getChromosome()), paired(record.isPaired()), nHits(record.getNHits()), size(0) {
			addPart(record);
			strand = strandednessFunction(record.getStrand(), record.isFirst(), record.isPaired());
		}
		void addPart(XamRecord &record) {
			Position s = record.getStart(), e = record.getStart();
			start      = std::min<Position>(start, record.getStart());
			first      = record.isFirst();
			for (auto &part: record.getCigar()) {
				char c = part.first;
				int  v = part.second;
				switch (c) {
					case 'M':
					case 'D':
					case '=':
					case 'X':
						e += v;
						break;
					case 'N':
						if (s != e) addExon(s, e-1);
						introns.push_back(Interval(e, e+v-1));
						e += v;
						s = e;
						break;
					case 'I':
					case 'S':
					case 'H':
					case 'P':
						break;
					default:
						MMERR << "Problem in the cigar: do not understand char " << c << std::endl;
				}
			}
			if (s != e) addExon(s, e-1);
			--e;
			if ((end == UNKNOWN) || (e > end)) end = e;
			nHits = std::min<unsigned int>(nHits, record.getNHits());
			size += record.getSize();
		}
		std::string &getChromosome()       { return chromosome; }
		bool         getStrand()     const { return strand; }
		unsigned int getNHits ()     const { return nHits; }
		size_t       getSize ()      const { return size; }
		bool         isPaired ()     const { return paired; }
		bool         hasFirst ()     const { return first; }
};

class Gene: public Interval {
	protected:
		std::string id, name;
		bool strand;
		std::vector <Transcript> transcripts;
		unsigned int chromosomeId;
	public:
		Gene (std::string i, std::string n, Position s, Position e, bool st, unsigned int c): Interval(s, e), id(i), name(n), strand(st), chromosomeId(c) {
			if (id.empty()) id = name;
			if (name.empty()) name = id;
		}		
		Gene (GtfLineParser &line, unsigned int c): Gene(line.getGeneId(), line.getGeneName(), line.getStart(), line.getEnd(), line.getStrand(), c) {
			if (name.empty()) name = id;
			if (id.empty()) id = name;
		}
		std::string      &getName ()               { return name; }
		std::string      &getId ()                 { return id; }
		bool         getStrand ()       const { return strand; }
		unsigned int getChromosomeId () const { return chromosomeId; }
		void correctName () { name += " (" + id + ")"; }
		void addTranscript (Transcript &t) {
			transcripts.push_back(t);
			start = std::min<Position>(start, t.getStart());
			end   = std::max<Position>(end,   t.getEnd());
		}
		void addTranscript (std::string n, Position s, Position e) {
			Transcript t(n, s, e);
			addTranscript(t);
		}
		void addExon (Interval &e, std::string &transcriptName) {
			auto pos = find_if(transcripts.begin(), transcripts.end(), [&transcriptName](Transcript &t){return (t.getName() == transcriptName);});
			if (pos == transcripts.end()) {
				Transcript t(transcriptName, e.getStart(), e.getEnd());
				t.addExon(e);
				transcripts.push_back(t);
			}
			else {
				pos->addExon(e);
			}
			start = std::min<Position>(start, e.getStart());
			end   = std::max<Position>(end,   e.getEnd());
		}
		void checkStructure () {
			if (transcripts.empty()) transcripts.push_back(Transcript(*this));
			for (Transcript &transcript: transcripts) transcript.checkStructure();
		}
		size_t overlaps (Read &read, bool isIntrons = false) {
			if (Interval::overlaps(read) == 0) return 0;
			size_t o = 0;
			for (Transcript &t: transcripts) {
				o = std::max<size_t>(o, t.overlaps(read, isIntrons));
			}
			return o;
		}
		bool includes (Read &read) {
			return (any_of(transcripts.begin(), transcripts.end(), [&read](Transcript &t){return t.Transcript::includes(read);}));
		}
		std::string getFeatureCountId () {
			std::vector <Interval> exons, mergedExons;
			for (Transcript &transcript: transcripts) {
				for (Interval &exon:transcript.getExons()) {
					exons.push_back(exon);
				}
			}
			sort(exons.begin(), exons.end());
			for (Interval &exon: exons) {
				if ((mergedExons.empty()) || (! exon.overlaps(mergedExons.back()))) {
					mergedExons.push_back(exon);
				}
				else {
					mergedExons.back().merge(exon);
				}
			}
			std::string output = getId() + "\t";
			std::vector <std::string> chrs, starts, ends, strands;
			Position size = 0;
			for (Interval &exon: mergedExons) {
				chrs.push_back(allChromosomes[chromosomeId]);
				starts.push_back(std::to_string(exon.getStart()));
				ends.push_back(std::to_string(exon.getEnd()));
				strands.push_back((strand)? "+": "-");
				size += exon.getSize();
			}
			join(chrs,    output, ";");
			output += "\t";
			join(starts,  output, ";");
			output += "\t";
			join(ends,    output, ";");
			output += "\t";
			join(strands, output, ";");
			output += "\t" + std::to_string(size);
			return output;
		}
		friend bool operator<(const Gene &g1, const Gene &g2) {
			return ((g1.getChromosomeId() < g2.getChromosomeId()) || ((g1.getChromosomeId() == g2.getChromosomeId()) && (g1.getStart() < g2.getStart())));
		}
};

inline bool geneInclusion (Parameters &parameters, Gene &g, Read &r) {
	return g.includes(r);
}
inline bool geneOverlapPc (Parameters &parameters, Gene &g, Read &r) {
	return (r.getSize() * parameters.overlap <= g.overlaps(r));
}
inline bool geneOverlap (Parameters &parameters, Gene &g, Read &r) {
	return (g.overlaps(r) >= parameters.overlap);
}

struct GeneListPosition {
	size_t chromosomeId, geneId;
	GeneListPosition (): chromosomeId(0), geneId(0) {}
	void reset () {
		chromosomeId = geneId = 0;
	}
};

class GeneList {
	protected:
		std::vector <std::string>                              unknownChromosomes;
		std::vector <Gene>                                     genes;
		std::vector <unsigned int>                             chrStarts;
		std::unordered_map <std::string, std::vector <size_t>> bins;
		Parameters                                            &parameters;

		void reduceOverlappingGeneList(Read &read, std::vector < unsigned int > &geneIdsIn, std::vector < unsigned int > &geneIdsOut, bool isIntron) {
			unsigned int maxOverlap = 0;
			std::vector < unsigned int > overlaps;
			for (unsigned int geneId: geneIdsIn) {
				unsigned int overlap = genes[geneId].overlaps(read, isIntron);
				maxOverlap = std::max<unsigned int>(maxOverlap, overlap);
				overlaps.push_back(overlap);
			}
			for (unsigned int i = 0; i < geneIdsIn.size(); i++) {
				if (((maxOverlap <= parameters.nOverlapDifference) || (overlaps[i] + parameters.nOverlapDifference >= maxOverlap)) || (maxOverlap <= overlaps[i] * parameters.pcOverlapDifference)) {
					geneIdsOut.push_back(geneIdsIn[i]);
				}
			}
		}
		void evaluateScan(Read &read, std::vector < unsigned int > &geneIds) {
			if (genes.size() <= 1) {
				return;
			}
			std::vector < unsigned int > keptGenes1, keptGenes2;
			reduceOverlappingGeneList(read, geneIds, keptGenes1, false);
			if (keptGenes1.size() == 1) {
				geneIds = keptGenes1;
				return;
			}
			reduceOverlappingGeneList(read, keptGenes1, keptGenes2, true);
			geneIds = keptGenes2;
		}

	public:
		GeneList(Parameters &p, std::string &fileName): parameters(p) {
			std::ifstream file (fileName.c_str());
			std::vector <std::unordered_map<std::string, unsigned int>> geneHash;
			std::vector <std::unordered_map<std::string, Gene>> unsortedGenes;
			std::string line, chromosome;
			unsigned long cpt;
			unsigned int chromosomeId = std::numeric_limits<unsigned int>::max();
			if (! parameters.quiet) MMERR << "Reading GTF file" << std::endl;
			for (cpt = 0; getline(file, line); cpt++) {
				if ((! line.empty()) && (line[0] != '#')) {
					GtfLineParser parsedLine(line);
					if (parsedLine.getChromosome() != chromosome) {
						geneHash.push_back(std::unordered_map<std::string, unsigned int>());
						chromosome = parsedLine.getChromosome();
						bool seen = false;
						for (unsigned int i = 0; (i < allChromosomes.size()) && (! seen); i++) {
							if (allChromosomes[i] == chromosome) {
								chromosomeId = i;
								seen         = true;
							}
						}
						if (! seen) {
							chromosomeId = allChromosomes.size();
							allChromosomes.push_back(chromosome);
						}
					}
					if (parsedLine.getType() == "gene") {
						Gene gene(parsedLine, chromosomeId);
						geneHash.back()[parsedLine.getGeneId()] = genes.size();
						genes.push_back(gene);
					}
					else if (parsedLine.getType() == "transcript") {
						std::string geneName = parsedLine.getGeneId();
						auto pos        = geneHash.back().find(geneName);
						Transcript t(parsedLine);
						if (pos == geneHash.back().end()) {
							Gene gene(parsedLine, chromosomeId);
							gene.addTranscript(t);
							geneHash.back()[parsedLine.getGeneId()] = genes.size();
							genes.push_back(gene);
						}
						else {
							genes[pos->second].addTranscript(t);
						}
					}
					else if (parsedLine.getType() == "exon") {
						std::string geneName = parsedLine.getGeneId();
						auto pos        = geneHash.back().find(geneName);
						Interval e(parsedLine);
						std::string transcriptName = parsedLine.getTranscriptId();
						if (pos == geneHash.back().end()) {
							Gene gene(parsedLine, chromosomeId);
							gene.addExon(e, transcriptName);
							geneHash.back()[parsedLine.getGeneId()] = genes.size();
							genes.push_back(gene);
						}
						else {
							genes[pos->second].addExon(e, transcriptName);
						}
					}
				}
				if (parameters.progress && (cpt % 100000 == 0)) MMERR << "\t" << cpt << " lines read.\r" << std::flush;
			}
			if (! parameters.quiet) MMERR << "\t" << cpt << " lines read, done.  " << genes.size() << " genes found." << std::endl;
			sort(genes.begin(), genes.end());
			chrStarts = std::vector<unsigned int>(allChromosomes.size());
			chromosomeId = std::numeric_limits<unsigned int>::max();
			for (unsigned int i = 0; i < genes.size(); i++) {
				genes[i].checkStructure();
				if (genes[i].getChromosomeId() != chromosomeId) {
					chrStarts[chromosomeId = genes[i].getChromosomeId()] = i;
				}
			}
			chromosomeId = 0;
			std::vector <std::string> names, duplicates;
			std::string previousName;
			for (Gene &gene: genes) {
				names.push_back(gene.getName());
			}
			sort(names.begin(), names.end());
			for (std::string &name: names) {
				if (name == previousName) {
					if ((duplicates.empty()) || (name != duplicates.back())) {
						duplicates.push_back(name);
					}
				}
				else {
					previousName = name;
				}
			}
			for (Gene &gene: genes) {
				if (binary_search(duplicates.begin(), duplicates.end(), gene.getName())) gene.correctName();
			}
			if (! parameters.allSorted) {
				for (unsigned int i = 0; i < genes.size(); i++) {
					std::vector<size_t> &chrBins    = bins[allChromosomes[genes[i].getChromosomeId()]];
					unsigned        bin        = genes[i].getEnd() / parameters.binSize;
					if (chrBins.size() <= bin) {
						chrBins.insert(chrBins.end(), bin-chrBins.size()+1, i);
					}
				}
			}
		}
		void scan(Read &read, std::vector <unsigned int> &matchingGenes, GeneListPosition &position, Strandedness strandedness, bool sorted) {
			if (sorted) {
				if (allChromosomes[position.chromosomeId] != read.getChromosome()) {
					if (find(unknownChromosomes.begin(), unknownChromosomes.end(), read.getChromosome()) != unknownChromosomes.end()) return;
					for (position.chromosomeId = 0; (position.chromosomeId < allChromosomes.size()) && (allChromosomes[position.chromosomeId] != read.getChromosome()); position.chromosomeId++) ;
					if (position.chromosomeId == allChromosomes.size()) {
						unknownChromosomes.push_back(read.getChromosome());
						position.reset();
						return;
					}
					position.geneId = chrStarts[position.chromosomeId];
				}
			}
			else {
				unsigned int bin = read.getStart() / parameters.binSize;
				auto         p   = bins.find(read.getChromosome());
				if (p == bins.end()) return;
				bin = std::min<unsigned int>(bin, p->second.size()-1);
				position.geneId       = p->second[bin];
				position.chromosomeId = genes[position.geneId].getChromosomeId();
			}
			while ((position.geneId < genes.size()) && (genes[position.geneId].getChromosomeId() == position.chromosomeId) && (genes[position.geneId].isBefore(read))) {
				position.geneId++;
			}
			size_t id = position.geneId;
			while ((id < genes.size()) && (genes[id].getChromosomeId() == position.chromosomeId) && (! genes[id].isAfter(read))) {
				Gene &gene = genes[id];
				if ((strandedness == Strandedness::U) || (read.getStrand() == gene.getStrand())) {
					bool m = parameters.geneOverlapFunction(parameters, gene, read);
					if (m) matchingGenes.push_back(id);
				}	
				id++;
			}
			evaluateScan(read, matchingGenes);
		}
		std::string getGeneName(unsigned int i) {
			if (parameters.featureCountStyle) return genes[i].getFeatureCountId();
			if (parameters.printGeneName)     return genes[i].getName();
			return genes[i].getId();
		}
};

class Counter {
	protected:
		GeneList &geneList;
		std::unordered_map<std::string, std::pair <unsigned int, std::vector <unsigned int>>> readCounts;
		std::unordered_map<std::vector<unsigned int>, unsigned int> geneCounts;
		std::vector<std::vector<unsigned int>> genes;
		unsigned int nHits, nUnique, nAmbiguous, nMultiple, nUnassigned;
		std::string fileName;
		Parameters &parameters;
		void addGeneCount (const std::vector <unsigned int> &g) {
			std::vector <unsigned int> s(g.begin(), g.end());
			sort(s.begin(), s.end());
			std::vector<unsigned int>::iterator it = unique(s.begin(), s.end());
			s.resize(distance(s.begin(), it)); 
			geneCounts[s]++;
		}
		void addCount(std::string &read, std::vector <unsigned int> &matchingGenes, unsigned int nHits) {
			if (matchingGenes.empty()) {
				nUnassigned++;
				return;
			}
			if      (matchingGenes.size() > 1) nAmbiguous++;
			else if (nHits == 1) nUnique++;
			if (nHits > 1) {
				nMultiple++;
				auto pos = readCounts.find(read);
				if (pos == readCounts.end()) {
					readCounts[read] = std::make_pair(nHits-1, matchingGenes);
				}
				else {
					pos->second.first--;
					pos->second.second.insert(pos->second.second.end(), matchingGenes.begin(), matchingGenes.end());
					if (pos->second.first == 0) {
						addGeneCount(pos->second.second);
						readCounts.erase(pos);
					}
				}
			}
			else {
				addGeneCount(matchingGenes);
			}
		}
	public:
		Counter (Parameters &p, GeneList &gl): geneList(gl), parameters(p) {}
		void clear () {
			readCounts.clear();
			geneCounts.clear();
			genes.clear();
			nHits = nUnique = nAmbiguous = nMultiple = nUnassigned = 0;
		}
		void read (std::string &f, Strandedness strandedness, StrandednessFunction strandednessFunction, bool sorted, ReadsFormat format) {
			Reader *reader;
			fileName = f;
			if      (format == ReadsFormat::bam) reader = new BamReader(parameters, fileName);
			else if (format == ReadsFormat::sam) reader = new SamReader(parameters, fileName);
			else {
				if (fileName.size() < 4) {
					MMERR << "Cannot deduce type from file name '" << fileName << "'.  Should be a .sam or .bam file.  Please specify it using the '-f' option.";
					MMEXIT;
				}
				std::string suffix = fileName.substr(fileName.size()-4);
				lower(suffix);
				if      (suffix == ".bam") reader = new BamReader(parameters, fileName);
				else if (suffix == ".sam") reader = new SamReader(parameters, fileName);
				else {
					MMERR << "Cannot deduce type from file name '" << fileName << "'.  Should be a .sam or .bam file.  Please specify it using the '-f' option.";
                    MMEXIT;					
				}
			}
			unsigned int cpt;
			GeneListPosition position;
			Position previousPos = 0;
			std::unordered_map < std::string, std::array < std::queue < Read >, 2 > > pendingReads;
			XamRecord &record = reader->getRecord();
			geneCounts.clear();
			for (cpt = 0; !record.isOver(); cpt++, reader->getNextRecord()) {
				if (record.isMapped()) {
					Read read;
					bool pending = true;
					if (record.isPaired()) {
						std::string readName = record.getName();
						size_t bucket   = record.isFirst()? 0: 1;
						auto   pos      = pendingReads.find(readName);
						if (pos != pendingReads.end()) {
							std::queue < Read > &reads = pos->second[1-bucket];
							if (! reads.empty()) {
								pending = false;
								read = reads.front();
								read.addPart(record);
								reads.pop();
								if ((reads.empty()) && (pos->second[bucket].empty())) pendingReads.erase(pos);
							}
						}
						if (pending) {
							read = Read(record, strandednessFunction);
							pendingReads[readName][bucket].push(read);
						}
					}
					else {
						pending = false;
						read = Read(record, strandednessFunction);
					}
					if (! pending) {
						nHits++;
						std::vector <unsigned int> matchingGenes;
						geneList.scan(read, matchingGenes, position, strandedness, sorted);
						addCount(read.getName(), matchingGenes, read.getNHits());
					}
					if (previousPos > record.getStart()) pendingReads.clear();
					previousPos = record.getStart();
				}
				if (parameters.progress && (parameters.nThreads == 1) && (cpt % 1000000 == 0)) MMERR << "\t" << cpt << " lines read.\r" << std::flush;
			}
			if (! parameters.quiet) MMERR << "\t" << cpt << " lines read, done." << std::endl;
			for (auto &e: readCounts) {
				addGeneCount(e.second.second);
			}
			pendingReads.clear();
			delete reader;
		}
		std::unordered_map<std::vector<unsigned int>, unsigned int> &getCounts () {
			return geneCounts;
		}
		void dump () {
		    if (! parameters.quiet) {
    			MMERR << "Results for " << fileName << ":" << std::endl;
    			if (nHits == 0) {
    				MMERR << "\tNo hit." << std::endl;
    			}
    			else {
    				MMERR << "\t# hits:                     " << nHits << "\n";
    				printStats(nUnique,     "# uniquely mapped reads:    ", nHits);
    				printStats(nAmbiguous,  "# ambiguous hits:           ", nHits);
    				printStats(nMultiple,   "# non-uniquely mapped hits: ", nHits);
    				printStats(nUnassigned, "# unassigned hits:          ", nHits);
    			}
		    }
		}
};
class TableCount {
	protected:
		GeneList &geneList;
		std::vector<std::vector<unsigned int>> genes;
		std::vector<std::vector<unsigned int>> table;
		std::vector<std::pair<std::string, std::vector<unsigned int>>> selectedTable;
		unsigned int nColumns;
		std::unordered_map<std::vector<unsigned int>, std::vector<unsigned int>> geneCounts;
		Parameters &parameters;
	public:
		TableCount(Parameters &p, GeneList &g): geneList(g), nColumns(0), parameters(p) {}
	    std::vector<std::pair<std::string, std::vector<unsigned int>>> &getTable () {
	        return selectedTable;
	    }
		void addCounter(Counter &counter) {
			auto &counts = counter.getCounts();
			for (auto &count: counts) {
				auto p = geneCounts.find(count.first);
				if (p == geneCounts.end()) {
					geneCounts[count.first] = std::vector <unsigned int> (parameters.nInputs, 0);
					std::vector <unsigned int> v (parameters.nInputs, 0);
					v[nColumns] = count.second;
					geneCounts[count.first] = v;
				}
				else {
					p->second[nColumns] = count.second;
				}
			}
			++nColumns;
		}
		void selectGenes() {
			unsigned int nGenes = geneCounts.size();
			genes.reserve(nGenes);
			table.reserve(nGenes);
			for (auto &e: geneCounts)
				genes.push_back(e.first);
			sort(genes.begin(), genes.end());
			for (std::vector<unsigned int> &gene: genes) {
				table.push_back(geneCounts[gene]);
			}
			unsigned int maxId = 0;
			for (std::vector<unsigned int> &gene: genes) {
				for (unsigned int g: gene) {
					maxId = std::max<unsigned int>(maxId, g);
				}
			}
			std::vector <std::vector <unsigned int>> groups (maxId+1);
			for (unsigned int i = 0; i < genes.size(); i++) {
				std::vector <unsigned int> &gene = genes[i];
				if (! table[i].empty()) {
					for (unsigned int g: gene) {
						groups[g].push_back(i);
					}
				}
			}
			std::vector <std::vector <unsigned int>> ancestors (genes.size());
			for (std::vector <unsigned int> &group: groups) {
				for (unsigned int i = 0; i < group.size(); i++) {
					for (unsigned int j = 0; j < group.size(); j++) {
						if (i != j) {
							if (includes(genes[group[i]].begin(), genes[group[i]].end(), genes[group[j]].begin(), genes[group[j]].end())) {
								ancestors[group[j]].push_back(group[i]);
							}
						}
					}
				}
			}
			for (unsigned int i = 0; i < ancestors.size(); i++) {
				sort(ancestors[i].begin(), ancestors[i].end());
				auto it = unique(ancestors[i].begin(), ancestors[i].end());
				ancestors[i].resize(distance(ancestors[i].begin(), it));
			}
			std::vector <std::vector <unsigned int>> parents (ancestors);
			for (unsigned int i = 0; i < ancestors.size(); i++) {
				std::vector <unsigned int> &theseAncestors = ancestors[i];
				std::vector <unsigned int> &theseParents   = parents[i];
				std::vector <unsigned int>::iterator parentsBegin = theseParents.begin(), parentsEnd = theseParents.end();
				for (unsigned int j: theseAncestors) {
					std::vector <unsigned int> &thoseAncestors = ancestors[j];
					parentsEnd = remove_if(parentsBegin, parentsEnd, [&thoseAncestors] (unsigned int k) {return (find(thoseAncestors.begin(), thoseAncestors.end(), k) != thoseAncestors.end());});
				}
				parents[i] = std::vector <unsigned int>(parentsBegin, parentsEnd);
			}
			bool changed;
			do {
				changed = false;
				for (unsigned int i = 0; i < genes.size(); i++) {
					if (! table[i].empty()) {
						bool goodParent;
						unsigned int j = i;
						do {
							goodParent = (parents[j].size() == 1);
							if (goodParent) j = parents[j].front();
						}
						while ((goodParent) && (table[j].empty()));
						if ((goodParent) && (! table[j].empty())) {
							bool underThreshold = true;
							for (unsigned int k = 0; (k < table[i].size()) && (underThreshold); k++) underThreshold = (static_cast<float>(table[i][k]) <= parameters.mergeThreshold * static_cast<float>(table[j][k]));
							if (underThreshold) {
								for (unsigned int k = 0; k < table[i].size(); k++) table[j][k] += table[i][k];
								table[i].clear();
								changed = true;
							}
						}
					}
				}
			} while (changed);
			for (unsigned int i = 0; i < table.size(); i++) {
				if (! table[i].empty()) {
					if (*max_element(table[i].begin(), table[i].end()) < parameters.countThreshold) table[i].clear();
				}
			}
		}
	    void prepareOutput() {
			for (unsigned int i = 0; i < genes.size(); i++) {
				if (! table[i].empty()) {
					std::stringstream ss;
					std::vector<unsigned int> &gene = genes[i];
					std::string s;
					std::vector<std::string> geneNames;
					geneNames.reserve(gene.size());
					for (unsigned int id: gene) {
						geneNames.push_back(geneList.getGeneName(id));
					}
					join(geneNames, s, "--");
					selectedTable.push_back(std::make_pair(s, table[i]));
				}
			}
			sort(selectedTable.begin(), selectedTable.end());
		}
	    void dump() {
			if (parameters.featureCountStyle) {
				parameters.getOS() << "# Program:mmquant v" << VERSION << "; Command:";
				for (std::string &arg: parameters.args) {
					parameters.getOS() << " \"" << arg << "\"";
				}
				parameters.getOS() << "\nGeneid\tChr\tStart\tEnd\tStrand\tLength";
			}
			else {
				parameters.getOS() << "Gene";
			}
			for (std::string &sample: parameters.names) {
				parameters.getOS() << "\t" << sample;
			}
			parameters.getOS() << "\n";
			for (auto &line: selectedTable) {
    			parameters.getOS() << line.first;
			    for (unsigned int i: line.second) {
			        parameters.getOS() << "\t" << i;
			    }
			    parameters.getOS() << "\n";
			}
		}
};

static void doWork (Parameters &parameters, GeneList &geneList, TableCount &table, std::atomic < unsigned int > &i, std::mutex &m1, std::mutex &m2) {
	Counter counter (parameters, geneList);
	while (i < parameters.nInputs) {
		unsigned int thisI;
		m1.lock();
		thisI = i;
		++i;
		m1.unlock();
		counter.clear();
		counter.read(parameters.readsFileNames[thisI], parameters.strandednesses[thisI], parameters.strandednessFunctions[thisI], parameters.sortednesses[thisI], parameters.formats[thisI]);
		m2.lock();
		counter.dump();
		m2.unlock();
		table.addCounter(counter);
	}
}

static int start (Parameters &parameters) {
    int code = parameters.check();
    if (code != 0) return code;
    if (! parameters.quiet) parameters.printState();
    std::locale comma_locale(std::locale(), new comma_numpunct());
    MMERR.imbue(comma_locale);
    GeneList geneList (parameters, parameters.gtfFileName);
    TableCount table (parameters, geneList);
    std::atomic < unsigned int > i(0);
    std::mutex m1, m2;
    std::vector < std::thread > threads;
    for (unsigned int it = 0; it < parameters.nThreads - 1; ++it) {
        threads.emplace_back(std::thread(doWork, std::ref(parameters), std::ref(geneList), std::ref(table), std::ref(i), std::ref(m1), std::ref(m2)));
    }
    doWork(parameters, geneList, table, i, m1, m2);
    for (std::thread &t: threads) {
        t.join();
    }
    table.selectGenes();
    table.prepareOutput();
    if (parameters.printStructure) {
        table.dump();
    }
    else {
        outputTable = table.getTable();
    }
    if (! parameters.quiet) MMERR << "Successfully done." << std::endl;
    return EXIT_SUCCESS;
}

static std::pair < std::vector< std::string >, std::vector < std::pair < std::string, std::vector < unsigned int > > > > getTable (Parameters &parameters) {
    return std::make_pair(parameters.names, outputTable);
}