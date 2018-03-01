#include <Rcpp.h>
#include "mmquant.h"

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]


/*
RCPP_MODULE(Parameters_module) {
  using namespace Rcpp;
 
  class_<Parameters>( "Parameters")
    .default_constructor("Default constructor")
    .method("setGtfFileName", &Parameters::setGtfFileName)
  ;
}
*/

Parameters parameters;

//' This is stuff
//'
//' @return nothing
//' @export
// [[Rcpp::export]]
void init () {
	parameters.geneOverlapFunction = geneInclusion;
    parameters.printStructure      = false;
    parameters.quiet               = true;
}

// [[Rcpp::export]]
void printUsage () {
	parameters.printUsage();
}

// [[Rcpp::export]]
void printState () {
	parameters.printState();
}

// [[Rcpp::export]]
void setGtfFileName(std::string &s) {
	parameters.setGtfFileName(s);
}

// [[Rcpp::export]]
void addReadsFileName(std::string &s) {
	parameters.addReadsFileName(s);
}

// [[Rcpp::export]]
void addName(std::string &s) {
	parameters.addName(s);
}

// [[Rcpp::export]]
void setOutputFileName(std::string &s) {
	parameters.setOutputFileName(s);
}

// [[Rcpp::export]]
void setOverlap(float f) {
	parameters.setOverlap(f);
}

// [[Rcpp::export]]
int addStrand(std::string &s) {
	return parameters.addStrand(s);
}

// [[Rcpp::export]]
int addSort(std::string &s) {
	return parameters.addSort(s);
}

// [[Rcpp::export]]
void setCountThreshold(unsigned int u) {
	parameters.setCountThreshold(u);
}

// [[Rcpp::export]]
void setMergeThreshold(float f) {
	parameters.setMergeThreshold(f);
}

// [[Rcpp::export]]
void setPrintGeneName(bool b) {
	parameters.setPrintGeneName(b);
}

// [[Rcpp::export]]
void setFeatureCountStyle(bool b) {
	parameters.setFeatureCountStyle(b);
}

// [[Rcpp::export]]
void setProgress(bool b) {
	parameters.setProgress(b);
}

// [[Rcpp::export]]
void setNThreads(int n) {
	parameters.setNThreads(n);
}

// [[Rcpp::export]]
int addFormat(std::string &s) {
	return parameters.addFormat(s);
}

// [[Rcpp::export]]
void setNOverlapDifference(int n) {
	parameters.setNOverlapDifference(n);
}

// [[Rcpp::export]]
void setPcOverlapDifference(float f) {
	parameters.setPcOverlapDifference(f);
}

// [[Rcpp::export]]
NumericMatrix start () {
	start(parameters);
    auto table        = getTable(parameters);
    auto &sampleNames = table.first;
    auto &rest        = table.second;
    NumericMatrix matrix(rest.size(), sampleNames.size());
    CharacterVector genes(rest.size()), samples(sampleNames.size());
    for (size_t i = 0; i < sampleNames.size(); ++i) {
        samples[i] = sampleNames[i];
    }
    for (size_t i = 0; i < rest.size(); ++i) {
        auto         &line = rest[i];
        NumericVector l    = wrap(line.second);
        genes[i]           = line.first;
        matrix.row(i)      = l;
    }
    colnames(matrix) = samples;
    rownames(matrix) = genes;
    return matrix;
}
