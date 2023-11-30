
/*
 * Copyright (C) 2003, 2004  Columbia Genome Center
 * All Rights Reserved.
 *
 * param.h --
 *
 * $Id: param.h,v 1.1.1.1 2008/10/08 19:30:44 manju Exp $
 */
#ifndef PARAM_H__
#define PARAM_H__

#include <string>
#include <vector>

using namespace std;

/**
  * Parameters for the algorithm
  */
struct Parameter{
	// Default parameters:
	static const double default_threshold;
        static const double default_pvalue;
        static const double default_eps;
	static const double default_sigma;
        static const int default_sample;
        static const double default_percent;
	static const double default_mean;
	static const double default_cv;
        static const double default_correction;

	double threshold; // mi threshold
        double pvalue;
        double eps;       // DPI tolerance
        double sigma;     // gaussian kernel width
        int sample;       // sample number
        double percent;   // high/low percentage for the conditional analysis
	double mean;      // filter mean
	double cv;        // filter coefficient of variance
        double correction;  // coorection for noise
	string verbose;
	string algorithm;
        string infile;
        string outfile;
        string adjfile;
        string hub;
        string subnetfile;
        string annotfile;
        string controlId;
        string condition;
        string home_dir;
        vector < string > subnet;
        vector < string > tf_list;

	Parameter():threshold(default_threshold), pvalue(default_pvalue), eps(default_eps),
                sigma(default_sigma),
                sample(default_sample), percent(default_percent),
                mean(default_mean), cv(default_cv), correction(default_correction),
                verbose("off"), algorithm("fixed_bandwidth"), infile(""), outfile(""),
                adjfile(""), hub(""), subnetfile(""), annotfile(""), controlId(""),
                condition(""), home_dir("./")
		{
		}
};

// service functions
string getFileName(std::string matrixName);
void checkParameter(Parameter &p);
void displayParameter(Parameter &p);
int readProbeList(std::string subnetfile, std::vector<std::string> & subnet);
int cmpIgnoreCase(const char* a, const char* b);
void createOutfileName(Parameter &p );
#endif
