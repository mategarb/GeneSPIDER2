
/*
 * Copyright (C) 2003, 2008 Center for Computational Biology and Bioinformatics, Columbia University 
 * All Rights Reserved.
 *
 * util.h --
 *
 */
#ifndef UTIL_H__
#define UTIL_H__

#define M_PI 3.14159265358979323846

using namespace std;

//Utility Math functions 

double median(double* sortedData, int n);

double interQuartileRange(double * sortedData, int size);

double normalPDF(double dx, double sigma);

double multinormalPDF(double dx, double dy, double sigmaX, double sigmaY);

#endif
