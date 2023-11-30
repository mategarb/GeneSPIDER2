/*
 * Copyright (C) 2003, 2008  Center for Computational Biology and Bioinformatics
 * All Rights Reserved.
 *
 * util.cpp
 *
 */

#include "util.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

double median(vector<double> sortedData, const int n){
  double median ;
  const int lhs = (n - 1) / 2 ;
  const int rhs = n / 2 ;

  if (n == 0)
    return 0.0 ;

  if (lhs == rhs)
    {
      median = sortedData.at(lhs) ;
    }
  else
    {
      median = (sortedData.at(lhs) + sortedData.at(rhs))/2.0 ;
    }

  return median;
}

double interQuartileRange(double * sortedData, const int size){

 // Q(1) = median(y(find(y<median(y))));

  sort ( sortedData , sortedData + size );
  
  int medianIndex = (size + 2) / 2;

  vector<double> subset;

  for (int i = 0; i < (medianIndex - 1); i++){
    subset.push_back( sortedData[i] );
  }

  double q1 = median( subset, subset.size() ); 

  // Q(3) = median(y(find(y>median(y))));

  medianIndex = (size + 1) / 2;

  subset.clear();

  for (int i = medianIndex; i < size; i++){
    subset.push_back( sortedData[i] );
  }

  double q3 = median( subset, subset.size() );

  subset.clear();

  return q3 - q1;
}

double normalPDF(double dx, double sigma){
  double y = exp(-0.5 * pow((dx/sigma), 2)) / (sqrt(2 * M_PI) * sigma);
  return y;
}

double multinormalPDF(double dx, double dy, double sigmaX, double sigmaY){
  double y = exp(-0.5 * (pow((dx / sigmaX), 2) + pow((dy / sigmaY), 2)));
  y = y /(2 * M_PI * sigmaX * sigmaY);
  return y;
}
