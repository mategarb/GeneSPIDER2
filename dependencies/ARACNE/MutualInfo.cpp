//
// Copyright (C) 2003 Columbia Genome Center
// All Rights Reserved.
//
// $Id: MutualInfo.cpp,v 1.2 2008/10/08 23:04:18 mb3113 Exp $
//

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <valarray>
#include "MutualInfo.h"
#include "util.h"

using namespace std;
typedef unsigned int uint;

Mutual_Info::Mutual_Info( int maNo, double sigma, int _type )
  : Microarray_Num( maNo ), type( _type ), variance2( 2*sigma*sigma )
{
  switch ( type )
  {
    case MI_ADAPTIVE_PARTITIONING:
      Set_Copula_Transform( false );
    break;

    case MI_GAUSSIAN:

      // set copula transform
      Set_Copula_Transform( true );

      // initialize the lookup table for computing Gaussian kernels
      prob_table = new double * [maNo];
      for ( int i = 0; i < maNo; i++ )
      {
        prob_table[i] = new double[maNo];
        for ( int j = 0; j < maNo; j++ )
        {
          prob_table[i][j] = -1.0;
        }
      }

      // initialize the lookup table for computing normalization factors
      norm_1D_table = new double[maNo];

      norm_2D_table = new double * [maNo];
      for ( int p = 0; p < maNo; p++ )
      {
        norm_2D_table[p] = new double[maNo];
      }

      Initialize_Norm_Table(sigma);

    break;

    case MI_GAUSSIAN_NO_COPULA:
      // do not use copula transform
      Set_Copula_Transform( false );
    break;

    default:
      throw std::string( "MI computation not supported" );
  }
}

void Mutual_Info::Initialize_Norm_Table( double sigma)
{
  // initialize the 1D table
  for ( int i=0; i<Microarray_Num; i++ )
  {
    // array id is from 0 to M-1, after copula transform, it becomes from 0 to 1
    // with 1/(2M) at both sides
    double x = Copula_Transform ? ( ( double(i+1)-0.5 )/Microarray_Num ) : double( i );
    norm_1D_table[i] = 0.5*( erf( (1-x)/(sigma*M_SQRT2) ) - erf( (0-x)/(sigma*M_SQRT2) ) );
  }

  // initialize the 2D table
  for (int j=0; j<Microarray_Num; j++)
  {
    for ( int k=j; k<Microarray_Num; k++ )
    {
      norm_2D_table[j][k] = norm_1D_table[j]*norm_1D_table[k];
      norm_2D_table[k][j] = norm_2D_table[j][k];
    }
  }
}

/** Get Score for MI_GAUSSIAN */
double Mutual_Info::Compute_Pairwise_MI( Pair_Vector & pairs , double sigmaX, double sigmaY)
{
  double mi;

  switch ( type )
  {
    case Mutual_Info::MI_ADAPTIVE_PARTITIONING:
      mi = Get_Mutual_Info_XY( pairs, Loki::Int2Type < Mutual_Info::MI_ADAPTIVE_PARTITIONING > () );
    break;

    case Mutual_Info::MI_GAUSSIAN:
      mi = Get_Mutual_Info_XY( pairs, Loki::Int2Type < Mutual_Info::MI_GAUSSIAN > () );
    break;

    case Mutual_Info::MI_GAUSSIAN_NO_COPULA:
      mi = Get_Mutual_Info_XY( pairs, sigmaX, sigmaY, Loki::Int2Type < Mutual_Info::MI_GAUSSIAN_NO_COPULA > () );
    break;

    default:
      throw std::string( "MI computation not supported" );
  }
  return ( mi );
}


Mutual_Info::~Mutual_Info()
{
  switch ( type )
  {
    case MI_ADAPTIVE_PARTITIONING:
    break;

    case MI_GAUSSIAN:
      for ( int i = 0; i < Microarray_Num; i++ )
      {
        delete[] prob_table[i];
        delete[] norm_2D_table[i];
      }

      delete[] prob_table;
      delete[] norm_2D_table;
      delete[] norm_1D_table;

    break;

    case MI_GAUSSIAN_NO_COPULA:
    break;

    default:
      throw std::string( "MI computation not supported" );
  }
}


double Mutual_Info::Get_Mutual_Info_XY( Pair_Vector & pairs, const Loki::Int2Type < Mutual_Info::MI_GAUSSIAN > & rankType )
{
  const uint size = pairs.size();

  // rank and copula transformation
  Sort_X X_Sorter;
  sort( pairs.begin(), pairs.end(), X_Sorter );
  for ( uint i = 0; i < size; i++ )
  {
    pairs[i].Set_XI(i);
    if (Copula_Transform)
    {
      pairs[i].Set_X( ( double(i+1)-0.5 )/size );
    }
  }

  Sort_Y Y_Sorter;
  sort( pairs.begin(), pairs.end(), Y_Sorter );
  for ( uint i = 0; i < size; i++ )
  {
    pairs[i].Set_YI( i );
    if (Copula_Transform)
    {
      pairs[i].Set_Y( ( double(i+1)-0.5 )/size );
    }
  }

  double sum = 0.0;
  for ( uint i=0; i<size; i++ )
  {
    double v = Get_Kernel( pairs, i , Loki::Int2Type < Mutual_Info::MI_GAUSSIAN > ());
    sum += log( v );
  }
  // double mi = sum / (static_cast <double>size);
  double mi = sum/(double)size;
  return ( std::max( mi, 0.0 ) );

}

double Mutual_Info::Get_Mutual_Info_XY( Pair_Vector & pairs, double sigmaX, double sigmaY, const Loki::Int2Type < Mutual_Info::MI_GAUSSIAN_NO_COPULA > & rankType )
{
  const uint size = pairs.size();

  double sum = 0.0;
  for ( uint i=0; i<size; i++ )
  {
    double v = Get_Kernel( pairs, i, sigmaX, sigmaY, Loki::Int2Type < Mutual_Info::MI_GAUSSIAN_NO_COPULA > ());
    sum += log( v );
  }
  double mi = sum/(double)size;
  return ( std::max( mi, 0.0 ) );
}

double Mutual_Info::Get_Kernel( Pair_Vector & pairs, unsigned int i, double sigmaX, double sigmaY, const Loki::Int2Type < Mutual_Info::MI_GAUSSIAN_NO_COPULA > & rankType )
{
  double fxy = 0.0;
  double fx = 0.0;
  double fy = 0.0;

  const uint size = pairs.size();

  for ( uint j = 0; j < size; j++ )
  {
    double mu_x = pairs[j].Get_X();
    double mu_y = pairs[j].Get_Y();

    double dx = std::abs( pairs[i].Get_X() - mu_x );
    double dy = std::abs( pairs[i].Get_Y() - mu_y );

    fx += normalPDF(dx, sigmaX);
    fy += normalPDF(dy, sigmaY);
    fxy += multinormalPDF(dx, dy, sigmaX, sigmaY);
  }

  return ( ( fxy * size ) / ( fx * fy ) );
}

double Mutual_Info::Get_Kernel( Pair_Vector & pairs, unsigned int i, const Loki::Int2Type < Mutual_Info::MI_GAUSSIAN > & rankType )
{
  // mu_x, mu_y - the center of each gaussian kernel
  // ix, iy     - distance to the gaussian center in unit of index.
  // dx, dy     - the actual distance used to compute the kernel denstiy.
  //              if compula transformed, they are equal to ix and iy rescaled
  //              between 0 and 1;

  double fxy = 0.0;
  double fx = 0.0;
  double fy = 0.0;

  const uint size = pairs.size();

  for ( uint j = 0; j < size; j++ )
  {
    int mu_x = pairs[j].Get_XI();
    int mu_y = pairs[j].Get_YI();

    int ix = std::abs( pairs[i].Get_XI() - mu_x );
    int iy = std::abs( pairs[i].Get_YI() - mu_y );

    double dx = pairs[i].Get_X() - pairs[j].Get_X();
    double dy = pairs[i].Get_Y() - pairs[j].Get_Y();

    // compute the kernel density as necessary
    if ( prob_table[ix][iy] == -1.0 )
    {
      // if (ix, iy) is not computed, at least one of ix and iy is not computed
      // otherwise, all three should have been computed
      // when either ix=0 or iy=0, the trctdGssnBi.getProbability(dx, dy) is
      // equal to trctdGssnUni.getProbability(dx) or trctdGssnUni.getProbability(dy)
      if ( prob_table[ix][0] == -1.0 )
      {
        prob_table[ix][0] = std::exp(-(dx*dx)/variance2);
        prob_table[0][ix] = prob_table[ix][0];
      }

      if ( prob_table[0][iy] == -1.0 )
      {
        prob_table[0][iy] = std::exp(-(dy*dy)/variance2);
        prob_table[iy][0] = prob_table[0][iy];
      }

      prob_table[ix][iy] = prob_table[ix][0]*prob_table[0][iy];
      prob_table[iy][ix] = prob_table[ix] [iy];
    }

    fx += prob_table[ix][0]/norm_1D_table[mu_x];
    fy += prob_table[0][iy]/norm_1D_table[mu_y];
    fxy += prob_table[ix][iy]/norm_2D_table[mu_x][mu_y];
  }

  return ( ( fxy * size ) / ( fx * fy ) );
}

double Mutual_Info::Get_Mutual_Info_XY( Pair_Vector & pairs, const Loki::Int2Type < Mutual_Info::MI_ADAPTIVE_PARTITIONING > & rankType )
{
  Sort_X X_Sorter;
  for ( uint i = 0; i < pairs.size(); i++ )
  {
    pairs.at(i).Set_XI(i);
    pairs.at(i).Set_YI(i);
  }
  sort( pairs.begin(), pairs.end(), X_Sorter );
  valarray<int> xranks(pairs.size());
  for (uint i = 0; i < pairs.size(); i++) {
    xranks[pairs.at(i).Get_XI()] = i + 1;
  }

  Sort_Y Y_Sorter;
  sort( pairs.begin(), pairs.end(), Y_Sorter );
  valarray<int> yranks(pairs.size());
  for (uint i = 0; i < pairs.size(); i++) {
    yranks[pairs.at(i).Get_YI()] = i + 1;
  }

  double xcor = 0; int npar = 1; int run = 0; int N = pairs.size();
  valarray<int> poc(1, 20); valarray<int> kon(N, 20); valarray<int> poradi(N);
  valarray<int> NN(4); valarray<int> marg(0, 80);
 
  for (int i = 1; i <= N; i++) {
    poradi[i-1] = i;
  }

  int t[] = {1, 1, N, N};
  marg[slice(0, 4, 20)] = valarray<int>(t, 4);

  while (npar > 0) {
    run++;
    int apoc = poc[npar - 1];
    int akon = kon[npar - 1];
    valarray<int> apor(akon - apoc + 1);
    for (int i = 0; i <= (akon - apoc); i++) {
      apor[i] = poradi[apoc + i - 1];
    }

    int Nex = apor.size();

    valarray<int> rownpar = marg[slice((npar - 1), 4, 20)];

    int ave1 = (int)floor((rownpar[0] + rownpar[2]) / 2);
    int ave2 = (int)floor((rownpar[1] + rownpar[3]) / 2);

    valarray<bool> J1(apor.size());

    valarray<size_t> apor1 = valarray<size_t>(apor.size());

    for (uint i = 0; i < apor.size(); i++) {
      apor1[i] = apor[i] - 1;
    }

    valarray<int> maskedxranks = valarray<int>(xranks[apor1]);

    J1[maskedxranks <= ave1] = true; J1[maskedxranks > ave1] = false;

    valarray<bool> J2(apor.size());

    valarray<int> maskedyranks = valarray<int>(yranks[apor1]);
    J2[maskedyranks <= ave2] = true; J2[maskedyranks > ave2] = false;

    valarray<bool> Ia(4 * apor.size());

    Ia[slice(0, apor.size(), 4)] = J1 && J2;
    Ia[slice(1, apor.size(), 4)] = J1 && (!J2);
    Ia[slice(2, apor.size(), 4)] = (!J1) && J2;
    Ia[slice(3, apor.size(), 4)] = (!J1) && (!J2);


    valarray<int> I( 4 * apor.size());

    I[Ia] = 1; I[!Ia] = 0;

    NN[0] = (valarray<int>(I[slice(0, apor.size(), 4)])).sum();
    NN[1] = (valarray<int>(I[slice(1, apor.size(), 4)])).sum();
    NN[2] = (valarray<int>(I[slice(2, apor.size(), 4)])).sum();
    NN[3] = (valarray<int>(I[slice(3, apor.size(), 4)])).sum();

    valarray<int> amarg(16);
 
    valarray<int> slicemarg = valarray<int>(marg[slice((npar - 1), 4, 20)]); 

    int t[] = { slicemarg[0], slicemarg[1] , ave1, ave2 };
    amarg[slice(0, 4, 4)] = valarray<int>(t, 4);

    int t1[] = { slicemarg[0], ave2 + 1, ave1, slicemarg[3] };
    amarg[slice(1, 4, 4)] = valarray<int>(t1, 4);

    int t2[] = { ave1 + 1, slicemarg[1], slicemarg[2], ave2 };
    amarg[slice(2, 4, 4)] = valarray<int>(t2, 4);

    int t3[] = { ave1 + 1, ave2 + 1, slicemarg[2], slicemarg[3] };
    amarg[slice(3, 4, 4)] = valarray<int>(t3, 4);

    valarray<double> NN2(4);
    NN2[0] = (double)NN[0];
    NN2[1] = (double)NN[1];
    NN2[2] = (double)NN[2];
    NN2[3] = (double)NN[3];

    double tst = 4 * valarray<double>(pow(NN2 - (double)(Nex/4), 2)).sum() / Nex;

    if (tst > 7.8 || run == 1) {
      --npar;
      for (uint i = 0; i < 4; i++) {
        if (NN[i] > 2) {
          ++npar;
          akon = apoc + NN[i] - 1;
          poc[npar - 1] = apoc;
          kon[npar - 1] = akon; 
          marg[slice((npar - 1), 4, 20)] = amarg[slice(i, 4, 4)];

          valarray<int> indices = valarray<int>(I[slice(i, apor.size(), 4)]);
          valarray<int> t = valarray<int>(apor[indices==1]);
          for (int kk = (apoc - 1); kk < akon; kk++){
            poradi[kk] = t[kk - apoc + 1];
	  }
          apoc = akon + 1;
        }
        else if (NN[i] > 0) {
          valarray<int> t = valarray<int>(amarg[slice(i, 4, 4)]);
          double Nx = t[2] - t[0] + 1;
          double Ny = t[3] - t[1] + 1;
          xcor += NN[i] * log(NN[i] / (Nx * Ny));
        }
      }
    }
    else {
      valarray<int> t = valarray<int>(marg[slice((npar - 1), 4, 20)]);
      double Nx = t[2] - t[0] + 1;
      double Ny = t[3] - t[1] + 1;
      xcor += Nex * log(Nex / (Nx * Ny));
      --npar;
    }
  }
  xcor = xcor / N + log(N);
  return xcor;
}
