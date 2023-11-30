//
// Copyright (C) 2003  Columbia Genome Center
// All Rights Reserved.
//
// MutualInfo.h -- Class and constant definitions for mutual
// information calculation progoram.
//
// $Id: MutualInfo.h,v 1.2 2008/10/23 21:09:22 manju Exp $
//

#ifndef MUTUALINFO_H__
  #define MUTUALINFO_H__
  #include "TypeManip.h"

#define M_SQRT2         1.41421356237309504880

using namespace std;


class Gene_Pair
{

  double x;
  double y;
  int xi; // index of x
  int yi; // index y
  int maId;

public:

  inline double Get_X() const { return x; }

  inline double Get_Y() const { return y; }

  inline int Get_XI() const { return xi; }

  inline int Get_YI() const { return yi; }

  inline int Get_MaID() const { return maId; }

  inline void Set_X( double X ) { x = X; }

  inline void Set_Y( double Y ) { y = Y; }

  inline void Set_XI( int XI ) { xi = XI; }

  inline void Set_YI( int YI ) { yi = YI; }

  inline void Set_MaID( int MaID ) { maId = MaID; }
};


typedef std::vector < Gene_Pair > Pair_Vector;


class Sort_X : std::binary_function < Gene_Pair, Gene_Pair, bool >
{

public:

  bool operator() ( const Gene_Pair & a, const Gene_Pair & b ) const
  {
    if ( a.Get_X() != b.Get_X() )
    {
      return ( a.Get_X() < b.Get_X() );
    }
    else
    {
      return ( a.Get_MaID() < b.Get_MaID() );
    }
  }
};


class Sort_Y : std::binary_function < Gene_Pair, Gene_Pair, bool >
{

public:

  bool operator() ( const Gene_Pair & a, const Gene_Pair & b ) const
  {
    if ( a.Get_Y() != b.Get_Y() )
    {
      return ( a.Get_Y() < b.Get_Y() );
    }
    else
    {
      return ( a.Get_MaID() < b.Get_MaID() );
    }
  }
};



/****************************************************************************
* Class MutualInfo
*
* The main MutualInfo class
*
**/
class Mutual_Info
{
  static const int MIBLOCKS = 2;
  static const int BINS = 1000;

  int     Microarray_Num;
  int     ** MI_Space;
  double  ** MI_Prob;
  double  MA_Per_MI_Step;
  double  * histogram;
  double  * background;
  int     Max_Histo_Bin;
  int     Max_Background_Bin;

  //type of caluclation
  int type;

  // variance2 = 2*sigma*sigma, doing this just same some computation
  double variance2;
  double ** prob_table;
  double *  norm_1D_table;
  double ** norm_2D_table;

  double * kernel_bandwidth;

  bool Copula_Transform;

public:

  static int const RANK_NONE = 0;
  static int const STD_REGRESSION = 1;
  static int const RANK_REGRESSION = 2;
  static int const MI_GAUSSIAN = 3;
  static int const MI_ADAPTIVE_PARTITIONING = 4;
  static int const MI_GAUSSIAN_NO_COPULA = 5;

  Mutual_Info( int, double, int );
  ~Mutual_Info();
  int Get_Type() const { return type; }
  void Set_Copula_Transform( bool t = true ) { Copula_Transform = t; }
  bool Is_Copula_Transform() const { return Copula_Transform; }
  double Compute_Pairwise_MI( Pair_Vector & , double sigmaX, double sigmaY);
  int Get_Microarray_Num() { return Microarray_Num; }

private:

  void Initialize_Norm_Table( double );
  double Get_Kernel( Pair_Vector & pairs, unsigned int i , const Loki::Int2Type < Mutual_Info::MI_GAUSSIAN > & );
  double Get_Kernel( Pair_Vector & pairs, unsigned int i , double sigmaX, double sigmaY, const Loki::Int2Type < Mutual_Info::MI_GAUSSIAN_NO_COPULA > & );
  double Get_Mutual_Info_XY( Pair_Vector &, const Loki::Int2Type < Mutual_Info::MI_ADAPTIVE_PARTITIONING > & );
  double Get_Mutual_Info_XY( Pair_Vector &, const Loki::Int2Type < Mutual_Info::MI_GAUSSIAN > & );
  double Get_Mutual_Info_XY( Pair_Vector &, double sigmaX, double sigmaY, const Loki::Int2Type < Mutual_Info::MI_GAUSSIAN_NO_COPULA > & );

};

#endif



