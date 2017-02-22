/***************************************************************************
  MFT.inl  - Math Template Foundation.
  -------------------------------------

  General porpose mathematical templates and algorithms bundled into the namespace: "mft".
  These utilities are ment for high-performance computation, so proper call and parameter check
  is left for the developer using them.

  begin     : Donnerstag, 6. September 2007
  copyright : (C) 2007 by Scienomics
  eMail     : Csaba.Kiss@scienomics.com

 ***************************************************************************
 *  This program and all subroutines,  data, and  files used by it are
 *  protected by copyright and hence may not be used, copied, modified
 *  transmitted, inspected, or executed by any means including the use
 *  of electronic data processing equipment,  xerography, or any other
 *  methods  without the express  written permission  of the copyright holder.
 *
 *  Copyright (C) 2007 Scienomics S. A.
 ***************************************************************************
*/

namespace mft {

#ifndef _WINDOWS
  #define _finite  finite
  #define _isnan   isnan
#endif


template< class Ty_> inline 
double Deg2Rad( const Ty_& deg)
  { 
    return static_cast<Ty_>( static_cast<double>( deg)*M_PI/180.0); 
  }


template< class Ty_> inline 
double Rad2Deg( const Ty_& rad) 
  { 
    return static_cast<Ty_>( static_cast<double>( rad)*180.0/M_PI); 
  }


//NOTE: Bohr radius (CODATA 2006): 0.52917720859(36)*1e-10. 
template< class Ty_> inline 
double ToAng( const Ty_& x)
  { 
    return static_cast<double>( 0.52917706*x);
  }

template< class Ty_> inline 
double ToAu( const Ty_& x)
  { 
    return static_cast<double>( 1.88972666*x);
  }


template< class Ty_> inline
Ty_ Log2( const Ty_& x)
  {
    static const double cd_InvLog2( 1.0/::log( 2.0));
    return static_cast<Ty_>( ::log( x)*cd_InvLog2);
  }

template< class Ty_> inline
valarray<Ty_> Log2( const valarray<Ty_>& v)
  {
    static const double cd_InvLog2( 1.0/::log( 2.0));
    return log( v)*cd_InvLog2;
  }


template< class Ty_> inline
Ty_ Abs( const Ty_& v)
  {
    return ( (v<(0)) ? (-v) : (v));
  }


//FIX: Use the MKL's vector functions vdFloor/vdCeil/vdTrunc/vdRound/vdNearbyInt/vdRint
inline 
double Round( const double& d) 
  { 
    double s( (d < 0) ? -1.0 : 1.0);
    double r( floor( d + 0.5 + s*numeric_limits<double>::epsilon()));

    return r;
  }


//FIX: replace with rounding from MKL
struct fuRound {
  template< class Ty_> Ty_ operator()( const Ty_& d) { 
      double s( ( (double)d < 0.0) ? -1.0 : 1.0);
      double r( floor( (double)d + 0.5 + s*numeric_limits<double>::epsilon()));
      return r;
    }
  };


inline 
const valarray<double>& Round( valarray<double>& m) 
  { 
    transform( &m[0], &m[m.size()], &m[0], fuRound());
    return m;
  }


inline 
const valarray<double>& Round( valarray<double>& mR, const valarray<double>& m) 
  { 
    SetM( mR, m);
    Round( mR);

    return mR;
  }


inline 
valarray<double> Round2( const valarray<double>& vA)
  { 
    valarray<double> vR( vA);
    Round( vR);

    return vR;
  }


//FIX: Use the MKL's vector functions vdFloor/vdCeil/vdTrunc/vdRound/vdNearbyInt/vdRint
template< class Ty_> inline
const valarray<Ty_>& Round( valarray<Ty_>& m) 
  { 
    m = ( m + static_cast<Ty_>(0.5) - numeric_limits<Ty_>::epsilon()).apply( floor);
    return m;
  }


//FIX: Use the MKL's vector functions vdFloor/vdCeil/vdTrunc/vdRound/vdNearbyInt/vdRint
template< class Ty_> inline
const valarray<Ty_>& Round( valarray<Ty_>& mR, const valarray<Ty_>& m) 
  { 
    mR.resize( m.size());
    mR = ( m + static_cast<Ty_>(0.5) - numeric_limits<Ty_>::epsilon()).apply( floor);

    return mR;
  }


//FIX: Use the MKL's vector functions vdFloor/vdCeil/vdTrunc/vdRound/vdNearbyInt/vdRint
template< class Ty_> inline
const valarray<Ty_>& RoundDec( valarray<Ty_>& m, const Ty_& prec)
  {
    if ( prec > numeric_limits<Ty_>::epsilon()) {
      // Precision as subunit epsylon
      if ( prec < 1.0) {
        valarray<Ty_> m_( m/(double)prec);
        Round( m_);
        m = m_*prec;
        }

      // Precision as decimal places
      else {
        Ty_ cP( pow( static_cast<Ty_>( 10), prec));
        valarray<Ty_> m_( m*cP);
        Round( m_);
        m = m_/cP;   //Assert: cP>0
        }
      }

    return m;    
  }


template< class Ty_> inline
valarray<Ty_> RoundDec( const valarray<Ty_>& m, const Ty_& prec)
  {
    valarray<Ty_> m2( m);
    RoundDec( m2, prec);

    return m2;
  }


inline 
valarray<double> Fix( valarray<double>& v)
  {
    valarray<double> vR( v.size());
    vdTrunc( (MKL_INT)v.size(), &v[0], &vR[0]);

    return vR;
  }


inline 
double Fix( const double& dV)
  {
    double dR;
    vdTrunc( 1, &dV, &dR);

    return dR;
  }


//DSG: returning 0 for underflow, is this ok?
template< class Ty_> inline
short Sign( const Ty_& v) 
  {
    if ( v < -numeric_limits<Ty_>::epsilon()) 
      return -1;

    else if ( v > numeric_limits<Ty_>::epsilon()) 
      return 1;

    return 0;
  }


template< class Ty_> inline
Ty_ gcd( const Ty_& u, const Ty_& v)
  {
    if ( !u || !v) {
      return max( abs(u), abs(v));
      }
    else {
      Ty_ r, u_( abs(u)), v_( abs(v));

      if ( u_ < v_) swap( u_, v_);
      while ( v_) {
        r = u_ % v_;
        u_ = v_;
        v_ = r;
        }

      return u_;
      }
  }


template< class Ty_> inline
Ty_ gcd_Stein( const Ty_& u, const Ty_& v)
  {
    if ( !u || !v) {
      return max( abs(u), abs(v));
      }
    else {
      Ty_ r, u_( abs(u)), v_( abs(v));
      size_t k(0);

      while( !(u_&1) && !(v_&1)) {
        u_ >>= 1;
        v_ >>= 1;
        k++;
        }

      Ty_ t( (u_&1) ? -v_ : u_);
      while ( t) {
        while( !(t&1)) t >>= 1;
        
        if ( t>0)
          u_ = t;
        else
          v_ = -t;

        t = u_ - v_;
        }
      return (u_<<k);
      }    
  }


template< class Ty_> inline 
Ty_ gcd( const valarray<Ty_>& vX)
  {
    const size_t nL( vX.size());

    if ( !nL || All( vX == Ty_(0))) 
      return 0;

    else if ( nL == 1)
      return abs( vX[0]);

    else {
      valarray<Ty_> vx( abs( vX));
      long k;
      Ty_ d;

      d = vx[ nL-1];
      k = nL - 2;

      while ( d != 1 && k >= 0)   //  Dirichlet, 1849
        d = gcd( d, vx[ k--]);
      
      return d;
      }          
  }


template< class Ty_> inline
const valarray<Ty_>& Gcd( valarray<Ty_>& vGCD, const Ty_& A, const Ty_& B)
  {
    vGCD.resize(3, Ty_(0));

    if ( A || B) {
      const Ty_ vU_[] = { 1, 0, A};
      const Ty_ vV_[] = { 0, 1, B};
      valarray<Ty_> vT(3), vV( vV_, 3), vU( vU_, 3);
      Ty_ q;

      while( vV[2]) {
        q = static_cast<Ty_>( floor( (double)vU[2]/vV[2]));
        vT = vU - q*vV;
        vU = vV;
        vV = vT;
        }
      
      SetM( vGCD, vU.cshift(2));
      if ( vGCD[0] < 0) vGCD *=-1;
      }

    return vGCD;
  }


template< class Ty_> inline
Ty_ Gcd( const Ty_& A, const Ty_& B)
  {
    valarray<Ty_> vGCD;
    return Gcd( vGCD, A, B)[0];
  }


template< class Ty_> inline
void Cross( Ty_* vC, const Ty_* vA, const Ty_* vB)
  {
    vC[0] = vA[1]*vB[2] - vA[2]*vB[1];
    vC[1] = vA[2]*vB[0] - vA[0]*vB[2];
    vC[2] = vA[0]*vB[1] - vA[1]*vB[0];
  }


/** cross( v12, v1, v2) */
template< class Ty_> inline
const valarray<Ty_>& Cross( valarray<Ty_>& vC, const valarray<Ty_>& vA, const valarray<Ty_>& vB)
  {
    vC.resize( 3);
    vC[0] = vA[1]*vB[2] - vA[2]*vB[1];
    vC[1] = vA[2]*vB[0] - vA[0]*vB[2];
    vC[2] = vA[0]*vB[1] - vA[1]*vB[0];

    return vC;
  }


/** v12 <-- cross( v1, v2) */      
//DSG: this is a superfluous special case of cross( mV1, mV2, nCols), with nCols=1 !
template< class Ty_> inline
valarray<Ty_> Cross( const valarray<Ty_>& vA, const valarray<Ty_>& vB)
  {
    valarray<Ty_> vC;
    Cross( vC, vA, vB);

    return vC;
  }


//DSG: does this really help on Linux?
template< class Ty_> inline
valarray<Ty_> Cross( const slice_array<Ty_>& vA, const slice_array<Ty_>& vB)
  {
    return Cross( valarray<Ty_>( vA), valarray<Ty_>( vB));
  }


template< class Ty_> inline
valarray<Ty_> Cross( const valarray<Ty_>& mA, const valarray<Ty_>& mB, const size_t& nc)
  {
    valarray<Ty_> mC;
    Cross( mC, mA, mB, nc);
    return mC;
  }


template< class Ty_> inline
const valarray<Ty_>& Cross( valarray<Ty_>& mC, const valarray<Ty_>& mA, const valarray<Ty_>& mB, const size_t& nc)
  {
    mC.resize( mA.size());
    for( size_t n = 0; n < nc; n++)
      mC[slice(n,3,nc)] = Cross( valarray<double>( mA[slice(n,3,nc)]), 
                                 valarray<double>( mB[slice(n,3,nc)]) );
    return mC;
  }


//DSG: replace with proper expression templates !!!
template< class Ty_> inline
Ty_ Mixed( const valarray<Ty_>& mV)
  {
    return Mixed( valarray<Ty_>( mV[slice(0,3,3)]), 
                  valarray<Ty_>( mV[slice(1,3,3)]), 
                  valarray<Ty_>( mV[slice(2,3,3)]) );
  }


template< class Ty_> inline
Ty_ Mixed( const valarray<Ty_>& v1, const valarray<Ty_>& v2, const valarray<Ty_>& v3)
  {
    valarray<Ty_> t;
    Cross( t, v1, v2);

    return Dot( t, v3);
  }


inline 
valarray<double> Mixed( const valarray<double>& mA, const valarray<double>& mB, const valarray<double>& mC, const size_t& nc)
  {
    valarray<double> mD( Cross( mB, mC, nc)); 
    return Dot( mA, mD, nc);
  }


template< class Ty_> inline
Ty_ MNorm( const valarray<Ty_>& mA, const size_t& nrA, const size_t& nNorm)
  {
    Ty_ norm( 0);

    switch( nNorm) {

      case 1: {
        valarray<Ty_> absV( abs( mA));
        norm = MaxVal( SumD1( absV, nrA));

        break;
        }

      case 2: {
        valarray<Ty_> vS;
        int iErr = Svd( vS, mA, nrA);
        if (!iErr) norm = MaxVal( vS);

        break;
        }

      case 3: {
        valarray<Ty_> absV( abs( mA));
        norm = MaxVal( SumD2( absV, nrA));

        break;
        }

      case 4: {
        size_t ncA( mA.size()/nrA);           // Asserting: nrA>0
        valarray<Ty_> mAA( MM( Trp( mA, nrA), mA, ncA, ncA));
        norm = sqrt( Trace( mAA, ncA));

        break;
        }
      }

    return norm;
  }



template< class Ty_> inline
Ty_ Norm( valarray<Ty_>& vX)
  {
    // convert to real, do real norm, convert back to Ty_
    return Norm( Conv<double,Ty_>(vX));
  }

template< class Ty_> inline
Ty_ Norm( const valarray<Ty_>& vX)
  {
    // convert to real, do real norm, convert back to Ty_
    return Norm( Conv<double,Ty_>(vX));
  }


template< class Ty_> inline
const valarray<Ty_>& Norm( valarray<Ty_>& vN, const valarray<Ty_>& mA, const size_t& nrA)
  {
    vN.resize(0);

    if ( nrA) {
      valarray<Ty_> mT( mA*mA);

      vN.resize( mA.size()/nrA);
      vN = sqrt( SumD1( mT, nrA));
      }

    return vN;
  }


template< class Ty_> inline
valarray<Ty_> Norm( const valarray<Ty_>& mA, const size_t& nrA)
  {
    valarray<Ty_> vN(0);
    
    if ( nrA) {
      valarray<Ty_> mT( mA*mA);

      //DSG: adapting to non-real Ty_!
      //vN.resize( mA.size()/nrA);
      //vN = sqrt( SumD1( mT, nrA));
      valarray<double> vdN(mA.size()/nrA);
      vdN = sqrt( Conv<double,Ty_>( SumD1( mT, nrA)));
      Conv<Ty_,double>( vN, vdN);
      }

    return vN;
  }


inline
const valarray<double>& Norm2( valarray<double>& vN2, const valarray<double>& mA, const size_t& nrA)
  {
    vN2.resize(0);

    if ( nrA) {
      valarray<double> mT( mA*mA);

      vN2.resize( mA.size()/nrA);
      vN2 = SumD1( mT, nrA);
      }

    return vN2;
  }


#if 1  // ERRORPRONE!
template< class Ty_> inline
const valarray<Ty_>& Unit( valarray<Ty_>& v)
  {
    Ty_ d( Norm(v));
    v = v/d;       // Asserting: d != 0

    return v;
  }
#endif

template< class Ty_> inline
valarray<Ty_> Unit( const valarray<Ty_>& v)
  {
    valarray<Ty_> vC( v);
    Unit( vC);

    return vC;
  }


template< class Ty_> inline
const valarray<Ty_>& Unit( valarray<Ty_>& mU, const valarray<Ty_>& mV, const size_t& nR)
  {
    SetM( mU, mV);
    mU /= RepCvec( Norm( mV, nR), nR, 1);     // Assert: all( mNrm > 0)

    return mU;
  }


//DSG: breaks on null norm!
template< class Ty_> inline
valarray<Ty_> Unit( const valarray<Ty_>& mV, const size_t& nR)
  {
    valarray<Ty_> mU;
    Unit( mU, mV, nR);

    return mU;
  }


/** Return a vector defined by a direction and length.*/
template< class Ty_> inline
valarray<Ty_> Vector( const valarray<Ty_>& vDir, const Ty_& dLen)
  {
    valarray<Ty_> vR( 0.0, vDir.size());
    Ty_ d( Norm(vDir));

    if ( abs(d) > numeric_limits<Ty_>::epsilon())      
      vR = vDir*(dLen/d);
    
    return vR;
  }


namespace obsolete {

//OBS: it is a little bit slower than indirect indexing
template< class Ty_> inline
valarray<Ty_> Trp( const valarray<Ty_>& m, const size_t& nR)
  {
    static size_t ix[] = { 
      /*2x2*/ 0, 2, 1, 3,
      /*3x3*/ 0, 3, 6, 1, 4, 7, 2, 5, 8,
      /*4x4*/ 0, 4, 8, 12,  1,  5,  9, 13,  2,  6, 10, 14,  3,  7, 11, 15};

    size_t r, c, nC;

    if ( !(nC = m.size()) || !nR)
      return valarray<Ty_>(0);
    else
      nC /=nR;

    // Transpose quadratic matrices: hard wired indirect indexing
    if ( nR == nC)
      switch(nR) {

        case 2: 
          return valarray<Ty_>( m[ valarray<size_t>( ix, 4)]);

        case 3: 
          return valarray<Ty_>( m[ valarray<size_t>( ix+4, 9)]);

        case 4: 
          return valarray<Ty_>( m[ valarray<size_t>( ix+4+9, 16)]);

        default: {
          valarray<Ty_> mt(m);
          for( r = 0; r < nR; r++) 
            for( c = r; c < nC; c++) swap( mt[r*nC+c], mt[c*nC+r]);

          return mt; 
          }
        }

    // Transpose general matrices: sequential slicing
    else {
      valarray<Ty_> mt(m.size());
      for( r = 0; r < nR; r++) 
        mt[ slice(r,nC,nR)] = valarray<Ty_>( m[ slice(r*nC,nC,1)]);

      return mt;
      }
  }
}  // namespace mft::obsolete


inline
valarray<size_t> Trp( const size_t& nR, const size_t& nC)
  {
    valarray<size_t> vT( nR*nC);
    size_t i, j, r;

    for ( r = 0, j = 0; j < nC; j++)
      for ( i = 0; i < nR; i++) 
        vT[r++] = i*nC + j;

    return vT;
  }


/*OBS:  
  Variants orderd by speed:
    - indirect indexing:  Ty_ = M[vI], vI index array
    - for with slice:        for Ty_[ slice()] = D[ slice()]
        valarray<Ty_> mt(m.size());
        for( r = 0; r < nR; r++) 
          mt[ slice(r,nC,nR)] = valarray<Ty_>( m[ slice(r*nC,nC,1)]);
        return mt;
    - gslice-ing:             Ty_ = M[ gslice()]
  Best performs the direct indexing, has double speed relative to gslice indexing.
  The second variant is just slightly worse  than direct indexing.
*/
template< class Ty_> inline
valarray<Ty_> Trp( const valarray<Ty_>& mA, const size_t& nR)
  {
    return mA[ Trp( nR, mA.size()/nR)];
  }


template< class Ty_> inline
void Trp( valarray<Ty_>& mT, const valarray<Ty_>& mA, const size_t& nR)
  {
    mT.resize( mA.size());
    mT = Trp( mA, nR);
  }


template< class Ty_> inline
valarray<Ty_>& SetM( valarray<Ty_>& mB, const valarray<Ty_>& mA)
  {
    mB.resize( mA.size());
    return mB = mA;
  }


template< class Ty_> inline
vector<Ty_>& SetM( vector<Ty_>& mB, const vector<Ty_>& mA)
  {
    vector<Ty_>& mA_( const_cast< vector<Ty_>&>( mA));
    mB.assign( &mA_[0], &mA_[ mA_.size()]);

    return mB;
  }


template< class Ty_> inline
size_t SetM( valarray<Ty_>& vTrg, const valarray<Ty_>& vSrc, const valarray<size_t>& vIdx)
  {
    const size_t nL( vIdx.size());
    
    vTrg.resize( nL);
    if ( nL && vSrc.size()) {
      for ( size_t n = 0; n < nL; n++)
        vTrg[n] = vSrc[ vIdx[ n]];
      }
      
    return nL;
  }


template< class Ty_> inline
valarray<Ty_>& SetM( valarray<Ty_>& vA, const Ty_* pA, const size_t na)
  {
    vA.resize( 0);

    if ( na) {
      vA.resize( na);
      copy( &pA[0], &pA[na], &vA[0]);
      }

    return vA;
  }


template< class Ty_> inline
valarray<Ty_>& SetM( valarray<Ty_>& vA, const set<Ty_>& setA)
  {
    vA.resize(0);

    if ( setA.size()) {
      vA.resize( setA.size());
      copy( setA.begin(), setA.end(), &vA[0]);
      }

    return vA;
  }


//DSG: transpose alternatives
#if 0

  template< class Ty_> inline
  valarray<Ty_> Trp3_1( const valarray<Ty_>& mD, size_t nR)
    {
      size_t nC = mD.size()/nR;
      size_t vL_[] = {nC,  nR},  vS_[] = {  1,  nC};
      valarray<size_t> vL(vL_,2), vS( vS_,2);

      return mD[ gslice( 0, vL, vS)];
    }


  template< class Ty_> inline
  valarray<Ty_> Trp3_2( const valarray<Ty_>& mD, size_t nR)
    {
      size_t nC = mD.size()/nR;
      valarray<size_t> vL(2), vS( 2);
      vL[0] = nC; vL[1] = nR;
      vS[0] = 1; vS[1] = nC;

      return mD[ gslice( 0, vL, vS)];
    }


  template< class Ty_, size_t nR, size_t nC> inline
  valarray<Ty_> Trp3_3( const valarray<Ty_>& mD)
    {
      size_t nC = mD.size()/nR;
      valarray<size_t> vL(2), vS( 2);
      vL[0] = nC; vL[1] = nR;
      vS[0] = 1; vS[1] = nC;

      return mD[ gslice( 0, vL, vS)];
    }

#endif


template< class Ty_> inline
valarray<Ty_> Skew( const valarray<Ty_>& v)
  {
    valarray<Ty_> mS( static_cast<Ty_>(0), 3*3);

    mS[1] = -v[2]; mS[2] = v[1]; mS[3] = v[2]; 
    mS[5] = -v[0]; mS[6] = -v[1]; mS[7] = v[0];

    return mS;
  }


template< class Ty_> inline
valarray<Ty_> Col( const valarray<Ty_>& M, const size_t& n1, const size_t& n2)
  {
    // 2D matrix coulmns
    if ( !n2) {
      valarray<Ty_> m( M.size());
      size_t n2( M.size()/n1);          // Asserting: n1>0

      for( size_t r = 0; r < n1; r++) 
        m[slice(r,n2,n1)] = valarray<Ty_>( M[slice(r*n2,n2,1)]);

      return m;
      }

    // 3D matrix columns
    else {
      valarray<Ty_> m( M.size());
      size_t n12( n1*n2), n3( M.size()/n12);  // Asserting: n1>0 && n2>0
      size_t p, r;

      for( p = 0; p < n3; p++)
        for( r = 0; r < n1; r++) 
          m[ slice( p*n12 + r, n2, n1)] = valarray<Ty_>( M[slice( p*n12 + r*n2, n2, 1)]);

      return m;
      }
  }


namespace obsolete {

  template< class Ty_> inline
  const valarray<Ty_>& MatVect( valarray<Ty_>& vR, const valarray<Ty_>& mA, const valarray<Ty_>& vB)
    {
      size_t nC( vB.size()), nR( mA.size()/nC);   // Assert C>0

      vR.resize( nR);
      for( size_t r = 0; r < nR; r++)
        vR[r] = Dot( valarray<Ty_>(mA[slice(r*nC,nC,1)]), vB);

      return vR;
    }

  template< class Ty_> inline
  valarray<Ty_> MatVect( const valarray<Ty_>& mA, const valarray<Ty_>& vB)
    {
      valarray<Ty_> vR;
      MatVect( vR, mA, vB);
      return vR;
    }

  template< class Ty_> inline
  const valarray<Ty_>& MatMat( valarray<Ty_>& mR, const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nrA, size_t ncB)
    {
      size_t ncA( mA.size()/nrA); // Assert nrA>0
      valarray<Ty_> vR, vT(ncA);

      mR.resize( nrA*ncB);
      for( size_t c = 0; c < ncB; c++) {
        vT = mB[ slice(c,ncA,ncB)];
        mR[ slice(c,nrA,ncB)] = MatVect( vR, mA, vT);
        }

      return mR;
    }

  template< class Ty_> inline
  valarray<Ty_> MatMat( const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nrA, size_t ncB)
    {
      valarray<Ty_> mR;
      MatMat( mR, mA, mB, nrA, ncB);
      return mR;
    }

} // namespace obsolete


template< class Ty_> inline
const valarray<Ty_>& Tensor( valarray<Ty_>& mR, const valarray<Ty_>& vA, const valarray<Ty_>& vB)
  {
    return MM( mR, vA, vB, vA.size(), vB.size());
  }

template< class Ty_> inline
valarray<Ty_> Tensor( const valarray<Ty_>& vA, const valarray<Ty_>& vB)
  {
    valarray<Ty_> mR;
    MM( mR, vA, vB, vA.size(), vB.size());

    return mR;
  }


template< class Ty_> inline
valarray<Ty_> Ones( size_t nR, size_t nC)
  {
    return valarray<Ty_>( 1, nR*nC);
  }


inline
valarray<size_t> ones( const size_t& nR, const size_t& nC)
  {
    return Ones<size_t>( nR, nC);
  }


template< class Ty_> inline
valarray<Ty_> Zeros( size_t nR, size_t nC)
  {
    return valarray<Ty_>( static_cast<Ty_>(0), nR*nC);
  }


inline 
valarray<size_t> zeros( size_t nR, size_t nC)
  {
    return Zeros<size_t>( nR, nC);
  }


template< class Ty_> inline
valarray<Ty_> Eye( const size_t& n)
  {
    valarray<Ty_> mE;
    Eye( mE, n);

    return mE;
  }


template< class Ty_> inline
void Eye( valarray<Ty_>& mE, const size_t& n)
  {
    mE.resize( n*n, Ty_(0));

    for ( size_t k = 1; k <= n; k++) 
      mE[(k-1)*n+k-1] = Ty_(1);
  }


inline
valarray<size_t> eye( const size_t& n)
  {
    return Eye<size_t>( n);
  }


inline
void eye( valarray<size_t>& mE, const size_t& n)
  {
    Eye( mE, n);
  }


template< class Ty_> inline
valarray<Ty_> Rep( const Ty_& V, size_t nR, size_t nC)
  {
    return valarray<Ty_>( V, nR*nC);
  }


template< class Ty_> inline 
void RepRvec( valarray<Ty_>& mR, const valarray<Ty_>& vT, size_t nR, size_t nC)
  {
    size_t nc( vT.size());

    mR.resize( nR*nC*nc);
    for ( size_t i = 0; i < nR*nC; i++)
      mR[slice(i*nc,nc,1)] = vT;
  }


template< class Ty_> inline 
valarray<Ty_> RepRvec( const valarray<Ty_>& vT, size_t nR, size_t nC)
  {
    valarray<Ty_> mR;
    RepRvec( mR, vT, nR, nC);

    return mR;
  }


template< class Ty_> inline 
void RepCvec( valarray<Ty_>& mR, const valarray<Ty_>& vT, size_t nR, size_t nC)
  {
    size_t i, nr( vT.size()), m( nr*nC);

    mR.resize( nR*m);
    if ( mR.size()) {
      for ( i = 0; i < nC; i++)
        mR[slice(i,nr,nC)] = vT;

      for ( i = 0; i < nR; i++)
        copy( &mR[0], &mR[m], &mR[i*m]);
      }
  }


template< class Ty_> inline 
valarray<Ty_> RepCvec( const valarray<Ty_>& vT, size_t nR, size_t nC)
  {
    valarray<Ty_> mR;
    RepCvec( mR, vT, nR, nC);

    return mR;
  }


//TODO: expand to 3D
//TODO: combine MXc & MXr into this!!!!
template< class Ty_> inline
valarray<Ty_> RepMat( const valarray<Ty_>& mA, const size_t nR, const size_t nC, const size_t nrA)
  {
    size_t ncA( mA.size()/nrA);      // Asserting: nrA>0

    valarray<size_t> ri, rj;
    RepCvec( ri, range(1,nrA), nR, 1);
    RepCvec( rj, range(1,ncA), nC, 1);
    
    valarray<Ty_> mR;
    MX( mR, mA, ri, rj, nrA);

    return mR;
  }


template< class Ty_> inline
size_t RepMat( valarray<Ty_>& mB, const valarray<Ty_>& mA, const size_t nR, const size_t nC, const size_t nrA)
  {
    mB.resize( nR*nC*mA.size());
    mB = RepMat( mA, nR, nC, nrA);

    return nR*nrA;
  }


inline
valarray<size_t> DiagIx( const long m, const long n, const long k)
  {
    valarray<size_t> ix( m+n);
    long r, i, j;

    for ( r = 0l, i = Max(0l,-k); i < Min( m, n-k); i++)
      for ( j = Max(i+k,0l); j <= Min(i+k,n-1); j++) 
        ix[ r++] = i*n + j;

    return ix[slice(0,r,1)];
  }


inline
size_t DiagLength( const long nrA, const long ncA, const long kDiag)
  {
    size_t lm, nL( 0);

    if ( -nrA < kDiag && kDiag < ncA) {
      lm = Min( nrA, ncA);
      nL = lm;

      if ( kDiag) {
        // m > n
        if ( nrA>ncA) { 
          if ( kDiag>0)
            nL = lm - kDiag;
          else if ( kDiag < (ncA-nrA))
            nL = lm + ( kDiag + nrA - ncA);
          }

        // m <= n
        else {          
          if ( kDiag<0)
            nL = lm + kDiag;
          else if ( kDiag > (ncA-nrA))
            nL = lm - ( kDiag + nrA - ncA);
          }
        }
      }

    return nL;
  }

template< class Ty_> inline 
size_t DiagV( valarray<Ty_>& vD, const valarray<Ty_>& mA, const size_t nrA, const long k)
  {
    vD.resize( 0);

    if ( mA.size() && nrA) {
      valarray<size_t> vxD( DiagIx( nrA, mA.size()/nrA, k));    //Asserting: nrA>0
      vD.resize( vxD.size());
      vD = mA[vxD];
      }

    return vD.size();
  }

template< class Ty_> inline 
valarray<Ty_> DiagV( const valarray<Ty_>& mA, const size_t nrA, const long k)
  {
    valarray<Ty_> vD;
    DiagV( vD, mA, nrA, k);

    return vD;
  }


template< class Ty_> inline
size_t Diag( valarray<Ty_>& mD, const valarray<Ty_>& vD, const long k)
  {
    size_t nR( 0);
    
    if ( vD.size()) {
      nR = vD.size() + abs(k);
      mD.resize( nR*nR, Ty_(0));
      mD[ DiagIx( nR, nR, k)] = vD;
      }

    else 
      mD.resize( nR);

    return nR;
  }


template< class Ty_> inline
valarray<Ty_> Diag( const valarray<Ty_>& vD, const long k)
  {
    valarray<Ty_> mD( 0);
    Diag( mD, vD, k);

    return mD;
  }


template<class Ty_> inline
void Diff( valarray<Ty_>& vD, const valarray<Ty_>& vX)
  {
    SetM( vD, vX);
    adjacent_difference( &vD[0], &vD[vD.size()], &vD[0]);
  }


template<class Ty_> inline
valarray<Ty_> Diff( const valarray<Ty_>& vX)
  {
    valarray<Ty_> vD( vX);
    adjacent_difference( &vD[0], &vD[vD.size()], &vD[0]);

    return vD;
  }


//TODO: make obsolete all functions based on Fortran convention!
//TODO: add the 3rd dimension too, ie. B = A( rI, rJ, rK)
template< class Ty_> inline
const valarray<Ty_>& MX( valarray<Ty_>& mB, const valarray<Ty_>& mA, const valarray<size_t>& mI, const valarray<size_t>& mJ, size_t nrA, size_t nrI, size_t nrJ)
  {
    size_t ncA( mA.size()/nrA);                   // Asserting: nrA>0
    size_t ni( mI.size()), nj( mJ.size());

    // Fortran ordering
    valarray<size_t> vI(ni), vJ(nj);
    if ( nrI > 1) vI = Trp( mI, nrI); else vI = mI;  //TODO: seems too expensive!
    if ( nrJ > 1) vJ = Trp( mJ, nrJ); else vJ = mJ;

    // Base-1 indexing
    vI -= 1; 
    vJ -= 1;

    // Replicate by linear indexing
    valarray<size_t> vL( ni*nj);
    for( size_t i = 0; i < ni; i++) fill_n( &vL[i*nj], nj, ncA*vI[i]);
    
    mB.resize( vL.size()); 
    mB = mA[ vL + RepCvec( vJ, ni, 1)];

    return mB;
  }


//TODO: make obsolete all functions based on Fortran convention!
template< class Ty_> inline
valarray<Ty_> MX( const valarray<Ty_>& mA, const valarray<size_t>& mI, const valarray<size_t>& mJ, size_t nrA, size_t nrI, size_t nrJ)
  {
    valarray<Ty_> mB;
    MX( mB, mA, mI, mJ, nrA, nrI, nrJ);

    return mB;
  }


//TODO: obsolete this, replace with gslice2() !
inline
gslice MGS( const size_t& naC, const size_t& nOff, const size_t& nR, const size_t& nC)
  {
    valarray<size_t> nL(2), nS(2);

    nL[0] = nR;   nL[1] = nC;
    nS[0] = naC; nS[1] = 1;
    
    return gslice( nOff, nL, nS); 
  }


inline
gslice gslice2( const size_t& nD1, const size_t& nD2, const size_t& nC, const size_t& nOff)
  {
    valarray<size_t> vnL(2), vnS(2);

    vnL[0] = nD1; vnL[1] = nD2;
    vnS[0] = nC; vnS[1] = 1;
    
    return gslice( nOff, vnL, vnS);
  }


template< class Ty_> inline
valarray<Ty_> MG( const valarray<Ty_>& mA, const size_t& naR, const size_t& nOff, const size_t& nR, const size_t& nC)
  {
    return valarray<Ty_>( mA[ MGS(mA.size()/naR, nOff, nR, nC)]);
  }

//TODO: make this obsolete, than rename MXc0 as MXc!
template< class Ty_> inline
const valarray<Ty_>& MXc( valarray<Ty_>& mB, const valarray<Ty_>& mA, const valarray<size_t>& mK, const size_t nrA)
  {
    size_t ncA( mA.size()/nrA);      // Asserting: nrA>0
    size_t nK( mK.size());

    valarray<size_t> ix( mK - static_cast<size_t>(1));
    mB.resize( nrA*nK);

    for( size_t k = 0; k < nrA; k++) {
      mB[ slice(k*nK,nK,1)] = mA[ix];
      ix += ncA;
      }

    return mB;
  }


//TODO: make this obsolete, than rename MXc0 as MXc!
template< class Ty_> inline
valarray<Ty_> MXc( const valarray<Ty_>& mA, const valarray<size_t>& mK, const size_t nrA)
  {
    valarray<Ty_> mB;
    MXc( mB, mA, mK, nrA);

    return mB;
  }


template< class Ty_> inline
valarray<Ty_> MXc0( const valarray<Ty_>& mA, const valarray<size_t>& vi, const size_t nrA)
  {
    valarray<Ty_> mB;
    MXc( mB, mA, vi + (size_t)1, nrA);  //TODO: remove this inefficiency!

    return mB;
  }


template< class Ty_> inline
const valarray<Ty_>& MXc0( valarray<Ty_>& mB, const valarray<Ty_>& mA, const valarray<size_t>& vi, const size_t nrA)
  {
    MXc( mB, mA, vi + (size_t)1, nrA);  //TODO: remove this inefficiency!
    return mB;
  }

template< class Ty_> inline
valarray<Ty_> MXc0( const valarray<Ty_>& mA, const size_t nK, const size_t nrA)
  {
//DSG: old code!
#if 0
    const size_t ncA( mA.size()/nrA);
    return valarray<Ty_>( mA[ slice(nK,nrA,ncA)]);

//DSG: new code!
#else
  valarray<Ty_> vR;
  MXc0( vR, mA, nK, nrA);
  return vR;

#endif
  }

template< class Ty_> inline 
const valarray<Ty_>& MXc0( valarray<Ty_>& vR, const valarray<Ty_>& mA, const size_t nK, const size_t nrA)
  {
    vR.resize( nrA);
    vR = mA[ slice(nK,nrA,mA.size()/nrA)];

    return vR;
  }


//TODO: make this obsolete! Than rename MXr0 as MXr!
template< class Ty_> inline
const valarray<Ty_>& MXr( valarray<Ty_>& mB, const valarray<Ty_>& mA, const valarray<size_t>& mK, size_t nrA)
  {
    size_t ncA( mA.size()/nrA),       // Asserting: nrA>0
           nK( mK.size());

    mB.resize( ncA*nK);
    for( size_t k = 0; k < nK; k++)
      mB[slice(k*ncA,ncA,1)] = mA[slice( (mK[k]-1)*ncA,ncA,1)];

    return mB;
  }


//TODO: make this obsolete! Than rename MXr0 as MXr!
template< class Ty_> inline
valarray<Ty_> MXr( const valarray<Ty_>& mA, const valarray<size_t>& mK, size_t nrA)
  {
    valarray<Ty_> mB;
    MXr( mB, mA, mK, nrA);

    return mB;
  }


template< class Ty_> inline
const valarray<Ty_>& MXr0( valarray<Ty_>& mB, const valarray<Ty_>& mA, const valarray<size_t>& vi, size_t nrA)
  {
    mB.resize(0);

    if ( mA.size() && vi.size() && nrA) {
      const size_t ncA( mA.size()/nrA);
      const size_t nK( vi.size());

      mB.resize( ncA*nK);
      for( size_t k = 0; k < nK; k++)
        mB[slice(k*ncA,ncA,1)] = mA[slice( (vi[k])*ncA,ncA,1)];
      }

    return mB;
  }


template< class Ty_> inline
valarray<Ty_> MXr0( const valarray<Ty_>& mA, const valarray<size_t>& vi, size_t nrA)
  {
    valarray<Ty_> mB;
    MXr0( mB, mA, vi, nrA);

    return mB;
  }


template< class Ty_> inline
const valarray<Ty_>& MXr0( valarray<Ty_>& vR, const valarray<Ty_>& mA, const size_t nR, size_t nrA)
  {
    const size_t ncA( mA.size()/nrA);

    vR.resize( ncA);
    vR = mA[slice(nR*ncA,ncA,1)];
    
    return vR;
  }


template< class Ty_> inline
valarray<Ty_> MXr0( const valarray<Ty_>& mA, const size_t nR, size_t nrA)
  {
    valarray<Ty_> vR;
    MXr0( vR, mA, nR, nrA);

    return vR;
  }


template< class Ty_> inline
const valarray<Ty_>& MSc( valarray<Ty_>& mA, const valarray<size_t>& vi, const valarray<Ty_>& mB)
  {
    if ( mA.size() && vi.size() && mB.size()) {
      const size_t ncB( vi.size());     // Assert: size(mB,2) == numel(vi)
      const size_t nrB( mB.size()/ncB);
      const size_t ncA( mA.size()/nrB); // Assert: size(mA,1) == size(mB,1)

      for ( size_t n = 0; n < ncB; n++)
        mA[ slice( vi[n], nrB, ncA)] = mB[slice( n, nrB, ncB)];
      }

    return mA;
  }


template< class Ty_> inline
const valarray<Ty_>& MSc( valarray<Ty_>& mA, const size_t ni, const valarray<Ty_>& vB)
  {
    if ( mA.size() && vB.size()) {
      const size_t nrB( vB.size());
      const size_t ncA( mA.size()/nrB); // Assert: size(mA,1) == size(mB,1)
      mA[ slice( ni, nrB, ncA)] = vB;
      }

    return mA;
  }


template< class Ty_> inline
const valarray<Ty_>& MSr( valarray<Ty_>& mA, const valarray<size_t>& vi, const valarray<Ty_>& mB)
  {
    if ( mA.size() && vi.size() && mB.size()) {
      const size_t nrB( vi.size());       // Assert: size(mB,1) == numel(vi)
      const size_t ncB( mB.size()/nrB);
      const size_t nrA( mA.size()/ncB);   // Assert: size(mA,2) == size(mB,2)

      for ( size_t n = 0; n < nrB; n++)
        mA[ slice( vi[n]*ncB, ncB, 1)] = mB[ slice( n*ncB, ncB, 1)];
      }

    return mA;
  }


template< class Ty_> inline
const valarray<Ty_>& MSr( valarray<Ty_>& mA, const size_t ni, const valarray<Ty_>& vB)
  {
    if ( mA.size() && vB.size()) {
      const size_t ncB( vB.size());
      mA[ slice( ni*ncB, ncB, 1)] = vB;
      }

    return mA;
  }


template< class Ty_> inline
const valarray<Ty_>& MDr( valarray<Ty_>& mA, const valarray<size_t>& vi, const size_t nrA)
  {
    if ( mA.size()>0 && nrA>0 && vi.size()) {
      const size_t ncA( mA.size()/nrA);

      valarray<size_t> viK( SetDiff( range(nrA), vi));
      const size_t nrK( viK.size());

      valarray<Ty_> mA_( nrK*ncA);
      for ( size_t n = 0; n < nrK; n++)
        mA_[ slice( n*ncA, ncA, 1)] = valarray<Ty_>( mA[ slice( viK[n]*ncA, ncA, 1)]);

      SetM( mA, mA_);
      }

    return mA;
  }


template< class Ty_> inline
const valarray<Ty_>& MDr( valarray<Ty_>& mA, const size_t ni, const size_t nrA)
  {
    if ( mA.size()>0 && nrA>0 && ni<nrA) {
      const size_t ncA( mA.size()/nrA);
      valarray<Ty_> mA_( (nrA-1)*ncA);

      size_t nr, nr2; 
      for ( nr = 0, nr2 = 0; nr < nrA; nr++)
        if (nr != ni) {
          mA_[ slice( nr2*ncA, ncA, 1)] = valarray<Ty_>( mA[ slice( nr*ncA, ncA, 1)]);
          nr2++;
          }

      SetM( mA, mA_);
      }

    return mA;
  }


template< class Ty_> inline
const valarray<Ty_>& MDc( valarray<Ty_>& mA, const valarray<size_t>& vi, const size_t nrA)
  {
    if ( mA.size()>0 && nrA>0 && vi.size()) {
      const size_t ncA( mA.size()/nrA);

      valarray<size_t> viK( SetDiff( range(ncA), vi));
      const size_t nrK( viK.size());

      valarray<Ty_> mA_( nrA*nrK);
      for ( size_t n = 0; n < nrK; n++)
        mA_[ slice( n, nrA, nrK)] = valarray<Ty_>( mA[ slice( viK[n], nrA, ncA)]);

      SetM( mA, mA_);
      }

    return mA;
  }


template< class Ty_> inline
const valarray<Ty_>& MDc( valarray<Ty_>& mA, const size_t ni, const size_t nrA)
  {
    if ( mA.size()>0 && nrA>0 && ni<nrA) {
      const size_t ncA( mA.size()/nrA);   // Assert: nrA > 0
      valarray<Ty_> mA_( nrA*(ncA-1));

      size_t nr, nr2; 
      for ( nr = 0, nr2 = 0; nr < ncA; nr++)
        if (nr != ni) {
          mA_[ slice( nr2, nrA, ncA-1)] = valarray<Ty_>( mA[ slice( nr, nrA, ncA)]);
          nr2++;
          }

      SetM( mA, mA_);
      }

    return mA;
  }


//TODO: fixed nc=0 case, when matrix was not changed, although all items should
//      have been removed. Test all implications!!!!!!!!!!!!!!!!!!!!!!!!!
template< class Ty_> inline
const valarray<Ty_>& MRc( valarray<Ty_>& mA, const size_t nc, const size_t nrA)
  {
    if ( mA.size() && nrA) {
      const size_t ncA( mA.size()/nrA);

      // Remove columns
      if ( !nc) {
        mA.resize(0);
        }

      if ( nc < ncA) {
        SetM( mA, MXc0( mA, range(nc), nrA));
        }

      // Add columns
      else if ( nc > ncA) {
        valarray<Ty_> mA_( nrA*nc);
        mA_[ gslice2(nrA,ncA,nc)] = mA;
        SetM( mA, mA_);
        }
      }

    return mA;
  }


//TODO: fixed error when nr==0. Matrix was not resized as expected to empty,
//      but kept the full size!
//      This was detected from MapStr, which returned full mapping with all zeros
//      when nothing was matched! 
//      Test the implications in other places!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
template< class Ty_> inline
const valarray<Ty_>& MRr( valarray<Ty_>& mA, const size_t nr, const size_t nrA)
  {
    if ( mA.size() && nrA) {
      const size_t ncA( mA.size()/nrA);

      // Remove rows
      if ( !nr) {
        mA.resize(0);
        }

      else if ( nr < nrA) {
        SetM( mA, MXr0( mA, range(nr), nrA));
        }

      // Add rows
      else if ( nr > nrA) {
        valarray<Ty_> mA_( nr*ncA);
        mA_[ gslice2(nrA,ncA,ncA)] = mA;
        SetM( mA, mA_);
        }
      }

    return mA;
  }


template< class ContR_, class ContA_, class ContB_> inline
void Join1( ContR_& mR, const ContA_& mA, const ContB_& mB, const size_t ncA)
  {
    size_t nrA( mA.size()/ncA),     // Asserting: ncA>0
           nrB( mB.size()/ncA);

    mR.resize(0);
    if ( nrA>0 && nrB>0) {
      size_t n, na( mA.size());
      mR.resize( (nrA+nrB)*ncA);
      for( n = 0; n < mA.size(); n++) mR[n] = mA[n];
      for( n = 0; n < mB.size(); n++) mR[na+n] = mB[n];
      }

    else if ( nrA){
      mR.resize( mA.size());
      mR = mA;
      }
      
    else if ( nrB) {
      mR.resize( mB.size());
      mR = mB;
      }

  }


template< class ContR_, class ContA_, class ContB_> inline
ContR_ Join1( const ContA_& mA, const ContB_& mB, const size_t ncA)
  {
    ContR_ mR;
    Join1( mR, mA, mB, ncA);

    return mR;
  }


template< class Ty_> inline
const valarray<Ty_>& join1( valarray<Ty_>& mR, const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nc)
  {
    size_t nrA( mA.size()/nc);    //Asserting: nc>0
    size_t nrB( mB.size()/nc);
    
    mR.resize(0);

    if ( nrA>0 && nrB>0) {
      mR.resize( (nrA+nrB)*nc);
      gslice gA( MGS( nc, 0, nrA, nc));
      gslice gB( MGS( nc, nc*nrA, nrB, nc));
      mR[gA] = mA;
      mR[gB] = mB;
      }

    else if ( nrA) {
      mR.resize( mA.size());
      mR = mA;
      }
      
    else if ( nrB) {
      mR.resize( mB.size());
      mR = mB;
      }

    return mR;
  }


template< class Ty_> inline
valarray<Ty_> join1( const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nc)
  {
    valarray<Ty_> mR;
    join1( mR, mA, mB, nc);

    return mR;
  }


template< class Ty_> inline
const valarray<Ty_>& join2( valarray<Ty_>& mR, const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nrA)
  {
    mR.resize(0);

    if ( nrA) {
      const size_t ncA( mA.size()/nrA);
      const size_t ncB( mB.size()/nrA);
      
      // Neither matrix is empty
      if ( ncA>0 && ncB>0) {
        mR.resize( (ncA+ncB)*nrA);
        gslice gA( MGS( ncA+ncB, 0, nrA, ncA));
        gslice gB( MGS( ncA+ncB, ncA, nrA, ncB));
        mR[gA] = mA;
        mR[gB] = mB;
        }

      // Matrix mB is empty
      else if ( ncA) {
        mR.resize( mA.size());
        mR = mA;
        }

      // Matrix mA is empty
      else if ( ncB) {
        mR.resize( mB.size());
        mR = mB;
        }
      }

    return mR;
  }


template< class Ty_> inline
valarray<Ty_> join2( const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nrA)
  {
    valarray<Ty_> mR;
    join2( mR, mA, mB, nrA);

    return mR;
  }


template< class ContR_, class ContA_, class ContB_> inline
void Join2( ContR_& mR, const ContA_& mA, const ContB_& mB, const size_t nrA)
  {
    mR.resize(0);

    if ( nrA) {
      const size_t ncA( mA.size()/nrA);
      const size_t ncB( mB.size()/nrA);

      if ( ncA>0 && ncB>0) {
        size_t m, n, ncR(ncA+ncB);
        mR.resize( nrA*ncR);

        for( m = 0; m < nrA; m++) {
          for( n = 0; n < ncA; n++) mR[ m*ncR + n] = mA[ m*ncA + n];
          for( n = 0; n < ncB; n++) mR[ m*ncR + ncA + n] = mB[ m*ncB + n];
          }
        }

      else if ( ncB) {
        mR.resize( mB.size());
        mR = mB;
        }

      else if ( ncA) {
        mR.resize( mA.size());
        mR = mA;
        }
      }

  }


template< class ContR_, class ContA_, class ContB_> inline
ContR_ Join2( const ContA_& mA, const ContB_& mB, const size_t nrA)
  {
    ContR_ mR;
    Join2( mR, mA, mB, nrA);

    return mR;
  }


//TODO: make this obsolete! Use dlaruv instead.
template< class Ty_> inline
Ty_ matRandUnit( void) 
  {
    return static_cast<Ty_>( rand()/(RAND_MAX+1.));
  }


//TODO: make this obsolete!
template< class Ty_> inline
const valarray<Ty_>& LoadRand( valarray<Ty_>& v, size_t N)
  {
    v.resize( N);
    generate( &v[0], &v[N], matRandUnit<Ty_>);

    return v;
  }


//TODO: make this obsolete!
template< class Ty_> inline
const valarray<Ty_>& LoadRand( valarray<Ty_>& v, size_t N, const Ty_& A, const Ty_& B)
  {
    v.resize(N);
    return v = (B-A)*LoadRand( v, N) + valarray<Ty_>( A, N);
  }


//----------------------------------------------
//DSG: make the LoadInc(), LoadRange(), LoadIncRange() functions obsolete ...
//     Replace them with the Range, range, Interval, interval family.

//TODO: make this obsolete! Use Range, range instead.
template< class Ty_> inline
const valarray<Ty_>& LoadInc( valarray<Ty_>& v, size_t N, const Ty_& start, const Ty_& stride)
  {
    if ( N) {
      v.resize( N, stride); v[0] = start;
      partial_sum( &v[0], &v[N], &v[0]);
      }
    else
      v.resize(0);

    return v;
  }

//TODO: make this obsolete! Use Range, range instead.
template< class Ty_> inline
valarray<Ty_> LoadInc( size_t N, const Ty_& start, const Ty_& stride)
  {
    valarray<Ty_> v;
    LoadInc( v, N, start, stride);

    return v;
  }


template< class Ty_> inline
const valarray<Ty_>& LoadRange( valarray<Ty_>& v, size_t N, const Ty_& A, const Ty_& B)
  {
    if ( N >= 2 ) {
      v.resize( N, (B-A)/(N-1));

      v[0] = A;
      partial_sum( &v[0], &v[N-1], &v[0]);
      v[N-1] = B;
      }

    return v;
  }

//TODO: make this obsolete, use Range/range instead.
template< class Ty_> inline
valarray<Ty_> LoadRange( size_t N, const Ty_& A, const Ty_& B)
  {
    valarray<Ty_> v;
    LoadRange( v, N, A, B);

    return v;
  }


template< class Ty_> inline
const valarray<Ty_>& Range( valarray<Ty_>& vR, const Ty_& start, const size_t& steps, const Ty_& inc)
  {
    vR.resize( steps, inc);
    if ( steps > 0) {
      vR[0] = start;
      partial_sum( &vR[0], &vR[steps], &vR[0]);
      }

    return vR;
  }


template< class Ty_> inline
valarray<Ty_> Range( const Ty_& start, const size_t& steps, const Ty_& inc)
  {
    valarray<Ty_> vR;
    Range( vR, start, steps, inc);

    return vR;
  }


template< class Ty_> inline
valarray<Ty_> Range( const size_t& nSteps)
  {
    valarray<Ty_> vR;
    Range( vR, Ty_(0), nSteps, Ty_(1));

    return vR;
  }


template< class Ty_> inline
size_t Range( valarray<Ty_>& vR, const size_t& nSteps)
  {
    Range( vR, Ty_(0), nSteps, Ty_(1));
    return vR.size();
  }


template< class Ty_> inline
valarray<size_t> range( const Ty_& start, const size_t& steps, const long& inc)
  {
    valarray<size_t> vR( inc, steps);
    
    if ( steps > 0) {
      vR[0] = start;
      partial_sum( &vR[0], &vR[steps], &vR[0]);
      }

    return vR;
  }


//TODO: make this better a template specialization !?
template< class Ty_> inline
size_t range( valarray<size_t>& vR, const Ty_& start, const size_t& steps, const long& inc)
  {
    vR.resize( steps);
    vR = range( start, steps, inc);

    return vR.size();
  }


inline
size_t range( valarray<size_t>& vR, const size_t& nSteps)
  { 
    range( vR, 0, nSteps, 1);
    return vR.size();
  }


inline
valarray<size_t> range( const size_t& nSteps)
  { 
    valarray<size_t> vR;
    range( vR, 0, nSteps, 1);

    return vR;
  }


template< class Ty_> inline
size_t countInRange( const valarray<Ty_>& vT, const pair<Ty_,Ty_>& prRange)
  {
    //valarray<bool> vb( prRange.first <= vT && vT < prRange.second); // not accepted by gcc!
    valarray<bool> vb( valarray<bool>( vT >= prRange.first) && 
                       valarray<bool>( vT < prRange.second));
    return count_if( &vb[0], &vb[vb.size()], bind2nd( equal_to<bool>(), true));
  }


//DSG: this is obsolete. Remove all references, inclusive this function.
template< class Ty_> inline 
const valarray<Ty_>& Interval( valarray<Ty_>& vR, const Ty_& begin, const Ty_& end, const Ty_& inc)
  {
    vR.resize( 0);

    #if defined(_WINDOWS)
      #pragma warning( push)
      #pragma warning( disable : 4146)
    #endif
    
    if (  numeric_limits<Ty_>::is_signed && 
          ( inc < -numeric_limits<Ty_>::epsilon() || numeric_limits<Ty_>::epsilon() < inc) ||
          numeric_limits<Ty_>::epsilon() < inc ) {
      long steps( long( double(end - begin)/inc));      //TODO: redesign this with typeid(). Double cast is needed(?) if typeid(Ty_)!={"doubl"|"float"}...
      if ( steps >= 0) Range( vR, begin, ++steps, inc);
      }

    #if defined(_WINDOWS)
      #pragma warning( pop)
    #endif

    return vR;
  }


//DSG: this is obsolete. Remove all references, inclusive this function.
template< class Ty_> inline
valarray<Ty_> Interval( const Ty_& begin, const Ty_& end, const Ty_& inc)
  {
    valarray<Ty_> vR;
    Interval( vR, begin, end, inc);

    return vR;
  }


template< class Ty_, class Ty2_> inline 
valarray<size_t> interval( const Ty_& begin, const Ty2_& end, const long& inc)
  {
    valarray<size_t> vR( 0);

    if ( inc) {
      long n( long(end - begin)/inc);
      if ( n >= 0) {
        vR.resize( ++n);
        vR = range( begin, n, inc); //Assert all values are in a valid range
        }
      }

    return vR;
  }


template< class Ty_> inline
Ty_ ClampToZero( const Ty_& v, const Ty_& eps)
  {
    return ( abs(v) < eps) ? 0 : v;
  }


//OBSOLETE: use MeshGrid2(). Replace all instances with MeshGrid2(), remove this version 
//  than rename all to keep the MeshGrid() name!
template< class Ty_> inline
valarray<size_t> MeshGrid( valarray<Ty_>& mX, valarray<Ty_>& mY, 
                           const valarray<Ty_>& vx, const valarray<Ty_>& vy)
  {
    size_t nx( vx.size()), ny( vy.size());

    mX.resize( nx*ny); 
    mX = MX( vx, ones(ny,1), LoadInc<size_t>(nx,1), 1, ny, 1);
    
    mY.resize( nx*ny); 
    mY = MX( vy, LoadInc<size_t>(ny,1), ones(nx,1), ny, 1, 1);

    valarray<size_t> vD(2);
    vD[0] = ny;   // rows
    vD[1] = nx;   // columns
    
    return vD;
  }


//OBSOLETE: use MeshGrid2(). Replace all instances with MeshGrid2(), remove this version 
//  than rename all to keep the MeshGrid() name!
template< class Ty_> inline
valarray<size_t> MeshGrid( valarray<Ty_>& mX, valarray<Ty_>& mY, valarray<Ty_>& mZ, 
                           const valarray<Ty_>& vx, const valarray<Ty_>& vy, const valarray<Ty_>& vz)
  {
    size_t nx( vx.size()), ny( vy.size()), nz( vz.size());
    size_t np( nx*ny);

		valarray<Ty_> mX0( MX( vx, ones(ny,1), LoadInc<size_t>(nx,1), 1, ny, 1));
		RepCvec( mX, mX0, nz, 1);

    valarray<Ty_> mY0( MX( vy, LoadInc<size_t>(ny,1), ones(nx,1), ny, 1, 1));
    RepCvec( mY, mY0, nz, 1);

    mZ.resize( mX.size());
    for( size_t i = 0; i < nz; i++) 
      mZ[ slice(i*np,np,1)] = vz[i];

    valarray<size_t> vD(3);
    vD[0] = ny;   // rows
    vD[1] = nx;   // columns
    vD[2] = nz;   // pages
  
    return vD;
  }


template< class Ty_> inline
void MeshGrid2( valarray<Ty_>& mX, valarray<Ty_>& mY, valarray<Ty_>& mZ, valarray<size_t>& vDim, 
                const valarray<Ty_>& vx, const valarray<Ty_>& vy, const valarray<Ty_>& vz)
  {
    const size_t nx( vx.size()), ny( vy.size()), nz( vz.size());
    const size_t np( nx*ny);

    valarray<Ty_> mX0, mY0;
    MX( mX0, vx, ones(ny,1), range(1,nx), 1, ny, 1);
    MX( mY0, vy, range(1,ny), ones(nx,1), ny, 1, 1);

    RepCvec( mX, mX0, nz, 1);
    RepCvec( mY, mY0, nz, 1);

    mZ.resize( mX.size());
    for( size_t i = 0; i < nz; i++) 
      mZ[ slice(i*np,np,1)] = vz[i];

    vDim.resize(3);
    vDim[0] = ny;   // rows
    vDim[1] = nx;   // columns
    vDim[2] = nz;   // pages
  }


template< class Ty_> inline
void MeshGrid2( valarray<Ty_>& mX, valarray<Ty_>& mY, valarray<size_t>& vDim, 
                const valarray<Ty_>& vx, const valarray<Ty_>& vy)
  {
    const size_t nx( vx.size()), ny( vy.size());

    MX( mX, vx, ones(ny,1), range(1,nx), 1, ny, 1);
    MX( mY, vy, range(1,ny), ones(nx,1), ny, 1, 1);

    vDim.resize(2);
    vDim[0] = ny;   // rows
    vDim[1] = nx;   // columns
  }


template< class Ty_> inline
const valarray<Ty_>& Kron( valarray<Ty_>& mK, const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nrA, size_t nrB)
  {
    mK.resize( mA.size()*mB.size());
    if ( mK.size() > 0) {
      size_t ncA( mA.size()/nrA), ncB( mB.size()/nrB);

      valarray<size_t> ia, ib, ja, jb;
      MeshGrid( ia, ib, range(1,nrA), range(1,nrB)); 
      MeshGrid( ja, jb, range(1,ncA), range(1,ncB));
      mK = MX(mA,ia,ja,nrA,nrB,ncB)*MX(mB,ib,jb,nrB,nrB,ncB);
      }

    return mK;
  }


template< class Ty_> inline
valarray<Ty_> Kron( const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nrA, size_t nrB)
  {
    valarray<Ty_> mK;
    Kron( mK, mA, mB, nrA, nrB);

    return mK;
  }


template< class TL, class TR> inline
const valarray<TL>& Conv( valarray<TL>& L, const valarray<TR>& R)
  {
    L.resize( R.size());

    for( size_t i = 0; i < R.size(); i++)
      L[i] = static_cast<TL>( R[i]);

    return L;
  }

template< class TL, class TR> inline
const valarray<TL>& Conv( valarray<TL>& L, const vector<TR>& R)
  {
    L.resize( R.size());

    for( size_t i = 0; i < R.size(); i++)
      L[i] = static_cast<TL>( R[i]);

    return L;
  }

template< class TL, class TR> inline
const vector<TL>& Conv( vector<TL>& L, const valarray<TR>& R)
  {
    L.resize( R.size());

    for( size_t i = 0; i < R.size(); i++)
      L[i] = static_cast<TL>( R[i]);

    return L;
  }

template< class TL, class TR> inline
valarray<TL> Conv( const valarray<TR>& R)
  {
    valarray<TL> L;
    return Conv( L, R);
  }


inline
valarray<double> Ceil( valarray<double>& v)
  {
    valarray<double> vR( v.size());
    vdCeil( (MKL_INT)v.size(), &v[0], &vR[0]);

    return vR;
  }


inline
valarray<double> Floor( valarray<double>& v)
  {
    valarray<double> vR( v.size());
    vdFloor( (MKL_INT)v.size(), &v[0], &vR[0]);

    return vR;
  }


template< class Ty_> inline
const valarray<Ty_>& DivCeil( valarray<Ty_>& vnR, const valarray<Ty_>& vnA, const Ty_& nD)
  {
    valarray<double> vdA;

    vnR.resize( vnA.size());
    Conv<double,Ty_>( vdA, vnA);
    
    vdA /= (double)nD;                // Asserting: nD>0
    Conv<Ty_,double>( vnR, Ceil(vdA));

    return vnR;
  }


inline void PeriodicIndex( valarray<size_t>& vI, const valarray<size_t>& vJ, const size_t nP)
  {
    valarray<size_t> vJ_( vJ + (size_t)1);
    DivCeil( vI, vJ_, nP);
    vI -= 1;
  }


inline valarray<size_t> PeriodicIndex( const valarray<size_t>& vJ, const size_t nP)
  {
    valarray<size_t> vI;
    PeriodicIndex( vI, vJ, nP);
    return vI;
  }


template< class Ty_> inline
Ty_ Mod( const Ty_& k, const Ty_& N)
  {
    Ty_ r(k);

    r %= N;
    if ( r < Ty_(0)) r += N;

    return r;
  }


template< class Ty_> inline
valarray<Ty_> Mod( const valarray<Ty_>& k, const Ty_& N)
  {
    valarray<Ty_> r(k), rp( Ty_(0), k.size());

    r %= N;
    valarray<bool> b( r<Ty_(0));
    rp[b] = Ty_(N);
    r = r + rp;

    return r;
  }


template< class Ty_> inline
valarray<Ty_> SumD1( const valarray<Ty_>& mA, const size_t nrA)
  {
    size_t ncA( mA.size()/nrA);       // Asserting: nrA>0
    valarray<Ty_> vs(ncA), va(nrA);

    for( size_t i = 0; i < ncA; i++) {
      va = mA[ slice(i,nrA,ncA)];
      vs[i] = va.sum();
      }

    return vs;
  }


template< class Ty_> inline
valarray<Ty_> SumD2( const valarray<Ty_>& mA, const size_t nrA)
  {
    size_t ncA( mA.size()/nrA);      // Asserting: nrA>0
    valarray<Ty_> vs(nrA), va(ncA);

    for( size_t i = 0; i < nrA; i++) {
      va = mA[slice(i*ncA,ncA,1)];
      vs[i] = va.sum();
      }

    return vs;      
  }


template< class Ty_> inline
valarray<Ty_> Sum( const valarray<Ty_>& mA, const size_t nrA, const size_t nDim)
  {
    switch( nDim) {

      case 1: {
        return SumD1( mA, nrA);
        break;
        }

      case 2: {
        return SumD2( mA, nrA);
        break;
        }

      //TODO: add more dimensions!!!
      //case 3: ...

      default: 
        return valarray<Ty_>(0);
      }
  }


template<class Ty_> inline
size_t Sum( valarray<Ty_>& vR, const valarray<Ty_>& mA, const size_t nrA, const size_t nDim)
  {
    SetM( vR, Sum( mA, nrA, nDim));
    return vR.size();
  }


template<class Ty_> inline
Ty_ ValRange( const valarray<Ty_>& vR)
  {
    valarray<Ty_>& vR_( const_cast< valarray<Ty_>&>( vR));

    if ( vR.size()) {
      double* pMin( min_element( &vR_[0], &vR_[vR_.size()]));
      double* pMax( max_element( &vR_[0], &vR_[vR_.size()]));

      return (*pMax - *pMin);
      }

    return Ty_(0);
  }


template< class Ty_> inline
Ty_ Trace( const valarray<Ty_>& mA, const size_t nrA)
  {
    return DiagV( mA, nrA).sum();  // sum of the main diagonal items
  }


//TODO: rewrite this to be more efficient!
inline
valarray<double> CycVector( const valarray<double>& vD, const valarray<long>& k)
  {
    valarray<long> v( k);
    v -= 1l;

    valarray<long> vI( Mod( v, static_cast<long>( vD.size())));
    return vD[Conv<size_t,long>(vI)];
  }


inline
double CycVector( const valarray<double>& vD, const long k)
  {
    long I = Mod( k-1l, static_cast<long>(vD.size()));
    return vD[ static_cast<size_t>(I)];
  }


template< class Ty_> inline
bool Find( valarray<size_t>& vI, const valarray<Ty_>& vA)
  {
    return SetM( vI, Find( vA)).size()>0;
  }


template< class Ty_> inline
valarray<size_t> Find( const valarray<Ty_>& vA)
  {
    valarray<size_t> ix( range( vA.size()));
    return valarray<size_t>( ix[ vA != Ty_(0)]);
  }


inline
valarray<size_t> Find( const valarray<double>& vA, const double dEps)
  {
    valarray<size_t> ix( range( vA.size()));
    return valarray<size_t>( ix[ abs(vA) > dEps]);
  }


//TODO: check this if ok
template< class Ty_> inline
size_t Find( valarray<size_t>& ixR, valarray<size_t>& ixC, const valarray<Ty_>& vA, const size_t nrA)
  {
    size_t nix, ncA( vA.size()/nrA);       // Asserting: nrA>0
    valarray<size_t> ix( Find( vA));

    if ( nix = ix.size()) { 
      ixC.resize( nix); 
      ixC = Mod( ix, ncA);

      ixR.resize( nix); 
      ixR = ix/static_cast<size_t>( ncA);
      }

    return nix;
  }


//FIX: this is badly conditioned on the machine precision boundary !!!
//  Never sort directly up to float/double machine precision, as it is platform 
//  and machine dependent.
inline 
ptrdiff_t FindNearest( const valarray<double>& vD, const double dT)
  {
    valarray<double> diff( abs( vD - dT));

    //FIX: make serach for minimum only within the finite double precision,
    //  otherways very small machine dependent differences might result
    //  in machine dependent results!    
    double* pMin = min_element( &diff[0], &diff[diff.size()]);
    
    return (pMin - &diff[0]);             
  }


inline 
ptrdiff_t FindReal( const valarray<double>& vD, const double dT, const double eps2)
  {
    //TODO: see note about machine precision in FindNearest
    ptrdiff_t ix( FindNearest( vD, dT));
    double d( vD[ix]);

    return ( abs(d-dT) > eps2 ) ? -1 : ix;
  }


template< class Ty_> inline
size_t MapVectors( valarray<size_t>& miM, valarray<size_t>& viAu, valarray<size_t>& viBu,
                   const valarray<Ty_>& mA, const valarray<Ty_>& mB, const size_t nR, const double eps2)
  {
    size_t nA( mA.size()/nR);   // ASSERT: nR == size(mA,1) == size( mB,1) > 0
    size_t nB( mB.size()/nR);
    double dEps2( eps2 < 0.0 ? Eps( 1.0) : eps2);

    miM.resize(0);
    viAu.resize(0);
    viBu.resize(0);
    
    // mA is empty
    if ( !nA && nB) {
      range( viBu, 0, nB);
      return 0;
      }

    // mB is empty
    else if ( nA && !nB) {
      range( viAu, 0, nA);
      return 0;
      }

    // Compare only nonempty sets
    else if ( nA && nB) {      
      //ASSERT: size(mA) == nR*nA, size(mB) == nR*nB!
      // All pairwise distances
      valarray<Ty_> mD( Kron( mA, Ones<Ty_>( 1, nB), nR, 1));
      mD -= RepMat( mB, 1, nA, nR);
      
      // Find the indices of matching ones
      valarray<Ty_> vN(mD.size()/nR);
      if ( nR>1)
        vN = Norm( mD, nR);
      else
        vN = abs( mD);

      valarray<size_t> viT;
      if ( Find( viT, valarray<bool>( vN <= Ty_(dEps2)))) { 
        size_t nT( viT.size());

        // Indices into mA, resp. mB
        valarray<size_t> vkA( viT / nB);
        valarray<size_t> vkB( viT % nB);

        miM.resize( 2*nT);
        miM[ slice(0,nT,1)] = vkA;
        miM[ slice(nT,nT,1)] = vkB;

        // Indices of unique vectors in mA, resp. mB
        SetDiff( viAu, range( nA), Unique( vkA));
        SetDiff( viBu, range( nB), Unique( vkB));
        }

      else {
        range( viAu, 0, nA);
        range( viBu, 0, nB);
        }
      }

    return miM.size();
  }


template< class Ty_> inline
size_t MapVectors( valarray<size_t>& miM, const valarray<Ty_>& mA, const valarray<Ty_>& mB, const size_t nR, const double eps2)
  {
    valarray<size_t> viAu, viBu;
    return MapVectors( miM, viAu, viBu, mA, mB, nR, eps2);
  }


template< class Ty_> inline
size_t MapVectors( valarray<size_t>& miM, const valarray<Ty_>& mA, const size_t nR, const double eps2)
  {
    valarray<size_t> vi;
    return MapVectors( miM, vi, vi, mA, mA, nR, eps2);
  }


template<class Ty_> inline
size_t MapColVectors( valarray<size_t>& miMA, valarray<size_t>& viU, valarray<size_t>& viR,
                      const valarray<Ty_>& mA, const size_t nrA)
  {
    miMA.resize(0);
    viU.resize(0);
    viR.resize(0);

    if ( nrA) {
      // Self mapping
      valarray<size_t> miM(0);
      MapVectors( miM, mA, nrA);
      const size_t ncM( miM.size()/2);

      // Column mappings, without selfies
      valarray<size_t> viD( miM[slice(0,ncM,1)]);
      viD -= miM[slice(ncM,ncM,1)];
      valarray<size_t> vid( Find( valarray<bool>( viD != size_t(0))));
      MXc0( miMA, miM, vid, 2);
      const size_t ncMA( miMA.size()/2);

      // Unique column indices of mA
      SetDiff<size_t>( viU, range(mA.size()/nrA), miMA[slice(0,ncMA,1)]);

      // Basis of the repeated columns
      set<size_t> setB, setC, setR;
      for( size_t k = 0; k < ncMA; k++) {
        if ( setC.find( miMA[k]) == setC.end()) {
          setB.insert( miMA[k]);
          setC.insert( miMA[ncMA+k]);
          }
        else
          setR.insert( k);
        }

      // Remove complementary mirror mappings
      viD.resize( setR.size());
      copy( setR.begin(), setR.end(), &viD[0]);
      MDc( miMA, viD, 2);
      
      viR.resize( setB.size());
      copy( setB.begin(), setB.end(), &viR[0]);
      }

    return ( viU.size() + viR.size());
  }


template<class Ty_> inline
size_t MapRowVectors( valarray<size_t>& miMA, valarray<size_t>& viUA, valarray<size_t>& viB,
                      const valarray<Ty_>& mA, const size_t nrA)
  {
    miMA.resize(0);
    viUA.resize(0);
    viB.resize(0);

    if ( nrA)
      return MapColVectors( miMA, viUA, viB, Trp(mA,nrA), mA.size()/nrA);
    else
      return 0;
  }


template<class Ty_> inline
size_t MatColBasis( valarray<Ty_>& mB, const valarray<Ty_>& mA, const size_t nrA)
  {
    mB.resize(0);

    valarray<size_t> miMA(0), viU(0), viR(0);
    size_t nB = MapColVectors( miMA, viU, viR, mA, nrA);
    if ( nB) {
      set<size_t> icB;
      icB.insert( &viU[0], &viU[viU.size()]);
      icB.insert( &viR[0], &viR[viR.size()]);

      viR.resize( icB.size());
      copy( icB.begin(), icB.end(), &viR[0]);
      MXc0( mB, mA, viR, nrA);
      }

    return viR.size();
  }


template<class Ty_> inline
size_t MatRowBasis( valarray<Ty_>& mB, const valarray<Ty_>& mA, const size_t nrA)
  {
    mB.resize(0);

    valarray<size_t> miMA(0), viU(0), viR(0);
    size_t nB = MapRowVectors( miMA, viU, viR, mA, nrA);
    if ( nB) {
      set<size_t> icB;
      icB.insert( &viU[0], &viU[viU.size()]);
      icB.insert( &viR[0], &viR[viR.size()]);

      viR.resize( icB.size());
      copy( icB.begin(), icB.end(), &viR[0]);
      MXr0( mB, mA, viR, nrA);
      }

    return viR.size();
  }


template< class Ty_> inline
bool All( const Ty_& vB)
  {
    for( size_t n = 0; n < vB.size(); n++) 
      if ( !vB[n]) return false;

    return true;
  }


template< class Ty_> inline
bool Any( const Ty_& vB)
  {
    for( size_t n = 0; n < vB.size(); n++) 
      if ( vB[n]) return true;

    return false;
  }


template< class Ty_> inline
bool Any( size_t& nI, const Ty_& vB)
  {
    for( nI = 0; nI < vB.size(); nI++) 
      if ( vB[nI]) return true;

    return false;
  }


template< class Ty_> inline
size_t SetDiff( valarray<Ty_>& vD, const valarray<Ty_>& vA, const valarray<Ty_>& vB)
  {
    const size_t nA( vA.size());
    vD.resize( 0);

    if ( nA) {
      valarray<Ty_> vA_( vA);
      sort( &vA_[0], &vA_[vA_.size()]);

      const size_t nB( vB.size());

      // Both mA and mB have values
      if ( nB) {
        valarray<Ty_> vB_( vB);
        sort( &vB_[0], &vB_[vB_.size()]);

        vA_ = Unique( vA_);
        vB_ = Unique( vB_);

        // vD <- ( vA \ vB)
        valarray<Ty_> vD_( nA); // most nA items
        Ty_* pD = set_difference( &vA_[0], &vA_[ vA_.size()], &vB_[0], &vB_[ vB_.size()], &vD_[0]);

        // Differences
        if ( pD > &vD_[0]) {
          size_t nD( static_cast<size_t>( pD - &vD_[0]));
          vD.resize( nD);
          vD = vD_[slice(0,nD,1)];
          }

        }

      // mB is empty
      else {
        Unique( vD, vA);  // all unique items of vA
        }
      }

    return vD.size();
  }


template< class Ty_> inline
valarray<Ty_> SetDiff( const valarray<Ty_>& vA, const valarray<Ty_>& vB)
  {
    valarray<Ty_> vD(0);
    SetDiff( vD, vA, vB);

    return vD;
  }


//TODO: simplify this, use sort & set_difference
inline 
bool SetDiff( valarray<size_t>& vDiff, const valarray<double>& vA, const valarray<double>& vB, const double eps2)
  {
    valarray<size_t> vU2A;
    valarray<double> dif;
    size_t k, n, nA( vA.size());

    if ( nA > 0 ) {
      vU2A.resize( nA);
      dif.resize( vB.size());
      for ( k = 0, n = 0; k < vA.size(); k++) {
        dif = vB - vA[k];
        if ( All( abs( dif)>eps2)) 
          vU2A[n++] = k;
        }

      vDiff.resize(n); vDiff = vU2A[slice(0,n,1)];
      }

    else 
      vDiff.resize(0);

    return vDiff.size()>0;
  }


//TODO: replace or reimplement this like SetDiff(). Source of error: overlapping in/out for set_difference.
inline 
bool SetDiffInt( valarray<size_t>& ixC, valarray<size_t>& ixU, const valarray<double>& vA, const valarray<double>& vB, const double eps2)
  {
    size_t na( vA.size());

    if ( na) {
      // A\B, 0-based, sorted
      SetDiff( ixU, vA, vB, eps2);  //TODO: fix SetDiff(...,eps) for machine precision!

      // A & B, 0-based, sorted
      valarray<size_t> ixc( LoadInc<size_t>(na,0,1));

      //FIX: Source and destination ranges for set_difference are not allowed to overlap !!!
      size_t* pDiffEnd = set_difference( &ixc[0], &ixc[ixc.size()], &ixU[0], &ixU[ixU.size()], &ixc[0]);
      if ( pDiffEnd > &ixc[0]) {
        size_t n = static_cast<size_t>(pDiffEnd - &ixc[0]);
        ixC.resize( n);
        ixC = ixc[ slice(0,n,1)];
        }

      else
        ixC.resize(0);

      return true;
      }

    else {
      ixC.resize(0);
      ixU.resize(0);

      return false;
      }
  }


//FIX: bad real sorting and const correctnes if called for real sets!
inline 
bool SetDiffInc( valarray<size_t>& vixDiff, valarray<size_t>& vixA, valarray<size_t>& vixB)
  {
    size_t na( vixA.size()), nb( vixB.size());
    size_t n;

    vixDiff.resize(0);

    if ( na > 0) {
      valarray<size_t> difAB( na);
      
      sort( &vixA[0], &vixA[na]);    //FIX: don't sort reals like this !!!
      sort( &vixB[0], &vixB[nb]);

      n = set_difference( &vixA[0], &vixA[vixA.size()], &vixB[0], &vixB[vixB.size()], &difAB[0]) - &difAB[0];
      if ( n) {
        vixDiff.resize( n);
        vixDiff = difAB[slice(0,n,1)];

        return true;
        }
      }

    return false;
  }


inline
pair< size_t, double> DiscreteIndex( const valarray<double>& vR, double r0)
  {
    size_t ix( MinIndex( valarray<double>( abs( vR - r0))));
    return make_pair( ix, vR[ix]);
  }


inline
void DiscreteRange( valarray<double>& vD, pair<size_t,size_t>& N, const valarray<double>& vR, double dS)
  {
    if ( vR[1] < vR[0]) {
      vD.resize(0);
      N = make_pair( 0, 0);
      }

    else {
      N.first = (size_t)floor( (vR[1] - vR[0])/abs(dS)) + 1;
      valarray<double> vr = LoadInc( N.first, vR[0], dS);
      
      N.second = MinIndex( valarray<double>( abs( vr)));
      vD.resize(2); vD = vR - vr[ N.second];
      }
  }


inline
size_t TriLen( size_t m, size_t n, size_t k)
  {
    if ( k<n) 
      return ( m >= n-k) ? (n-k)*(n-k+1)/2 : 
                           ( (n-k)*(n-k+1) - (n-m-k)*(n-m-k+1))/2;
    return 0;
  }


inline
valarray<size_t> TriU( long m, long n, long k)
  {
    valarray<size_t> ix( m*n);
    long r, i, j;

    for ( r = 0l, i = 0l; i < Min( m, n-k); i++)
      for ( j = Max(i+k,0l); j < n; j++) ix[ r++] = i*n + j;

    return ix[slice(0,r,1)];
  }


template< class Ty_> inline
valarray<Ty_>& TriU( valarray<Ty_>& mR, const valarray<Ty_>& mD, size_t nR, long k)
  {
    SetM( mR, mD);
    mR[ TriL( nR, mD.size()/nR, k-1)] = Ty_( 0);  // Asserting: nR>0

    return mR;
  }


template< class Ty_> inline
valarray<Ty_> TriU( const valarray<Ty_>& mD, size_t nR, long k)
  {
    valarray<Ty_> mR;
    return TriU( mR, mD, nR, k);
  }


inline
valarray<size_t> TriL( long m, long n, long k)
  {
    valarray<size_t> ix( m*n);
    long r, i, j;

    for ( r = 0l, i = Max(0l,-k); i < m; i++)
      for ( j = 0l; j <= Min(i+k,n-1); j++) 
        ix[ r++] = i*n + j;

    return ix[ slice(0,r,1)];
  }


template< class Ty_> inline
valarray<Ty_>& TriL( valarray<Ty_>& mR, const valarray<Ty_>& mD, size_t nR, long k)
  {
    SetM( mR, mD);
    mR[ TriU( nR, mD.size()/nR, k+1)] = Ty_( 0);    // Asserting: nR>0

    return mR;
  }


template< class Ty_> inline
valarray<Ty_> TriL( const valarray<Ty_>& mD, size_t nR, long k)
  {
    valarray<Ty_> mR;
    return TriL( mR, mD, nR, k);
  }


template< class Ty_> inline
void Mirroru( valarray<Ty_>& mAU, const size_t nR)
  {
    //Assert: mAU quadratic
    for( size_t n = 0; n < nR; n++)
      mAU[ slice( n*(nR+1), nR-n, nR)] = valarray<Ty_>( mAU[ slice( n*(nR+1), nR-n, 1)]);
  }


template< class Ty_> inline
size_t MirrorU( valarray<Ty_>& mAU, const valarray<Ty_>& mA, const size_t nR)
  {
    SetM( mAU, mA);
    Mirroru( mAU, nR);

    return mAU.size();
  }


template< class Ty_> inline
valarray<Ty_> MirrorU( const valarray<Ty_>& mA, const size_t nR)
  {
    valarray<Ty_> mAU( mA);
    Mirroru( mAU, nR);

    return mAU;
  }


template< class Ty_>
void Mirrorl( valarray<Ty_>& mAL, const size_t nR)
  {
    //Assert: mAL quadratic
    for( size_t n = 0; n < nR; n++)
      mAL[ slice( n*(nR+1), nR-n, 1)] = valarray<Ty_>( mAL[ slice( n*(nR+1), nR-n, nR)]);
  }


template< class Ty_> inline
size_t MirrorL( valarray<Ty_>& mAL, const valarray<Ty_>& mA, const size_t nR)
  {
    SetM( mAL, mA);
    Mirrorl( mAL, nR);

    return mAL.size();
  }


template< class Ty_> inline
valarray<Ty_> MirrorL( const valarray<Ty_>& mA, const size_t nR)
  {
    valarray<Ty_> mAL( mA);
    Mirrorl( mAL, nR);

    return mAL;
  }


template< class Ty_> inline
const valarray<Ty_>& CumSum( valarray<Ty_>& mS, const valarray<Ty_>& mT, size_t nR, size_t nDim)
  {
    size_t nC( mT.size()/nR), n;     // Asserting: nR>0

    mS.resize( mT.size());
    mS = mT;

    if ( nDim == 1 && nR>1)
      for ( n = 1; n < nR; n++)
        mS[slice(n*nC,nC,1)] += mS[slice((n-1)*nC,nC,1)];

    else // assert nDim=2
      for ( n = 0; n < nR; n++)
        partial_sum( &mS[n*nC], &mS[(n+1)*nC], &mS[n*nC]);

    return mS;
  }


template< class Ty_> inline
valarray<Ty_> CumSum( const valarray<Ty_>& mT, size_t nR, size_t nDim)
  {
    valarray<Ty_> mS;
    CumSum( mS, mT, nR, nDim);

    return mS;
  }


template< class Ty_> inline
size_t RemoveSorted( valarray<Ty_>& vD, const valarray<size_t>& vI)
  {
    size_t nR( 0);

    if ( vI.size() && vD.size()) {
      size_t nD( vD.size()), nI( vI.size());
      size_t n, i;

      for( i = 0, n = 0; n < nD; n++)
        if ( n < vI[i] || i >= nI)
          vD[n-nR] = vD[n];
        else {
          ++nR;
          if ( i < nI) ++i;
          }

      } // if

    return nR;   // Return the number of removed items
  }


template< class Ty_> inline
size_t Remove( valarray<Ty_>& vD, const valarray<size_t>& vI)
  {
    if ( vI.size() && vD.size()) {
      valarray<Ty_> vT( vD);
      valarray<size_t> vIqs( vI);
      ptrdiff_t nUq( 0);

      sort( &vIqs[0], &vIqs[vIqs.size()]);
      nUq = unique( &vIqs[0], &vIqs[vIqs.size()]) - &vIqs[0];
      
      size_t nR( 0);
      if ( nUq == (ptrdiff_t)vIqs.size())
        nR = RemoveSorted( vT, vIqs);
      else
        nR = RemoveSorted( vT, valarray<size_t>( &vIqs[0], nUq));
       
      size_t n( vD.size() - nR);
      vD.resize( n);
      vD = valarray<Ty_>( &vT[0], n);

      return nR;     //Nr. of items removed
      }

    return 0;
  }


template< class Ty_> inline
valarray<Ty_> Remove( const valarray<Ty_>& vD, const valarray<size_t>& vI)
  {
    valarray<Ty_> vR( vD);
    Remove( vR, vI);

    return vR;
  }


//NOTE: this function is platform dependent, when applied to reals!
//      Use Uniquev() instead! 
template< class Ty_> inline 
size_t Unique( valarray<Ty_>& vU, const valarray<Ty_>& vX)
  {
    vU.resize( 0);

    if ( vX.size()) {
      valarray<Ty_> vx( vX);
      ptrdiff_t nUq =0;

      sort( &vx[0], &vx[vx.size()]);
      nUq = unique( &vx[0], &vx[vx.size()]) - &vx[0];

      vU.resize( nUq);
      vU = valarray<Ty_>( &vx[0], nUq);
      }
    
    return vU.size();
  }


//NOTE: this function is platform dependent, when applied to reals!
//      Use Uniquev() instead!
template< class Ty_> inline
valarray<Ty_> Unique( const valarray<Ty_>& vX)
  {
    valarray<Ty_> vU;
    Unique( vU, vX);

    return vU;
  }


template< class Ty_> inline
const valarray<size_t>& Unique3D( valarray<size_t>& viU, const valarray<Ty_>& mP, const size_t nR, const double eps2)
  {
    viU.resize( 0);

    if ( mP.size() && nR) {
      size_t nM( mP.size()/nR), nU(0);
      viU.resize( nM);  // of maximum size

      valarray<size_t> viT( range(nM)), vt, miM;
      while ( viT.size() > 1) {
        // Reference vector and leftover ones
        viU[ nU++] = viT[0];
        MXc0( vt, viT, range(1,viT.size()-1), 1);

        // Compare reference vector to the left over ones,
        // than remove matches from the running index list
        MapVectors( miM, MXc0(mP,viT[0],3), MXc0(mP,vt,3), 3, eps2);  //TODO: this makes Unique3D work only for nR=3!
        if ( miM.size())                                              // Alternative1: fix this and rename back Unique3D
          SetDiff( viT, vt, MXc0(vt,MXr0(miM,1,2),1));                // Alternative1: fix this and rename back Unique3D
        else
          MDc( viT, 0, 1);  // delete first item
        }

        if ( viT.size())      // left over 
          viU[ nU++] = viT[0];

        MRc( viU, nU, 1);  // remove unused items
        }

    return viU;
  }


template< class Ty_> inline
valarray<size_t> Unique3D( const valarray<Ty_>& mP, const size_t nR, const double eps2)
  {
    valarray<size_t> viU;
    Unique3D( viU, mP, nR, eps2);

    return viU;
  }


template< class Ty_> inline
size_t IsNan( valarray<bool>& vB, const valarray<Ty_>& vX)
  {
    size_t nN( vX.size()), nNans(0);

    vB.resize( nN);
    if ( nN) {
      valarray<Ty_>& vX_( const_cast<valarray<Ty_>&>(vX));
      valarray<int> vN( nN);

      transform( &vX_[0], &vX_[nN], &vN[0], ptr_fun<double,int>( _isnan));
      vB = ( vN != 0);

      nNans = count( &vB[0], &vB[nN], true);
      }

    return nNans;
  }


template< class Ty_> inline
valarray<bool> IsNan( const valarray<Ty_>& vX)
  {
    valarray<bool> vB;
    IsNan( vB, vX);

    return vB;
  }


template< class Ty_> inline
size_t IsFinite( valarray<bool>& vB, const valarray<Ty_>& vX)
  {
    size_t nN( vX.size()), nFins(0);

    vB.resize( nN);
    if ( nN) {
      valarray<Ty_>& vX_( const_cast<valarray<Ty_>&>(vX));
      valarray<int> vN( nN);

      transform( &vX_[0], &vX_[nN], &vN[0], ptr_fun<double,int>( _finite));
      vB = ( vN != 0);

      nFins = count( &vB[0], &vB[nN], true);
      }

    return nFins;
  }


template< class Ty_> inline
valarray<bool> IsFinite( const valarray<Ty_>& vX)
  {
    valarray<bool> vB;
    IsFinite( vB, vX);

    return vB;
  }


/** Test if a real number vector is identic zero.*/
template< class Ty_> inline
bool IsZero( const valarray<Ty_>& vA, const double dEpsScale)
  {
    return All( abs( vA) <= Eps( 1.0)*dEpsScale);
  }


template< class Ty_> inline
bool IsEqual( const valarray<Ty_>& vA, const valarray<Ty_>& vB, const double dEpsScale)
  {
    if ( vA.size() == vB.size()) {
      valarray<Ty_> vMax;
      valarray<double> vEps;
      
      MaxAbs( vMax, vA, vB);
      Eps( vEps, vMax);
      vEps *= dEpsScale;

      valarray<bool> bL( vB-vEps <= vA);
      valarray<bool> bH( vA <= vB+vEps);

      return All( bL && bH);
      }

    return false;
  }


template< class Ty_> inline
bool IsSymmetric( const valarray<Ty_>& mX, const size_t n, const double dEpsScale)
  {
    size_t N( n);

    if ( N==0)
      N = static_cast<size_t>( sqrt( (float)mX.size()));

    if ( mX.size() == N*N)
      return IsEqual( mX, Trp( mX, N), dEpsScale);

    return false;
  }


template< class Ty_> inline
size_t Count( const valarray<Ty_>& vX, const Ty_& x)
  {
    valarray<Ty_>& vX_( const_cast<valarray<Ty_>&>(vX));
    return count( &vX_[0], &vX_[ vX_.size()], x);
  }


template< class Ty_> inline
size_t Countf( const valarray<Ty_>& vX, const Ty_& x, const double& dEpsScale)
  {
    valarray<bool> vbEq( abs( vX - x) <= Eps( 1.0)*dEpsScale);
    return Count( vbEq, true);
  }



namespace {

template< class Ty_>
class ClampToRange : public unary_function<Ty_, Ty_>
  {
    public:
      ClampToRange( const Ty_& rMin, const Ty_& rMax)
        : m_rMin( rMin), m_rMax( rMax) {}

      Ty_ operator()( const Ty_& v) {
        return Clamp( v, m_rMin, m_rMax);
        }

    private:
      const Ty_& m_rMin;
      const Ty_& m_rMax;
  };

} // noname


template< class Ty_> inline
size_t MapToRange( valarray<double>& vT, const valarray<Ty_>& vS, 
                   const valarray<double>& vSRange, const valarray<double>& vTRange, const bool bClamp)
  {
    double dSR( vSRange[1]-vSRange[0]),   // source range
           dTR( vTRange[1]-vTRange[0]);   // target range

    vT.resize( 0);
    if ( abs( dSR) > numeric_limits<double>::epsilon()) {
      double dU( dTR/dSR);
      Conv<double>( vT, vS);

      if ( bClamp) 
        transform( &vT[0], &vT[vT.size()], &vT[0], ClampToRange<double>( vSRange[0], vSRange[1]));

      vT *= dU;
      vT += vTRange[0] - dU*vSRange[0];
      }

    return vT.size();
  }


template< class Ty_> inline
size_t MapToRangeAuto( valarray<double>& vT, const valarray<Ty_>& vS, 
                       const valarray<double>& vTRange)
  {
    valarray<double> vSRange(2);
    vSRange[0] = MinVal( vS);
    vSRange[1] = MaxVal( vS);

    return MapToRange( vT, vS, vSRange, vTRange, false);
  }


template< class Ty_> inline
size_t MapToSigRange( valarray<double>& vT, const valarray<Ty_>& vS, 
                      const valarray<double>& vSRange, const bool bClamp)
  {
    const double vSigRange_[] = { -1, 1};
    const valarray<double> vSigRange( vSigRange_, SizeOfArray( vSigRange_));

    return MapToRange( vT, vS, vSRange, vSigRange, bClamp);
  }

template< class Ty_> inline
valarray<Ty_> MapToSigRange( const valarray<Ty_>& vS, const valarray<double>& vSRange, const bool bClamp)
  {
    valarray<double> vT;
    MapToSigRange( vT, vS, vSRange, bClamp);

    return vT;
  }


inline
void Sigmoid( valarray<double>& vR, const valarray<double>& vX, const valarray<double>& vRange, const double dA)
  {
    vR.resize( vX.size(), vRange[1] - vRange[0]);
    vR /= 1.0 + exp( -dA*vX);
    vR += vRange[0];
  }

inline
valarray<double> Sigmoid( const valarray<double>& vX, const valarray<double>& vRange, const double dA)
  {
    valarray<double> vR;
    Sigmoid( vR, vX, vRange, dA);

    return vR;
  }


template< class Ty_> inline
Ty_ Det3( const valarray<Ty_>& mA)
  {
    return   mA[0*3+0]*( mA[1*3+1]*mA[2*3+2] - mA[1*3+2]*mA[2*3+1]) 
           - mA[0*3+1]*( mA[1*3+0]*mA[2*3+2] - mA[2*3+0]*mA[1*3+2])
           + mA[0*3+2]*( mA[1*3+0]*mA[2*3+1] - mA[2*3+0]*mA[1*3+1]);
  }



}    //namespace mft

