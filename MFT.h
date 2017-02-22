/***************************************************************************
  MFT.h  - Math Foundation Templates
  -------------------------------------

  General porpose mathematical templates and algorithms bundled into the namespace: "mft".
  These utilities are meant for high-performance computation, so proper call and parameter check
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

#ifndef __MFT__H__
#define __MFT__H__

#include "extreme.h"
#include "streamprint.h"


namespace mft {

  using namespace std;


/** Size of a static C-array. */
#ifndef SizeOfArray
#define SizeOfArray( a)      ( sizeof(a)/sizeof((a)[0]))
#endif

// NOTE: the SetM() functions are mainly required for STL standards 
//   compliance under Linux.

/** Assign a valarray vA to a valarray vB. */
template< class Ty_> 
valarray<Ty_>& SetM( valarray<Ty_>& vB, const valarray<Ty_>& vA);

/** Assign STL vector vA to a vector vB. */
template< class Ty_>
vector<Ty_>& SetM( vector<Ty_>& vB, const vector<Ty_>& vA);


/** 
  Assign a subset vi of a valarray vA to another valarray vB.

  Params:
    - vA: source vector, [na,1]
    - vi: indices in the range 0:(na-1).

    - vB: on return the target vector with the vA[vi] values.
*/
template< class Ty_>
size_t SetM( valarray<Ty_>& vB, const valarray<Ty_>& vA, const valarray<size_t>& vi);

template< class Ty_>
valarray<Ty_> SetM( const valarray<Ty_>& vA, const valarray<size_t>& vi);


/** 
  Assign a C-vector to a valarray.

  Params:
    - pA: pointer to a vector, [na,1]
    - na: size of the C-vector
    
    - vA: on return contains the pA[0:(na-1)] values
*/
template< class Ty_>
valarray<Ty_>& SetM( valarray<Ty_>& vA, const Ty_* pA, const size_t na);

/**
  Assign the content of an STL set to a valarray.

  Params:
    - setA: source set

    - vA: on return the values of the setA set.

  Return value is a reference to vA.
*/
template< class Ty_>
valarray<Ty_>& SetM( valarray<Ty_>& vA, const set<Ty_>& setA);


/** Angle conversions between degrees and radians. */
template< class Ty_> 
double Deg2Rad( const Ty_& deg);

template< class Ty_> 
double Rad2Deg( const Ty_& rad);


/* Unit conversions. */
template<class Ty_> 
double ToAng( const Ty_& x);

template<class Ty_> 
double ToAu( const Ty_& x);


/** Compute log2(). */
template< class Ty_> 
Ty_ Log2( const Ty_& x);

template< class Ty_> 
valarray<Ty_> Log2( const valarray<Ty_>& v);


/** Test if all items are non-zero. */
template< class Ty_> 
bool All( const Ty_& vB);


/** Test if any item is non-zero. */
template< class Ty_> 
bool Any( const Ty_& vB);

/** 
  Test if any item is non-zero and return its index.
  Note: if all items are non-zero, nI==numel(vB).
*/
template< class Ty_> 
bool Any( size_t& nI, const Ty_& vB);

/** Absolute values. */
template< class Ty_> 
Ty_ Abs( const Ty_& v);


/** Round a scalar to the nearest integer. */ 
double Round( const double& d);
const valarray<double>& Round( valarray<double>& mR, const valarray<double>& m);

template< class Ty_>
const valarray<Ty_>& Round( valarray<Ty_>& m);

template< class Ty_>
const valarray<Ty_>& Round( valarray<Ty_>& mR, const valarray<Ty_>& m);


/** 
  Round real values to a given number of decimals or epsylon precision.

  Params:
    - m: values to be rounded ( in place)
    - prec: decimal precision. It can be an integer number of decimal positions or a 
            real subunit epsylon. If prec is an epsylon, it should be only negative powers of 10!

  Return value is a convenience reference to m.
*/
template< class Ty_> 
const valarray<Ty_>& RoundDec( valarray<Ty_>& m, const Ty_& prec);

template< class Ty_> 
valarray<Ty_> RoundDec( const valarray<Ty_>& m, const Ty_& prec);


/** Ceil on a vector.*/
MATH_API valarray<double> Ceil( valarray<double>& v);


/** Floor on a vector.*/
MATH_API valarray<double> Floor( valarray<double>& v);


/** 
  Real divison and ceil for integral type vectors.

  Params:
    - vA: integral type numeric vector
    - nD: integral divisor

  Return value is the vR vector with ceil( vA/nD).
*/
template< class Ty_>
const valarray<Ty_>& DivCeil( valarray<Ty_>& vR, const valarray<Ty_>& vA, const Ty_& nD);


/** 
  Periodic indices.

  Params:
    - vJ: (0-based) index vector
    - nP: division (period)

  Return value is the vI (0-based) index vector: ceil( vJ/nP).
*/
void PeriodicIndex( valarray<size_t>& vI, const valarray<size_t>& vJ, const size_t nP);
valarray<size_t> PeriodicIndex( const valarray<size_t>& vJ, const size_t nP);


/** Round toward 0, ie. return the integer part. */
MATH_API valarray<double> Fix( valarray<double>& v);
MATH_API double Fix( const double& dV);

template< class Ty_>
short Sign( const Ty_& v); 


/** Type conversion. */
template< class _TyL,class _TyR>
const valarray<_TyL>& Conv( valarray<_TyL>& L, const valarray<_TyR>& R);

template< class TL, class TR> inline
const valarray<TL>& Conv( valarray<TL>& L, const vector<TR>& R);

template< class TL, class TR> inline
const vector<TL>& Conv( vector<TL>& L, const valarray<TR>& R);

template< class _TyL,class _TyR>
valarray<_TyL> Conv( const valarray<_TyR>& R);


/** 
  Modulo operator: |a|.

  Params:
    - k: numerical value/vector
    - N: modulo value

  Return the value/vector of (k mod N).
*/
template< class Ty_>
Ty_ Mod( const Ty_& k, const Ty_& N);

template< class Ty_>
valarray<Ty_> Mod( const valarray<Ty_>& k, const Ty_& N);


/** 
  Greatest common divisor of two values: g = gcd( a, b).
  See also: Gcd.
*/
template< class Ty_>
Ty_ gcd( const Ty_& a, const Ty_& b);


/** 
  Greatest common divisor of n values: g = gcd( v1, v2, ...).
  See also: Gcd.
*/
template< class Ty_>
Ty_ gcd( const valarray<Ty_>& vX);

/**
  Binary implementation of gcd with the method of Stein.
*/
template< class Ty_>
Ty_ gcd_Stein( const Ty_& u, const Ty_& v);


/** 
  Extended greatest common divisor of two values: [g,c,d] = gcd( a, b).

  Params:
    - a, b: integral typed numerical values.
    - vR: solution vector [g,c,d] of the Diophantine equation a*c + b*d = g.

  See also: Extended Euclidean in D.E. Knuth, "The Art of Computer Programming", Vol.2
*/
template< class Ty_> 
const valarray<Ty_>& Gcd( valarray<Ty_>& vR, const Ty_& a, const Ty_& b);

template< class Ty_> 
Ty_ Gcd( const Ty_& a, const Ty_& b);


/** 
  Matrix and vector norms.

  Params:
    - mA: matrix, [nrA,N]
    - nrA: number of rows of mA
    - nNorm: norm type to be computed:
        1 = 1-norm, largest column sum
        2 = largest singular value ( default)
        3 = infinity norm, largest row sum
        4 = Euclidean or Frobenius norm

  Return value is the requested norm of the mA matrix.
*/
template< class Ty_> 
Ty_ MNorm( const valarray<Ty_>& mA, const size_t& nrA, const size_t& nNorm =2);


/** 
  Euclidean norm of column vectors.

  Params:
    - vX: column/row vector of size n.
    - mA: matrix of column vectors, [nrA,ncA]
    - nrA: number of rows of mA

    - vN: return vector of size ncA, formed by the Euclidean norms
          for each mA[:,i] vector, i=1..ncA.
    
  Return value for the scalar version is the Euclidean norm of vX.
*/

template< class Ty_> 
Ty_ Norm( valarray<Ty_>& vX);

MATH_API double Norm( valarray<double>& vX);        // based on MKL cblas_dnrm2


template< class Ty_> 
Ty_ Norm( const valarray<Ty_>& vX);

MATH_API double Norm( const valarray<double>& vX);  // based on MKL cblas_dnrm2


template< class Ty_> 
valarray<Ty_> Norm( const valarray<Ty_>& mA, const size_t& nrA);

template< class Ty_> 
const valarray<Ty_>& Norm( valarray<Ty_>& vN, const valarray<Ty_>& mA, const size_t& nrA);


/** 
  Squared Euclidean norm of a set of column vectors.

  Params:
    - mA: column vectors, side by side as a matrix, [nrA,ncA]
    - nrA: number of rows of mA
    - vN2: return vector of size ncA, formed by the squared Euclidean norms
           for each mA[:,i] vector, i=1..ncA.
*/
const valarray<double>& Norm2( valarray<double>& vN2, const valarray<double>& mA, const size_t& nrA);


/**
  Compute the angle between two vectors.

  Params:
    - vA: first vector, [N,1]
    - vB: second vector, [N,1]
    - bDeg: if true the returned angle is in degrees, otherways radians.

  Return value is the angle between too N dimensional vectors, in the range [0..pi].
*/
MATH_API double AbsAngle( const valarray<double>& vA, const valarray<double>& vB, const bool bDeg=false);


/** Clamp near-zero values.*/
MATH_API void ClampToZero( valarray<double>& vOut, const valarray<double>& vIn, const double dEps);

template< class Ty_>
Ty_ ClampToZero( const Ty_& v, const Ty_& eps);


/**
  Unit length scaling of vectors.

  Params:
    - vX: non-zero length vector, [N,1]
    - mV: set of non-zero length column vectors, [N,M]
    - nr: number of rows in mV, ie: nr=M
*/
template< class Ty_>
const valarray<Ty_>& Unit( valarray<Ty_>& vX);

template< class Ty_>
valarray<Ty_> Unit( const valarray<Ty_>& vX);

template< class Ty_>
valarray<Ty_> Unit( const valarray<Ty_>& mV, const size_t& nR);

template< class Ty_>
const valarray<Ty_>& Unit( valarray<Ty_>& mU, const valarray<Ty_>& mV, const size_t& nR);


/** Dot product of vectors. */
MATH_API double Dot( double* va, double* vb, const size_t n);
MATH_API double Dot( const valarray<double>& va, const valarray<double>& vb);
MATH_API double Dot( const valarray<double>& v);
MATH_API double Dot( valarray<double>& va, valarray<double>& vb);
MATH_API void Dot( valarray<double>& vD, const valarray<double>& mA, const valarray<double>& mB, const size_t nc);
MATH_API valarray<double> Dot( const valarray<double>& mA, const valarray<double>& mB, const size_t nc);


/** 
  Cross product of 2 vectors: vC = cross( vA, vB).

  Params:
    - pA: address of the left vector, [3]
    - pB: address of the right vector, [3]
    - pC: address of the output vector, [3]. It ha s to be preallocated by the caller.
*/
template< class Ty_>
void Cross( Ty_* pC, const Ty_* pA, const Ty_* pB);

template< class Ty_>
valarray<Ty_> Cross( const valarray<Ty_>& vA, const valarray<Ty_>& vB);

template< class Ty_>
const valarray<Ty_>& Cross( valarray<Ty_>& vC, const valarray<Ty_>& vA, const valarray<Ty_>& vB);

template< class Ty_>
valarray<Ty_> Cross( const slice_array<Ty_>& vA, const slice_array<Ty_>& vB);


/** 
  Cross product for N column vectors:  mC(:,i) = mA(:,i)*mB(:,i) with i=[1:N].
  
  Params:
    - mA: left-set of 3D column vectors, [3,N]
    - mB: right-set of 3D column vectors, [3,N]
    
  Return value is the [3,N] matrix, formed by column vectors ci = cross( ai, bi), 1=[1,N].
*/
template< class Ty_>
valarray<Ty_> Cross( const valarray<Ty_>& mA, const valarray<Ty_>& mB, const size_t& nc);

template< class Ty_>
const valarray<Ty_>& Cross( valarray<Ty_>& mC, const valarray<Ty_>& mA, const valarray<Ty_>& mB, const size_t& nc);


/** Mixed product of 3D vectors, c = (v1xv2)*v3. */
template< class Ty_>
Ty_ Mixed( const valarray<Ty_>& mV);

template< class Ty_>
Ty_ Mixed( const valarray<Ty_>& v1, const valarray<Ty_>& v2, const valarray<Ty_>& v3);

valarray<double> Mixed( const valarray<double>& mA, const valarray<double>& mB, const valarray<double>& mC, const size_t& nc);


/** Tensor product of vectors. */
template< class Ty_>
const valarray<Ty_>& Tensor( valarray<Ty_>& mR, const valarray<Ty_>& vA, const valarray<Ty_>& vB);

template< class Ty_>
valarray<Ty_> Tensor( const valarray<Ty_>& vA, const valarray<Ty_>& vB);


/** Transpose rectangular matrices. */
valarray<size_t> Trp( const size_t& nR, const size_t& nC);

template< class Ty_> 
valarray<Ty_> Trp( const valarray<Ty_>& mA, const size_t& nR);

template< class Ty_> 
void Trp( valarray<Ty_>& mT, const valarray<Ty_>& mA, const size_t& nR);


/** Skew symmetric matrix of a vector. */
template< class Ty_> 
valarray<Ty_> Skew( const valarray<Ty_>& v);


/** 
  Colon operator: A(:).
  Transform a 2D/3D matrix with Fortran ordering, into a column vector.

  Params:
    - A: matrix of [n1,N] or [n1,n2,P] dimensions
    - n1, n2: first and second dimension. For 2D matrices use only n1, but for 3D matrices
          specify both (n1,n2)=(rows,cols) parameters.

  Return the elements of A in a "column" vector format, of size [n1*N], respective [n1*n2*P].
*/
template< class Ty_> 
valarray<Ty_> Col( const valarray<Ty_>& A, const size_t& n1, const size_t& n2 =0);


/** 
  General matrix-vector multiplication: y = A*x
  
  Params:
    - mM: matrix, [M,N]
    - vX: vector, [N,1]
    - vY: vector, [M,1]

  See also: BLAS DGEMV.
*/
MATH_API const valarray<double>& MV( valarray<double>& vY, const valarray<double>& mM, const valarray<double>& vX);
MATH_API valarray<double> MV( const valarray<double>& mM, const valarray<double>& vX);


/** 
  General matrix-matrix multiplication: C = A*B.
  
  Params:
    - mA: left matrix, [nrA,*]
    - mB: right matrix, [*,ncB]
    - mC: matrix product mA*mB, computed with BLAS DGEMM. Size [nrA,ncB].

  Return value is a reference, respective the product matrix mC.
*/
MATH_API const valarray<double>& MM( valarray<double>& mC, const valarray<double>& mA, const valarray<double>& mB, const size_t& nrA, const size_t& ncB);
MATH_API valarray<double> MM( const valarray<double>& mA, const valarray<double>& mB, const size_t& nrA, const size_t& ncB);


/** 
  Generate a matrix of ones: ones(nr,nc).

  Params:
    - nr: number of row replicas
    - nc: number of column replicas

  Return value is a matrix of [nr,nc] ones.
*/
template< class Ty_>
valarray<Ty_> Ones( size_t nr =1, size_t nc =1);


/** Create a matrix of 1s (of type size_t). */
valarray<size_t> ones( const size_t& nr =1, const size_t& nc =1);


/** 
  Generate a matrix of zeros: zeros(nr,nc).

  Params:
    - nR: number of row replicas
    - nC: number of column replicas

  Return value is a matrix of [nR,nC] zeros.
*/
template< class Ty_>
valarray<Ty_> Zeros( size_t nR =1, size_t nC =1);

valarray<size_t> zeros( size_t nR =1, size_t nC =1);


/** Generate an identity matrix: eye(n).

  Params:
    - nL: length of the diagonal.

  Return value is an identity matrix, [nL,nL].
*/
template< class Ty_>
void Eye( valarray<Ty_>& mE, const size_t& nL =1);

template< class Ty_>
valarray<Ty_> Eye( const size_t& nL =1);

valarray<size_t> eye( const size_t& nL =1);
void eye( valarray<size_t>& mE, const size_t& nL =1);


/** Determinant of a 3x3 matrix: det(A). */
template< class Ty_>
Ty_ Det3( const valarray<Ty_>& mA);


/** 
  Replicate a scalar value.

  Params:
    - a: scalar value
    - nR, nC: number of rows and columns

  Return value is a matrix of size [nR,nC].
*/
template< class Ty_>
valarray<Ty_> Rep( const Ty_& a, size_t nR =1, size_t nC =1);


/** 
  Replicate row and column vectors: [ v,v; v,v...].

  Params:
    - vT: row [1,n] or column [n,1] vector
    - nR, nC: number of row and column replicas of vT

  Return value of RepRvec is a matrix of replicated row vectors, [nR, n*nC].
  Return value of RepCvec is a matrix of replicated column vectors, [n*nR,nC].
*/
template< class Ty_>
void RepRvec( valarray<Ty_>& mR, const valarray<Ty_>& vT, size_t nR =1, size_t nC =1);

template< class Ty_>
valarray<Ty_> RepRvec( const valarray<Ty_>& vT, size_t nR =1, size_t nC =1);

template< class Ty_>
void RepCvec( valarray<Ty_>& mR, const valarray<Ty_>& vT, size_t nR =1, size_t nC =1);

template< class Ty_>
valarray<Ty_> RepCvec( const valarray<Ty_>& vT, size_t nR =1, size_t nC =1);


/** 
  Replicate a 2D matrix: [ A,A; A,A...].

  Params:
    - mA: sample matrix, [nrA,ncA]
    - nR: number of rows to generate with the mA submatrix
    - nC: number of columns to generate with the mA submatrix
    - nrA: number of rows in  mA

  Return value is a matrix of [nR,nC] submatrices formed with mA, size [nR*].
*/
template< class Ty_>
valarray<Ty_> RepMat( const valarray<Ty_>& mA, const size_t nR, const size_t nC, size_t const nrA);


/** 
  Replicate a 2D matrix.

  Params:
    - mA: sample matrix to be replicated, [nrA,ncA]
    - nR: number of rows to generate with the mA submatrix
    - nC: number of columns to generate with the mA submatrix
    - nrA: rows of the mA matrix

    - mB: Generated output matrix

  Return value is the number of rows [nR*nrA] of mB.
*/
template< class Ty_>
size_t RepMat( valarray<Ty_>& mB, const valarray<Ty_>& mA, const size_t nR, const size_t nC, const size_t nrA);


/** 
  Diagonal element indices of a matrix.

  Params:
    - m, n: matrix dimensions, [m,n]
    - k: diagonal order. The main diagonal has order zero ( k=0). Subdiagonals have 
         negative order (k<0), supra diagonals have positive order ( k>0).

  Return value is the vector of linear indices of the k-order diagonal elements.
*/
MATH_API valarray<size_t> DiagIx( const long m, const long n, const long k);


/** 
  Length of the k-th diagonal of a matrix: numel( diag(A,k)).

  Params:
    - nR, nC: row and column dimensions of the matrix.
    - kDiag: the diagonal index. The main diagonal is kDiag=0,
          subdiagonals are kDiag<0, superdiagonals kDiag>0.

  Return value of DiagLength( nR, nC, kDiag) is the length of the kDiag
  diagonal of a matrix of size [nR,nC].
*/
MATH_API size_t DiagLength( const long nR, const long nC, const long kDiag);


/** 
  Diagonal elements of a matrix: diag(A).
  
  Params:
    - mX: data matrix, [nR,nC].
    - nR: rows of mX.
    - kDiag: number of diagonal. Main diagonal is kDiag=0, subdiagonals 
          are kDiag<0, superdiagonal kDiag>0.

    - vD: vector of the elements of mX on the k-th diagonal.

  Return value of DiagV( vD, mX, nR, kDiag) is the size of vD.
  Return value of DiagV( mX, nR, kDiag) is the vector of k-th diagonal elements in mX.
*/
template< class Ty_>
size_t DiagV( valarray<Ty_>& vD, const valarray<Ty_>& mX, const size_t nR, const long kDiag =0l);

template< class Ty_>
valarray<Ty_> DiagV( const valarray<Ty_>& mX, const size_t nR, const long kDiag =0l);


/** 
  Return an element from one of the diagonals of a matrix.

  Params:
    - mA: numeric matrix, [nrA,ncA]
    - nrA: number of rows of mA
    - n: element index of the k-th diagonal of mA
    - kDiag: diagonal order. See DiagIx() for more.

  If n is out of range, the return value is NaN.
*/
template< class Ty_>
Ty_ DiagV( const valarray<Ty_>& mA, const size_t nrA, const size_t n, const long kDiag =0l);


/** 
  Create a diagonal matrix from a vector: diag(v).

  Params:
    - vD: diagonal entries, [N,1].
    - k: the order of the diagonal. The main diagonal has an order of zero.

    - mD: diagonal matrix of size [M,M], with the k-order diagonal elements 
          taken from vD.

  The first function returns M or zero if vD is empty. 
  The second function returns the diagonal matrix or an empty matrix if fails.
*/
template< class Ty_>
size_t Diag( valarray<Ty_>& mD, const valarray<Ty_>& vD, const long k =0l);

template< class Ty_>
valarray<Ty_> Diag( const valarray<Ty_>& vD, const long k =0l);


/** 
  Compute discrete differences on consecutive vector elements.

  Params:
    - vX: vector of values to differenciate, [N,1]
    - vD: differences of adjacent values in vX, [N-1,1]
*/
template< class Ty_>
void Diff( valarray<Ty_>& vD, const valarray<Ty_>& vX);

template< class Ty_>
valarray<Ty_> Diff( const valarray<Ty_>& vX);


//TODO: make this obsolete!
/**
  Select matrix items specified by matrices of row and colum indices 
  with Fortran convention.
  
  Params:
    - mA: matrix A, [nrA,ncA].
    - mI: matrix of row indices, [nrI,ncI].
    - mJ: matrix of column indices, [nrJ,ncJ].

    - mR: matrix of the selected mA values of size [?]
*/
template< class Ty_>
const valarray<Ty_>& MX( valarray<Ty_>& mR, const valarray<Ty_>& mA, const valarray<size_t>& mI, 
                         const valarray<size_t>& mJ, size_t nrA, size_t nrI =1, size_t nrJ =1);

template< class Ty_>
valarray<Ty_> MX( const valarray<Ty_>& mA, const valarray<size_t>& mI, const valarray<size_t>& mJ, 
                  size_t nrA, size_t nrI =1, size_t nrJ =1);


//TODO: make this obsolete! Always use MXc0 with 0-based indices.
/** 
  Select columns from a matrix: mA(:,vi).

  Params:
    - mA: source matrix, [nrA,ncA]

    - vi: vector of column indices ( 1-based), [1,N].
        Valid values for vi are in the range of [1:ncA].
        Note: indices are 1-based.

    - nrA: number of rows of mA

    - mR: columns from mA selected by the vi indices, [naR,N].

  Return value is a reference to mR or a copy of the resulting matrix.
*/
template< class Ty_>
const valarray<Ty_>& MXc( valarray<Ty_>& mR, const valarray<Ty_>& mA, const valarray<size_t>& vi, const size_t nrA);

template< class Ty_>
valarray<Ty_> MXc( const valarray<Ty_>& mA, const valarray<size_t>& vi, const size_t nrA);


/** 
  Select a set of columns from a matrix: A(:,vi).

  Params:
    - mA: matrix of values, [nrA,ncA].
    - vi: vector of column indices into mA, [1,n]. 
        Valid values for vi are in the range [0:ncA-1].

    - mR: selected column values as a matrix of size [nrA,n].
  
  See also: MXc.
*/
template< class Ty_>
const valarray<Ty_>& MXc0( valarray<Ty_>& mR, const valarray<Ty_>& mA, const valarray<size_t>& vi, const size_t nrA);

template< class Ty_>
valarray<Ty_> MXc0( const valarray<Ty_>& mA, const valarray<size_t>& vi, const size_t nrA);


/** Return a column of a matrix: A(:,n). */
template< class Ty_>
valarray<Ty_> MXc0( const valarray<Ty_>& mA, const size_t ni, const size_t nrA);

template< class Ty_>
const valarray<Ty_>& MXc0( valarray<Ty_>& mR, const valarray<Ty_>& mA, const size_t ni, const size_t nrA);


//TODO: make this obsolete in favor of MXr0!
/** 
  Select a set of rows from a matrix: mA(vi,:).

  Params:
    - mA: matrix of values, [nrA,ncA].
    - vi: vector of row indices into mA, [1,n]. Valid values are in the
        range [1:nrA] by Fortran convention.

    - mR: matrix of the vi rows selected from mA, [n,ncA].
*/
template< class Ty_>
const valarray<Ty_>& MXr( valarray<Ty_>& mR, const valarray<Ty_>& mA, const valarray<size_t>& vi, size_t nrA);

template< class Ty_>
valarray<Ty_> MXr( const valarray<Ty_>& mA, const valarray<size_t>& vi, size_t nrA);


//TODO: after removing the function MXr reuse the name 'MXr' for 'MXr0'.
/** 
  Select a set of rows from a matrix: mA(vi,:).

  Params:
    - mA: matrix of values, [nrA,ncA].
    - viR: vector of row indices into mA, [1,n]. Valid values are in the
        range [0:nrA-1].
    - nR: one row index.
    - mR: matrix of the vi rows selected from mA, [n,ncA].
*/
template< class Ty_>
const valarray<Ty_>& MXr0( valarray<Ty_>& mR, const valarray<Ty_>& mA, const valarray<size_t>& viR, size_t nrA);

template< class Ty_>
valarray<Ty_> MXr0( const valarray<Ty_>& mA, const valarray<size_t>& viR, size_t nrA);

/** Select one row from a matrix: mA(nr,:). */
template< class Ty_>
const valarray<Ty_>& MXr0( valarray<Ty_>& vR, const valarray<Ty_>& mA, const size_t nR, size_t nrA);

template< class Ty_>
valarray<Ty_> MXr0( const valarray<Ty_>& mA, const size_t nR, size_t nrA);


/** 
  Column-wise assignement:  A(:,vi) = B.

  Params:
    - mA: matrix to be modified, [nrA,ncA]. On return the mA( :,vi)
        columns are overwritten with the columns of mB.

    - mB: matrix of values, [nrA,n]. Each column of mB overwrites 
        a column of mA, indicated by the vi index vector.

    - vi: vector of column indices into mA, [1,n]. Valid index values 
        are in the range [0:ncA-1].
*/
template< class Ty_>
const valarray<Ty_>& MSc( valarray<Ty_>& mA, const valarray<size_t>& vi, const valarray<Ty_>& mB);

template< class Ty_>
const valarray<Ty_>& MSc( valarray<Ty_>& mA, const size_t ni, const valarray<Ty_>& vB);


/** 
  Row-wise assignement:  A(vi,:) = B.

  Params:
    - mA: matrix to be modified, [nrA,ncA]. On return the
        mA( vi,:) rows are overwritten with the rows of mB.

    - mB: matrix of values, [n,ncA]. Each row of mB is assigned 
        to a row of mA indicated by the vi index vector.

    - vi: vector of row indices into mA, [1,n]. Valid index values 
        are in the range [0:nrA-1].
*/
template< class Ty_>
const valarray<Ty_>& MSr( valarray<Ty_>& mA, const valarray<size_t>& vi, const valarray<Ty_>& mB);

template< class Ty_>
const valarray<Ty_>& MSr( valarray<Ty_>& mA, const size_t ni, const valarray<Ty_>& vB);


/** 
  Delete a set of rows from a matrix: A(vi,:) = []. 
  
  Params:
    - mA: matrix from which rows are deleted, [nrA,ncA].
       Final size of mA is [nrA-n, ncA].

    - vi: vector of row indices, [1,n]. Valid index values 
        are in the range [0:nrA-1].
  */
template< class Ty_>
const valarray<Ty_>& MDr( valarray<Ty_>& mA, const valarray<size_t>& vi, const size_t nrA);

/** 
  Delete one row from a matrix: A(ni,:) = []. 
  
  Params:
    - mA: matrix from which the ni row is deleted, [nrA,ncA].
       Final size of mA is [nrA-1,ncA].

    - ni: index of the row to be deleted. Valid range is [0:nrA-1].
*/
template< class Ty_>
const valarray<Ty_>& MDr( valarray<Ty_>& mA, const size_t ni, const size_t nrA);


/** 
  Delete a set of columns from a matrix: A(:,vi) = [].

  Params:
    - mA: matrix to be modified, [nrA,ncA]. On return the
        colums indicated by the vi indices are removed from 
        mA, suchat thet the rest is unchanged.

    - vi: column indices into mA, [1,n]. Valid index values 
        are in the range [0:ncA-1].
*/
template< class Ty_>
const valarray<Ty_>& MDc( valarray<Ty_>& mA, const valarray<size_t>& vi, const size_t nrA);

/** 
  Delete one column from a matrix: A(:,ni) = []. 
  
  Params:
    - mA: matrix from which the ni column is deleted, [nrA,ncA].
       Final size of mA is [nrA,ncA-1].

    - ni: index of the column to be deleted. Valid range is [0:ncA-1].
*/
template< class Ty_>
const valarray<Ty_>& MDc( valarray<Ty_>& mA, const size_t ni, const size_t nrA);


/**
  Reshape a matrix by changing the number of rows.

  Params:
    - mA: matrix to be reshaped, [nrA,ncA]. On return mA has the size 
        [nr,ncA] with rows added ( nr>nrA) or removed ( nr<nrA). The 
        [0:min(nr,nrA)] rows keep the original content, while the rest
        is uninitialized.

    - nr: number of rows to resize mA to. If nr > nrA, mA is padded
        with (nr-nrA) additional, uninitialized rows. If nr < nrA 
        the last (nrA-nr) rows of mA are removed.
*/
template< class Ty_>
const valarray<Ty_>& MRr( valarray<Ty_>& mA, const size_t nr, const size_t nrA);


/**
  Reshape a matrix by changing the number of columns.

  Params:
    - mA: matrix to be reshaped, [nrA,ncA]. On return mA has the size 
        [nrA,nc] with columns added ( nc>ncA) or removed ( nc<ncA). The 
        [0:min(nc,ncA)] columns keep the original content, while the rest
        is uninitialized.

    - nc: number of columns to resize mA to. If nc > ncA, mA is padded
        with (nc-ncA) additional, uninitialized columns. If nc < ncA 
        the last (ncA-nc) columns of mA are removed.
*/
template< class Ty_> 
const valarray<Ty_>& MRc( valarray<Ty_>& mA, const size_t nc, const size_t nrA);


/** 
  Rectangular submatrix indexing: M( [c1:c2], [r1:r2])

  Params:
    - mA: two dimensional matrix, [naR,M]
    - naR: matrix rows
    - nOff: offset of the first element in the submatrix, ( 0-based, C-order)
    - nR: number of submatrix rows
    - nC: number of submatrix columns

  MGS returns a submatrix index object.
  MG returns a submatrix of size [nR,nC], with elements from mA.
*/
MATH_API gslice MGS( const size_t& naC, const size_t& nOff, const size_t& nR, const size_t& nC);

template< class Ty_>
valarray<Ty_> MG( const valarray<Ty_>& mA, const size_t& naR, const size_t& nOff, const size_t& nR, const size_t& nC);


/** 
  Create a rectangular sub-matrix indexing object
  for a matrix A of size [nrA,ncA].

  Params:
    - nr: number of rows of the sub-matrix. Valid values are
        in the range [1:nrA].

    - nc: number of colums of the sub-matrix. Valid values are
        in the range [1:ncA].

    - ncA: number of columns of the matrix A.

    - nOff: offset of the first element of the sub-matrix.
        Note: offset is the index of the element in C-matrix
        convention.

  Return value is a gslice object for indexing a [nr,nc] 
  rectangular submatrix in the matrix A, with the first
  element (upper right corner) at offset nOff.
*/
gslice gslice2( const size_t& nr, const size_t& nc, const size_t& ncA =1, const size_t& nOff =0);


/** 
  Stack matrices vertically, along the 1nd dimension: R = [ A; B]. 

  Params:
    - mA: first matrix, [nrA,ncA].
    - mB: second matrix, [nrB,ncB].
    - nc: columns of the stacked matrices. By default nc == ncA == ncB.
    
    - mR: matrix obtained by vertically stacking mA and mB, [ nrA+nrB, nc].
        If mA is empty, mR == mB, respective if mB is empty, mR == mA. If both
        matrices are empty mR is empty. If the number of columns of mA and mB are 
        different ( except being 0), the results are unpredictable.

  Note: empty matrices are ignored.
*/
template< class Ty_>
const valarray<Ty_>& join1( valarray<Ty_>& mR, const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t ncA);

template< class Ty_>
valarray<Ty_> join1( const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nc);


/** 
  Concatenate two matrices vertically, along the 1st dimension: R = [ A; B].
  Allow dynamic type coercion. For more details see join1.
*/
template< class ContR_, class ContA_, class ContB_>
void Join1( ContR_& mR, const ContA_& mA, const ContB_& mB, const size_t ncA);

template< class ContR_, class ContA_, class ContB_>
ContR_ Join1( const ContA_& mA, const ContB_& mB, const size_t ncA);


/** 
  Concatenate matrices horizontally, along the 2nd dimension: R = [ A, B]. 

  Params:
    - mA: first matrix, [nrA,ncA].
    - mB: second matrix, [nrB,ncB].
    - nr: rows of the concatenated matrices. By default nr == nrA == nrB.
    
    - mR: matrix obtained by horizontally concatenating mA and mB, [ nr, ncA+ncB].
        If mA is empty, mR == mB, respective if mB is empty, mR == mA. If both
        matrices are empty mR is empty. If the number of rows of mA and mB are 
        different ( except being 0), the results are unpredictable.

  Note: empty matrices are ignored.
*/
template< class Ty_>
const valarray<Ty_>& join2( valarray<Ty_>& mR, const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nr);

template< class Ty_>
valarray<Ty_> join2( const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nrA);


/** 
  Concatenate matrices horizontally (2nd dimension): R = [ A, B]. 
  Allow dynamic type coercion. For more details see join2.
*/
template< class ContR_, class ContA_, class ContB_>
void Join2( ContR_& mR, const ContA_& mA, const ContB_& mB, const size_t nrA);

template< class ContR_, class ContA_, class ContB_>
ContR_ Join2( const ContA_& mA, const ContB_& mB, const size_t nrA);


//OBSOLETE! TODO: remove all references to the following functions!
/** Return a random number in [0..1).*/
template< class Ty_> Ty_ matRandUnit( void);


//OBSOLETE! TODO: remove all references to the following functions!
/** Generate random vectors.*/
template< class Ty_> const valarray<Ty_>& LoadRand( valarray<Ty_>& v, size_t N);
template< class Ty_> const valarray<Ty_>& LoadRand( valarray<Ty_>& v, size_t N, const Ty_& A, const Ty_& B);


/** 
  Arithmetic progression: [ a, a+r, a+2r, ... a+(n-1)r].

  Params:
    - nA: start value.
    - nSteps: number of steps. If nSteps=0 an empty sequence is generated.
    - nInc: step increment (by default 1).

    - vR: vector of the generated value sequence [ a, a+nD ... a+(nN-1)*nD].
*/
template<class Ty_> 
const valarray<Ty_>& Range( valarray<Ty_>& vR, const Ty_& nA, const size_t& nSteps, const Ty_& nInc =Ty_(1));

template<class Ty_> 
valarray<Ty_> Range( const Ty_& nA, const size_t& nSteps, const Ty_& nInc =Ty_(1));

/** 
  Generate an "indexing" range: [0,1,...].
  See range(nSteps).
*/
template< class Ty_>
size_t Range( valarray<Ty_>& vR, const size_t& nSteps);

template< class Ty_>
valarray<Ty_> Range( const size_t& nSteps);


/** 
  Arithmetic progression for positive integral value sequences.

  Params:
    - nA: integral start value. The Ty_ type should be an integral numeric type. 
    - nSteps: number of steps.
    - nInc: increment value (default 1). Integral positive or negative.

    - vR: resulting vector of size [nSteps,1] with the positive values 
        [nA, nA+nInc, nA+2*nInc, ...]. 

  Note: the return value is expected to be positive for any combination
    of the nA, nSteps and nInc values. If this condition is not satisfied,
    vR contains huge size_t values, instead of the negative values.
*/
template<class Ty_>
size_t range( valarray<size_t>& vR, const Ty_& nA, const size_t& nSteps, const long& nInc =1l);

template<class Ty_>
valarray<size_t> range( const Ty_& nA, const size_t& nSteps, const long& nInc =1l);

/** 
  Generate an "indexing" range: [0,1,...].

  Params:
    - nSteps: nr. of steps

  Return value is a vector of 0-based indices of length nSteps
  and increment 1.
*/
size_t range( valarray<size_t>& vR, const size_t& nSteps);
valarray<size_t> range( const size_t& nSteps);


//TODO: review design
/** 
  Generates a sequence of N values, starting from an initial value and 
  using the given stride.
  
  Params:
    - N: number of values (>=1) to generate.
    - Start: start value. Default value 0.
    - nStride: stride. Default value 1.
    - vL: vector to be loaded with the incremental sequence, [1,N]

  Return value is a reference to vL.
*/
template< class Ty_>
const valarray<Ty_>& LoadInc( valarray<Ty_>& vL, size_t N, const Ty_& Start, const Ty_& nStride);

template< class Ty_>
valarray<Ty_> LoadInc( size_t N = 100, const Ty_& Start  =Ty_(0), const Ty_& nStride =Ty_(1));


//TODO: review design
/** 
  Generate an nL element subdivision of the closed range [nA,nB].

  Params:
    - nL: number of elements (default 100) to be generated in the range [nA..nB].
    - nA, nB: limits of the colsed range.

    - vR: generated closed range vector, with exactly N elemnts. The first and
        last items in are vR[0]=nA, respective vR[N-1]=nB.

  Return value is a reference to vR.
*/
template< class Ty_>
const valarray<Ty_>& LoadRange( valarray<Ty_>& vR, size_t nL, const Ty_& nA, const Ty_& nB);

template< class Ty_>
valarray<Ty_> LoadRange( size_t N =100, const Ty_& nA =Ty_(0), const Ty_& nB =Ty_(1));


template< class Ty_>
size_t countInRange( const valarray<Ty_>& vT, const pair<Ty_,Ty_>& prRange);


//TODO: review design.
/** 
  Generate an interval: [ (t0=a), t1, t2 ... (tn<=b) )

  Params:
    - nA: first element
    - nB: last element
    - nD: incremental resolution, by default 1.

    - vR: generated sequence

  The Ty_ type can be any numerical type.

  Sequence generated for Interval( nA, nB, nD) is [ nA, nA+nD ... nB'), with step
  size nS=(nB-nA)/nD and th alst element nB'<=nB. If nD=0 or (nA+nD) is not in the 
  closed range [nA,nB] an empty vR=[] sequence is generated.
*/
template< class Ty_> 
const valarray<Ty_>& Interval( valarray<Ty_>& vR, const Ty_& nA, const Ty_& nB, const Ty_& nD =Ty_(1));

//OBSOLETE!  TODO: remove all references.
template< class Ty_> 
valarray<Ty_> Interval( const Ty_& begin, const Ty_& end, const Ty_& inc =Ty_(1));


/** 
  Bounded interval for size_t typed numerical values.

  Params:
    - nA: left boundary ( nA>=0). The type Ty1_ can be any integral type.
    - nB: right boundary ( nB>=0). The type Ty2_ can be any integral type.
    - nD: incremental value of type long ( default nD=1). This parameter can be positive or nagative, 
        as long as the generated values are positive or zero. 

  Both parameter types, Ty_ and _Ty2 are expected to be integral types.
  The return value is an unsigned integral range of positive semidefinite values.

  NOTE: if the resulting values vould be negative huge values are generated as
  the values are wrapped over through max<size_t>!
*/
template< class Ty1_, class Ty2_> 
valarray<size_t> interval( const Ty1_& nA, const Ty2_& nB, const long& nD =1l);


/** 
  Generate 2D/3D matrices for plane and volum tesselation.

  Params:
    - vx, vy, vz: range vectors, of size nX, nY, respective nZ.
        For the 2D case only vx and vz are used.

    - mX: rows of mX are copies of vx. For the 2D case mX is of size [nY,nX].
        For the 3D case mX is a 3D matrix of size [nY,nX,nZ], each plane has
        the same values as mX in the 2D case. 

    - mY: columns of mY are copies of vy. For the 2D case mY is of size [nY,nX].
        For the 3D case mY is a 3D matrix of size [nY,nX,nZ], each plane has
        the same values as mY in the 2D case. 

    - mZ: used only for the 3D case, of size [nX,nY,nZ]. Each nZ plane of mZ has 
        a constant value from vz.

  2D case: 
    MeshGrid( mX, mY, vx, vy) with nX=|vx|, nY=|vy| produces 2D matrices mX, mY.

  3D case:
    MeshGrid( mX, mY, mZ, vx, vy, vz) with nX=|vx|, nY=|vy| and nZ=|vz| produces 
    3D matrices mX, mY, mZ as described above.      

  Return value is a vector with the dimension values.
*/
//OBSOLETE: use MeshGrid2()
template< class Ty_> 
valarray<size_t> MeshGrid( valarray<Ty_>& mX, valarray<Ty_>& mY, const valarray<Ty_>& vx, const valarray<Ty_>& vy);

//OBSOLETE: use MeshGrid2()
template< class Ty_> 
valarray<size_t> MeshGrid( valarray<Ty_>& mX, valarray<Ty_>& mY, valarray<Ty_>& mZ, const valarray<Ty_>& vx, const valarray<Ty_>& vy, const valarray<Ty_>& vz);


template< class Ty_> 
void MeshGrid2( valarray<Ty_>& mX, valarray<Ty_>& mY, valarray<size_t>& vDim, const valarray<Ty_>& vx, const valarray<Ty_>& vy);

template< class Ty_> 
void MeshGrid2( valarray<Ty_>& mX, valarray<Ty_>& mY, valarray<Ty_>& mZ, valarray<size_t>& vDim, const valarray<Ty_>& vx, const valarray<Ty_>& vy, const valarray<Ty_>& vz);


/** 
  Kronecker product of two matirces: kron( A, B).

  Params:
    - mA: first matrix, [nrA,ncA].
    - mB: second matrix, [nrB,ncB].
    - nrA: number of rows of mA.
    - nrB: number of rows of mB.

    - mK: the resulting kron(A,B) matrix, [ nrA*nrB, ncA*ncB].
*/
template< class Ty_>
const valarray<Ty_>& Kron( valarray<Ty_>& mK, const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nrA, size_t nrB);

template< class Ty_>
valarray<Ty_> Kron( const valarray<Ty_>& mA, const valarray<Ty_>& mB, size_t nrA, size_t nrB);


/** 
  General matrix inversion: inv(A).

  Params:
    - mA: general quadratic matrix, [na,na]. 
    - na: order of mA.

    - mInv: generated inverse matrix, inv(mA) of size [na,na]. 
        If mA is singular or badly conditionned, mInv will be a [1,1] 
        scalar error value.

  Note: expensive and should never be used for solving a linear system.

  See also: dgetrf, dgetri.
*/
MATH_API const valarray<double>& Inv( valarray<double>& mInv, const valarray<double>& mA, const int& na);
MATH_API valarray<double> Inv( const valarray<double>& mA, const int& na);


/** 
  Invers of a quadratic symmetric matrix: inv(A).

  Params:
    - mA: symmetric matrix, [na,na]
    - na: order of mA

  Return value is zero if the inverse of mA could be computed.
  Return value is k, the order of the Bunch-Kaufman factoriizations 
  diagonal which causes mA to be singular.

  See also: dsytrf, dsytri.
*/
MATH_API int InvSym( valarray<double>& mA, const int& na);


/** 
  Invert a symmetric positive definite matrix.

  Params:
    - mA: symmetric positive-definite quadratic matrix, [nR,nR]
    - nR: number of rows
    - mInv: inverse of mA. If A is singular or badly conditionned 
        mInv is a scalar with the error value.

  See also: dpotrf, dpotri.
*/
MATH_API int InvSymPd( valarray<double>& mA, const int& nR);

MATH_API const valarray<double>& InvSymPd( valarray<double>& mInv, const valarray<double>& mA, const int& nR);

MATH_API valarray<double> InvSymPd( const valarray<double>& mA, const int& nR);


/** 
  Solve a lienar system, with multiple right hand sides and Fortran 
  matrix convention: AX = B.
    
  Params:
    - mAt: transposed version of the mA matrix, [na,na]
    - mBt: transpose version of the mB matrix, [na,nb]
    - mXt: solution matrix mX, [na,nb]. If an error occured mXt[1,1] 
        is the error code.

  See also: dgetrf, dgetrs.
*/
MATH_API const valarray<double>& Solve_AXBt( valarray<double>& mXt, const valarray<double>& mAt, 
                                             const valarray<double>& mBt, const int& na, const int& nb);

MATH_API valarray<double> Solve_AXBt( const valarray<double>& mAt, const valarray<double>& mBt, 
                                      const int& na, const int& nb);


/** 
  Solve a linear system with multiple right hand sides: AX = B.

  Params:
    - mA: matrix mA, [na,na]
    - mB: matrix mB, [na,nb]
    - mX: solution matrix mX, [na.nb]. If an error occured mX[1,1] 
        is the error code.

  See also: dgetrf, dgetrs.
*/
MATH_API const valarray<double>& Solve_AXB( valarray<double>& mX, const valarray<double>& mA, 
                                           const valarray<double>& mB, const int& na, const int& nb);

MATH_API valarray<double> Solve_AXB( const valarray<double>& mA, const valarray<double>& mB, 
                                     const int& na, const int& nb);


/** 
  Return cyclic values from a vector.
*/
MATH_API valarray<double> CycVector( const valarray<double>& vD, const valarray<long>& k);

MATH_API double CycVector( const valarray<double>& vD, const long k);


/**
  Find the indices of the non-zero elements of an integral typed
  numeric vector. 
  
  Return value is true if found non-zero elements.
*/
template< class Ty_>
bool Find( valarray<size_t>& vI, const valarray<Ty_>& vA);

/**
  Return the vector of indices of the non-zero elements of 
  an integral typed numeric vector.
*/
template< class Ty_> 
valarray<size_t> Find( const valarray<Ty_>& vA);

/** 
  Return the vector of indices of the non-zero elements of 
  a real valued numeric vector. Precision of matching is dEps.
*/
valarray<size_t> Find( const valarray<double>& vA, const double dEps);

/** 
  Find the row and colum indices of the non-zero elements for
  an integral typed numeric matrix.

  Params:
    - mA: integral typed numeric matrix, [nrA,ncA].

    - viR: vector of row indices of the zero elements, [1,n].
    - viC: vector of column indices of the zero elements, [1,n].

  Return value is the number of non-zero elements in mA.
*/
template< class Ty_> 
size_t Find( valarray<size_t>& viR, valarray<size_t>& viC, const valarray<Ty_>& mA, const size_t nrA);


//FIX: results might be machine dependent for very small differences in 
//     hitting dT. Introduce an eps mathching precision!
/** 
  Find the index of an item of a real valued vector, nearest 
  to a reference real value.

  Params:
    - vD: list of real values, [N,1]
    - dT: real value to be approximated

  The return value is an index into the vD vector in the range [0..N-1].
*/
ptrdiff_t FindNearest( const valarray<double>& vD, const double dT);


//TODO: see machine precision note from FindNearest.
/** 
  Find the index of first occurance of a real value in a real 
  numerical vector.

  Params:
    - vD: real valuesd vector, [1,n].
    - dT: real reference value
    - eps2: precision of matching.

  Return value is the address of the first matching element up to a precision eps2.
*/
ptrdiff_t FindReal( const valarray<double>& vD, const double dT, const double eps2);


/**
  Mapping of the matching columns from the numerical matrix mA to matrix mB.
  Find also the indices of the unique columns in both matrices.
  
  Params:
    - mA: reference matrix, [nr,ncA]. The mA(:,i), i=[0,ncA) vectors 
        are matched against the vectors mB(:,j), j=[0,ncB).
    - mB: matrix of reference vector set, [nr,ncB].
    - nr: number of rows in mA and mB.
    - eps2: matching precision. If eps2<0, the highest precision( Eps(1)) is used.
        For non-real numerical valued matrices, eps2 is converted to the Ty_ type.

    - miM: mapping of matching column indices from mA(:,*) to mB(:,*). 
        Each column miM(:,k) = [i;j] represents the pair of indices of
        the matching columns mA(:,i) and mB(:,j), within precision eps2.
        Size [2,nM], where nM is the number of matches.
    - viAu: indices of the mA column vectors which are unique to mA, relative to mB.
    - viBu: indices of the mB column vectors which are unique to mB, realtive to mA.

  Return value is the number of matches, nM. 
*/
template< class Ty_>
size_t MapVectors( valarray<size_t>& miM, valarray<size_t>& viAu, valarray<size_t>& viBu,
                   const valarray<Ty_>& mA, const valarray<Ty_>& mB, const size_t nr, const double eps2 =-1);

/**
  Return the mappings from the columns of a numerical matrix mA to another matrix mB.
  Note: mirror mappings are kept.
*/
template< class Ty_>
size_t MapVectors( valarray<size_t>& miM, 
                   const valarray<Ty_>& mA, const valarray<Ty_>& mB, const size_t nr, const double eps2 =-1);

/**
  Return the mappings of the columns in a numerical matrix mA.
  Note: mirror mappings are kept.
*/
template< class Ty_>
size_t MapVectors( valarray<size_t>& miM, 
                   const valarray<Ty_>& mA, const size_t nr, const double eps2 =-1);


/**
  Find the column mappings of a matrix.

  Params:
    - mA: matrix of size [nrA,ncA]
    - nrA: number of rows of mA

    - miM : K column mappings of mA without mirrors, [2,K].
    - viU : indices of the unique columns in mA.
    - viR : indices of the repeated columns.

  The column basis of mA is viB = [viU,viR] and ncA == size(viB) + K.
  Return value is the size of the column basis (viB) of matrix mA.

  NOTE: for returning the basis vectors of mA use MatColBasis() or MatRowBasis().
*/
template<class Ty_>
size_t MapColVectors( valarray<size_t>& miM, valarray<size_t>& viU, valarray<size_t>& viR,
                      const valarray<Ty_>& mA, const size_t nrA);

/**
  Find the row mappings of a matrix.
  See MapColVectors().
*/
template<class Ty_>
size_t MapRowVectors( valarray<size_t>& miM, valarray<size_t>& viU, valarray<size_t>& viR,
                      const valarray<Ty_>& mA, const size_t nrA);


/** Find the column basis of matrix mA.
    Return values is the number of column basis vectors. */
template<class Ty_>
size_t MatColBasis( valarray<Ty_>& mB, const valarray<Ty_>& mA, const size_t nrA);


/** Find the row basis of matrix mA.
    Return values is the number of row basis vectors. */
template<class Ty_>
size_t MatRowBasis( valarray<Ty_>& mB, const valarray<Ty_>& mA, const size_t nrA);


/** Return the range of values in a vector.*/
template<class Ty_>
Ty_ ValRange( const valarray<Ty_>& vR);


/** 
  Set difference for integral valued vectors: C = A\B.

  Params:
    - vA: integral typed number vector, [1,n].
    - vB: integral typed number vector, [1,n].

    - vD: the resulting set difference (A \ B).

  Return value is the size of the D.
*/
template< class Ty_>
size_t SetDiff( valarray<Ty_>& vD, const valarray<Ty_>& vA, const valarray<Ty_>& vB);

template< class Ty_>
valarray<Ty_> SetDiff( const valarray<Ty_>& vA, const valarray<Ty_>& vB);


/** 
  Set-difference for real valued vectors: C = A\B.

  Params:
    - vA: real number vector, set A
    - vB: real number vector, set B
    - eps2: precision
    - vDiff: set A\B, 0-based and sorted indices
*/
bool SetDiff( valarray<size_t>& vDiff, const valarray<double>& vA, const valarray<double>& vB, const double eps2);


/** 
  Set-difference and set-intersection of two real valued vectors: C = A\B.

  Params:
    - vA: real number vector, set A
    - vB: real number vector, set B
    - eps2: precision

    - ixC: set-intersection= (A & B), 0-based, sorted
    - ixU: 
*/
//FIX: overlapping in/out for set_difference!
bool SetDiffInt( valarray<size_t>& ixC, valarray<size_t>& ixU, const valarray<double>& vA, const valarray<double>& vB, const double eps2);


/** Set-difference of index ranges. */
bool SetDiffInc( valarray<size_t>& vixDiff, valarray<size_t>& vIxA, valarray<size_t>& vIxB);


/** 
  Sum the elements of a matrix along the 1st dimension ( ie. column sums).

  Params:
    - mA: matrix, [nrA,ncA]
    - nrA: nr. of rows of mA

  Return value is a row-vector of the column sums, size [1,ncA].
*/
template< class Ty_> 
valarray<Ty_> SumD1( const valarray<Ty_>& mA, const size_t nrA);


/** 
  Sum the elements of a matrix along the 2nd dimension ( ie. row sums).

  Params:
    - mA: matrix, [nrA,ncA]
    - nrA: nr. of rows of mA

  Return value is a column-vector of the row sums, size [nrA,1].
*/
template< class Ty_> 
valarray<Ty_> SumD2( const valarray<Ty_>& mA, const size_t nrA);


/** 
  Sum the elements of a matrix along a dimension.

  Params:
    - mA: matrix, [nrA,ncA]
    - nrA: nr. of rows of mA
    - nDim: dimension to be summed. Use 1 for column sums ( default), 2 for row sums.

  Return a vector of row sums ( size [nrA]) or column sums ( size [ncA]).
*/
template< class Ty_> 
valarray<Ty_> Sum( const valarray<Ty_>& mA, const size_t nrA, const size_t nDim =1);

template< class Ty_> 
size_t Sum( valarray<Ty_>& vR, const valarray<Ty_>& mA, const size_t nrA, const size_t nDim =1);


/** Trace of a matrix. */
template< class Ty_> 
Ty_ Trace( const valarray<Ty_>& mA, const size_t nrA);



pair<size_t,double> DiscreteIndex( const valarray<double>& vR, double r0);
void DiscreteRange( valarray<double>& vD, pair<size_t,size_t>& N, const valarray<double>& vR, double dS);


/** 
  Number of elements of a triangular (2D ) submatrix above the kth diagonal.

  Params:
    - m,n : dimensions of the matrix
    - k: order of the diagonal. Use only k>=0, where k=0 is the main diagonal.

  Return: number of upper( left) triangular elements.
*/
size_t TriLen( size_t m, size_t n, size_t k);


/** 
  Upper-triangular elements of a 2D matrix.

  Params:
    - nR, nC : matrix dimensions, [nR,nC]
    - k: standard diagonal order. k=0 represents the main diagonal, k<0 lower 
          subdiagonals, k>0 upper subdiagonals, inclusive the k-th subdiagonal.

  The first function returns the linear indices of the k-th subdiagonal and 
  all elements above it. The second function returns a copy of these elements.
*/
valarray<size_t> TriU( long nR, long nC, long k =0l);

template< class Ty_> 
valarray<Ty_>& TriU( valarray<Ty_>& mR, const valarray<Ty_>& mD, size_t nR, long k =0l);

template< class Ty_> 
valarray<Ty_> TriU( const valarray<Ty_>& mD, size_t nR, long k =0l);


/** 
  Lower-triangular elements of a 2D matrix.

  Params:
    - nR, nC : matrix dimensions, [nR,nC]
    - k: standard diagonal order. k=0 represents the main diagonal, k<0 lower 
         subdiagonals, k>0 upper subdiagonals.

  The first function returns the linear indices of the k-th subdiagonal and 
  all elements below it. The second function returns a copy of these elements.
*/
valarray<size_t> TriL( long m, long n, long k =0l);

template< class Ty_> 
valarray<Ty_>& TriL( valarray<Ty_>& mR, const valarray<Ty_>& mD, size_t nR, long k =0l);

template< class Ty_> 
valarray<Ty_> TriL( const valarray<Ty_>& mD, size_t nR, long k =0l);


/** 
  Mirror ( inplace) the upper triangular part of a quadratic matrix.

  Params:
    - mAU: quadratic matrix, [nR,nR]
    - nR: order of mAU

  On return mAU contains the original matrix, but with the lower
  triangle part overwritten by the mirror image of the upper triangular part.
*/
template< class Ty_>
void Mirroru( valarray<Ty_>& mAU, const size_t nR);


/** Mirror the upper triangular part of a quadratic matrix. */
template< class Ty_>
size_t MirrorU( valarray<Ty_>& mAU, const valarray<Ty_>& mA, const size_t nR);

template< class Ty_>
valarray<Ty_> MirrorU( const valarray<Ty_>& mA, const size_t nR);


/** 
  Mirror (inplace) the lower triangular part of a quadratic matrix.

  Params:
    - mAL: quadratic matrix, [nR,nR]
    - nR: order of mAL

  On return mAL contains the original matrix, but with the upper
  triangle part overwritten by the mirror image of the lower triangular part.
*/
template< class Ty_>
void Mirrorl( valarray<Ty_>& mAL, const size_t nR);

/** Mirror the lower triangular part of a quadratic matrix.*/
template< class Ty_>
size_t MirrorL( valarray<Ty_>& mAL, const valarray<Ty_>& mA, const size_t nR);

template< class Ty_>
valarray<Ty_> MirrorL( const valarray<Ty_>& mA, const size_t nR);


/** 
  Cumulative sums.

  Params:
    - mT: matrix of [nR,nC]
    - nR: number of rows of mT.
    - nDim: cumulative sum is computed along this diemension, nDim= {1,2}. If nR=1 and (nC>1), 
            mT is a row vector and a cumulative sum is computed always along nDim=2.
    - mS: matrix of cumulative sums
*/
template< class Ty_> 
const valarray<Ty_>& CumSum( valarray<Ty_>& mS, const valarray<Ty_>& mT, size_t nR =1, size_t nDim =1);

template< class Ty_> 
valarray<Ty_> CumSum( const valarray<Ty_>& mT, size_t nR =1, size_t nDim =1);


/** 
  Remove a set of elements from a vector, specified by a list of unique, sorted indices.
    
  Params:
    - vD: vector of size N, to be trimmed in place. Items remaining in vD are stored
            at the beginning N-M positions. The content of the rest of vD is unpredictable.
    - vI: vector of M, zero based indices. The elements of vI have to be sorted in
          ascending order, without duplicates. Indices greater than length(vD) are ignored.

  Return value is the number of elements removed from vD. If vI is empty vD is not changed.
*/
template< class Ty_> 
size_t RemoveSorted( valarray<Ty_>& vD, const valarray<size_t>& vI);


/** 
  Remove a set of elements from a vector, specified by a list of indices.

  Params:
    - vI: list of indices. Duplicates or out of range values are ignored.
    - vD: vector operated on.

  The return value of the first function is the number of elements removed from vD
  and vD is modified in-place. If vI=[] no changes are done and the return value is zero.

  The return value of the second function is the vector vD without the elements
  with indices in vI. If vI=[] the return value is a full copy of vD.
*/
template< class Ty_> 
size_t Remove( valarray<Ty_>& vD, const valarray<size_t>& vI);

template< class Ty_> 
valarray<Ty_> Remove( const valarray<Ty_>& vD, const valarray<size_t>& vI);


/** 
  Find the unique values in a vector for integral and string typed items.

  Note: avoid using this for real values, with full precision as the 
        behavior is platform dependent!

  Params:
    - vX: vector of items, [N,1]
    - vU: vector with unique items from vX

  Return value of the first function is the size of vU.
  Return value of the second function is the unique values vector (vU).
*/
//TODO: make these obsolete in favour of Uniquev()!
template< class Ty_> 
size_t Unique( valarray<Ty_>& vU, const valarray<Ty_>& vX);

template< class Ty_> 
valarray<Ty_> Unique( const valarray<Ty_>& vX);


/**
  Find the subset of unique coordinates in a list of 3D points.

  Params:
    - mP: set of column vectors, [M,N]
    - eps2: matching precision
    
  Return value is the vector of column indices of the 
  unique subset in mP.
*/
//TODO: make this obsolete!!!
template< class Ty_>
const valarray<size_t>& Unique3D( valarray<size_t>& viU, const valarray<Ty_>& mP, 
                                const size_t nR, const double eps2);

//TODO: make this obsolete!!!
template< class Ty_>
valarray<size_t> Unique3D( const valarray<Ty_>& mP, const size_t nR, const double eps2);


/** 
  Test if the elements of a vector are NaN.

  Params:
    - vX: vector of double or float numbers.

    - vB: vector of booleans, of the same size as vX. For each NaN element
        of vX, the corresponding element of vB is set to true, otherways false.

  Return value of IsNan( vB, vX) is the number of NaN's in vX.
  Return value of IsNan( vX) is the vector of booleans, with true for each
  NaN element of vX.
*/
template< class Ty_>
size_t IsNan( valarray<bool>& vB, const valarray<Ty_>& vX);

template< class Ty_>
valarray<bool> IsNan( const valarray<Ty_>& vX);


/** 
  Test if the elements of a vector are finite.
  A number x is considered finite if it is not infinite, -Inf < x < Inf,
  respective is not a NaN.

  Params:
    - vX: vector of double or float numbers.

    - vB: vector of booleans, of the same size as vX. For each finite element
        of vX, the corresponding element of vB is set to true, otherways false.

  Return value of IsFinite( vB, vX) is the number of finite elements of vX.
  Return value of IsFinite( vX) is the vector of booleans, with true for each
  finite element of vX.
*/
template< class Ty_>
size_t IsFinite( valarray<bool>& vB, const valarray<Ty_>& vX);

template< class Ty_>
valarray<bool> IsFinite( const valarray<Ty_>& vX);


/** 
  Test if a real number vector is identic to a zero vector.

  Params:
    - vA: real number vector
    - dEpsScale: scale factor for the eps(0) uncertainty limit.
*/
template< class Ty_> 
bool IsZero( const valarray<Ty_>& vA, const double dEpsScale =1.0);


/** 
  Test if two real number vectors/matrices are equal, within a numerical uncertainty.

  Params:
    - vA: first vector
    - vB: second vector
    - dEpsScale: common scale factor for the eps( max(|A|,|B|)) uncertainty limit.
*/
template< class Ty_> 
bool IsEqual( const valarray<Ty_>& vA, const valarray<Ty_>& vB, const double dEpsScale =1.0);


/** 
  Test if a real number matrix is symmetric.

  Params:
    - mX: real number matrix, assumed of size [n,n]
    - n: order of the mX matrix. If n=0, order is detected runtime.
    - dEpsScale: common scale factor for the eps( max(|A|,|B|)) uncertainty limit.
*/
template< class Ty_> 
bool IsSymmetric( const valarray<Ty_>& mX, const size_t n =0, const double dEpsScale =1.0);


/** Count the number of occurances of a value in a vector.*/
template< class Ty_> 
size_t Count( const valarray<Ty_>& vX, const Ty_& x);

/** Count the number of occurances of a real value in a real vector/matrix.*/
template< class Ty_> 
size_t Countf( const valarray<Ty_>& vX, const Ty_& x, const double& dEpsScale =1.0);


/** 
  Rank of a matrix.

  Params:
    - mX: matrix of size [nR,nC]
    - nR: nr. of rows in mX
    - dEps: tolerance of singular values

  Return value is the rank of mX.
*/
MATH_API int Rank( const valarray<double>& mX, const size_t& nR, const double& dEps =0.0);


/** 
  Compute the condition number of a matrix using different norms.

  Params:
    - mA: matrix, [nrA,N]
    - nrA: number of rows of mA
    - nNorm: type of norm, values in the range {1,2,3,4}:
          1 = 1-norm, largest column sum
          2 = largest singular value. This is the default matrix/vector norm.
          3 = infinity-norm, largest row sum
          4 = Euclidean or Frobenius norm
          See also MNorm() for more description.

  The 2-norm works with rectangular matrices of any shape and
  the condition number is computed as the ratio of largest 
  singular value to the smallest one.

  The 1-norm, infinity-norm and Frobenius-norm are only for
  quadratic matrices, as the condition number is based on the 
  inverse of mA. If mA is singular, the condition number is 'inf'.

  On success the return value is the condition number of mA, otherways
  infinity if:
    - 2-norm is used and mA has a near-zero singular value
    - 1-, infinity- or Frobenius-norm is used and mA is singular
*/
MATH_API double Cond( const valarray<double>& mA, const size_t& nrA, const size_t& nNorm =2);


/** 
  Estimate the reciprocal of the condition number of a matrix in infinte- or 1-norm.
    K_1 = norm1(X)*norm1(inv(X))
    K_inf = norminf(X)*norminf(inv(X))

  Params:
    - mA: general square matrix, [N,N]
    - nrA: order of mA, ( nrA==N)
    - nNorm: 1=1-norm or 3=infinity-norm. Default is 1-norm.

  Obs: 1-norm of a matrix is the maximum column sum of absolute values,
  while infinity-norm is the similar maximum, but along the rows.
*/
MATH_API double Rcond( const valarray<double>& mA, const size_t& nrA, const size_t& nNorm =1);


/** 
  Singular Value Decomposition of a general matrix: [U,S,V] = Svd(X) 
  such that: X = USV'.

  Params:
    - X: matrix to be decomposed, [M,N]
    - nR: number of rows of X, (notation: M=nR)
    - bEco: if true the economic size SVD is computed, otherwayz the full one

    - U: left singular vectors, for full SVD size of U is [M,M], otherways [M,min(M,N)]
    - S: singular values as a diagonal matrix, for full SVD size of S is [M,N], otherways [ min(M,N), min(M,N)]
    - V: right singular vectors, for full SVD size of V is [N,N], otherways [N,min(M,N)]                               

  On error the return value of Svd() is nonzero and U, S, V are empty ( zero size).
  On success the return value is zero and the U, S, V have the following sizes:
    1) full SVD :  U[M,M], S[M,N], V[N,N]
    2) econ. SVD:  U[M,p], S[p,p], V[N,p], with p=min(M,N)

  Obs1: the singular values in S are unique and most min(M,N) are nonzero. But the U, V singular 
  vectors are not unique. Matlab might use beside dgesvd() also other means of SVD decomposition, 
  as for certain matrices A ( with cond(A) small), the U and V vectors might differ considerably.

  Obs2: Negative return values indicate wrong parameters to dgesvd(), positive values indicate 
  the number of superdiagonals which could not be reduced to zero in the intermediate dbdsqr()
  bidiagonalization step.
*/
//!TODO: probably there is a fine tuning in Matlab, which treats M<N and M>N cases different,
//          something similar to row balancing in eigenvalue problems.
MATH_API int Svd( valarray<double>& U, valarray<double>& S, valarray<double>& V, const valarray<double>& X, const size_t& nR, const bool& bEco =false);


/** 
  Compute the singular values of a matrix: S = svd(X).

  Params:
    - X: data matrix, [M,N]
    - nR: number of rows of X, (notation: nR==M)
    - S: singular values vector, [ min(M,N)], where X = USV'

  Return value is zero on success. On failure S is reset to zero size.
*/
MATH_API int Svd( valarray<double>& S, const valarray<double>& X, const size_t& nR);


/** 
  Compute the K-th singular value and left/right vectors of a matrix: [u,s,v] = SvdK( X, k).

  Params:
    - mX: data matrix, [M,N]
    - nK: number of the singular value set to be returned, range [1 .. min(M,N)]
    - vU: nK-th left singular vector, [M]
    - dS: nK-th singular value
    - vV: nK-th  right singular vector, [N]

  Return value is zero on success.
  For more details see Svd(U,S,V,X,...).
*/
MATH_API int SvdK( valarray<double>& vU, double& dS, valarray<double>& vV, const valarray<double>& mX, const size_t& nR, const size_t& nK =1);


/** 
  Compute eigenvalues and eigenvectors for a real symmetric matrix.
  For an eigendecomposition guilt: XV = VD.
 
  Params:
    - mX: real symmetric matrix, [n,n]
    - bVec: if true both eigen vectors and values are computed, otherways only the eigenvalues

    - mV: eigenvectors,[n,n], computed only when bVec=true. Otherways mV is emptyed.
    - mD: eigenvalues, always computed. If bVec=true, mD is diagonal of size [n,n] with the
          eigenvalues on the diagonal. If bVec=false, mD is a vector [n].
*/  
MATH_API int EigSym( valarray<double>& mV, valarray<double>& mD, const valarray<double>& mX, const bool& bVec =true);


/** Driver for real asymmetric eigenvalue problems. */
MATH_API int EigAsym( valarray<double>& mVL, valarray<double>& mVR, valarray<double>& vreS, valarray<double>& vimS,
                      const valarray<double>& mX, const bitset<8>& bsOpt);

/** 
  Compute eigenvalues and eigenvectors for real asymmetric matrices.

  Params:
    - mX: quadratic data matrix, [n,n]
    - bVect: if true, compute eigenvectors and values, otherways only eigenvalues
    - bRight: if true, compute the right eigenvectors, otherways the right ones

    - mV: left or right eigenvector matrix
    - mD: eigenvalues. If the eigenvalues are real size(mD)=[n,n], else size(mD)=[2*n,n]
         The real part of the complex eigenvalues are in [1:n,1:n], the imaginary part
         are stored in the lower part, [n:2n,n]
*/
MATH_API int EigAsym( valarray<double>& mV, valarray<double>& mD, const valarray<double>& mX, 
         const bool& bVect =true, const bool& bRight =true, const bool& bNoBalance =false);


/** 
  Compute eigenvalues and eigenvectors for real matrices.

  Params:
    - mX: real quadratic data matrix, [N,N]. If mX is symmetric bRight and bBalance are ignored.
    - bVect: compute eigenvectors and eigenvalues if bVect=true, otherways compute only eigenvalues.
    - bRight: if mX is asymmetric and bRight=true, compute the right eigenvectors of mX, 
        otherways if bRight=false compute the left eigenvectors. This options is ignored
        if mX is symmetric.
    - bNoBalance: skip rowbalancing in asymmetric matrices with small entries. By default
        row balancing is used in asymmetric eigenvalue problems to improve conditioning.
        A sideffect of balancing is the amplification of noise, which might overrun small
        matrix element. Set bNoBalance=true to skip default balancing and reduce noise amplification.

    - mV: right ( bRight=true) or left ( bRight=false) eigenvector matrix. Eigenvectors are
        computed only if bVect is true. If the eigenvectors are real the size of mV is [N,N]
        otherways [2N,N]. For complex eigenvectors the lower half of mV, mV[N+1:2N,:] contains
        the imaginary parts of the complex eigenvector coefficients.
        
    - mD: diagonal matrix of eigenvalues. If the eigenvalues are real mD is of size [N,N],
         with the eigenvalues on the main diagonal. If the eigenvalues are complex mD is of
         size [2*N,N] and the imaginary parts are stored in the diagonal elements of the lower 
         half submatrix of mD, mD[N+1:2N,:].
*/
MATH_API int Eig( valarray<double>& mV, valarray<double>& mD, const valarray<double>& mX, 
         const bool& bVect =true, const bool& bRight =true, const bool& bNoBalance =false);


/** 
  Principal components analysis.

  Params:
    - mX: data matrix, [n,p]. Rows are the samples, columns the variables.
    - nR: number of samples, ( =n)
    - bEcon: if true return only nonzero eigenvalues and the 
             corresponding  scores and loadings

    - mCoeff: loadings, [p,p] matrix of PCA coefficients in decreasing order
          of component variance. Each colum is the coefficients for one
          principal component.

    - mScore: scores matrix of size [n,p], the projections of the samples into
          the PCA space. Rows of mScore are the observations, columns are the 
          components. If n<=p, mScore(:,n:p)=0 and vLatent(n:p)=0.

    - vLatent: eigenvalues of the covariance matrix of mX or variances of 
          the columns of mScore.

    - vTSq: Hotellins's T^2 statistic for each data point. Represents the measure of
          multivariate distance of each observation from the center of the dataset.
*/
MATH_API pair<long,string> Princomp( valarray<double>& mCoeff, valarray<double>& mScore, valarray<double>& vLatent, valarray<double>& vTSq, 
                                     const valarray<double>& mX, const size_t& nR, const bool& bEcon =false);


/** 
  Solve a system of linear equations with upper/lower triangular matrix 
  and multiple right-hand sides vectors: A*X = B.

  Params:
    - A: quadratic triangular matrix, [n,n]
    - B: right hand side vectors, [n,k]. If k=1 use SolveAxB_Tri instead
        of this function.
    - nrA: order of A, ( nrA=n)
    - bUpper: A is upper triangular, otherways lowertriangular
    - bNoTrp: solve AX=B, otherways A'X=B
    - bUnitTri: A is unit triangular (ie.: has unit main diagonal)

    - X: solution matrix of AX=B, [n,r]

  Return value is zero on success.
*/
MATH_API int SolveAxB_TriMrhs( valarray<double>& X, const valarray<double>& A, const valarray<double>& B, const size_t& nrA, 
                               const bool& bUpper =true, const bool& bNoTrp =true, const bool& bUnitTri =true);

/** 
  Solve a system of linear eq. with triangular matrix and a right-hand vetor: A*x = b.

  Params:
    - A: quadratic triangular matrix, [n,n]
    - B: right hand side vector, [n]
    - bUpper: A is upper triangular, otherways lowertriangular
    - bNoTrp: solve AX=B, otherways A'X=B
    - bUnitTri: A is unit triangular ( main diagonal ia all 1)

  Return value is zero on success.
*/
MATH_API void SolveAxB_Tri( valarray<double>& x, const valarray<double>& A, const valarray<double>& b,
                            const bool& bUpper =true, const bool& bNoTrp =true, const bool& bUnitTri =true);

/** 
  QR orthogonal decomposition: X*E = Q*R.

  Params:
    - mX: data atrix, [nR,nC]
    - mQ: unitary orthogonal matrix. If nR>nC, size(mQ)=[nR,nR]
    - mR: upper diagonal matrix, [nR,nC] 
    - vE: column permutation vector, [nC].

    If nR >= nC, [nR,nC] = size(Q), [nC,nC] = size(R).
    If  nR < nC, [nR,nR] = size(Q), [nR,nC] = size(R).
    Permuted decomposition:  mX(:,vE) = Q*R.

  Return value is zero on success.
*/
MATH_API int QR( valarray<double>& mQ, valarray<double>& mR, valarray<size_t>& mE,
                 const valarray<double>& mX, const size_t& nR);

/** 
  Linear mapping of values to a range.

  Params:
    - vS: values in the input range
    - vSRange: source range [ srcMin, srcMax]
    - vTRange: target range [ trgMin, trgMax]
    - bClamp: clamp values in vS to the vSRange range

    - vT: the vS values mapped to the target range
*/
template< class Ty_>
size_t MapToRange( valarray<double>& vT, const valarray<Ty_>& vS, const valarray<double>& vSRange, 
                   const valarray<double>& vTRange, const bool bClamp =false);

/** Map the input range to the sigmoid range [-1,+1].*/
template< class Ty_>
size_t MapToSigRange( valarray<double>& vT, const valarray<Ty_>& vS, const valarray<double>& vSRange, 
                      const bool bClamp);

template< class Ty_>
valarray<Ty_> MapToSigRange( const valarray<Ty_>& vS, const valarray<double>& vSRange, const bool bClamp);


/** 
  Linear mapping of values to a range. Source range is automatically detected. 

  Params:
    - vS: source values
    - vRanges: mapping range, [ trgMin, trgMax]

    - vT: the vS values mapped to the target range

  The source value range is automatically detected from vS.
*/
template< class Ty_>
size_t MapToRangeAuto( valarray<double>& vT, const valarray<Ty_>& vS, const valarray<double>& vTRange);


/** Compute a sigmoid mapping of a linear value range.*/
void Sigmoid( valarray<double>& vR, const valarray<double>& vX, const valarray<double>& vRange, const double dA);

valarray<double> Sigmoid( const valarray<double>& vX, const valarray<double>& vRange, const double dA);


/** 
  Return a vector defined by a direction and length.
*/
template< class Ty_>
valarray<Ty_> Vector( const valarray<Ty_>& vDir, const Ty_& dLen);

} //namespace mft

// Inline implementations
#include "MFT.inl"

#endif // __MFT__H__
