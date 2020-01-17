//=================================================================================================
/*!
//  Copyright (C] 2012-2019 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised] BSD License. Redistribution and use in source and binary
//  forms][ with or without modification][ are permitted provided that the following conditions
//  are mstd::complex<double>:
//
//  1. Redistributions of source code must rstd::complex<double>ain the above copyright notice][ this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice][ this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES][ INCLUDING][ BUT NOT LIMITED TO][ THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT][ INDIRECT][
//  INCIDENTAL][ SPECIAL][ EXEMPLARY][ OR CONSEQUENTIAL DAMAGES (INCLUDING][ BUT NOT LIMITED
//  TO][ PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE][ DATA][ OR PROFITS; OR
//  BUSINESS INTERRUPTION] HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY][ WHstd::complex<double>HER IN
//  CONTRACT][ STRICT LIABILITY][ OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE] ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE][ EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _H_INVERSION_CPLX_
#define _H_INVERSION_CPLX_
/*
lapacke uses C99 _Complex by default, which Visual C++ doesn't support.
You can define those macros as suggested to use std::complex instead, which are supported by Visual C++:
https://stackoverflow.com/questions/24853450/errors-using-lapack-c-header-in-c-with-visual-studio-2010
*/
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include "lapacke.h"

inline std::complex<double> DetCplx(std::complex<double>** mat, int dim) {
  switch (dim) {
  case 1: return DetCplx1x1(mat);
    break;
  case 2:return DetCplx2x2(mat);
    break;
  case 3:return DetCplx3x3(mat);
    break;
  case 4:return DetCplx4x4(mat);
    break;
  case 5:return DetCplx5x5(mat);
    break;
  case 6:return DetCplx6x6(mat);
    break;
  default:return DetCplxNxN(mat,dim);
    break;
  }
}

inline std::complex<double> InvertCplx(std::complex<double>** A, std::complex<double>** B, int dim) {
  switch (dim) {
  case 1: InvertCplx1x1(A,B);
    break;
  case 2: InvertCplx2x2(A, B);
    break;
  case 3:InvertCplx3x3(A, B);
    break;
  case 4:InvertCplx4x4(A, B);
    break;
  case 5:InvertCplx5x5(A, B);
    break;
  case 6:InvertCplx6x6(A, B);
    break;
  default:InvertCplxNxN(A, B, dim);
    break;
  }
}

inline std::complex<double> DetCplxNxN(std::complex<double>** mat, int dim) {
  int i;
  std::complex<double> det;
  double t;
  int n;
  int* idx;
  std::complex<double>* tmat;
  n = dim;

  idx = new int[n];
  tmat = new std::complex<double>[dim * dim];
  //TODO
  //memcpy : mat -> tmat;
  for (i = 0; i < n; i++) {
    memcpy(tmat + i * n * sizeof(std::complex<double>), mat[i], sizeof(std::complex<double>) * n);
  }

  //TODO
  //row ? col ?
  LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, tmat, n, idx);

  det = { 1.0,0.0 };
  for (i = 0; i < n; i++) {
    if (i + 1 != idx[i])
      det *= -tmat[i * n + i];
    else
      det *= tmat[i * n + i];
  }
  delete[] idx;
  delete[] tmat;
  return det;
}

inline void InvertCplxNxN(std::complex<double>**A, std::complex<double>**B, int n) {
  int* idx;
  std::complex<double>* tmat;
  tmat = new std::complex<double>[n * n];

  idx = new int[n];
  for (int i = 0; i < n; i++) {
    memcpy(tmat + i * n * sizeof(std::complex<double>), A[i], sizeof(std::complex<double>) * n);
  }

  LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n, tmat, n, idx);
  LAPACKE_zgetri(LAPACK_COL_MAJOR, n, tmat, n, idx);
  
  for (int i = 0; i < n; i++) {
    memcpy(B[i],tmat + i * n * sizeof(std::complex<double>), sizeof(std::complex<double>) * n);
  }

  delete[] idx;
}

inline void InvertCplx1x1(std::complex<double>** A, std::complex<double>** B) {
  B[0][0] = 1.0 / A[0][0];
}

inline std::complex<double> DetCplx1x1(std::complex<double>** A) {
  std::complex<double> det = (1.0 / A[0][0]);
  return det;
}

inline void InvertCplx2x2(std::complex<double>** A, std::complex<double>** B)
{

  const std::complex<double> det(A[0][0] * A[1][1] - A[0][1] * A[1][0]);

  const std::complex<double> idet(1.0 / det);
  const std::complex<double> a11(A[0][0] * idet);

  A[0][0] = A[1][1] * idet;
  A[1][0] = -A[1][0] * idet;
  A[0][1] = -A[0][1] * idet;
  A[1][1] = a11;
}

inline std::complex<double> DetCplx2x2(std::complex<double>** A)
{
  const std::complex<double> det(A[0][0] * A[1][1] - A[0][1] * A[1][0]);
  return det;
}


inline void InvertCplx3x3(std::complex<double>** A, std::complex<double>** B)
{

  B[0][0] = A[1][1] * A[2][2] - A[1][2] * A[2][1];
  B[1][0] = A[1][2] * A[2][0] - A[1][0] * A[2][2];
  B[2][0] = A[1][0] * A[2][1] - A[1][1] * A[2][0];

  const std::complex<double> det(A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0]);

  B[0][1] = A[0][2] * A[2][1] - A[0][1] * A[2][2];
  B[1][1] = A[0][0] * A[2][2] - A[0][2] * A[2][0];
  B[2][1] = A[0][1] * A[2][0] - A[0][0] * A[2][1];
  B[0][2] = A[0][1] * A[1][2] - A[0][2] * A[1][1];
  B[1][2] = A[0][2] * A[1][0] - A[0][0] * A[1][2];
  B[2][2] = A[0][0] * A[1][1] - A[0][1] * A[1][0];

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      B[i][j] /= det;
}

inline std::complex<double> DetCplx3x3(std::complex<double>** A)
{

  std::complex<double>b0 = A[1][1] * A[2][2] - A[1][2] * A[2][1];
  std::complex<double>b1 = A[1][2] * A[2][0] - A[1][0] * A[2][2];
  std::complex<double>b2 = A[1][0] * A[2][1] - A[1][1] * A[2][0];
  std::complex<double> det(A[0][0] * b0 + A[0][1] * b1 + A[0][2] * b2);
  return det;
}


inline void InvertCplx4x4(std::complex<double>** A, std::complex<double>** B)
{
  std::complex<double> tmp1(A[2][2] * A[3][3] - A[2][3] * A[3][2]);
  std::complex<double> tmp2(A[2][1] * A[3][3] - A[2][3] * A[3][1]);
  std::complex<double> tmp3(A[2][1] * A[3][2] - A[2][2] * A[3][1]);

  B[0][0] = A[1][1] * tmp1 - A[1][2] * tmp2 + A[1][3] * tmp3;
  B[0][1] = A[0][2] * tmp2 - A[0][1] * tmp1 - A[0][3] * tmp3;

  std::complex<double> tmp4(A[2][0] * A[3][3] - A[2][3] * A[3][0]);
  std::complex<double> tmp5(A[2][0] * A[3][2] - A[2][2] * A[3][0]);

  B[1][0] = A[1][2] * tmp4 - A[1][0] * tmp1 - A[1][3] * tmp5;
  B[1][1] = A[0][0] * tmp1 - A[0][2] * tmp4 + A[0][3] * tmp5;

  tmp1 = A[2][0] * A[3][1] - A[2][1] * A[3][0];

  B[2][0] = A[1][0] * tmp2 - A[1][1] * tmp4 + A[1][3] * tmp1;
  B[2][1] = A[0][1] * tmp4 - A[0][0] * tmp2 - A[0][3] * tmp1;
  B[3][0] = A[1][1] * tmp5 - A[1][0] * tmp3 - A[1][2] * tmp1;
  B[3][1] = A[0][0] * tmp3 - A[0][1] * tmp5 + A[0][2] * tmp1;

  tmp1 = A[0][2] * A[1][3] - A[0][3] * A[1][2];
  tmp2 = A[0][1] * A[1][3] - A[0][3] * A[1][1];
  tmp3 = A[0][1] * A[1][2] - A[0][2] * A[1][1];

  B[0][2] = A[3][1] * tmp1 - A[3][2] * tmp2 + A[3][3] * tmp3;
  B[0][3] = A[2][2] * tmp2 - A[2][1] * tmp1 - A[2][3] * tmp3;

  tmp4 = A[0][0] * A[1][3] - A[0][3] * A[1][0];
  tmp5 = A[0][0] * A[1][2] - A[0][2] * A[1][0];

  B[1][2] = A[3][2] * tmp4 - A[3][0] * tmp1 - A[3][3] * tmp5;
  B[1][3] = A[2][0] * tmp1 - A[2][2] * tmp4 + A[2][3] * tmp5;

  tmp1 = A[0][0] * A[1][1] - A[0][1] * A[1][0];

  B[2][2] = A[3][0] * tmp2 - A[3][1] * tmp4 + A[3][3] * tmp1;
  B[2][3] = A[2][1] * tmp4 - A[2][0] * tmp2 - A[2][3] * tmp1;
  B[3][2] = A[3][1] * tmp5 - A[3][0] * tmp3 - A[3][2] * tmp1;
  B[3][3] = A[2][0] * tmp3 - A[2][1] * tmp5 + A[2][2] * tmp1;

  const std::complex<double> det(A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0] + A[0][3] * B[3][0]);

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      B[i][j] /= det;
}

inline std::complex<double> DetCplx4x4(std::complex<double>** A)
{
  std::complex<double> tmp1(A[2][2] * A[3][3] - A[2][3] * A[3][2]);
  std::complex<double> tmp2(A[2][1] * A[3][3] - A[2][3] * A[3][1]);
  std::complex<double> tmp3(A[2][1] * A[3][2] - A[2][2] * A[3][1]);

  std::complex<double>b0 = A[1][1] * tmp1 - A[1][2] * tmp2 + A[1][3] * tmp3;

  std::complex<double> tmp4(A[2][0] * A[3][3] - A[2][3] * A[3][0]);
  std::complex<double> tmp5(A[2][0] * A[3][2] - A[2][2] * A[3][0]);

  std::complex<double>b1 = A[1][2] * tmp4 - A[1][0] * tmp1 - A[1][3] * tmp5;

  tmp1 = A[2][0] * A[3][1] - A[2][1] * A[3][0];

  std::complex<double>b2 = A[1][0] * tmp2 - A[1][1] * tmp4 + A[1][3] * tmp1;
  std::complex<double>b3 = A[1][1] * tmp5 - A[1][0] * tmp3 - A[1][2] * tmp1;


  std::complex<double> det(A[0][0] * b0 + A[0][1] * b1 + A[0][2] * b2 + A[0][3] * b3);
  return  det;
}

inline void InvertCplx5x5(std::complex<double>** A, std::complex<double>** B)
{

  std::complex<double> tmp1(A[3][3] * A[4][4] - A[3][4] * A[4][3]);
  std::complex<double> tmp2(A[3][2] * A[4][4] - A[3][4] * A[4][2]);
  std::complex<double> tmp3(A[3][2] * A[4][3] - A[3][3] * A[4][2]);
  std::complex<double> tmp4(A[3][1] * A[4][4] - A[3][4] * A[4][1]);
  std::complex<double> tmp5(A[3][1] * A[4][3] - A[3][3] * A[4][1]);
  std::complex<double> tmp6(A[3][1] * A[4][2] - A[3][2] * A[4][1]);
  std::complex<double> tmp7(A[3][0] * A[4][4] - A[3][4] * A[4][0]);
  std::complex<double> tmp8(A[3][0] * A[4][3] - A[3][3] * A[4][0]);
  std::complex<double> tmp9(A[3][0] * A[4][2] - A[3][2] * A[4][0]);
  std::complex<double> tmp10(A[3][0] * A[4][1] - A[3][1] * A[4][0]);

  std::complex<double> tmp11(A[2][2] * tmp1 - A[2][3] * tmp2 + A[2][4] * tmp3);
  std::complex<double> tmp12(A[2][1] * tmp1 - A[2][3] * tmp4 + A[2][4] * tmp5);
  std::complex<double> tmp13(A[2][1] * tmp2 - A[2][2] * tmp4 + A[2][4] * tmp6);
  std::complex<double> tmp14(A[2][1] * tmp3 - A[2][2] * tmp5 + A[2][3] * tmp6);
  std::complex<double> tmp15(A[2][0] * tmp1 - A[2][3] * tmp7 + A[2][4] * tmp8);
  std::complex<double> tmp16(A[2][0] * tmp2 - A[2][2] * tmp7 + A[2][4] * tmp9);
  std::complex<double> tmp17(A[2][0] * tmp3 - A[2][2] * tmp8 + A[2][3] * tmp9);

  B[0][0] = A[1][1] * tmp11 - A[1][2] * tmp12 + A[1][3] * tmp13 - A[1][4] * tmp14;
  B[0][1] = -A[0][1] * tmp11 + A[0][2] * tmp12 - A[0][3] * tmp13 + A[0][4] * tmp14;
  B[1][0] = -A[1][0] * tmp11 + A[1][2] * tmp15 - A[1][3] * tmp16 + A[1][4] * tmp17;
  B[1][1] = A[0][0] * tmp11 - A[0][2] * tmp15 + A[0][3] * tmp16 - A[0][4] * tmp17;

  std::complex<double> tmp18(A[2][0] * tmp4 - A[2][1] * tmp7 + A[2][4] * tmp10);
  std::complex<double> tmp19(A[2][0] * tmp5 - A[2][1] * tmp8 + A[2][3] * tmp10);
  std::complex<double> tmp20(A[2][0] * tmp6 - A[2][1] * tmp9 + A[2][2] * tmp10);

  B[2][0] = A[1][0] * tmp12 - A[1][1] * tmp15 + A[1][3] * tmp18 - A[1][4] * tmp19;
  B[2][1] = -A[0][0] * tmp12 + A[0][1] * tmp15 - A[0][3] * tmp18 + A[0][4] * tmp19;
  B[3][0] = -A[1][0] * tmp13 + A[1][1] * tmp16 - A[1][2] * tmp18 + A[1][4] * tmp20;
  B[3][1] = A[0][0] * tmp13 - A[0][1] * tmp16 + A[0][2] * tmp18 - A[0][4] * tmp20;
  B[4][0] = A[1][0] * tmp14 - A[1][1] * tmp17 + A[1][2] * tmp19 - A[1][3] * tmp20;
  B[4][1] = -A[0][0] * tmp14 + A[0][1] * tmp17 - A[0][2] * tmp19 + A[0][3] * tmp20;

  tmp11 = A[1][2] * tmp1 - A[1][3] * tmp2 + A[1][4] * tmp3;
  tmp12 = A[1][1] * tmp1 - A[1][3] * tmp4 + A[1][4] * tmp5;
  tmp13 = A[1][1] * tmp2 - A[1][2] * tmp4 + A[1][4] * tmp6;
  tmp14 = A[1][1] * tmp3 - A[1][2] * tmp5 + A[1][3] * tmp6;
  tmp15 = A[1][0] * tmp1 - A[1][3] * tmp7 + A[1][4] * tmp8;
  tmp16 = A[1][0] * tmp2 - A[1][2] * tmp7 + A[1][4] * tmp9;
  tmp17 = A[1][0] * tmp3 - A[1][2] * tmp8 + A[1][3] * tmp9;
  tmp18 = A[1][0] * tmp4 - A[1][1] * tmp7 + A[1][4] * tmp10;
  tmp19 = A[1][0] * tmp5 - A[1][1] * tmp8 + A[1][3] * tmp10;

  B[0][2] = A[0][1] * tmp11 - A[0][2] * tmp12 + A[0][3] * tmp13 - A[0][4] * tmp14;
  B[1][2] = -A[0][0] * tmp11 + A[0][2] * tmp15 - A[0][3] * tmp16 + A[0][4] * tmp17;
  B[2][2] = A[0][0] * tmp12 - A[0][1] * tmp15 + A[0][3] * tmp18 - A[0][4] * tmp19;

  tmp1 = A[0][2] * A[1][3] - A[0][3] * A[1][2];
  tmp2 = A[0][1] * A[1][3] - A[0][3] * A[1][1];
  tmp3 = A[0][1] * A[1][2] - A[0][2] * A[1][1];
  tmp4 = A[0][0] * A[1][3] - A[0][3] * A[1][0];
  tmp5 = A[0][0] * A[1][2] - A[0][2] * A[1][0];
  tmp6 = A[0][0] * A[1][1] - A[0][1] * A[1][0];
  tmp7 = A[0][2] * A[1][4] - A[0][4] * A[1][2];
  tmp8 = A[0][1] * A[1][4] - A[0][4] * A[1][1];
  tmp9 = A[0][0] * A[1][4] - A[0][4] * A[1][0];
  tmp10 = A[0][3] * A[1][4] - A[0][4] * A[1][3];

  tmp11 = A[2][2] * tmp10 - A[2][3] * tmp7 + A[2][4] * tmp1;
  tmp12 = A[2][1] * tmp10 - A[2][3] * tmp8 + A[2][4] * tmp2;
  tmp13 = A[2][1] * tmp7 - A[2][2] * tmp8 + A[2][4] * tmp3;
  tmp14 = A[2][1] * tmp1 - A[2][2] * tmp2 + A[2][3] * tmp3;
  tmp15 = A[2][0] * tmp10 - A[2][3] * tmp9 + A[2][4] * tmp4;
  tmp16 = A[2][0] * tmp7 - A[2][2] * tmp9 + A[2][4] * tmp5;
  tmp17 = A[2][0] * tmp1 - A[2][2] * tmp4 + A[2][3] * tmp5;

  B[0][3] = A[4][1] * tmp11 - A[4][2] * tmp12 + A[4][3] * tmp13 - A[4][4] * tmp14;
  B[0][4] = -A[3][1] * tmp11 + A[3][2] * tmp12 - A[3][3] * tmp13 + A[3][4] * tmp14;
  B[1][3] = -A[4][0] * tmp11 + A[4][2] * tmp15 - A[4][3] * tmp16 + A[4][4] * tmp17;
  B[1][4] = A[3][0] * tmp11 - A[3][2] * tmp15 + A[3][3] * tmp16 - A[3][4] * tmp17;

  tmp18 = A[2][0] * tmp8 - A[2][1] * tmp9 + A[2][4] * tmp6;
  tmp19 = A[2][0] * tmp2 - A[2][1] * tmp4 + A[2][3] * tmp6;
  tmp20 = A[2][0] * tmp3 - A[2][1] * tmp5 + A[2][2] * tmp6;

  B[2][3] = A[4][0] * tmp12 - A[4][1] * tmp15 + A[4][3] * tmp18 - A[4][4] * tmp19;
  B[2][4] = -A[3][0] * tmp12 + A[3][1] * tmp15 - A[3][3] * tmp18 + A[3][4] * tmp19;
  B[3][3] = -A[4][0] * tmp13 + A[4][1] * tmp16 - A[4][2] * tmp18 + A[4][4] * tmp20;
  B[3][4] = A[3][0] * tmp13 - A[3][1] * tmp16 + A[3][2] * tmp18 - A[3][4] * tmp20;
  B[4][3] = A[4][0] * tmp14 - A[4][1] * tmp17 + A[4][2] * tmp19 - A[4][3] * tmp20;
  B[4][4] = -A[3][0] * tmp14 + A[3][1] * tmp17 - A[3][2] * tmp19 + A[3][3] * tmp20;

  tmp11 = A[3][1] * tmp7 - A[3][2] * tmp8 + A[3][4] * tmp3;
  tmp12 = A[3][0] * tmp7 - A[3][2] * tmp9 + A[3][4] * tmp5;
  tmp13 = A[3][0] * tmp8 - A[3][1] * tmp9 + A[3][4] * tmp6;
  tmp14 = A[3][0] * tmp3 - A[3][1] * tmp5 + A[3][2] * tmp6;

  tmp15 = A[3][1] * tmp1 - A[3][2] * tmp2 + A[3][3] * tmp3;
  tmp16 = A[3][0] * tmp1 - A[3][2] * tmp4 + A[3][3] * tmp5;
  tmp17 = A[3][0] * tmp2 - A[3][1] * tmp4 + A[3][3] * tmp6;

  B[3][2] = A[4][0] * tmp11 - A[4][1] * tmp12 + A[4][2] * tmp13 - A[4][4] * tmp14;
  B[4][2] = -A[4][0] * tmp15 + A[4][1] * tmp16 - A[4][2] * tmp17 + A[4][3] * tmp14;

  const std::complex<double> det(A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0] + A[0][3] * B[3][0] + A[0][4] * B[4][0]);

  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
      B[i][j] /= det;
}

inline std::complex<double> DetCplx5x5(std::complex<double>** A)
{

  std::complex<double> tmp1(A[3][3] * A[4][4] - A[3][4] * A[4][3]);
  std::complex<double> tmp2(A[3][2] * A[4][4] - A[3][4] * A[4][2]);
  std::complex<double> tmp3(A[3][2] * A[4][3] - A[3][3] * A[4][2]);
  std::complex<double> tmp4(A[3][1] * A[4][4] - A[3][4] * A[4][1]);
  std::complex<double> tmp5(A[3][1] * A[4][3] - A[3][3] * A[4][1]);
  std::complex<double> tmp6(A[3][1] * A[4][2] - A[3][2] * A[4][1]);
  std::complex<double> tmp7(A[3][0] * A[4][4] - A[3][4] * A[4][0]);
  std::complex<double> tmp8(A[3][0] * A[4][3] - A[3][3] * A[4][0]);
  std::complex<double> tmp9(A[3][0] * A[4][2] - A[3][2] * A[4][0]);
  std::complex<double> tmp10(A[3][0] * A[4][1] - A[3][1] * A[4][0]);

  std::complex<double> tmp11(A[2][2] * tmp1 - A[2][3] * tmp2 + A[2][4] * tmp3);
  std::complex<double> tmp12(A[2][1] * tmp1 - A[2][3] * tmp4 + A[2][4] * tmp5);
  std::complex<double> tmp13(A[2][1] * tmp2 - A[2][2] * tmp4 + A[2][4] * tmp6);
  std::complex<double> tmp14(A[2][1] * tmp3 - A[2][2] * tmp5 + A[2][3] * tmp6);
  std::complex<double> tmp15(A[2][0] * tmp1 - A[2][3] * tmp7 + A[2][4] * tmp8);
  std::complex<double> tmp16(A[2][0] * tmp2 - A[2][2] * tmp7 + A[2][4] * tmp9);
  std::complex<double> tmp17(A[2][0] * tmp3 - A[2][2] * tmp8 + A[2][3] * tmp9);

  std::complex<double>b0 = A[1][1] * tmp11 - A[1][2] * tmp12 + A[1][3] * tmp13 - A[1][4] * tmp14;
  std::complex<double>b1 = -A[1][0] * tmp11 + A[1][2] * tmp15 - A[1][3] * tmp16 + A[1][4] * tmp17;

  std::complex<double> tmp18(A[2][0] * tmp4 - A[2][1] * tmp7 + A[2][4] * tmp10);
  std::complex<double> tmp19(A[2][0] * tmp5 - A[2][1] * tmp8 + A[2][3] * tmp10);
  std::complex<double> tmp20(A[2][0] * tmp6 - A[2][1] * tmp9 + A[2][2] * tmp10);

  std::complex<double>b2 = A[1][0] * tmp12 - A[1][1] * tmp15 + A[1][3] * tmp18 - A[1][4] * tmp19;
  std::complex<double>b3 = -A[1][0] * tmp13 + A[1][1] * tmp16 - A[1][2] * tmp18 + A[1][4] * tmp20;
  std::complex<double>b4 = A[1][0] * tmp14 - A[1][1] * tmp17 + A[1][2] * tmp19 - A[1][3] * tmp20;

  std::complex<double> det(A[0][0] * b0 + A[0][1] * b1 + A[0][2] * b2 + A[0][3] * b3 + A[0][4] * b4);
  return det;


}

inline void InvertCplx6x6(std::complex<double>** A, std::complex<double>** B)
{
  std::complex<double> tmp1(A[4][4] * A[5][5] - A[4][5] * A[5][4]);
  std::complex<double> tmp2(A[4][3] * A[5][5] - A[4][5] * A[5][3]);
  std::complex<double> tmp3(A[4][3] * A[5][4] - A[4][4] * A[5][3]);
  std::complex<double> tmp4(A[4][2] * A[5][5] - A[4][5] * A[5][2]);
  std::complex<double> tmp5(A[4][2] * A[5][4] - A[4][4] * A[5][2]);
  std::complex<double> tmp6(A[4][2] * A[5][3] - A[4][3] * A[5][2]);
  std::complex<double> tmp7(A[4][1] * A[5][5] - A[4][5] * A[5][1]);
  std::complex<double> tmp8(A[4][1] * A[5][4] - A[4][4] * A[5][1]);
  std::complex<double> tmp9(A[4][1] * A[5][3] - A[4][3] * A[5][1]);
  std::complex<double> tmp10(A[4][1] * A[5][2] - A[4][2] * A[5][1]);
  std::complex<double> tmp11(A[4][0] * A[5][5] - A[4][5] * A[5][0]);
  std::complex<double> tmp12(A[4][0] * A[5][4] - A[4][4] * A[5][0]);
  std::complex<double> tmp13(A[4][0] * A[5][3] - A[4][3] * A[5][0]);
  std::complex<double> tmp14(A[4][0] * A[5][2] - A[4][2] * A[5][0]);
  std::complex<double> tmp15(A[4][0] * A[5][1] - A[4][1] * A[5][0]);

  std::complex<double> tmp16(A[3][3] * tmp1 - A[3][4] * tmp2 + A[3][5] * tmp3);
  std::complex<double> tmp17(A[3][2] * tmp1 - A[3][4] * tmp4 + A[3][5] * tmp5);
  std::complex<double> tmp18(A[3][2] * tmp2 - A[3][3] * tmp4 + A[3][5] * tmp6);
  std::complex<double> tmp19(A[3][2] * tmp3 - A[3][3] * tmp5 + A[3][4] * tmp6);
  std::complex<double> tmp20(A[3][1] * tmp1 - A[3][4] * tmp7 + A[3][5] * tmp8);
  std::complex<double> tmp21(A[3][1] * tmp2 - A[3][3] * tmp7 + A[3][5] * tmp9);
  std::complex<double> tmp22(A[3][1] * tmp3 - A[3][3] * tmp8 + A[3][4] * tmp9);
  std::complex<double> tmp23(A[3][1] * tmp4 - A[3][2] * tmp7 + A[3][5] * tmp10);
  std::complex<double> tmp24(A[3][1] * tmp5 - A[3][2] * tmp8 + A[3][4] * tmp10);
  std::complex<double> tmp25(A[3][1] * tmp6 - A[3][2] * tmp9 + A[3][3] * tmp10);
  std::complex<double> tmp26(A[3][0] * tmp1 - A[3][4] * tmp11 + A[3][5] * tmp12);
  std::complex<double> tmp27(A[3][0] * tmp2 - A[3][3] * tmp11 + A[3][5] * tmp13);
  std::complex<double> tmp28(A[3][0] * tmp3 - A[3][3] * tmp12 + A[3][4] * tmp13);
  std::complex<double> tmp29(A[3][0] * tmp4 - A[3][2] * tmp11 + A[3][5] * tmp14);
  std::complex<double> tmp30(A[3][0] * tmp5 - A[3][2] * tmp12 + A[3][4] * tmp14);
  std::complex<double> tmp31(A[3][0] * tmp6 - A[3][2] * tmp13 + A[3][3] * tmp14);
  std::complex<double> tmp32(A[3][0] * tmp7 - A[3][1] * tmp11 + A[3][5] * tmp15);
  std::complex<double> tmp33(A[3][0] * tmp8 - A[3][1] * tmp12 + A[3][4] * tmp15);
  std::complex<double> tmp34(A[3][0] * tmp9 - A[3][1] * tmp13 + A[3][3] * tmp15);
  std::complex<double> tmp35(A[3][0] * tmp10 - A[3][1] * tmp14 + A[3][2] * tmp15);

  std::complex<double> tmp36(A[2][2] * tmp16 - A[2][3] * tmp17 + A[2][4] * tmp18 - A[2][5] * tmp19);
  std::complex<double> tmp37(A[2][1] * tmp16 - A[2][3] * tmp20 + A[2][4] * tmp21 - A[2][5] * tmp22);
  std::complex<double> tmp38(A[2][1] * tmp17 - A[2][2] * tmp20 + A[2][4] * tmp23 - A[2][5] * tmp24);
  std::complex<double> tmp39(A[2][1] * tmp18 - A[2][2] * tmp21 + A[2][3] * tmp23 - A[2][5] * tmp25);
  std::complex<double> tmp40(A[2][1] * tmp19 - A[2][2] * tmp22 + A[2][3] * tmp24 - A[2][4] * tmp25);
  std::complex<double> tmp41(A[2][0] * tmp16 - A[2][3] * tmp26 + A[2][4] * tmp27 - A[2][5] * tmp28);
  std::complex<double> tmp42(A[2][0] * tmp17 - A[2][2] * tmp26 + A[2][4] * tmp29 - A[2][5] * tmp30);
  std::complex<double> tmp43(A[2][0] * tmp18 - A[2][2] * tmp27 + A[2][3] * tmp29 - A[2][5] * tmp31);
  std::complex<double> tmp44(A[2][0] * tmp19 - A[2][2] * tmp28 + A[2][3] * tmp30 - A[2][4] * tmp31);

  B[0][0] = A[1][1] * tmp36 - A[1][2] * tmp37 + A[1][3] * tmp38 - A[1][4] * tmp39 + A[1][5] * tmp40;
  B[0][1] = -A[0][1] * tmp36 + A[0][2] * tmp37 - A[0][3] * tmp38 + A[0][4] * tmp39 - A[0][5] * tmp40;
  B[1][0] = -A[1][0] * tmp36 + A[1][2] * tmp41 - A[1][3] * tmp42 + A[1][4] * tmp43 - A[1][5] * tmp44;
  B[1][1] = A[0][0] * tmp36 - A[0][2] * tmp41 + A[0][3] * tmp42 - A[0][4] * tmp43 + A[0][5] * tmp44;

  std::complex<double> tmp45(A[2][0] * tmp20 - A[2][1] * tmp26 + A[2][4] * tmp32 - A[2][5] * tmp33);
  std::complex<double> tmp46(A[2][0] * tmp21 - A[2][1] * tmp27 + A[2][3] * tmp32 - A[2][5] * tmp34);
  std::complex<double> tmp47(A[2][0] * tmp22 - A[2][1] * tmp28 + A[2][3] * tmp33 - A[2][4] * tmp34);
  std::complex<double> tmp48(A[2][0] * tmp23 - A[2][1] * tmp29 + A[2][2] * tmp32 - A[2][5] * tmp35);
  std::complex<double> tmp49(A[2][0] * tmp24 - A[2][1] * tmp30 + A[2][2] * tmp33 - A[2][4] * tmp35);

  B[2][0] = A[1][0] * tmp37 - A[1][1] * tmp41 + A[1][3] * tmp45 - A[1][4] * tmp46 + A[1][5] * tmp47;
  B[2][1] = -A[0][0] * tmp37 + A[0][1] * tmp41 - A[0][3] * tmp45 + A[0][4] * tmp46 - A[0][5] * tmp47;
  B[3][0] = -A[1][0] * tmp38 + A[1][1] * tmp42 - A[1][2] * tmp45 + A[1][4] * tmp48 - A[1][5] * tmp49;
  B[3][1] = A[0][0] * tmp38 - A[0][1] * tmp42 + A[0][2] * tmp45 - A[0][4] * tmp48 + A[0][5] * tmp49;

  std::complex<double> tmp50(A[2][0] * tmp25 - A[2][1] * tmp31 + A[2][2] * tmp34 - A[2][3] * tmp35);

  B[4][0] = A[1][0] * tmp39 - A[1][1] * tmp43 + A[1][2] * tmp46 - A[1][3] * tmp48 + A[1][5] * tmp50;
  B[4][1] = -A[0][0] * tmp39 + A[0][1] * tmp43 - A[0][2] * tmp46 + A[0][3] * tmp48 - A[0][5] * tmp50;
  B[5][0] = -A[1][0] * tmp40 + A[1][1] * tmp44 - A[1][2] * tmp47 + A[1][3] * tmp49 - A[1][4] * tmp50;
  B[5][1] = A[0][0] * tmp40 - A[0][1] * tmp44 + A[0][2] * tmp47 - A[0][3] * tmp49 + A[0][4] * tmp50;

  tmp36 = A[1][2] * tmp16 - A[1][3] * tmp17 + A[1][4] * tmp18 - A[1][5] * tmp19;
  tmp37 = A[1][1] * tmp16 - A[1][3] * tmp20 + A[1][4] * tmp21 - A[1][5] * tmp22;
  tmp38 = A[1][1] * tmp17 - A[1][2] * tmp20 + A[1][4] * tmp23 - A[1][5] * tmp24;
  tmp39 = A[1][1] * tmp18 - A[1][2] * tmp21 + A[1][3] * tmp23 - A[1][5] * tmp25;
  tmp40 = A[1][1] * tmp19 - A[1][2] * tmp22 + A[1][3] * tmp24 - A[1][4] * tmp25;
  tmp41 = A[1][0] * tmp16 - A[1][3] * tmp26 + A[1][4] * tmp27 - A[1][5] * tmp28;
  tmp42 = A[1][0] * tmp17 - A[1][2] * tmp26 + A[1][4] * tmp29 - A[1][5] * tmp30;
  tmp43 = A[1][0] * tmp18 - A[1][2] * tmp27 + A[1][3] * tmp29 - A[1][5] * tmp31;
  tmp44 = A[1][0] * tmp19 - A[1][2] * tmp28 + A[1][3] * tmp30 - A[1][4] * tmp31;
  tmp45 = A[1][0] * tmp20 - A[1][1] * tmp26 + A[1][4] * tmp32 - A[1][5] * tmp33;
  tmp46 = A[1][0] * tmp21 - A[1][1] * tmp27 + A[1][3] * tmp32 - A[1][5] * tmp34;
  tmp47 = A[1][0] * tmp22 - A[1][1] * tmp28 + A[1][3] * tmp33 - A[1][4] * tmp34;
  tmp48 = A[1][0] * tmp23 - A[1][1] * tmp29 + A[1][2] * tmp32 - A[1][5] * tmp35;
  tmp49 = A[1][0] * tmp24 - A[1][1] * tmp30 + A[1][2] * tmp33 - A[1][4] * tmp35;
  tmp50 = A[1][0] * tmp25 - A[1][1] * tmp31 + A[1][2] * tmp34 - A[1][3] * tmp35;

  B[0][2] = A[0][1] * tmp36 - A[0][2] * tmp37 + A[0][3] * tmp38 - A[0][4] * tmp39 + A[0][5] * tmp40;
  B[1][2] = -A[0][0] * tmp36 + A[0][2] * tmp41 - A[0][3] * tmp42 + A[0][4] * tmp43 - A[0][5] * tmp44;
  B[2][2] = A[0][0] * tmp37 - A[0][1] * tmp41 + A[0][3] * tmp45 - A[0][4] * tmp46 + A[0][5] * tmp47;
  B[3][2] = -A[0][0] * tmp38 + A[0][1] * tmp42 - A[0][2] * tmp45 + A[0][4] * tmp48 - A[0][5] * tmp49;
  B[4][2] = A[0][0] * tmp39 - A[0][1] * tmp43 + A[0][2] * tmp46 - A[0][3] * tmp48 + A[0][5] * tmp50;
  B[5][2] = -A[0][0] * tmp40 + A[0][1] * tmp44 - A[0][2] * tmp47 + A[0][3] * tmp49 - A[0][4] * tmp50;

  tmp1 = A[0][3] * A[1][4] - A[0][4] * A[1][3];
  tmp2 = A[0][2] * A[1][4] - A[0][4] * A[1][2];
  tmp3 = A[0][2] * A[1][3] - A[0][3] * A[1][2];
  tmp4 = A[0][1] * A[1][4] - A[0][4] * A[1][1];
  tmp5 = A[0][1] * A[1][3] - A[0][3] * A[1][1];
  tmp6 = A[0][1] * A[1][2] - A[0][2] * A[1][1];
  tmp7 = A[0][0] * A[1][4] - A[0][4] * A[1][0];
  tmp8 = A[0][0] * A[1][3] - A[0][3] * A[1][0];
  tmp9 = A[0][0] * A[1][2] - A[0][2] * A[1][0];
  tmp10 = A[0][0] * A[1][1] - A[0][1] * A[1][0];
  tmp11 = A[0][3] * A[1][5] - A[0][5] * A[1][3];
  tmp12 = A[0][2] * A[1][5] - A[0][5] * A[1][2];
  tmp13 = A[0][1] * A[1][5] - A[0][5] * A[1][1];
  tmp14 = A[0][0] * A[1][5] - A[0][5] * A[1][0];
  tmp15 = A[0][4] * A[1][5] - A[0][5] * A[1][4];

  tmp16 = A[2][3] * tmp15 - A[2][4] * tmp11 + A[2][5] * tmp1;
  tmp17 = A[2][2] * tmp15 - A[2][4] * tmp12 + A[2][5] * tmp2;
  tmp18 = A[2][2] * tmp11 - A[2][3] * tmp12 + A[2][5] * tmp3;
  tmp19 = A[2][2] * tmp1 - A[2][3] * tmp2 + A[2][4] * tmp3;
  tmp20 = A[2][1] * tmp15 - A[2][4] * tmp13 + A[2][5] * tmp4;
  tmp21 = A[2][1] * tmp11 - A[2][3] * tmp13 + A[2][5] * tmp5;
  tmp22 = A[2][1] * tmp1 - A[2][3] * tmp4 + A[2][4] * tmp5;
  tmp23 = A[2][1] * tmp12 - A[2][2] * tmp13 + A[2][5] * tmp6;
  tmp24 = A[2][1] * tmp2 - A[2][2] * tmp4 + A[2][4] * tmp6;
  tmp25 = A[2][1] * tmp3 - A[2][2] * tmp5 + A[2][3] * tmp6;
  tmp26 = A[2][0] * tmp15 - A[2][4] * tmp14 + A[2][5] * tmp7;
  tmp27 = A[2][0] * tmp11 - A[2][3] * tmp14 + A[2][5] * tmp8;
  tmp28 = A[2][0] * tmp1 - A[2][3] * tmp7 + A[2][4] * tmp8;
  tmp29 = A[2][0] * tmp12 - A[2][2] * tmp14 + A[2][5] * tmp9;
  tmp30 = A[2][0] * tmp2 - A[2][2] * tmp7 + A[2][4] * tmp9;
  tmp31 = A[2][0] * tmp3 - A[2][2] * tmp8 + A[2][3] * tmp9;
  tmp32 = A[2][0] * tmp13 - A[2][1] * tmp14 + A[2][5] * tmp10;
  tmp33 = A[2][0] * tmp4 - A[2][1] * tmp7 + A[2][4] * tmp10;
  tmp34 = A[2][0] * tmp5 - A[2][1] * tmp8 + A[2][3] * tmp10;
  tmp35 = A[2][0] * tmp6 - A[2][1] * tmp9 + A[2][2] * tmp10;

  tmp36 = A[3][2] * tmp16 - A[3][3] * tmp17 + A[3][4] * tmp18 - A[3][5] * tmp19;
  tmp37 = A[3][1] * tmp16 - A[3][3] * tmp20 + A[3][4] * tmp21 - A[3][5] * tmp22;
  tmp38 = A[3][1] * tmp17 - A[3][2] * tmp20 + A[3][4] * tmp23 - A[3][5] * tmp24;
  tmp39 = A[3][1] * tmp18 - A[3][2] * tmp21 + A[3][3] * tmp23 - A[3][5] * tmp25;
  tmp40 = A[3][1] * tmp19 - A[3][2] * tmp22 + A[3][3] * tmp24 - A[3][4] * tmp25;
  tmp41 = A[3][0] * tmp16 - A[3][3] * tmp26 + A[3][4] * tmp27 - A[3][5] * tmp28;
  tmp42 = A[3][0] * tmp17 - A[3][2] * tmp26 + A[3][4] * tmp29 - A[3][5] * tmp30;
  tmp43 = A[3][0] * tmp18 - A[3][2] * tmp27 + A[3][3] * tmp29 - A[3][5] * tmp31;
  tmp44 = A[3][0] * tmp19 - A[3][2] * tmp28 + A[3][3] * tmp30 - A[3][4] * tmp31;

  B[0][4] = -A[5][1] * tmp36 + A[5][2] * tmp37 - A[5][3] * tmp38 + A[5][4] * tmp39 - A[5][5] * tmp40;
  B[0][5] = A[4][1] * tmp36 - A[4][2] * tmp37 + A[4][3] * tmp38 - A[4][4] * tmp39 + A[4][5] * tmp40;
  B[1][4] = A[5][0] * tmp36 - A[5][2] * tmp41 + A[5][3] * tmp42 - A[5][4] * tmp43 + A[5][5] * tmp44;
  B[1][5] = -A[4][0] * tmp36 + A[4][2] * tmp41 - A[4][3] * tmp42 + A[4][4] * tmp43 - A[4][5] * tmp44;

  tmp45 = A[3][0] * tmp20 - A[3][1] * tmp26 + A[3][4] * tmp32 - A[3][5] * tmp33;
  tmp46 = A[3][0] * tmp21 - A[3][1] * tmp27 + A[3][3] * tmp32 - A[3][5] * tmp34;
  tmp47 = A[3][0] * tmp22 - A[3][1] * tmp28 + A[3][3] * tmp33 - A[3][4] * tmp34;
  tmp48 = A[3][0] * tmp23 - A[3][1] * tmp29 + A[3][2] * tmp32 - A[3][5] * tmp35;
  tmp49 = A[3][0] * tmp24 - A[3][1] * tmp30 + A[3][2] * tmp33 - A[3][4] * tmp35;

  B[2][4] = -A[5][0] * tmp37 + A[5][1] * tmp41 - A[5][3] * tmp45 + A[5][4] * tmp46 - A[5][5] * tmp47;
  B[2][5] = A[4][0] * tmp37 - A[4][1] * tmp41 + A[4][3] * tmp45 - A[4][4] * tmp46 + A[4][5] * tmp47;
  B[3][4] = A[5][0] * tmp38 - A[5][1] * tmp42 + A[5][2] * tmp45 - A[5][4] * tmp48 + A[5][5] * tmp49;
  B[3][5] = -A[4][0] * tmp38 + A[4][1] * tmp42 - A[4][2] * tmp45 + A[4][4] * tmp48 - A[4][5] * tmp49;

  tmp50 = A[3][0] * tmp25 - A[3][1] * tmp31 + A[3][2] * tmp34 - A[3][3] * tmp35;

  B[4][4] = -A[5][0] * tmp39 + A[5][1] * tmp43 - A[5][2] * tmp46 + A[5][3] * tmp48 - A[5][5] * tmp50;
  B[4][5] = A[4][0] * tmp39 - A[4][1] * tmp43 + A[4][2] * tmp46 - A[4][3] * tmp48 + A[4][5] * tmp50;
  B[5][4] = A[5][0] * tmp40 - A[5][1] * tmp44 + A[5][2] * tmp47 - A[5][3] * tmp49 + A[5][4] * tmp50;
  B[5][5] = -A[4][0] * tmp40 + A[4][1] * tmp44 - A[4][2] * tmp47 + A[4][3] * tmp49 - A[4][4] * tmp50;

  tmp36 = A[4][2] * tmp16 - A[4][3] * tmp17 + A[4][4] * tmp18 - A[4][5] * tmp19;
  tmp37 = A[4][1] * tmp16 - A[4][3] * tmp20 + A[4][4] * tmp21 - A[4][5] * tmp22;
  tmp38 = A[4][1] * tmp17 - A[4][2] * tmp20 + A[4][4] * tmp23 - A[4][5] * tmp24;
  tmp39 = A[4][1] * tmp18 - A[4][2] * tmp21 + A[4][3] * tmp23 - A[4][5] * tmp25;
  tmp40 = A[4][1] * tmp19 - A[4][2] * tmp22 + A[4][3] * tmp24 - A[4][4] * tmp25;
  tmp41 = A[4][0] * tmp16 - A[4][3] * tmp26 + A[4][4] * tmp27 - A[4][5] * tmp28;
  tmp42 = A[4][0] * tmp17 - A[4][2] * tmp26 + A[4][4] * tmp29 - A[4][5] * tmp30;
  tmp43 = A[4][0] * tmp18 - A[4][2] * tmp27 + A[4][3] * tmp29 - A[4][5] * tmp31;
  tmp44 = A[4][0] * tmp19 - A[4][2] * tmp28 + A[4][3] * tmp30 - A[4][4] * tmp31;
  tmp45 = A[4][0] * tmp20 - A[4][1] * tmp26 + A[4][4] * tmp32 - A[4][5] * tmp33;
  tmp46 = A[4][0] * tmp21 - A[4][1] * tmp27 + A[4][3] * tmp32 - A[4][5] * tmp34;
  tmp47 = A[4][0] * tmp22 - A[4][1] * tmp28 + A[4][3] * tmp33 - A[4][4] * tmp34;
  tmp48 = A[4][0] * tmp23 - A[4][1] * tmp29 + A[4][2] * tmp32 - A[4][5] * tmp35;
  tmp49 = A[4][0] * tmp24 - A[4][1] * tmp30 + A[4][2] * tmp33 - A[4][4] * tmp35;
  tmp50 = A[4][0] * tmp25 - A[4][1] * tmp31 + A[4][2] * tmp34 - A[4][3] * tmp35;

  B[0][3] = A[5][1] * tmp36 - A[5][2] * tmp37 + A[5][3] * tmp38 - A[5][4] * tmp39 + A[5][5] * tmp40;
  B[1][3] = -A[5][0] * tmp36 + A[5][2] * tmp41 - A[5][3] * tmp42 + A[5][4] * tmp43 - A[5][5] * tmp44;
  B[2][3] = A[5][0] * tmp37 - A[5][1] * tmp41 + A[5][3] * tmp45 - A[5][4] * tmp46 + A[5][5] * tmp47;
  B[3][3] = -A[5][0] * tmp38 + A[5][1] * tmp42 - A[5][2] * tmp45 + A[5][4] * tmp48 - A[5][5] * tmp49;
  B[4][3] = A[5][0] * tmp39 - A[5][1] * tmp43 + A[5][2] * tmp46 - A[5][3] * tmp48 + A[5][5] * tmp50;
  B[5][3] = -A[5][0] * tmp40 + A[5][1] * tmp44 - A[5][2] * tmp47 + A[5][3] * tmp49 - A[5][4] * tmp50;

  const std::complex<double> det(A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0] +
    A[0][3] * B[3][0] + A[0][4] * B[4][0] + A[0][5] * B[5][0]);

  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      B[i][j] /= det;
}

inline std::complex<double> DetCplx6x6(std::complex<double>** A)
{

  std::complex<double> tmp1(A[4][4] * A[5][5] - A[4][5] * A[5][4]);
  std::complex<double> tmp2(A[4][3] * A[5][5] - A[4][5] * A[5][3]);
  std::complex<double> tmp3(A[4][3] * A[5][4] - A[4][4] * A[5][3]);
  std::complex<double> tmp4(A[4][2] * A[5][5] - A[4][5] * A[5][2]);
  std::complex<double> tmp5(A[4][2] * A[5][4] - A[4][4] * A[5][2]);
  std::complex<double> tmp6(A[4][2] * A[5][3] - A[4][3] * A[5][2]);
  std::complex<double> tmp7(A[4][1] * A[5][5] - A[4][5] * A[5][1]);
  std::complex<double> tmp8(A[4][1] * A[5][4] - A[4][4] * A[5][1]);
  std::complex<double> tmp9(A[4][1] * A[5][3] - A[4][3] * A[5][1]);
  std::complex<double> tmp10(A[4][1] * A[5][2] - A[4][2] * A[5][1]);
  std::complex<double> tmp11(A[4][0] * A[5][5] - A[4][5] * A[5][0]);
  std::complex<double> tmp12(A[4][0] * A[5][4] - A[4][4] * A[5][0]);
  std::complex<double> tmp13(A[4][0] * A[5][3] - A[4][3] * A[5][0]);
  std::complex<double> tmp14(A[4][0] * A[5][2] - A[4][2] * A[5][0]);
  std::complex<double> tmp15(A[4][0] * A[5][1] - A[4][1] * A[5][0]);

  std::complex<double> tmp16(A[3][3] * tmp1 - A[3][4] * tmp2 + A[3][5] * tmp3);
  std::complex<double> tmp17(A[3][2] * tmp1 - A[3][4] * tmp4 + A[3][5] * tmp5);
  std::complex<double> tmp18(A[3][2] * tmp2 - A[3][3] * tmp4 + A[3][5] * tmp6);
  std::complex<double> tmp19(A[3][2] * tmp3 - A[3][3] * tmp5 + A[3][4] * tmp6);
  std::complex<double> tmp20(A[3][1] * tmp1 - A[3][4] * tmp7 + A[3][5] * tmp8);
  std::complex<double> tmp21(A[3][1] * tmp2 - A[3][3] * tmp7 + A[3][5] * tmp9);
  std::complex<double> tmp22(A[3][1] * tmp3 - A[3][3] * tmp8 + A[3][4] * tmp9);
  std::complex<double> tmp23(A[3][1] * tmp4 - A[3][2] * tmp7 + A[3][5] * tmp10);
  std::complex<double> tmp24(A[3][1] * tmp5 - A[3][2] * tmp8 + A[3][4] * tmp10);
  std::complex<double> tmp25(A[3][1] * tmp6 - A[3][2] * tmp9 + A[3][3] * tmp10);
  std::complex<double> tmp26(A[3][0] * tmp1 - A[3][4] * tmp11 + A[3][5] * tmp12);
  std::complex<double> tmp27(A[3][0] * tmp2 - A[3][3] * tmp11 + A[3][5] * tmp13);
  std::complex<double> tmp28(A[3][0] * tmp3 - A[3][3] * tmp12 + A[3][4] * tmp13);
  std::complex<double> tmp29(A[3][0] * tmp4 - A[3][2] * tmp11 + A[3][5] * tmp14);
  std::complex<double> tmp30(A[3][0] * tmp5 - A[3][2] * tmp12 + A[3][4] * tmp14);
  std::complex<double> tmp31(A[3][0] * tmp6 - A[3][2] * tmp13 + A[3][3] * tmp14);
  std::complex<double> tmp32(A[3][0] * tmp7 - A[3][1] * tmp11 + A[3][5] * tmp15);
  std::complex<double> tmp33(A[3][0] * tmp8 - A[3][1] * tmp12 + A[3][4] * tmp15);
  std::complex<double> tmp34(A[3][0] * tmp9 - A[3][1] * tmp13 + A[3][3] * tmp15);
  std::complex<double> tmp35(A[3][0] * tmp10 - A[3][1] * tmp14 + A[3][2] * tmp15);

  std::complex<double> tmp36(A[2][2] * tmp16 - A[2][3] * tmp17 + A[2][4] * tmp18 - A[2][5] * tmp19);
  std::complex<double> tmp37(A[2][1] * tmp16 - A[2][3] * tmp20 + A[2][4] * tmp21 - A[2][5] * tmp22);
  std::complex<double> tmp38(A[2][1] * tmp17 - A[2][2] * tmp20 + A[2][4] * tmp23 - A[2][5] * tmp24);
  std::complex<double> tmp39(A[2][1] * tmp18 - A[2][2] * tmp21 + A[2][3] * tmp23 - A[2][5] * tmp25);
  std::complex<double> tmp40(A[2][1] * tmp19 - A[2][2] * tmp22 + A[2][3] * tmp24 - A[2][4] * tmp25);
  std::complex<double> tmp41(A[2][0] * tmp16 - A[2][3] * tmp26 + A[2][4] * tmp27 - A[2][5] * tmp28);
  std::complex<double> tmp42(A[2][0] * tmp17 - A[2][2] * tmp26 + A[2][4] * tmp29 - A[2][5] * tmp30);
  std::complex<double> tmp43(A[2][0] * tmp18 - A[2][2] * tmp27 + A[2][3] * tmp29 - A[2][5] * tmp31);
  std::complex<double> tmp44(A[2][0] * tmp19 - A[2][2] * tmp28 + A[2][3] * tmp30 - A[2][4] * tmp31);

  std::complex<double>b0 = A[1][1] * tmp36 - A[1][2] * tmp37 + A[1][3] * tmp38 - A[1][4] * tmp39 + A[1][5] * tmp40;
  std::complex<double>b1 = -A[1][0] * tmp36 + A[1][2] * tmp41 - A[1][3] * tmp42 + A[1][4] * tmp43 - A[1][5] * tmp44;

  std::complex<double> tmp45(A[2][0] * tmp20 - A[2][1] * tmp26 + A[2][4] * tmp32 - A[2][5] * tmp33);
  std::complex<double> tmp46(A[2][0] * tmp21 - A[2][1] * tmp27 + A[2][3] * tmp32 - A[2][5] * tmp34);
  std::complex<double> tmp47(A[2][0] * tmp22 - A[2][1] * tmp28 + A[2][3] * tmp33 - A[2][4] * tmp34);
  std::complex<double> tmp48(A[2][0] * tmp23 - A[2][1] * tmp29 + A[2][2] * tmp32 - A[2][5] * tmp35);
  std::complex<double> tmp49(A[2][0] * tmp24 - A[2][1] * tmp30 + A[2][2] * tmp33 - A[2][4] * tmp35);

  std::complex<double>b2 = A[1][0] * tmp37 - A[1][1] * tmp41 + A[1][3] * tmp45 - A[1][4] * tmp46 + A[1][5] * tmp47;
  std::complex<double>b3 = -A[1][0] * tmp38 + A[1][1] * tmp42 - A[1][2] * tmp45 + A[1][4] * tmp48 - A[1][5] * tmp49;
  std::complex<double> tmp50(A[2][0] * tmp25 - A[2][1] * tmp31 + A[2][2] * tmp34 - A[2][3] * tmp35);

  std::complex<double>b4 = A[1][0] * tmp39 - A[1][1] * tmp43 + A[1][2] * tmp46 - A[1][3] * tmp48 + A[1][5] * tmp50;
  std::complex<double>b5 = -A[1][0] * tmp40 + A[1][1] * tmp44 - A[1][2] * tmp47 + A[1][3] * tmp49 - A[1][4] * tmp50;

  std::complex<double> det(A[0][0] * b0 + A[0][1] * b1 + A[0][2] * b2 + A[0][3] * b3 + A[0][4] * b4 + A[0][5] * b5);
  return det;

}


#endif