// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2019 Thomas Davies <thomasryantrd@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include "biharmonic_distance.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include "eigs.h"
#include <iostream>

template <typename DerivedV,typename DerivedF,typename DerivedD>
IGL_INLINE void igl::biharmonic_distance(
  const Eigen::MatrixBase<DerivedV> &V,
  const Eigen::MatrixBase<DerivedF> &F,
  const int i,
  const double p,
  Eigen::PlainObjectBase<DerivedD> &D)
{
  
  Eigen::SparseMatrix<double> L,M;
  Eigen::MatrixXd eU; //eigen vectors 
  Eigen::VectorXd eS; //eigen values 

  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
  const size_t k = 10;

  if(!eigs(L,M,k+1,igl::EIGS_TYPE_SM,eU,eS))
  {
    std::cout<<"eigen decomp failed."<<std::endl;
  }

  Eigen::MatrixXd eSD = eS.asDiagonal();
  eSD = eSD.block(1,1,eSD.rows()-1,eSD.rows()-1);
  eU = eU.middleCols(1,eU.cols()-1);
  
  //biharmonic embedding
  Eigen::MatrixXd B = eU * eSD.array().abs().matrix().inverse();
  //create query array from source point
  D = (B.row(i).replicate(B.rows(),1) - B)
                          .array().square()
                          .matrix().rowwise().sum()
                          .cwiseSqrt();

}