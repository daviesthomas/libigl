// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2019 Thomas Davies <thomasryantrd@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_BIHARMONIC_DISTANCE_H
#define IGL_BIHARMONIC_DISTANCE_H

#include "igl_inline.h"
#include <Eigen/Core>

namespace igl 
{
  // Biharmonic Distance algorithm
  // https://gfx.cs.princeton.edu/pubs/Lipman_2010_BD/index.php
  //
  // Inputs:
  //   V  #V by 3 list of 3D vertex positions
  //   F  #F by 3 list of mesh faces
  //   VS #VS by 1 vector specifying indice of source vertices, vertice to calcuate distance to
  //   p  exponent above eigen values (default to 2 for biharmonic)

  // Output:
  //   D  #V by 1 vector of biharmonic distances of each target w.r.t source
  //
  template <typename DerivedV,typename DerivedF,typename DerivedD>
  IGL_INLINE void biharmonic_distance(
    const Eigen::MatrixBase<DerivedV> &V,
    const Eigen::MatrixBase<DerivedF> &F,
    const int i,
    const double p,
    Eigen::PlainObjectBase<DerivedD> &D);
};

#ifndef IGL_STATIC_LIBRARY
#  include "biharmonic_distance.cpp"
#endif

#endif