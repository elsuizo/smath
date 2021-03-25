//-------------------------------------------------------------------------
// @file smatrix.rs
//
// @date 03/25/21 12:32:43
// @author Martin Noblia
// @email mnoblia@disroot.org
//
// @brief
//
// @detail
//
//  Licence:
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
//--------------------------------------------------------------------------
use std::ops::{Deref, DerefMut, Index, IndexMut};
use std::ops::{Add, Mul, Sub, AddAssign, SubAssign};

use num::{Num, Zero, One};
use crate::traits::LinearAlgebra;

//-------------------------------------------------------------------------
//                        code
//-------------------------------------------------------------------------
/// A generic static matrix of MxN shape
#[derive(Debug, PartialEq)]
pub struct SMatrix<T, const M: usize, const N: usize>([[T; N]; M]);

impl<T, const M: usize, const N: usize> SMatrix<T, M, N> {
    pub fn new(data_input: [[T; N]; M]) -> Self {
        Self(data_input)
    }

    pub const fn shape(&self) -> (usize, usize) {
        (M, N)
    }
}

//-------------------------------------------------------------------------
//                        indexing
//-------------------------------------------------------------------------
impl<T, const M: usize, const N: usize> Index<(usize, usize)> for SMatrix<T, M, N> {

    type Output = T;

    fn index(&self, index: (usize, usize)) -> &T {
        &self.0[index.0][index.1]
    }
}

impl<T, const M: usize, const N: usize> IndexMut<(usize, usize)> for SMatrix<T, M, N> {

    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        &mut self.0[index.0][index.1]
    }
}

//-------------------------------------------------------------------------
//              rectangular numerics matrix family
//-------------------------------------------------------------------------
impl<T: Num + Copy + AddAssign, const N: usize> SMatrix<T, N, N> {

    pub fn zeros() -> Self {
        <Self as Zero>::zero()
    }

    pub fn identity() -> Self {
        <Self as One>::one()
    }
}

// Zero
impl<T: Num + Copy + AddAssign, const N: usize> Zero for SMatrix<T, N, N> {

    fn zero() -> SMatrix<T, N, N> {
        Self::new([[T::zero(); N]; N])
    }

    fn is_zero(&self) -> bool {
        *self == Self::zero()
    }
}

// One
impl<T: Num + Copy + AddAssign, const N: usize> One for SMatrix<T, N, N> {
    /// create a identity Matrix
    fn one() -> SMatrix<T, N, N> {
        let one = T::one();
        let mut result = Self::zeros();
        for i in 0..N {
            for j in 0..N {
                if i == j {
                    result[(i, j)] = one;
                }
            }
        }
        result
    }
}

// Add
impl<T: Num + Copy + AddAssign, const N: usize> Add for SMatrix<T, N, N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut result = Self::zero();
        for i in 0..N {
            for j in 0..N {
                result[(i, j)] = self[(i, j)] + rhs[(i, j)];
            }
        }
        result
    }
}

// Sub
impl<T: Num + Copy + AddAssign, const N: usize> Sub for SMatrix<T, N, N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let mut result = Self::zeros();
        for i in 0..N {
            for j in 0..N {
                result[(i, j)] = self[(i, j)] - rhs[(i, j)];
            }
        }
        result
    }
}

// TODO(elsuizo:2020-09-04): this is a naive impl of Mul
// Mul
impl<T: Num + Copy + AddAssign, const N: usize> Mul for SMatrix<T, N, N> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut result = Self::zeros();
        for i in 0..N {
            for j in 0..N {
                for k in 0..N {
                    result[(i, j)] += self[(i, k)] * rhs[(k, j)];
                }
            }
        }
        result
    }
}

impl<T: Num + Copy + AddAssign, const N: usize> LinearAlgebra<T> for SMatrix<T, N, N> {

    fn det(&self) -> T {
        match self.shape() {
            (2, 2) => {
                (self[(0, 0)] * self[(1, 1)]) - (self[(1, 0)] * self[(0, 1)])
            }
            (3, 3) => {
                let a_00 = self[(0, 0)];
                let a_01 = self[(0, 1)];
                let a_02 = self[(0, 2)];
                let a_10 = self[(1, 0)];
                let a_11 = self[(1, 1)];
                let a_12 = self[(1, 2)];
                let a_20 = self[(2, 0)];
                let a_21 = self[(2, 1)];
                let a_22 = self[(2, 2)];

                a_00 * (a_11 * a_22 - a_21 * a_12)
                - a_01 * (a_10 * a_22 - a_12 * a_20)
                + a_02 * (a_10 * a_21 - a_11 * a_20)
            }
            (4, 4) => {
                let a1 = self[(0, 0)];
                let a2 = self[(0, 1)];
                let a3 = self[(0, 2)];
                let a4 = self[(0, 3)];
                let a5 = self[(1, 0)];
                let a6 = self[(1, 1)];
                let a7 = self[(1, 2)];
                let a8 = self[(1, 3)];
                let a9 = self[(2, 0)];
                let a10 = self[(2, 1)];
                let a11 = self[(2, 2)];
                let a12 = self[(2, 3)];
                let a13 = self[(3, 0)];
                let a14 = self[(3, 1)];
                let a15 = self[(3, 2)];
                let a16 = self[(3, 3)];

                a1 * a10 * a15 * a8
                - a1 * a10 * a16 * a7
                - a1 * a11 * a14 * a8
                + a1 * a11 * a16 * a6
                + a1 * a12 * a14 * a7
                - a1 * a12 * a15 * a6
                - a10 * a13 * a3 * a8
                + a10 * a13 * a4 * a7
                - a10 * a15 * a4 * a5
                + a10 * a16 * a3 * a5
                + a11 * a13 * a2 * a8
                - a11 * a13 * a4 * a6
                + a11 * a14 * a4 * a5
                - a11 * a16 * a2 * a5
                - a12 * a13 * a2 * a7
                + a12 * a13 * a3 * a6
                - a12 * a14 * a3 * a5
                + a12 * a15 * a2 * a5
                + a14 * a3 * a8 * a9
                - a14 * a4 * a7 * a9
                - a15 * a2 * a8 * a9
                + a15 * a4 * a6 * a9
                + a16 * a2 * a7 * a9
                - a16 * a3 * a6 * a9
            }
            (_, _) => panic!("under construction...sorry")
        }
    }
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod matrix_tests {
    use crate::traits::LinearAlgebra;
    use crate::smatrix::SMatrix;

    #[test]
    fn test_matrix_zeros() {
        let m: SMatrix<f64, 260, 260> = SMatrix::zeros();

        for i in 0..m.shape().0 {
            for j in 0..m.shape().1 {
                assert_eq!(m[(i, j)], 0.0);
            }
        }
    }

    #[test]
    fn m22_det_test() {
        let m = SMatrix::new([[1, 2], [3, 4]]);
        let result = m.det();
        let expected = -2;
        assert_eq!(result, expected);
    }

    #[test]
    fn m33_det_test() {
        let m = SMatrix::new([[0, 1, 2], [3, 4, 5], [6, 7, 9]]);
        let result = m.det();
        let expected = -3;
        assert_eq!(result, expected);
    }

    #[test]
    fn m44_det_test() {
        let m = SMatrix::new([[1, 2, 3, 1],
                             [5, 6, 7, 8],
                             [9, 0, 11, 12],
                             [13,1, 15, 16]]);
        let result = m.det();
        let expected = 168;
        assert_eq!(result, expected);
    }
}


