//-------------------------------------------------------------------------
// @file svector.rs
//
// @date 03/25/21 12:33:46
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
use num::{Float, Zero, Num, Signed};
use std::ops::{Add, Sub, Div, Mul, SubAssign, AddAssign, Neg};
use std::ops::{Deref, DerefMut};
use std::iter::Sum;
use crate::utils::find_max_min;

/// generic static array
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct SVector<T, const N: usize>([T; N]);

impl<T, const N: usize> SVector<T, N> {
    pub fn new(input: [T; N]) -> Self {
        Self(input)
    }

    pub const fn shape(&self) -> usize {
        N
    }
}

// norm2 this impl requiere a Float because sqrt
impl<T: Float, const N: usize> SVector<T, N> {
    pub fn norm2(&self) -> T {
        self.into_iter().fold(T::zero(), |n, &i| (i * i) + n).sqrt()
    }
}

// norm infinity
impl<T: Num + Copy + std::cmp::PartialOrd, const N: usize> SVector<T, N> {
    pub fn norm_inf(&self) -> T {
        let max_min = find_max_min(&self.0);
        max_min.max.0
    }
}

// norm l
impl<T: Num + Copy + Signed + std::iter::Sum, const N: usize> SVector<T, N> {
    pub fn norm_l(&self) -> T {
        self.into_iter().map(|element| element.abs()).sum()
    }
}

// projection
impl<T: Num + Copy + Sum, const N: usize> SVector<T, N> {
    pub fn project_over_y(&self, y: &Self) -> T {
        (*self) * (*y) / ((*y) * (*y))
    }
}
// Zero trait
impl<T: Num + Copy, const N: usize> Zero for SVector<T, N> {
    fn zero() -> SVector<T, N> {
        SVector::new([T::zero(); N])
    }

    fn is_zero(&self) -> bool {
        *self == SVector::zero()
    }
}

// SVector<T, N> + SVector<T, N>
impl<T: Num + Copy, const N: usize> Add for SVector<T, N> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut result = Self::zeros();
        for index in 0..N {
            result[index] = self[index] + rhs[index];
        }
        result
    }
}

// SVector<T, N> += SVector<T, N>
impl<T: Num + Copy, const N: usize> AddAssign for SVector<T, N> {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other
    }
}

// NOTE(elsuizo:2021-01-03): capaz que seria mejor hacerlo con iterators
// Vector<T, N> - Vector<T, N>
impl<T: Num + Copy, const N: usize> Sub for SVector<T, N> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        let mut result = Self::zeros();
        for index in 0..N {
            result[index] = self[index] - rhs[index];
        }
        result
    }
}

// SVector<T, N> -= SVector<T, N>
impl<T: Num + Copy, const N: usize> SubAssign for SVector<T, N> {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other
    }
}

impl<T: Num + Copy + Signed, const N: usize> Neg for SVector<T, N> {
    type Output = Self;

    fn neg(self) -> Self {
        let mut result = Self::zeros();
            for index in 0..N {
                result[index] = -self[index];
            }
        result
    }
}

// SVector * SVector (dot product)
impl<T: Num + Copy + Sum, const N: usize> Mul for SVector<T, N> {
    type Output = T;

    fn mul(self, rhs: Self) -> T {
        self.into_iter().zip(rhs.into_iter()).map(|(&a, &b)| a * b).sum()
    }
}

// SVector * constant
impl<T: Num + Copy, const N: usize> Mul<T> for SVector<T, N> {
    type Output = Self;

    fn mul(self, rhs: T) -> SVector<T, N> {
        let mut result = Self::zeros();
        for index in 0..N {
            result[index] = self[index] * rhs;
        }
        result
    }
}

// f32 * SVector<f32, N>
impl<const N: usize> Mul<SVector<f32, N>> for f32 {
    type Output = SVector<f32, N>;

    fn mul(self, rhs: SVector<f32, N>) -> SVector<f32, N> {
        let mut result = SVector::zeros();
        for index in 0..N {
            result[index] = rhs[index] * self;
        }
        result
    }
}

// i32 * SVector<i32, N>
impl<const N: usize> Mul<SVector<i32, N>> for i32 {
    type Output = SVector<i32, N>;

    fn mul(self, rhs: SVector<i32, N>) -> SVector<i32, N> {
        let mut result = SVector::zeros();
        for index in 0..N {
            result[index] = rhs[index] * self;
        }
        result
    }
}

// u32 * SVector<u32, N>
impl<const N: usize> Mul<SVector<u32, N>> for u32 {
    type Output = SVector<u32, N>;

    fn mul(self, rhs: SVector<u32, N>) -> SVector<u32, N> {
        let mut result = SVector::zeros();
        for index in 0..N {
            result[index] = rhs[index] * self;
        }
        result
    }
}

impl<T, const N: usize> Deref for SVector<T, N> {
    type Target = [T; N];
    #[inline]
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T, const N: usize> DerefMut for SVector<T, N> {
    #[inline]
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T, const N: usize> From<[T; N]> for SVector<T, N> {
    fn from(data: [T; N]) -> Self {
        Self(data)
    }
}

// NOTE(elsuizo:2021-01-01): nice!!!, podemos definir metodos para una familia de
// parametros dejando uno fijo, en este caso la dimension del vector
/// Cross product
impl<T: Num + Copy> SVector<T, 3> {
    pub fn cross(&self, rhs: Self) -> Self {
        let u_x = rhs[0];
        let u_y = rhs[1];
        let u_z = rhs[2];

        Self::new([u_y * self[2] - u_z * self[1],
                   u_z * self[0] - u_x * self[2],
                   u_x * self[1] - u_y * self[0]])
    }
}

impl<T: Num + Copy, const N: usize> SVector<T, N> {
    /// create a Vector with all elements zero
    pub fn zeros() -> Self {
        <Self as Zero>::zero()
    }

    /// create a Vector with all elements one
    pub fn ones() -> Self {
        let one = T::one();
        Self::new([one; N])
    }
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod test_vector {
    use super::SVector;
    #[test]
    fn vector_creation_test() {
        let v = SVector::new([1, 2, 3]);
        let v2 = SVector::new([true, false, false, true]);

        assert_eq!(v[0], 1);
        assert_eq!(v[1], 2);
        assert_eq!(v[2], 3);

        assert_eq!(v2[0], true);
        assert_eq!(v2[1], false);
        assert_eq!(v2[2], false);
        assert_eq!(v2[3], true);
    }

    #[test]
    fn negation_test() {
        let v = SVector::new([1, 2, 3, 4]);
        let result = -v;
        let expected = SVector::new([-1, -2, -3, -4]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn add_assign_test() {
        let mut result = SVector::new([1.0, 3.0]);
        let v2 = SVector::new([3.0, 3.0]);
        let expected = SVector::new([4.0, 6.0]);
        result += v2;
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn sub_assign_test() {
        let mut result = SVector::new([1.0, 3.0]);
        let v2 = SVector::new([3.0, 3.0]);
        let expected = SVector::new([-2.0, 0.0]);
        result -= v2;
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn cross_test() {
        let x = SVector::new([1.0, 0.0, 0.0]);
        let y = SVector::new([0.0, 1.0, 0.0]);

        let result = y.cross(x);
        // z
        let expected = SVector::new([0.0, 0.0, 1.0]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn product_test() {
        let x = SVector::new([1, 2, 3]);
        let y = SVector::new([1, 1, 1]);

        let result = x * y;
        let expected = 6;
        assert_eq!(result, expected);
    }

    #[test]
    fn norm2_test() {
        let x = SVector::new([2.0, 2.0, 2.0, 2.0]);
        let result = x.norm2();
        let expected = 4.0;
        assert_eq!(result, expected);
    }

    #[test]
    fn product_constant_test() {
        let x = SVector::new([1, 2, 3, 4, 5]);
        let result = x * 2;
        let expected = SVector::new([2, 4, 6, 8, 10]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }

    #[test]
    fn product_constant_svector() {
        let x = SVector::new([1u32, 2, 3, 4, 5]);
        let result = 2 * x;
        let expected = SVector::new([2, 4, 6, 8, 10]);
        assert_eq!(
            &result[..],
            &expected[..],
            "\nExpected\n{:?}\nfound\n{:?}",
            &result[..],
            &expected[..]
        );
    }
}


