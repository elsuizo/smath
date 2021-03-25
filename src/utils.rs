//-------------------------------------------------------------------------
// @file utils.rs
//
// @date 03/25/21 12:34:16
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
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MaxMin<T> {
    pub max: (T, usize),
    pub min: (T, usize),
}

/// generic function to fin min, max values and the position in a slice
pub fn find_max_min<T: std::cmp::PartialOrd + Copy>(slice: &[T]) -> MaxMin<T> {
    let mut max = &slice[0];
    let mut min = &slice[0];

    let mut max_pos: usize = 0;
    let mut min_pos: usize = 0;

    for index in 1..slice.len() {
        if slice[index] < *min { min = &slice[index]; min_pos = index;}
        if slice[index] > *max { max = &slice[index]; max_pos = index;}
    }

    MaxMin{max: (*max, max_pos), min: (*min, min_pos)}
}


