extern crate nalgebra as na;

use na::{Complex, SMatrix};
use std::{f64::consts::PI, fs::File, io::prelude::*};

const NX: usize = 101;
const NY: usize = 101;
const XRANGE: f64 = 15.0;
const YRANGE: f64 = 15.0;

type GridMatrix = SMatrix<Complex<f64>, NY, NX>;

fn solve_mean_flow(u: f64, alpha: f64) -> GridMatrix {
    let velocity: Complex<f64> = Complex::new(u * alpha.cos(), u * alpha.sin());
    let new_velocity = GridMatrix::from_element(velocity);
    return new_velocity;
}

fn _solve_vortex(grid: GridMatrix, m: f64, z0: Complex<f64>) -> GridMatrix {
    let mut new_velocity = GridMatrix::zeros();
    for i in 0..NX {
        for j in 0..NY {
            new_velocity[(j, i)] = Complex::i() * (m / (2. * PI)) * (1. / (grid[(j, i)] - z0));
        }
    }
    return new_velocity;
}

fn solve_source_sink(grid: GridMatrix, m: f64, z0: Complex<f64>) -> GridMatrix {
    let mut new_velocity = GridMatrix::zeros();
    for i in 0..NX {
        for j in 0..NY {
            new_velocity[(j, i)] = (-m / (2. * PI)) * (1. / (grid[(j, i)] - z0));
        }
    }
    return new_velocity;
}

fn gen_grid(nx: usize, ny: usize, xrange: f64, yrange: f64) -> GridMatrix {
    let mut output_grid = GridMatrix::zeros();
    let dx: f64 = (2. * xrange) / (nx as f64 - 1.);
    let dy: f64 = (2. * yrange) / (ny as f64 - 1.);
    for i in 0..nx {
        for j in 0..ny {
            let xp: f64 = -xrange + dx * i as f64;
            let yp: f64 = -yrange + dy * j as f64;
            output_grid[(j, i)] = Complex::new(yp, xp);
        }
    }

    return output_grid;
}

fn main() {
    let mut velocity_field = GridMatrix::zeros();
    let zgrid: GridMatrix = gen_grid(NX, NY, XRANGE, YRANGE);

    velocity_field +=
        solve_source_sink(zgrid, 5.0, Complex::new(1.0, 0.0)) + solve_mean_flow(1.50, 0.0);
    let mut vfile = File::create("./velfield.txt").expect("Failed to create file!");
    for i in 0..NX {
        for j in 0..NY {
            write!(vfile, "{}\t", velocity_field[(j, i)]).expect("Failed to write")
        }
        write!(vfile, "\n").expect("Failed to create new line");
    }
}
