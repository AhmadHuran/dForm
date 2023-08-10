use anyhow::{Result, anyhow};
use ndarray::{Array2, s};

/*
 *                      Code generated with SymPy 1.10.1
 *
 *              See http://www.sympy.org/ for more information.
 *
 */
fn a_inv_det(
    a11: f64,
    a12: f64,
    a13: f64,
    a22: f64,
    a23: f64,
    a33: f64,
    ) -> Vec<f64> {

    let out1 = a22*a33 - a23.powi(2);
    let out2 = -a12*a33 + a13*a23;
    let out3 = a12*a23 - a13*a22;
    let out4 = -a12*a33 + a13*a23;
    let out5 = a11*a33 - a13.powi(2);
    let out6 = -a11*a23 + a12*a13;
    let out7 = a12*a23 - a13*a22;
    let out8 = -a11*a23 + a12*a13;
    let out9 = a11*a22 - a12.powi(2);

    vec![out1, out2, out3, out4, out5, out6, out7, out8, out9]
}

/*
 *                      Code generated with SymPy 1.10.1
 *
 *              See http://www.sympy.org/ for more information.
 */
fn a_det(
    a11: f64,
    a12: f64,
    a13: f64,
    a22: f64,
    a23: f64,
    a33: f64,
    ) -> f64 {
    let out1 = a11 * a22 * a33
        - a11 * a23.powi(2)
        - a12.powi(2) * a33
        + 2.0 * a12*a13*a23
        - a13.powi(2) * a22;
    out1
}

pub fn inverse_3x3(arr: &Array2<f64>) -> Result<Array2<f64>> {
    let dim1 = arr.len_of(ndarray::Axis(0));
    let dim2 = arr.len_of(ndarray::Axis(1));
    if dim1 != dim2 || dim1 > 3 {
        return Err(anyhow!("inverse_3x3: supports at most 3x3"))
    }

    let det: f64;
    let inv_det: Array2<f64>;

    match dim1 {
        1 => {
            det = arr[[0, 0]];
            inv_det = Array2::ones((1,1));
        },

        2 => {
            let a11 = arr[[0, 0]];
            let a12 = arr[[0, 1]];
            let a13 = 0.0;

            let a22 = arr[[1, 1]];
            let a23 = 0.0;

            let a33 = 1.0;

            det = a_det(a11, a12, a13, a22, a23, a33);
            inv_det = Array2::from_shape_vec(
                (3,3),
                a_inv_det(a11, a12, a13, a22, a23, a33)
                )?.slice_move(s![..2, ..2]);
        },


        3 => {
            let a11 = arr[[0, 0]];
            let a12 = arr[[0, 1]];
            let a13 = arr[[0, 2]];

            let a22 = arr[[1, 1]];
            let a23 = arr[[1, 2]];

            let a33 = arr[[2, 2]];

            det = a_det(a11, a12, a13, a22, a23, a33);
            inv_det = Array2::from_shape_vec(
                (3,3),
                a_inv_det(a11, a12, a13, a22, a23, a33)
                )?;
        },

        0 => return Err(anyhow!("inverse_3x3: Cannot recover!")),
        _ => return Err(anyhow!("inverse_3x3: something went really bad. Cannot recover!"))

    };

    if det.abs() < f64::EPSILON {
        return Err(anyhow!("inverse_3x3: matrix is sigular"))
    };

    Ok(inv_det/det)
}
