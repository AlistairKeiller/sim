use nalgebra::{SMatrix,Complex};

const SIZE: usize = 50;
const DT: f64 = 1.0;
const DX: f64 = 1.0;

type MatrixSIZEc = SMatrix<Complex<f64>, SIZE, SIZE>;

fn taylor_exp(x: MatrixSIZEc, y: MatrixSIZEc, v: MatrixSIZEc, n: usize) -> (MatrixSIZEc, MatrixSIZEc, MatrixSIZEc) {
    let mut current_x_in_taylor = x.clone();
    let mut current_y_in_taylor = y.clone();
    let mut current_v_in_taylor = v.clone();
    let mut final_x: MatrixSIZEc = x;
    let mut final_y: MatrixSIZEc = SMatrix::identity()+y;
    let mut final_v: MatrixSIZEc = v;

    for i in 2..n {
        current_x_in_taylor = (current_x_in_taylor*x+y*current_x_in_taylor+v.component_mul(&current_x_in_taylor)).unscale(i as f64);
        current_y_in_taylor = (current_y_in_taylor*x+y*current_y_in_taylor+v.component_mul(&current_y_in_taylor)).unscale(i as f64);
        current_v_in_taylor = (current_v_in_taylor*x+y*current_v_in_taylor+v.component_mul(&current_v_in_taylor)).unscale(i as f64);
        final_x += current_x_in_taylor;
        final_y += current_y_in_taylor;
        final_v += current_v_in_taylor;
    }

    (final_x,final_y,final_v)
}

fn main() {
    // let x = MatrixSIZEc::from_fn(|i,j| if i == j { Complex::new(0.0,-DX*DX) } else if i.abs_diff(j) == 1 { Complex::new(0.0,DX*DX/2.0) } else { Complex::new(0.0,0.0) });
    let (x,y,v) = taylor_exp(
        MatrixSIZEc::from_fn(|i,j| if i == j { Complex::new(0.0,-DX*DX) } else if i.abs_diff(j) == 1 { Complex::new(0.0,DX*DX/2.0) } else { Complex::new(0.0,0.0) }),
        MatrixSIZEc::from_fn(|i,j| if i == j { Complex::new(0.0,-DX*DX) } else if i.abs_diff(j) == 1 { Complex::new(0.0,DX*DX/2.0) } else { Complex::new(0.0,0.0) }),
        MatrixSIZEc::from_fn(|i,j| Complex::new(0.0,0.0)),
        20
    );
    println!("{}",x);
}