extern crate feanor_la;
extern crate mosek;

use feanor_la::combinatorics::iters::*;
use feanor_la::combinatorics::bitset::*;
use feanor_la::la::mat::*;
use feanor_la::algebra::rat::*;

fn find_feasible_solution<M, V>(A: Matrix<M, r64>, b: Vector<V, r64>) -> Option<Vector<VectorOwned<f64>, f64>>
    where M: MatrixView<r64>, V: VectorView<r64>
{
    let cast = |x: r64| (x.num() as f64) / (x.den() as f64);

    let mut env = mosek::Env::new().unwrap();
    let mut task = env.task().unwrap();
    task.append_vars(A.col_count() as i32).unwrap();
    task.append_cons(A.row_count() as i32).unwrap();
    for j in 0..A.col_count() {
        task.put_c_j(j as i32, 1.).unwrap();
        task.put_var_bound(j as i32, mosek::MSK_BK_LO, 0.0, f64::INFINITY).unwrap();
    }
    for i in 0..A.row_count() {
        let indices = A.row(i).iter().enumerate().filter(|(_, x)| **x != r64::ZERO).map(|(i, _)| i as i32).collect::<Vec<_>>();
        let values = A.row(i).iter().filter(|x| **x != r64::ZERO).map(|x| cast(*x)).collect::<Vec<_>>();
        task.put_a_row(i as i32, &indices, &values).unwrap();
        task.put_con_bound(i as i32, mosek::MSK_BK_FX, cast(*b.at(i)), cast(*b.at(i))).unwrap();
    }
    task.optimize().unwrap();
    if task.get_sol_sta(0).unwrap() == mosek::MSK_SOL_STA_PRIM_INFEAS_CER {
        return None;
    };
    let mut result = (0..A.col_count()).map(|_| 0.0).collect::<Vec<_>>();
    task.get_xx(0, &mut result).unwrap();
    return Some(Vector::from_fn(result.len(), |i| result[i]));
}

fn main() {
    let edges = |L: &[Box<[usize]>], U: &[Box<[usize]>]| cartesian_product(0..L.len(), 0..U.len())
        .filter(|&(a, b)| {
            (0..L[a].len()).all(|j| L[a][j] <= U[b][j])
        }).collect::<Vec<_>>();

    let k = 4;
    let d = 7;
    let superset = (0..d).map(|_| k).collect::<Vec<_>>();
    let layers: Vec<_> = (0..(k * d)).map(|i| {
        multiset_combinations(&superset, i, clone_slice).collect::<Vec<_>>()
    }).collect();

    for i in 1..(k * d) {
        let edges = edges(&layers[i - 1], &layers[i]);
        let e = edges.len();
        let mut A: Matrix<_, r64> = Matrix::zero(layers[i - 1].len() + layers[i].len() - 1, e + 1).into_owned();
        {
            let mut lower_conds = A.submatrix_mut(..(layers[i - 1].len() - 1), 1..);
            for (j, (u, _)) in edges.iter().enumerate() {
                if *u == 0 {
                    lower_conds.col_mut(j).assign(Vector::constant(layers[i - 1].len() - 1, -r64::ONE));
                } else {
                    *lower_conds.at_mut(u - 1, j) = r64::ONE;
                }
            }
        }
        {
            let mut upper_conds = A.submatrix_mut((layers[i - 1].len() - 1)..(A.row_count() - 1), 1..);
            for (j, (_, v)) in edges.iter().enumerate() {
                if *v == 0 {
                    upper_conds.col_mut(j).assign(Vector::constant(layers[i].len() - 1, -r64::ONE));
                } else {
                    *upper_conds.at_mut(v - 1, j) = r64::ONE;
                }
            }
        }
        A.row_mut(A.row_count() - 1).assign(Vector::constant(e + 1, r64::ONE));
        let solution = find_feasible_solution(A.submatrix(.., 1..), A.col(0));
        println!("(k, d) = ({}, {}), i = {}:        sol: {}", k, d, i, solution.is_some());
    }
}