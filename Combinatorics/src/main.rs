extern crate feanor_la;
extern crate mosek;

use feanor_la::combinatorics::iters::*;
use feanor_la::combinatorics::bitset::*;
use feanor_la::la::mat::*;
use feanor_la::algebra::rat::*;

///
/// Uses the LP solver to find x >= 0
/// with Ax = b or states that this does
/// not exist.
/// 
fn find_feasible_solution<M, V>(
    A: Matrix<M, r64>, 
    b: Vector<V, r64>
) -> Option<Vector<VectorOwned<f64>, f64>>
    where M: MatrixView<r64>, V: VectorView<r64>
{
    let cast = |x: r64| (x.num() as f64) / (x.den() as f64);

    let mut env = mosek::Env::new().unwrap();
    let mut task = env.task().unwrap();
    task.append_vars(A.col_count() as i32).unwrap();
    task.append_cons(A.row_count() as i32).unwrap();

    // initialize variable bounds and target function
    for j in 0..A.col_count() {
        task.put_c_j(j as i32, 1.).unwrap();
        task.put_var_bound(
            j as i32, mosek::MSK_BK_LO, 0.0, f64::INFINITY
        ).unwrap();
    }
    // initialize conditions
    for i in 0..A.row_count() {
        // we need sparse format
        let indices = A.row(i).iter().enumerate()
            .filter(|(_, x)| **x != r64::ZERO)
            .map(|(i, _)| i as i32)
            .collect::<Vec<_>>();

        let values = A.row(i).iter()
            .filter(|x| **x != r64::ZERO)
            .map(|x| cast(*x)).collect::<Vec<_>>();

        task.put_a_row(i as i32, &indices, &values).unwrap();
        task.put_con_bound(
            i as i32, mosek::MSK_BK_FX, 
            cast(*b.at(i)), 
            cast(*b.at(i))
        ).unwrap();
    }
    task.optimize().unwrap();
    let state = task.get_sol_sta(0).unwrap();
    if state == mosek::MSK_SOL_STA_PRIM_INFEAS_CER {
        // problem is infeasible
        return None;
    };
    // parse the result
    let mut result = (0..A.col_count())
        .map(|_| 0.0)
        .collect::<Vec<_>>();
    task.get_xx(0, &mut result).unwrap();
    return Some(Vector::from_fn(result.len(), |i| result[i]));
}

///
/// Returns a matrix defining conditions on a
/// vector indexed by u < v, u in L, v in U
/// where U, V are two consecutive layers in 
/// a graded poset.
/// 
/// Builds the matrix (b | A) with an mxn matrix A
/// (m is the number of pairs u < v) representing
/// the conditions:
///  - sum over all u < v is 1
///  - for u in L the sum over all u < v is
///    constant
///  - for v in U the sum over all u < v is
///    constant
/// 
fn build_matrix<T>(
    L: &Vec<T>, 
    U: &Vec<T>,
    edges: Vec<(usize, usize)>
) -> Matrix<MatrixOwned<r64>, r64> {
    let e = edges.len();
    let mut A: Matrix<_, r64> = Matrix::zero(
        L.len() + U.len() - 1, e + 1
    ).into_owned();
    {
        // add the conditions for the lower layer L
        let mut L_conds = A.submatrix_mut(..(L.len() - 1), 1..);
        for (j, (u, _)) in edges.iter().enumerate() {
            if *u == 0 {
                L_conds.col_mut(j).assign(
                    Vector::constant(L.len() - 1, -r64::ONE)
                );
            } else {
                *L_conds.at_mut(u - 1, j) = r64::ONE;
            }
        }
    }
    {
        // add the conditions for the upper layer U
        let mut U_conds = A.submatrix_mut(
            (L.len() - 1)..(A.row_count() - 1), 1..
        );
        for (j, (_, v)) in edges.iter().enumerate() {
            if *v == 0 {
                U_conds.col_mut(j).assign(
                    Vector::constant(U.len() - 1, -r64::ONE)
                );
            } else {
                *U_conds.at_mut(v - 1, j) = r64::ONE;
            }
        }
    }
    A.row_mut(A.row_count() - 1).assign(
        Vector::constant(e + 1, r64::ONE)
    );
    return A;
}

fn simple_lp_alg(k: usize, d: usize) {

    // a function to find the edges between two layers
    let edges = |L: &[Box<[usize]>], U: &[Box<[usize]>]| 
        cartesian_product(0..L.len(), 0..U.len())
        .filter(|&(a, b)| {
            (0..L[a].len()).all(|j| L[a][j] <= U[b][j])
        }).collect::<Vec<_>>();

    let superset = (0..d).map(|_| k).collect::<Vec<_>>();
    let layers: Vec<_> = (0..(k * d)).map(|i|
        multiset_combinations(&superset, i, clone_slice)
        .collect::<Vec<_>>()
    ).collect();

    for i in 1..=(k * d) {
        let edges = edges(&layers[i - 1], &layers[i]);
        let A = build_matrix(&layers[i - 1], &layers[i], edges);
        let solution = find_feasible_solution(
            A.submatrix(.., 1..), 
            A.col(0)
        );
        println!(
            "(k, d) = ({}, {}), i = {}:        sol: {}", 
            k, d, i, solution.is_some()
        );
    }
}

///
/// Returns a matrix defining conditions on a
/// vector indexed by u < v, u in L, v in U
/// where U, V are two consecutive layers in 
/// a graded poset.
/// 
/// Builds the matrix (b | A) with an mxn matrix A
/// (m is the number of pairs u < v) representing 
/// the conditions:
///  - sum over all pairs u < v is 1
///  - for u in L, the sum over all u < v
///    is fraction(u)
///  - for v in U, the sum over all u < v
///    is fraction(v)
/// 
fn build_symmetry_using_matrix<T, F>(
    L: &Vec<T>, 
    U: &Vec<T>, 
    edges: Vec<(usize, usize)>, 
    mut fraction: F
) -> Matrix<MatrixOwned<r64>, r64> 
    where F: FnMut(&T) -> r64
{
    let e = edges.len();
    let mut A: Matrix<_, r64> = Matrix::zero(
        L.len() + U.len() + 1, e + 1
    ).into_owned();
    {
        for u in 0..L.len() {
            *A.at_mut(u, 0) = fraction(&L[u]);
        }
        let mut L_conds = A.submatrix_mut(..L.len(), 1..);
        for (j, (u, _)) in edges.iter().enumerate() {
            *L_conds.at_mut(*u, j) = r64::ONE;
        }
    }
    {
        for v in 0..U.len() {
            *A.at_mut(L.len() + v, 0) = fraction(&U[v]);
        }
        let mut U_conds = A.submatrix_mut(
            L.len()..(A.row_count() - 1), 1..
        );
        for (j, (_, v)) in edges.iter().enumerate() {
            *U_conds.at_mut(*v, j) = r64::ONE;
        }
    }
    A.row_mut(A.row_count() - 1).assign(
        Vector::constant(e + 1, r64::ONE)
    );
    return A;
}

fn factorial(x: u128) -> u128 {
    if x == 0 { 1 } else { x * factorial(x - 1) }
}

fn multinomial_coefficient(d: usize, n: &[usize]) -> usize {
    (
        factorial(d as u128) / 
        n.iter()
            .map(|x| factorial(*x as u128))
            .product::<u128>()
    ) as usize
}

fn symmetry_using_lp_alg(k: usize, d: usize) {
    let bar_positions = (0..=d).map(|_| k).collect::<Vec<_>>();

    // this generates all (n_0, ..., n_k) with sum(n_i) = d
    let core_poset = multiset_combinations(
        &bar_positions[..], k, 
        |pos: &[usize]| {
            let mut result = pos.iter().enumerate()
                .flat_map(|(i, c)| (0..*c).map(move |_| i))
                .chain(std::iter::once(d))
                .collect::<Vec<_>>();
            for i in (1..=k).rev() {
                result[i] -= result[i - 1];
            }
            return result;
        }
    );

    // group the elements of C_{k, d} into layers
    let mut layers = (0..=(k * d))
        .map(|_| Vec::new())
        .collect::<Vec<_>>();

    for p in core_poset {
        let layer = p.iter().enumerate()
            .map(|(i, x)| i * *x)
            .sum::<usize>();

        layers[layer].push(p);
    }

    // a function to find the edges between two layers
    let edges = |L: &Vec<Vec<usize>>, U: &Vec<Vec<usize>>|
        cartesian_product(0..L.len(), 0..U.len())
        .filter(|&(a, b)| {
            (0..k).any(|delta| 
                U[b][delta] + 1 == L[a][delta] && 
                U[b][delta + 1] == L[a][delta + 1] + 1
            )
        }).collect::<Vec<_>>();
    
    for i in 1..=(k * d) {
        let edges = edges(&layers[i - 1], &layers[i]);
        let A = build_symmetry_using_matrix(
            &layers[i - 1], &layers[i], edges, |element| {
                let layer = element.iter().enumerate()
                    .map(|(i, x)| i * *x)
                    .sum::<usize>();
                let total_size = layers[layer].iter()
                    .map(|e| multinomial_coefficient(d, e))
                    .sum::<usize>();
                let size = multinomial_coefficient(d, element);
                return r64::new(size as i64, total_size as i64);
            }
        );
        let solution = find_feasible_solution(
            A.submatrix(.., 1..), A.col(0)
        );
        println!(
            "(k, d) = ({}, {}), i = {}:        sol: {}", 
            k, d, i, solution.is_some()
        );
    }
}

fn main() {
    symmetry_using_lp_alg(8, 9);
}