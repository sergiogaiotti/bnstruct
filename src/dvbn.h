#ifndef DVBN_H
#define DVBN_H

#include <R.h>
#include <Rinternals.h>

/* Algorithm 1: discretize one continuous variable (cuts only).
 * Arguments:
 *  data_cont:   numeric matrix (n_cases x n_cont) of continuous variables (column-major)
 *  data_disc:   integer matrix (n_cases x n_nodes) of all variables in discrete code (1..ns[j], NA_INTEGER for NA)
 *  n_nodes:     integer scalar (#columns in data_disc)
 *  n_cases:     integer scalar (#rows)
 *  ns:          integer vector length n_nodes (cardinalities; values for continuous vars can be 0 or ignored)
 *  cont_col:    integer scalar (0-based column in data_cont to discretize)
 *  target:      integer scalar (0-based original index in data_disc corresponding to cont_col)
 *  parents_list:list length n_nodes; each element is an integer vector (0-based parent indices for that var)
 *  children_list:list length n_nodes; each element is an integer vector (0-based children indices for that var)
 *  approx_parents: integer scalar (0/1) — if 1, approximate parents term as sum over single-parents (as in Julia option)
 *  l_card:       integer scalar (lambda parameter; if <=0, set to max cardinality of markov blanket of the continuous variable)
 *
 * Returns: REALSXP vector of cut points (midpoints), excluding min/max; length = #cuts (can be 0).
 */
SEXP bnstruct_dvbn_discretize_one(SEXP data_cont, SEXP data_disc,
                                  SEXP n_nodes, SEXP n_cases,
                                  SEXP ns, SEXP cont_col, SEXP target,
                                  SEXP parents_list, SEXP children_list,
                                  SEXP approx_parents, SEXP l_card);

/* Algorithm 2: discretize all continuous variables iteratively until convergence or n_cycles (columns of X in reverse topological order for better performance).
 * Arguments:
 *  data_cont, data_disc, n_nodes, n_cases, ns as above
 *  n_cont:       integer scalar (#columns in data_cont)
 *  cont_index:   integer vector length n_cont mapping data_cont col j -> original index in data_disc
 *  parents_list, children_list: lists as above
 *  n_cycles:     integer scalar (max iterations)
 *  approx_parents: integer scalar (0/1)
 *  l_card:       integer scalar (lambda parameter for Algorithm 1; if <=0, set to max cardinality of markov blanket of the continuous variable)
 *
 * Returns: list length n_cont, each element is REALSXP vector of cut points (min and max excluded) for that continuous variable.
 */
SEXP bnstruct_dvbn_discretize_all(SEXP data_cont, SEXP data_disc,
                                  SEXP n_nodes, SEXP n_cases,
                                  SEXP n_cont, SEXP cont_index,
                                  SEXP ns, SEXP parents_list, SEXP children_list,
                                  SEXP n_cycles, SEXP approx_parents, SEXP l_card);

#endif /* DVBN_H */
