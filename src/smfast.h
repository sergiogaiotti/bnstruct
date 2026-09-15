#ifndef _SMFAST_H_
#define _SMFAST_H_

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include "util.h"

/* find best parents */
SEXP bnstruct_fbp(SEXP aflml);

/* find best sinks */
SEXP bnstruct_fbs(SEXP bps, SEXP aflml);

/* find best sinks and return subset scores */
SEXP bnstruct_fbs_with_scores(SEXP bps, SEXP aflml);

/* fast upwards moebius transform */
SEXP bnstruct_fumt_mask(SEXP n_elements, SEXP pattern);

/* computes log marginal likelihood for all possible families */
SEXP bnstruct_all_fam_log_marg_lik(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP iss, SEXP func);

/* compute BDeu leave-one-out bounds for all families */
SEXP bnstruct_bdeu_loo_bounds(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP ess);
SEXP bnstruct_score_loo_bounds(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP ess, SEXP func);

/* compute conservative Bg bounds */
SEXP bnstruct_compute_bg_cons(SEXP lower, SEXP upper, SEXP imp_fam_mask);

/* compute tight Bg bounds */
SEXP bnstruct_compute_bg_tight(SEXP aflml, SEXP bps, SEXP lower, SEXP upper, SEXP imp_fam_mask);

/* compute asymmetric Lg/Ug bounds (conservative) */
SEXP bnstruct_compute_lg_ug_cons(SEXP bps, SEXP lower, SEXP upper, SEXP imp_fam_mask);

/* compute asymmetric Lg/Ug bounds (tight) */
SEXP bnstruct_compute_lg_ug_tight(SEXP aflml, SEXP bps, SEXP lower, SEXP upper, SEXP imp_fam_mask);

/* compute asymmetric LF/UF bounds over subsets */
SEXP bnstruct_compute_lf_uf(SEXP F, SEXP sinks, SEXP bps, SEXP aflml, SEXP Lg, SEXP Ug);

/* incremental exact leave-one-out: optimal parent sets per removed observation */
SEXP bnstruct_loo_parents(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP ess, SEXP aflml);
SEXP bnstruct_loo_parents_score(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP ess, SEXP aflml, SEXP func);

/* coupled parent margins for every (node, candidate set) and exact per-family coupling for the sink condition */
SEXP bnstruct_coupled_parent_margins(SEXP aflml, SEXP bps, SEXP imp_fam_mask, SEXP node_sizes, SEXP lower_k, SEXP upper_k, SEXP lower, SEXP upper);
SEXP bnstruct_coupling_rows(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP ess, SEXP func, SEXP refs);
SEXP bnstruct_parent_margins_rows(SEXP aflml, SEXP bps, SEXP imp_fam_mask, SEXP values, SEXP refs, SEXP ref_index);

/* compute BF bounds over subsets */
SEXP bnstruct_compute_bf(SEXP Bg);

/*
	very internal function
*/

/* fill _bits_ with the indices of ones in _word_ (which can be at most _size_)
	and return the number of ones */
unsigned int get_bits(unsigned int word, unsigned int *bits, unsigned int size);

#endif
