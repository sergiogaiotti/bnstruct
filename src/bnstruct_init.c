#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP bnstruct_all_fam_log_marg_lik(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_compute_counts(SEXP, SEXP);
extern SEXP bnstruct_compute_counts_nas(SEXP, SEXP);
extern SEXP bnstruct_fbp(SEXP);
extern SEXP bnstruct_fbs(SEXP, SEXP);
extern SEXP bnstruct_fbs_with_scores(SEXP, SEXP);
extern SEXP bnstruct_bdeu_loo_bounds(SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_score_loo_bounds(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_all_fam_scores_bounds(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_compute_bg_cons(SEXP, SEXP, SEXP);
extern SEXP bnstruct_compute_bg_tight(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_compute_lg_ug_cons(SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_compute_lg_ug_tight(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_compute_lf_uf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_loo_parents(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_loo_parents_score(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_coupled_parent_margins(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_coupling_rows(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_parent_margins_rows(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_compute_bf(SEXP);
extern SEXP bnstruct_fumt_mask(SEXP, SEXP);
extern SEXP bnstruct_g2_stat(SEXP, SEXP);
extern SEXP bnstruct_heom_dist(SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_in_tabu(SEXP, SEXP);
extern SEXP bnstruct_is_acyclic(SEXP);
extern SEXP bnstruct_next_comb(SEXP, SEXP);
extern SEXP bnstruct_score_node(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP bnstruct_dvbn_discretize_one(SEXP, SEXP,
                                         SEXP, SEXP,
                                         SEXP, SEXP, SEXP,
                                         SEXP, SEXP,
                                         SEXP, SEXP);

extern SEXP bnstruct_dvbn_discretize_all(SEXP, SEXP,
                                         SEXP, SEXP,
                                         SEXP, SEXP,
                                         SEXP, SEXP, SEXP,
                                         SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"bnstruct_all_fam_log_marg_lik", (DL_FUNC)&bnstruct_all_fam_log_marg_lik, 5},
    {"bnstruct_compute_counts", (DL_FUNC)&bnstruct_compute_counts, 2},
    {"bnstruct_compute_counts_nas", (DL_FUNC)&bnstruct_compute_counts_nas, 2},
    {"bnstruct_fbp", (DL_FUNC)&bnstruct_fbp, 1},
    {"bnstruct_fbs", (DL_FUNC)&bnstruct_fbs, 2},
    {"bnstruct_fbs_with_scores", (DL_FUNC)&bnstruct_fbs_with_scores, 2},
    {"bnstruct_bdeu_loo_bounds", (DL_FUNC)&bnstruct_bdeu_loo_bounds, 4},
    {"bnstruct_score_loo_bounds", (DL_FUNC)&bnstruct_score_loo_bounds, 5},
    {"bnstruct_all_fam_scores_bounds", (DL_FUNC)&bnstruct_all_fam_scores_bounds, 5},
    {"bnstruct_compute_bg_cons", (DL_FUNC)&bnstruct_compute_bg_cons, 3},
    {"bnstruct_compute_bg_tight", (DL_FUNC)&bnstruct_compute_bg_tight, 5},
    {"bnstruct_compute_lg_ug_cons", (DL_FUNC)&bnstruct_compute_lg_ug_cons, 4},
    {"bnstruct_compute_lg_ug_tight", (DL_FUNC)&bnstruct_compute_lg_ug_tight, 5},
    {"bnstruct_compute_lf_uf", (DL_FUNC)&bnstruct_compute_lf_uf, 6},
    {"bnstruct_loo_parents", (DL_FUNC)&bnstruct_loo_parents, 5},
    {"bnstruct_loo_parents_score", (DL_FUNC)&bnstruct_loo_parents_score, 6},
    {"bnstruct_coupled_parent_margins", (DL_FUNC)&bnstruct_coupled_parent_margins, 8},
    {"bnstruct_coupling_rows", (DL_FUNC)&bnstruct_coupling_rows, 6},
    {"bnstruct_parent_margins_rows", (DL_FUNC)&bnstruct_parent_margins_rows, 6},
    {"bnstruct_compute_bf", (DL_FUNC)&bnstruct_compute_bf, 1},
    {"bnstruct_fumt_mask", (DL_FUNC)&bnstruct_fumt_mask, 2},
    {"bnstruct_g2_stat", (DL_FUNC)&bnstruct_g2_stat, 2},
    {"bnstruct_heom_dist", (DL_FUNC)&bnstruct_heom_dist, 4},
    {"bnstruct_in_tabu", (DL_FUNC)&bnstruct_in_tabu, 2},
    {"bnstruct_is_acyclic", (DL_FUNC)&bnstruct_is_acyclic, 1},
    {"bnstruct_next_comb", (DL_FUNC)&bnstruct_next_comb, 2},
    {"bnstruct_score_node", (DL_FUNC)&bnstruct_score_node, 6},
    {"bnstruct_dvbn_discretize_one", (DL_FUNC)&bnstruct_dvbn_discretize_one, 11},
    {"bnstruct_dvbn_discretize_all", (DL_FUNC)&bnstruct_dvbn_discretize_all, 15},
    {NULL, NULL, 0}};

void R_init_bnstruct(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
