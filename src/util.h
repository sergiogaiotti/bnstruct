#ifndef _UTIL_H_
#define _UTIL_H_

#include <R.h>
#include <Rinternals.h>
#include <math.h>

/* compute the BDEu score of node _ni_, with parents _pars_, from dataset
   _data_ and with equivalent sample size _ess_
*/
SEXP bnstruct_score_node( SEXP data, SEXP node_sizes, SEXP ni, SEXP pars, SEXP func, SEXP ess );

// compute score for a node as above, but called from c, not from R
double score_node_1( unsigned int* data, int ncols_data, int nrows_data, unsigned int* node_sizes, unsigned int ni, unsigned int* pars, int length_pars, int func, double ess );

/* test wether the graph _g_ is acyclic
*/
SEXP bnstruct_is_acyclic( SEXP g );

/* computes contingency table, excluding NAs: 
*/
SEXP bnstruct_compute_counts_nas( SEXP data, SEXP node_sizes );

/* computes contingency table, no NAs assumed */
SEXP bnstruct_compute_counts( SEXP data, SEXP node_sizes );

/* computes log likelihood for the family (ni,pa), accounting for NAs
*/
double bdeu_score( unsigned int * d, unsigned int n_nodes, unsigned int n_cases, unsigned int * ns, 
	unsigned int ni, unsigned int * pa, unsigned int n_pa, double alpha );
  
/**
 * penalized log likelihood; func = 1 -> AIC penalty (dim), func = 2 -> BIC
 * penalty (dim/2 * log N, with N the family's complete-case row count)
 */
double log_likelihood( unsigned int * d, unsigned int n_nodes, unsigned int n_cases, unsigned int * ns,
  unsigned int ni, unsigned int * pa, unsigned int n_pa, int func );

/* the two family scores above from a contingency table counts[pc * r_i + k] */
double bdeu_score_counts( const double * counts, int prod_sizes, int prod_sizes_pa, int r_i, double alpha );
double log_likelihood_counts( const double * counts, int prod_sizes_pa, unsigned int r_i, int func );
	
#endif
