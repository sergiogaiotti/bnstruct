#include "smfast.h"

static void update_best_pair(double cand_val, int cand_idx,
							 double *best_val, int *best_idx,
							 double *second_val, int *second_idx)
{
	if (cand_idx < 0 || !R_finite(cand_val))
		return;

	if (cand_idx == *best_idx)
	{
		if (cand_val > *best_val)
			*best_val = cand_val;
		return;
	}

	if (cand_val > *best_val)
	{
		*second_val = *best_val;
		*second_idx = *best_idx;
		*best_val = cand_val;
		*best_idx = cand_idx;
		return;
	}

	if (cand_val > *second_val)
	{
		*second_val = cand_val;
		*second_idx = cand_idx;
	}
}

SEXP bnstruct_fbp(SEXP aflml)
{
	R_len_t ni, si, bi, pos, pos2;
	R_len_t n_nodes = nrows(aflml);
	R_len_t pow_nodes = ncols(aflml);
	R_len_t mask = 1;
	int *bps; // force to 32bit integer
	double *ml = REAL(aflml), *bss;
	SEXP res;

	// allocate output
	PROTECT(res = allocMatrix(INTSXP, n_nodes, pow_nodes));
	bps = INTEGER(res);
	memset(bps, 0, sizeof(int) * n_nodes * pow_nodes);

	// allocate internal structures
	bss = (double *)R_alloc(n_nodes * pow_nodes, sizeof(double));
	memset(bss, 0, sizeof(double) * n_nodes * pow_nodes);

	for (ni = 0; ni < n_nodes; ni++)
		for (si = 0; si < pow_nodes; si++)
		{
			/* a node can not be parent of itself */
			if (si & (mask << ni))
				continue;

			pos = si * n_nodes + ni;
			bps[pos] = si + 1;
			bss[pos] = ml[pos];

			/* cycle through all parents and remove one in turn */
			for (bi = 0; bi < n_nodes; bi++)
				if (si & (mask << bi))
				{
					pos2 = (si ^ (mask << bi)) * n_nodes + ni; /* remove bi-th parent */
					if (bss[pos2] > bss[pos])
					{
						bps[pos] = bps[pos2];
						bss[pos] = bss[pos2];
					}
				}
		}

	UNPROTECT(1);
	return (res);
}

SEXP bnstruct_fbs(SEXP bps, SEXP aflml)
{
	R_len_t si, sink, nop;
	R_len_t n_nodes = nrows(aflml);
	R_len_t pow_nodes = ncols(aflml);
	R_len_t mask = 1;
	double skore;
	int *sinks; // force to 32bit integer
	int *ps = INTEGER(bps);
	double *ml = REAL(aflml), *scores;
	SEXP res;

	// allocate output
	PROTECT(res = allocVector(INTSXP, pow_nodes));
	sinks = INTEGER(res);
	for (si = 0; si < pow_nodes; si++)
		sinks[si] = -1;

	// allocate internal structures
	scores = (double *)R_alloc(pow_nodes, sizeof(double));
	memset(scores, 0, sizeof(double) * pow_nodes);

	for (si = 0; si < pow_nodes; si++)
		for (sink = 0; sink < n_nodes; sink++)
			if (si & (mask << sink))
			{
				nop = si ^ (mask << sink);
				skore = scores[nop] + ml[(ps[nop * n_nodes + sink] - 1) * n_nodes + sink];

				if (sinks[si] == -1 || skore > scores[si])
				{
					scores[si] = skore;
					sinks[si] = sink + 1;
				}
			}

	UNPROTECT(1);
	return (res);
}

SEXP bnstruct_fbs_with_scores(SEXP bps, SEXP aflml)
{
	R_len_t si, sink, nop;
	R_len_t n_nodes = nrows(aflml);
	R_len_t pow_nodes = ncols(aflml);
	R_len_t mask = 1;
	double skore;
	int *sinks; // force to 32bit integer
	int *ps = INTEGER(bps);
	double *ml = REAL(aflml), *scores;
	SEXP res, sinks_sexp, scores_sexp, names;

	PROTECT(res = allocVector(VECSXP, 2));
	PROTECT(sinks_sexp = allocVector(INTSXP, pow_nodes));
	PROTECT(scores_sexp = allocVector(REALSXP, pow_nodes));
	PROTECT(names = allocVector(STRSXP, 2));

	SET_STRING_ELT(names, 0, mkChar("sinks"));
	SET_STRING_ELT(names, 1, mkChar("scores"));
	setAttrib(res, R_NamesSymbol, names);

	sinks = INTEGER(sinks_sexp);
	scores = REAL(scores_sexp);
	for (si = 0; si < pow_nodes; si++)
	{
		sinks[si] = -1;
		scores[si] = 0.0;
	}

	for (si = 0; si < pow_nodes; si++)
		for (sink = 0; sink < n_nodes; sink++)
			if (si & (mask << sink))
			{
				nop = si ^ (mask << sink);
				skore = scores[nop] + ml[(ps[nop * n_nodes + sink] - 1) * n_nodes + sink];

				if (sinks[si] == -1 || skore > scores[si])
				{
					scores[si] = skore;
					sinks[si] = sink + 1;
				}
			}

	SET_VECTOR_ELT(res, 0, sinks_sexp);
	SET_VECTOR_ELT(res, 1, scores_sexp);

	UNPROTECT(4);
	return (res);
}

SEXP bnstruct_bdeu_loo_bounds(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP ess)
{
	unsigned int i, j, n_pa, pos;

	// get inputs
	unsigned int *d = (unsigned int *)INTEGER(data);
	unsigned int n_nodes = ncols(data);
	unsigned int n_ex = nrows(data);
	unsigned int *ns = (unsigned int *)INTEGER(node_sizes);
	unsigned int *ifm = (unsigned int *)INTEGER(imp_fam_mask);
	unsigned int pow_nodes = ncols(imp_fam_mask);
	double alpha = *(REAL(ess));
	unsigned int *pa = (unsigned int *)R_alloc(n_nodes, sizeof(unsigned int));

	// allocate outputs
	SEXP res, lower_sexp, upper_sexp, names;
	PROTECT(res = allocVector(VECSXP, 2));
	PROTECT(lower_sexp = allocMatrix(REALSXP, n_nodes, pow_nodes));
	PROTECT(upper_sexp = allocMatrix(REALSXP, n_nodes, pow_nodes));
	PROTECT(names = allocVector(STRSXP, 2));

	SET_STRING_ELT(names, 0, mkChar("lower"));
	SET_STRING_ELT(names, 1, mkChar("upper"));
	setAttrib(res, R_NamesSymbol, names);

	double *lower = REAL(lower_sexp);
	double *upper = REAL(upper_sexp);
	for (i = 0; i < n_nodes * pow_nodes; i++)
	{
		lower[i] = 0.0;
		upper[i] = 0.0;
	}

	for (i = 0; i < n_nodes; i++)
		for (j = 0; j < pow_nodes; j++)
		{
			pos = j * n_nodes + i;
			if (!ifm[pos])
				continue;

			n_pa = get_bits(j, pa, n_nodes);

			int cum_prod_sizes[n_pa + 2];
			int strides[n_pa + 1];
			int prod_sizes_pa = 1;
			int prod_sizes;
			int k;
			int pc;
			unsigned int ex;

			cum_prod_sizes[0] = 1;
			cum_prod_sizes[1] = ns[i];
			for (k = 0; k < (int)n_pa; k++)
			{
				cum_prod_sizes[k + 2] = cum_prod_sizes[k + 1] * ns[pa[k]];
				prod_sizes_pa *= ns[pa[k]];
			}

			prod_sizes = prod_sizes_pa * ns[i];
			strides[0] = i * n_ex;
			for (k = 0; k < (int)n_pa; k++)
				strides[k + 1] = pa[k] * n_ex;

			double *counts = (double *)calloc(prod_sizes, sizeof(double));
			if (counts == NULL)
				error("memory allocation failed in bnstruct_bdeu_loo_bounds");

			unsigned int r_i = ns[i];
			double alpha_ij = alpha / (double)prod_sizes_pa;
			double alpha_ijk = alpha / (double)prod_sizes;
			double min_delta = R_PosInf;
			double max_delta = R_NegInf;
			int seen = 0;

			for (ex = 0; ex < n_ex; ex++)
			{
				int index = 0;
				int elmt;
				unsigned int jj;
				for (jj = 0; jj < n_pa + 1; jj++)
				{
					elmt = d[ex + strides[jj]];
					if (elmt == NA_INTEGER)
						break;
					index += (elmt - 1) * cum_prod_sizes[jj];
				}
				if (jj < n_pa + 1)
					continue;

				counts[index] += 1;
			}

			for (pc = 0; pc < prod_sizes_pa; pc++)
			{
				int base = pc * r_i;
				double n_ij = 0.0;
				for (k = 0; k < (int)r_i; k++)
					n_ij += counts[base + k];
				if (n_ij < 1.0)
					continue;
				for (k = 0; k < (int)r_i; k++)
				{
					double n_ijk = counts[base + k];
					if (n_ijk < 1.0)
						continue;
					double delta = log(alpha_ij + n_ij - 1.0) - log(alpha_ijk + n_ijk - 1.0);
					if (delta < min_delta)
						min_delta = delta;
					if (delta > max_delta)
						max_delta = delta;
					seen = 1;
				}
			}

			if (!seen)
			{
				min_delta = 0.0;
				max_delta = 0.0;
			}

			lower[pos] = min_delta;
			upper[pos] = max_delta;

			free(counts);
		}

	SET_VECTOR_ELT(res, 0, lower_sexp);
	SET_VECTOR_ELT(res, 1, upper_sexp);

	UNPROTECT(4);
	return res;
}

SEXP bnstruct_compute_bg_cons(SEXP lower, SEXP upper, SEXP imp_fam_mask)
{
	R_len_t ni, si, bi, pos, pos2;
	R_len_t n_nodes = nrows(lower);
	R_len_t pow_nodes = ncols(lower);
	R_len_t mask = 1;
	int *ifm = INTEGER(imp_fam_mask);
	double *low = REAL(lower);
	double *up = REAL(upper);
	SEXP res;

	PROTECT(res = allocMatrix(REALSXP, n_nodes, pow_nodes));
	double *bg = REAL(res);

	for (pos = 0; pos < n_nodes * pow_nodes; pos++)
	{
		if (ifm[pos])
			bg[pos] = fmax(fabs(low[pos]), fabs(up[pos]));
		else
			bg[pos] = R_NegInf;
	}

	for (ni = 0; ni < n_nodes; ni++)
		for (si = 0; si < pow_nodes; si++)
		{
			pos = si * n_nodes + ni;
			for (bi = 0; bi < n_nodes; bi++)
				if (si & (mask << bi))
				{
					pos2 = (si ^ (mask << bi)) * n_nodes + ni;
					if (bg[pos2] > bg[pos])
						bg[pos] = bg[pos2];
				}
		}

	UNPROTECT(1);
	return res;
}

SEXP bnstruct_compute_bg_tight(SEXP aflml, SEXP bps, SEXP lower, SEXP upper, SEXP imp_fam_mask)
{
	R_len_t ni, si, bi, pos, pos2;
	R_len_t n_nodes = nrows(aflml);
	R_len_t pow_nodes = ncols(aflml);
	R_len_t mask = 1;
	int *ifm = INTEGER(imp_fam_mask);
	int *ps = INTEGER(bps);
	double *ml = REAL(aflml);
	double *low = REAL(lower);
	double *up = REAL(upper);
	SEXP res;

	PROTECT(res = allocMatrix(REALSXP, n_nodes, pow_nodes));
	double *bg = REAL(res);

	for (ni = 0; ni < n_nodes; ni++)
	{
		double *best_val = (double *)R_alloc(pow_nodes, sizeof(double));
		double *second_val = (double *)R_alloc(pow_nodes, sizeof(double));
		int *best_idx = (int *)R_alloc(pow_nodes, sizeof(int));
		int *second_idx = (int *)R_alloc(pow_nodes, sizeof(int));

		for (si = 0; si < pow_nodes; si++)
		{
			pos = si * n_nodes + ni;
			if (ifm[pos])
			{
				best_val[si] = ml[pos] + up[pos];
				best_idx[si] = (int)si + 1;
			}
			else
			{
				best_val[si] = R_NegInf;
				best_idx[si] = -1;
			}
			second_val[si] = R_NegInf;
			second_idx[si] = -1;
		}

		for (si = 0; si < pow_nodes; si++)
			for (bi = 0; bi < n_nodes; bi++)
				if (si & (mask << bi))
				{
					pos2 = (si ^ (mask << bi));
					update_best_pair(best_val[pos2], best_idx[pos2],
									 &best_val[si], &best_idx[si],
									 &second_val[si], &second_idx[si]);
					update_best_pair(second_val[pos2], second_idx[pos2],
									 &best_val[si], &best_idx[si],
									 &second_val[si], &second_idx[si]);
				}

		for (si = 0; si < pow_nodes; si++)
		{
			pos = si * n_nodes + ni;
			if (!ifm[pos])
			{
				bg[pos] = R_NegInf;
				continue;
			}

			int p_star = ps[pos];
			int p_star_pos = (p_star - 1) * n_nodes + ni;
			double s_star = ml[p_star_pos];
			double lower_bound = low[p_star_pos];
			double upper_bound = up[p_star_pos];
			double upper_g = upper_bound;

			if (best_idx[si] != -1)
			{
				double cand_best = R_NegInf;
				if (best_idx[si] == p_star)
				{
					if (second_idx[si] != -1)
						cand_best = second_val[si];
				}
				else
				{
					cand_best = best_val[si];
				}

				if (R_finite(cand_best))
					upper_g = fmax(upper_g, cand_best - s_star);
			}

			bg[pos] = fmax(fabs(lower_bound), fabs(upper_g));
		}
	}

	UNPROTECT(1);
	return res;
}

SEXP bnstruct_compute_bf(SEXP Bg)
{
	R_len_t si, sink, nop;
	R_len_t n_nodes = nrows(Bg);
	R_len_t pow_nodes = ncols(Bg);
	R_len_t mask = 1;
	double skore;
	SEXP res;

	PROTECT(res = allocVector(REALSXP, pow_nodes));
	double *bf = REAL(res);
	double *bg = REAL(Bg);

	for (si = 0; si < pow_nodes; si++)
		bf[si] = R_NegInf;
	bf[0] = 0.0;

	for (si = 1; si < pow_nodes; si++)
		for (sink = 0; sink < n_nodes; sink++)
			if (si & (mask << sink))
			{
				nop = si ^ (mask << sink);
				skore = bf[nop] + bg[nop * n_nodes + sink];
				if (skore > bf[si])
					bf[si] = skore;
			}

	UNPROTECT(1);
	return res;
}

SEXP bnstruct_na_rows_int(SEXP mat)
{
	int i, j, stride, *res;
	int n_rows = nrows(mat);
	int n_cols = ncols(mat);
	int *m = INTEGER(mat);
	SEXP(result);

	// allocate oputput
	PROTECT(result = allocVector(INTSXP, n_rows));
	res = INTEGER(result);
	memset(res, 0, sizeof(int) * n_rows);

	for (j = 0; j < n_cols; j++)
	{
		stride = n_rows * j;
		for (i = 0; i < n_rows; i++)
			res[i] |= (int)(NA_INTEGER == (m[i + stride]));
	}

	UNPROTECT(1);
	return (result);
}

SEXP bnstruct_fumt_mask(SEXP n_elements, SEXP pattern)
{
	int i, times, ind, stride, start;
	int bitmask = 1;
	int n_el = asInteger(n_elements);
	PROTECT(pattern = coerceVector(pattern, INTSXP));
	int *pat = INTEGER(pattern);
	int l_pat = LENGTH(pattern);
	int pow_el = bitmask << n_el; // faster than power for int
	SEXP result;

	// allocate output
	PROTECT(result = allocVector(INTSXP, pow_el));
	int *mask = INTEGER(result);
	memset(mask, 0, sizeof(int) * pow_el);

	// assign pattern
	for (i = 0; i < l_pat; i++)
		mask[bitmask << ((int)pat[i] - 1)] = 1;

	// compute fumt
	for (i = 0; i < n_el; i++)
	{
		stride = bitmask << i;
		for (times = 0; times < (bitmask << (n_el - i - 1)); times++)
		{
			start = times * (bitmask << (i + 1));
			for (ind = 0; ind < stride; ind++)
				mask[start + ind + stride] += mask[start + ind];
		}
	}

	UNPROTECT(2);
	return (result);
}

SEXP bnstruct_all_fam_log_marg_lik(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP iss, SEXP func)
{
	unsigned int i, j, n_pa, pos;

	// get inputs
	unsigned int *d = (unsigned int *)INTEGER(data);
	unsigned int n_nodes = ncols(data);
	unsigned int n_ex = nrows(data);
	unsigned int *ns = (unsigned int *)INTEGER(node_sizes);
	unsigned int *ifm = (unsigned int *)INTEGER(imp_fam_mask);
	unsigned int pow_nodes = ncols(imp_fam_mask);
	double alpha = *(REAL(iss));
	unsigned int scoring_func = *INTEGER(func);
	unsigned int *pa = (unsigned int *)R_alloc(n_nodes, sizeof(unsigned int));

	// allocate and initialize output
	SEXP result;
	PROTECT(result = allocMatrix(REALSXP, n_nodes, pow_nodes));
	double *aflml = REAL(result);
	for (i = 0; i < n_nodes * pow_nodes; i++)
		aflml[i] = R_NegInf;

	// compute log likelihood
	for (i = 0; i < n_nodes; i++)
		for (j = 0; j < pow_nodes; j++)
		{
			pos = j * n_nodes + i;
			if (ifm[pos])
			{
				// Rprintf("get bits\n");
				n_pa = get_bits(j, pa, n_nodes);
				// Rprintf("log lik, node %d, n parents %d\n",i,n_pa);
				aflml[pos] = score_node_1((unsigned int *)d, n_nodes, n_ex, (unsigned int *)ns, i, (unsigned int *)pa, n_pa, scoring_func, alpha);
				// bdeu_score( d, n_nodes, n_ex, ns, i, pa, n_pa, alpha );
				// Rprintf("end\n");
			}
		}

	UNPROTECT(1);
	return result;
}

unsigned int get_bits(unsigned int word, unsigned int *bits, unsigned int size)
{
	unsigned int i, count = 0, bitmask = 1;
	for (i = 0; i < size; i++)
		if (word & bitmask << i)
			bits[count++] = i;

	return count;
}

/*
SEXP llna( SEXP parents, int node, SEXP node_sizes, double iss, SEXP data )
{
	int i, j, index, flag, elmt;
	int * pa = INTEGER(parents);
	int * ns = INTEGER(node_sizes);
	int * d = INTEGER(data);
	int n_cases = nrows(data);
	int n_pa = length(pa);

	// convenience variables
	long cum_size[ n_pa+1 ];
	long sum_cum_size_pa = 1;
	cum_size[0] = 1;
	for( i = 1; i < n_pa; i++ )
	{
		cum_size[i] = cum_size[i-1] * node_sizes[pa[i-1]-1];
		sum_cum_size_pa += cum_size[i];
	}
	cum_size[n_pa] = cum_size[n_pa-1] * node_sizes[pa[n_pa]-1];
	sum_cum_size = sum_cum_size_pa + cum_size[n_pa];
	long prod_size_pa = cum_size[n_pa];
	long prod_size = prod_size_pa * node_sizes[node-1];

	// allocate space for storing counts
	counts = (int *) R_alloc( prod_size, sizeof(int) );
	memset( counts, 0, sizeof(int) * prod_size );
	counts_pa = (int *) R_alloc( prod_size_pa, sizeof(int) );
	memset( counts_pa, 0, sizeof(int) * prod_size_pa );

	// compute counts
	for( i = 0; i < n_cases; i++ )
	{
		index = 0;
		flag = 0;
		for( j = 0; j < n_pa; j++ )
		{
			elmt = d[ i + n_cases*(pa[j]-1) ];
			if( NA_INTEGER(elmt) ) // if one elmt is NA, skip case
			{
				flag = 1;
				break;
			}
			index += elmt * cum_size[i];
		}
		if( flag ) continue;

		elmt = d[ i + n_cases*(node-1) ];
		if( !NA_INTEGER(elmt) )
		{
			counts_pa[ index - sum_cum_size_pa ]++;
			counts[ index*elmt - sum_cum_size ]++;
		}
	}

}
*/
