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

/* Lexicographic order of the rows of an n_ex x n_nodes integer matrix stored by
 * column (NA_INTEGER compares as an ordinary value), ties broken by row index, so
 * equal rows form consecutive runs with the smallest index first. qsort has no
 * context argument, hence the file-scope pointers; they are set and used only
 * inside loo_sort_rows. */
static const unsigned int *loo_sort_d = NULL;
static unsigned int loo_sort_n = 0, loo_sort_p = 0;

static int loo_row_cmp(const void *a, const void *b)
{
	unsigned int ra = *(const unsigned int *)a, rb = *(const unsigned int *)b;
	unsigned int j;
	for (j = 0; j < loo_sort_p; j++)
	{
		int va = (int)loo_sort_d[ra + (R_xlen_t)j * loo_sort_n];
		int vb = (int)loo_sort_d[rb + (R_xlen_t)j * loo_sort_n];
		if (va != vb)
			return (va < vb) ? -1 : 1;
	}
	return (ra < rb) ? -1 : (ra > rb);
}

static void loo_sort_rows(const unsigned int *d, unsigned int n_ex, unsigned int n_nodes, unsigned int *ord)
{
	loo_sort_d = d;
	loo_sort_n = n_ex;
	loo_sort_p = n_nodes;
	qsort(ord, n_ex, sizeof(unsigned int), loo_row_cmp);
	loo_sort_d = NULL;
}

static int loo_rows_equal(const unsigned int *d, unsigned int n_ex, unsigned int n_nodes, unsigned int ra, unsigned int rb)
{
	unsigned int j;
	for (j = 0; j < n_nodes; j++)
		if (d[ra + (R_xlen_t)j * n_ex] != d[rb + (R_xlen_t)j * n_ex])
			return 0;
	return 1;
}

/* ------------------------------------------------------------------------------------
 * Family counting over distinct row patterns, shared by every routine that counts all
 * families: the learner's family scores, the removal-delta bounds, exact leave-one-out
 * and the per-row coupling tables. Outputs are identical to counting row by row.
 *
 * (a) Rows enter a family's counts only through their values. The N rows are sorted
 *     once and replaced by their K distinct patterns with multiplicities (missing
 *     values included); counts are sums of integers, hence identical.
 * (b) The parent-configuration index of a pattern in parent set P,
 *         pc(P) = sum_k (x_{pa[k]} - 1) prod_{l<k} r_{pa[l]}   (parents increasing),
 *     satisfies pc(P u {e}) = pc(P) + (x_e - 1) prod_{l in P} r_l when e > max(P).
 *     A depth-first walk over parent sets that adds elements in increasing order
 *     computes it in O(1) per (pattern, set), keeping one index array per depth
 *     (memory (p + 1) K). Counting one family then costs O(K) instead of O(N |P|).
 *     The cell of a pattern in family (i, P) is pc(P) * r_i + (x_i - 1), the layout
 *     of bdeu_score; pc = -1 marks a pattern with a missing parent value.
 * ------------------------------------------------------------------------------------ */
typedef struct
{
	unsigned int n;		 /* number of distinct patterns K */
	unsigned int *ord;	 /* row indices in lexicographic order (length N) */
	unsigned int *start; /* rows of pattern u: ord[start[u]] .. ord[start[u + 1] - 1] */
	double *w;			 /* multiplicity of pattern u */
	int *x;				 /* K x p values of the patterns, column-major, NA_INTEGER kept */
} row_patterns;

static void build_row_patterns(const unsigned int *d, unsigned int n_ex, unsigned int n_nodes, row_patterns *rp)
{
	unsigned int g, h, u, j, K = 0;
	rp->ord = (unsigned int *)R_alloc(n_ex > 0 ? n_ex : 1, sizeof(unsigned int));
	for (g = 0; g < n_ex; g++)
		rp->ord[g] = g;
	if (n_ex > 1)
		loo_sort_rows(d, n_ex, n_nodes, rp->ord);
	for (g = 0; g < n_ex; g = h, K++)
		for (h = g + 1; h < n_ex && loo_rows_equal(d, n_ex, n_nodes, rp->ord[g], rp->ord[h]); h++)
			;
	rp->n = K;
	rp->start = (unsigned int *)R_alloc(K + 1, sizeof(unsigned int));
	rp->w = (double *)R_alloc(K > 0 ? K : 1, sizeof(double));
	rp->x = (int *)R_alloc((R_xlen_t)(K > 0 ? K : 1) * n_nodes, sizeof(int));
	for (g = 0, u = 0; g < n_ex; g = h, u++)
	{
		for (h = g + 1; h < n_ex && loo_rows_equal(d, n_ex, n_nodes, rp->ord[g], rp->ord[h]); h++)
			;
		rp->start[u] = g;
		rp->w[u] = (double)(h - g);
		for (j = 0; j < n_nodes; j++)
			rp->x[u + (R_xlen_t)j * K] = (int)d[rp->ord[g] + (R_xlen_t)j * n_ex];
	}
	rp->start[K] = n_ex;
}

typedef struct
{
	unsigned int p, K, maxsize, depth;
	int started;
	const int *x, *ns;
	unsigned int *set;	/* parent set at each depth */
	unsigned int *next; /* next element to try below each depth */
	R_xlen_t *q;		/* number of parent configurations at each depth */
	R_xlen_t *pc;		/* (p + 1) x K parent-configuration indices */
} subset_walk;

static void walk_init(subset_walk *w, const row_patterns *rp, unsigned int p, const int *ns, unsigned int maxsize)
{
	R_xlen_t u, K = rp->n;
	w->p = p;
	w->K = rp->n;
	w->maxsize = maxsize;
	w->depth = 0;
	w->started = 0;
	w->x = rp->x;
	w->ns = ns;
	w->set = (unsigned int *)R_alloc(p + 1, sizeof(unsigned int));
	w->next = (unsigned int *)R_alloc(p + 1, sizeof(unsigned int));
	w->q = (R_xlen_t *)R_alloc(p + 1, sizeof(R_xlen_t));
	w->pc = (R_xlen_t *)R_alloc((R_xlen_t)(p + 1) * (K > 0 ? K : 1), sizeof(R_xlen_t));
	w->set[0] = 0;
	w->next[0] = 0;
	w->q[0] = 1;
	for (u = 0; u < K; u++)
		w->pc[u] = 0;
}

/* Move to the next parent set in pre-order (the empty set first); 0 when every set
 * of at most maxsize elements has been visited. */
static int walk_next(subset_walk *w)
{
	if (!w->started)
	{
		w->started = 1;
		return 1;
	}
	for (;;)
	{
		unsigned int dd = w->depth;
		if (dd < w->maxsize && w->next[dd] < w->p)
		{
			unsigned int e = w->next[dd]++;
			R_xlen_t u, K = w->K, qd = w->q[dd];
			const R_xlen_t *src = w->pc + (R_xlen_t)dd * K;
			R_xlen_t *dst = w->pc + (R_xlen_t)(dd + 1) * K;
			const int *xe = w->x + (R_xlen_t)e * K;
			for (u = 0; u < K; u++)
				dst[u] = (src[u] < 0 || xe[u] == NA_INTEGER) ? -1 : src[u] + (R_xlen_t)(xe[u] - 1) * qd;
			w->set[dd + 1] = w->set[dd] | (1U << e);
			w->q[dd + 1] = qd * w->ns[e];
			w->next[dd + 1] = e + 1;
			w->depth = dd + 1;
			return 1;
		}
		if (dd == 0)
			return 0;
		w->depth = dd - 1;
	}
}

#define WALK_SET(w) ((w)->set[(w)->depth])
#define WALK_Q(w) ((w)->q[(w)->depth])
#define WALK_PC(w) ((w)->pc + (R_xlen_t)(w)->depth * (w)->K)

/* Largest admissible parent set, and the largest admissible count table (cells). */
static void admissible_extent(const int *ifm, const int *ns, unsigned int n_nodes, unsigned int pow_nodes,
							  unsigned int *maxsize, R_xlen_t *maxcells, R_xlen_t *maxq)
{
	unsigned int s, i, k;
	*maxsize = 0;
	*maxcells = 1;
	*maxq = 1;
	for (s = 0; s < pow_nodes; s++)
	{
		unsigned int size = 0;
		R_xlen_t q = 1;
		for (k = 0; k < n_nodes; k++)
			if (s & (1U << k))
			{
				size++;
				q *= ns[k];
			}
		for (i = 0; i < n_nodes; i++)
			if (ifm[(R_xlen_t)s * n_nodes + i])
			{
				if (size > *maxsize)
					*maxsize = size;
				if (q * ns[i] > *maxcells)
					*maxcells = q * ns[i];
				if (q > *maxq)
					*maxq = q;
			}
	}
}

/* f(x) = x log x - (x-1) log(x-1), the per-cell term of the BIC/AIC removal
 * delta (0 log 0 = 0, so f(1) = 0). The log-likelihood part of removing one
 * row from cell (j,k) is f(n_ij) - f(n_ijk). */
static double loo_f_cell(double x)
{
	double v = x * log(x);
	if (x > 1.0)
		v -= (x - 1.0) * log(x - 1.0);
	return v;
}

/* One pass over all admissible families on the distinct row patterns, computing any of:
 *   - the family scores (score_func 0 = BDeu, 1 = AIC, 2 = BIC), bit-identical to
 *     bdeu_score / log_likelihood (same count tables, same scoring code);
 *   - the per-family min/max of the one-row removal delta and its per-child-value
 *     min/max (lower, upper, lower_k, upper_k).
 * The BIC penalty uses the family's complete-case row count, matching log_likelihood();
 * its N -> N-1 relief term is row-independent within the family, and is zero for AIC
 * (no N in the penalty) and for BDeu. Removing a row that is incomplete in the family
 * shifts its score by 0, so 0 enters the bounds whenever such rows exist.
 * Returns list(aflml, lower, upper, lower_k, upper_k) with NULL for components not
 * requested. Cost O(N p log N) to find the patterns, O(K) per family for counting and
 * bounds, plus O(cells) per family for the scores (the score sums over all cells). */
static SEXP family_scores_bounds(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, double alpha,
								 int scoring_func, int want_scores, int want_bounds)
{
	const unsigned int *d = (const unsigned int *)INTEGER(data);
	unsigned int n_nodes = ncols(data), n_ex = nrows(data);
	const int *ns = INTEGER(node_sizes), *ifm = INTEGER(imp_fam_mask);
	unsigned int pow_nodes = ncols(imp_fam_mask);
	unsigned int i, k, max_r = 1, maxsize;
	R_xlen_t t, u, maxcells, maxq, npk = (R_xlen_t)n_nodes * pow_nodes;
	for (i = 0; i < n_nodes; i++)
		if ((unsigned int)ns[i] > max_r)
			max_r = ns[i];

	SEXP res, names, aflml_s = R_NilValue, lower_s = R_NilValue, upper_s = R_NilValue;
	SEXP lower_k_s = R_NilValue, upper_k_s = R_NilValue, dim_k;
	int nprot = 2;
	PROTECT(res = allocVector(VECSXP, 5));
	PROTECT(names = allocVector(STRSXP, 5));
	SET_STRING_ELT(names, 0, mkChar("aflml"));
	SET_STRING_ELT(names, 1, mkChar("lower"));
	SET_STRING_ELT(names, 2, mkChar("upper"));
	SET_STRING_ELT(names, 3, mkChar("lower_k"));
	SET_STRING_ELT(names, 4, mkChar("upper_k"));
	setAttrib(res, R_NamesSymbol, names);
	double *aflml = NULL, *lower = NULL, *upper = NULL, *lower_k = NULL, *upper_k = NULL;
	if (want_scores)
	{
		PROTECT(aflml_s = allocMatrix(REALSXP, n_nodes, pow_nodes));
		nprot++;
		aflml = REAL(aflml_s);
		for (t = 0; t < npk; t++)
			aflml[t] = R_NegInf;
		SET_VECTOR_ELT(res, 0, aflml_s);
	}
	if (want_bounds)
	{
		PROTECT(lower_s = allocMatrix(REALSXP, n_nodes, pow_nodes));
		PROTECT(upper_s = allocMatrix(REALSXP, n_nodes, pow_nodes));
		PROTECT(lower_k_s = allocVector(REALSXP, npk * max_r));
		PROTECT(upper_k_s = allocVector(REALSXP, npk * max_r));
		PROTECT(dim_k = allocVector(INTSXP, 3));
		nprot += 5;
		INTEGER(dim_k)[0] = n_nodes;
		INTEGER(dim_k)[1] = pow_nodes;
		INTEGER(dim_k)[2] = max_r;
		setAttrib(lower_k_s, R_DimSymbol, dim_k);
		setAttrib(upper_k_s, R_DimSymbol, duplicate(dim_k));
		lower = REAL(lower_s);
		upper = REAL(upper_s);
		lower_k = REAL(lower_k_s);
		upper_k = REAL(upper_k_s);
		for (t = 0; t < npk; t++)
		{
			lower[t] = 0.0;
			upper[t] = 0.0;
		}
		for (t = 0; t < npk * max_r; t++)
		{
			lower_k[t] = R_PosInf;
			upper_k[t] = R_NegInf;
		}
		SET_VECTOR_ELT(res, 1, lower_s);
		SET_VECTOR_ELT(res, 2, upper_s);
		SET_VECTOR_ELT(res, 3, lower_k_s);
		SET_VECTOR_ELT(res, 4, upper_k_s);
	}

	admissible_extent(ifm, ns, n_nodes, pow_nodes, &maxsize, &maxcells, &maxq);
	row_patterns rp;
	build_row_patterns(d, n_ex, n_nodes, &rp);
	R_xlen_t K = rp.n;
	/* count tables stay all-zero between families: only touched cells are reset */
	double *counts = (double *)R_alloc(maxcells, sizeof(double));
	double *nij = (double *)R_alloc(maxq, sizeof(double));
	for (t = 0; t < maxcells; t++)
		counts[t] = 0.0;
	for (t = 0; t < maxq; t++)
		nij[t] = 0.0;
	R_xlen_t *touched = (R_xlen_t *)R_alloc(K > 0 ? K : 1, sizeof(R_xlen_t)), nt;
	double *low_k_local = (double *)R_alloc(max_r, sizeof(double));
	double *up_k_local = (double *)R_alloc(max_r, sizeof(double));
	double *n_k_incomplete = (double *)R_alloc(max_r, sizeof(double));

	subset_walk w;
	walk_init(&w, &rp, n_nodes, ns, maxsize);
	while (walk_next(&w))
	{
		unsigned int P = WALK_SET(&w);
		R_xlen_t q = WALK_Q(&w);
		const R_xlen_t *pc = WALK_PC(&w);
		for (i = 0; i < n_nodes; i++)
		{
			R_xlen_t pos = (R_xlen_t)P * n_nodes + i;
			if (!ifm[pos])
				continue;
			int r_i = ns[i];
			const int *xi = rp.x + (R_xlen_t)i * K;
			double n_child_na = 0.0, n_fam = 0.0;
			for (k = 0; k < (unsigned int)r_i; k++)
				n_k_incomplete[k] = 0.0;
			nt = 0;
			for (u = 0; u < K; u++)
			{
				if (xi[u] == NA_INTEGER)
				{
					n_child_na += rp.w[u];
					continue;
				}
				if (pc[u] < 0)
				{
					n_k_incomplete[xi[u] - 1] += rp.w[u];
					continue;
				}
				R_xlen_t cell = pc[u] * r_i + (xi[u] - 1);
				if (counts[cell] == 0.0)
					touched[nt++] = cell;
				counts[cell] += rp.w[u];
				nij[pc[u]] += rp.w[u];
				n_fam += rp.w[u];
			}

			if (want_scores)
				aflml[pos] = (scoring_func == 0)
								 ? bdeu_score_counts(counts, (int)(q * r_i), (int)q, r_i, alpha)
								 : log_likelihood_counts(counts, (int)q, (unsigned int)r_i, scoring_func);

			if (want_bounds)
			{
				double alpha_ij = alpha / (double)q;
				double alpha_ijk = alpha / (double)(q * r_i);
				double min_delta = R_PosInf, max_delta = R_NegInf;
				int seen = 0;
				for (k = 0; k < (unsigned int)r_i; k++)
				{
					low_k_local[k] = R_PosInf;
					up_k_local[k] = R_NegInf;
				}
				double pen_relief = 0.0;
				if (scoring_func == 2 && n_fam >= 2.0)
					pen_relief = 0.5 * (double)q * (double)(r_i - 1) * log(n_fam / (n_fam - 1.0));
				/* one delta per occupied cell (min/max do not depend on the order) */
				for (t = 0; t < nt; t++)
				{
					int kk = (int)(touched[t] % r_i);
					double n_ij = nij[touched[t] / r_i], n_ijk = counts[touched[t]];
					double delta = (scoring_func == 0)
									   ? log(alpha_ij + n_ij - 1.0) - log(alpha_ijk + n_ijk - 1.0)
									   : loo_f_cell(n_ij) - loo_f_cell(n_ijk) + pen_relief;
					if (delta < min_delta)
						min_delta = delta;
					if (delta > max_delta)
						max_delta = delta;
					if (delta < low_k_local[kk])
						low_k_local[kk] = delta;
					if (delta > up_k_local[kk])
						up_k_local[kk] = delta;
					seen = 1;
				}
				if (!seen)
				{
					min_delta = 0.0;
					max_delta = 0.0;
				}
				if (n_child_na > 0.0)
				{
					if (min_delta > 0.0)
						min_delta = 0.0;
					if (max_delta < 0.0)
						max_delta = 0.0;
				}
				for (k = 0; k < (unsigned int)r_i; k++)
					if (n_k_incomplete[k] > 0.0)
					{
						if (low_k_local[k] > 0.0)
							low_k_local[k] = 0.0;
						if (up_k_local[k] < 0.0)
							up_k_local[k] = 0.0;
						if (min_delta > 0.0)
							min_delta = 0.0;
						if (max_delta < 0.0)
							max_delta = 0.0;
					}
				lower[pos] = min_delta;
				upper[pos] = max_delta;
				for (k = 0; k < (unsigned int)r_i; k++)
				{
					lower_k[pos + npk * k] = low_k_local[k];
					upper_k[pos + npk * k] = up_k_local[k];
				}
			}

			for (t = 0; t < nt; t++)
			{
				counts[touched[t]] = 0.0;
				nij[touched[t] / r_i] = 0.0;
			}
		}
	}
	UNPROTECT(nprot);
	return res;
}

/* Per-family min/max (and per-child-value min/max) of the one-row removal delta, for
 * the scoring function selected by func: 0 = BDeu, 1 = AIC, 2 = BIC. */
SEXP bnstruct_score_loo_bounds(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP ess, SEXP func)
{
	SEXP all, res, names;
	PROTECT(all = family_scores_bounds(data, node_sizes, imp_fam_mask, *REAL(ess), *INTEGER(func), 0, 1));
	PROTECT(res = allocVector(VECSXP, 4));
	PROTECT(names = allocVector(STRSXP, 4));
	int k;
	for (k = 0; k < 4; k++)
	{
		SET_VECTOR_ELT(res, k, VECTOR_ELT(all, k + 1));
		SET_STRING_ELT(names, k, STRING_ELT(getAttrib(all, R_NamesSymbol), k + 1));
	}
	setAttrib(res, R_NamesSymbol, names);
	UNPROTECT(3);
	return res;
}

/* Family scores and removal-delta bounds from one counting pass:
 * list(aflml, lower, upper, lower_k, upper_k). */
SEXP bnstruct_all_fam_scores_bounds(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP ess, SEXP func)
{
	return family_scores_bounds(data, node_sizes, imp_fam_mask, *REAL(ess), *INTEGER(func), 1, 1);
}

/* backward-compatible BDeu entry point */
SEXP bnstruct_bdeu_loo_bounds(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP ess)
{
	SEXP func, res;
	PROTECT(func = ScalarInteger(0));
	res = bnstruct_score_loo_bounds(data, node_sizes, imp_fam_mask, ess, func);
	UNPROTECT(1);
	return res;
}

/*
 * Incremental exact leave-one-out: recover, for every removed observation, the
 * optimal parent set of every node WITHOUT recomputing the family scores from
 * scratch. The full-data family scores (aflml) are supplied; a removal shifts
 * each family score by the closed-form delta of the single contingency cell the
 * removed row occupies (precomputed once per cell), after which only the SM
 * dynamic program is re-run.
 *
 * Phase 1 counts every family once on the K distinct row patterns (O(K) per family,
 * see subset_walk) and stores the removal delta of every occupied cell (memory: all
 * admissible count tables, independent of N). Phase 2 runs once per distinct pattern:
 * the score shift costs O(2^p + families) (the parent-configuration index of a
 * candidate set extends that of the set without its highest element) and the DP
 * O(p^2 2^p); rows with the same pattern share the result. Total
 * O(N p log N + K p 2^p + cells + K p^2 2^p).
 *
 * Returns an integer matrix [n_ex, n_nodes]; entry [o, i] is the 0-indexed
 * parent-set bitmask of node i in the optimal DAG for D \ {row o}.
 */
SEXP bnstruct_loo_parents_score(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP ess, SEXP aflml_sexp, SEXP func)
{
	unsigned int i, k;
	R_len_t s;
	const unsigned int *d = (const unsigned int *)INTEGER(data);
	unsigned int n_nodes = ncols(data);
	unsigned int n_ex = nrows(data);
	const int *ns = INTEGER(node_sizes);
	const int *ifm = INTEGER(imp_fam_mask);
	unsigned int pow_nodes = ncols(imp_fam_mask);
	double alpha = *(REAL(ess));
	int scoring_func = *INTEGER(func);
	double *aflml = REAL(aflml_sexp);
	unsigned int *pa = (unsigned int *)R_alloc(n_nodes, sizeof(unsigned int));
	unsigned int mask1 = 1;
	R_xlen_t nfam = (R_xlen_t)n_nodes * pow_nodes, t, u;

	/* ---- Phase 1: per-cell removal deltas for every feasible family ---- */
	R_xlen_t *offset = (R_xlen_t *)R_alloc(nfam, sizeof(R_xlen_t));
	R_xlen_t total = 0;
	for (s = 0; s < (R_len_t)pow_nodes; s++)
		for (i = 0; i < n_nodes; i++)
		{
			R_xlen_t pos = (R_xlen_t)s * n_nodes + i;
			if (!ifm[pos])
			{
				offset[pos] = -1;
				continue;
			}
			unsigned int npa = get_bits(s, pa, n_nodes);
			R_xlen_t q = 1;
			for (k = 0; k < npa; k++)
				q *= ns[pa[k]];
			offset[pos] = total;
			total += q * ns[i];
		}
	double *cell_delta = (double *)R_alloc(total > 0 ? total : 1, sizeof(double));
	for (t = 0; t < total; t++)
		cell_delta[t] = 0.0;

	unsigned int maxsize;
	R_xlen_t maxcells, maxq;
	admissible_extent(ifm, ns, n_nodes, pow_nodes, &maxsize, &maxcells, &maxq);
	row_patterns rp;
	build_row_patterns(d, n_ex, n_nodes, &rp);
	R_xlen_t K = rp.n;
	double *counts = (double *)R_alloc(maxcells, sizeof(double));
	double *nij = (double *)R_alloc(maxq, sizeof(double));
	for (t = 0; t < maxcells; t++)
		counts[t] = 0.0;
	for (t = 0; t < maxq; t++)
		nij[t] = 0.0;
	R_xlen_t *touched = (R_xlen_t *)R_alloc(K > 0 ? K : 1, sizeof(R_xlen_t));

	subset_walk w;
	walk_init(&w, &rp, n_nodes, ns, maxsize);
	while (walk_next(&w))
	{
		unsigned int P = WALK_SET(&w);
		R_xlen_t q = WALK_Q(&w);
		const R_xlen_t *pc = WALK_PC(&w);
		for (i = 0; i < n_nodes; i++)
		{
			R_xlen_t pos = (R_xlen_t)P * n_nodes + i;
			if (offset[pos] < 0)
				continue;
			int r_i = ns[i];
			const int *xi = rp.x + (R_xlen_t)i * K;
			double n_fam = 0.0;
			R_xlen_t nt = 0, tt;
			for (u = 0; u < K; u++)
				if (xi[u] != NA_INTEGER && pc[u] >= 0)
				{
					R_xlen_t cell = pc[u] * r_i + (xi[u] - 1);
					if (counts[cell] == 0.0)
						touched[nt++] = cell;
					counts[cell] += rp.w[u];
					nij[pc[u]] += rp.w[u];
					n_fam += rp.w[u];
				}
			double a_ij = alpha / (double)q;
			double a_ijk = alpha / (double)(q * r_i);
			double pen_relief = 0.0;
			if (scoring_func == 2 && n_fam >= 2.0)
				pen_relief = 0.5 * (double)q * (double)(r_i - 1) *
							 log(n_fam / (n_fam - 1.0));
			for (tt = 0; tt < nt; tt++)
			{
				R_xlen_t cell = touched[tt];
				double n_ij = nij[cell / r_i], n_ijk = counts[cell];
				cell_delta[offset[pos] + cell] =
					(scoring_func == 0)
						? log(a_ij + n_ij - 1.0) - log(a_ijk + n_ijk - 1.0)
						: loo_f_cell(n_ij) - loo_f_cell(n_ijk) + pen_relief;
			}
			for (tt = 0; tt < nt; tt++)
			{
				counts[touched[tt]] = 0.0;
				nij[touched[tt] / r_i] = 0.0;
			}
		}
	}

	/* ---- Phase 2: per distinct row pattern, shift scores and re-run the SM DP ---- */
	SEXP res;
	PROTECT(res = allocMatrix(INTSXP, n_ex, n_nodes));
	int *out = INTEGER(res);
	for (t = 0; t < (R_xlen_t)n_ex * n_nodes; t++)
		out[t] = 0;
	double *A = (double *)R_alloc(nfam, sizeof(double));
	int *bps = (int *)R_alloc(nfam, sizeof(int));
	double *bss = (double *)R_alloc(nfam, sizeof(double));
	int *sinks = (int *)R_alloc(pow_nodes, sizeof(int));
	double *fsc = (double *)R_alloc(pow_nodes, sizeof(double));

	/* parent-configuration strides and highest element of every candidate set */
	R_xlen_t *qs = (R_xlen_t *)R_alloc(pow_nodes, sizeof(R_xlen_t));
	int *hi = (int *)R_alloc(pow_nodes, sizeof(int));
	R_xlen_t *pcs = (R_xlen_t *)R_alloc(pow_nodes, sizeof(R_xlen_t));
	qs[0] = 1;
	hi[0] = -1;
	for (s = 1; s < (R_len_t)pow_nodes; s++)
	{
		hi[s] = (s == 1) ? 0 : hi[s >> 1] + 1;
		qs[s] = qs[s ^ (mask1 << hi[s])] * ns[hi[s]];
	}

	for (u = 0; u < K; u++)
	{
		unsigned int o = rp.ord[rp.start[u]];
		for (t = 0; t < nfam; t++)
			A[t] = aflml[t];

		/* Shift each feasible family (i, s) by the delta of the pattern's cell,
		 *     pc(s) = sum_k (x_{pa[k]} - 1) * prod_{l<k} r_{pa[l]},
		 * cells laid out by phase 1 as [offset + pc * r_i + (child - 1)]; pc(s) is
		 * obtained from pc(s without its highest element), -1 if a parent is missing
		 * (every family of s then shifts by 0). */
		pcs[0] = 0;
		for (s = 1; s < (R_len_t)pow_nodes; s++)
		{
			unsigned int rest = (unsigned int)s ^ (mask1 << hi[s]);
			int ev = rp.x[u + (R_xlen_t)hi[s] * K];
			pcs[s] = (pcs[rest] < 0 || ev == NA_INTEGER) ? -1 : pcs[rest] + (R_xlen_t)(ev - 1) * qs[rest];
		}
		for (s = 0; s < (R_len_t)pow_nodes; s++)
		{
			if (pcs[s] < 0)
				continue;
			for (i = 0; i < n_nodes; i++)
			{
				R_xlen_t pos = (R_xlen_t)s * n_nodes + i;
				if (offset[pos] < 0)
					continue;
				int ei = rp.x[u + (R_xlen_t)i * K];
				if (ei == NA_INTEGER)
					continue;
				A[pos] += cell_delta[offset[pos] + pcs[s] * ns[i] + (ei - 1)];
			}
		}

		/* find best parents: bss[i,s] = max_{P subset s} A[i,P] */
		for (i = 0; i < n_nodes; i++)
			for (s = 0; s < (R_len_t)pow_nodes; s++)
			{
				R_xlen_t pos = (R_xlen_t)s * n_nodes + i;
				if (s & (mask1 << i))
					continue;
				bps[pos] = (int)s + 1;
				bss[pos] = A[pos];
				unsigned int bi;
				for (bi = 0; bi < n_nodes; bi++)
					if (s & (mask1 << bi))
					{
						R_xlen_t pos2 = (R_xlen_t)(s ^ (mask1 << bi)) * n_nodes + i;
						if (bss[pos2] > bss[pos])
						{
							bps[pos] = bps[pos2];
							bss[pos] = bss[pos2];
						}
					}
			}

		/* find best sinks */
		for (s = 0; s < (R_len_t)pow_nodes; s++)
			sinks[s] = -1;
		fsc[0] = 0.0;
		for (s = 1; s < (R_len_t)pow_nodes; s++)
		{
			double best = R_NegInf;
			for (i = 0; i < n_nodes; i++)
				if (s & (mask1 << i))
				{
					unsigned int nop = s ^ (mask1 << i);
					double sk = fsc[nop] + bss[(R_xlen_t)nop * n_nodes + i];
					if (sinks[s] == -1 || sk > best)
					{
						best = sk;
						sinks[s] = (int)i;
					}
				}
			fsc[s] = best;
		}

		/* backtrack to parent sets */
		unsigned int S = pow_nodes - 1;
		while (S > 0)
		{
			int sink = sinks[S];
			unsigned int nop = S ^ (mask1 << sink);
			int ps = bps[(R_xlen_t)nop * n_nodes + sink];
			out[o + (R_xlen_t)sink * n_ex] = ps - 1;
			S = nop;
		}

		/* same pattern, same optimum */
		unsigned int g;
		for (g = rp.start[u] + 1; g < rp.start[u + 1]; g++)
			for (i = 0; i < n_nodes; i++)
				out[rp.ord[g] + (R_xlen_t)i * n_ex] = out[o + (R_xlen_t)i * n_ex];
	}

	UNPROTECT(1);
	return res;
}

/* backward-compatible BDeu entry point */
SEXP bnstruct_loo_parents(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP ess, SEXP aflml_sexp)
{
	SEXP func, res;
	PROTECT(func = ScalarInteger(0));
	res = bnstruct_loo_parents_score(data, node_sizes, imp_fam_mask, ess, aflml_sexp, func);
	UNPROTECT(1);
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

SEXP bnstruct_compute_lg_ug_cons(SEXP bps, SEXP lower, SEXP upper, SEXP imp_fam_mask)
{
	R_len_t ni, si, bi, pos, pos2;
	R_len_t n_nodes = nrows(lower);
	R_len_t pow_nodes = ncols(lower);
	R_len_t mask = 1;
	int *ifm = INTEGER(imp_fam_mask);
	int *ps = INTEGER(bps);
	double *low = REAL(lower);
	double *up = REAL(upper);
	SEXP res, lower_sexp, upper_sexp, names;

	PROTECT(res = allocVector(VECSXP, 2));
	PROTECT(lower_sexp = allocMatrix(REALSXP, n_nodes, pow_nodes));
	PROTECT(upper_sexp = allocMatrix(REALSXP, n_nodes, pow_nodes));
	PROTECT(names = allocVector(STRSXP, 2));

	SET_STRING_ELT(names, 0, mkChar("lower"));
	SET_STRING_ELT(names, 1, mkChar("upper"));
	setAttrib(res, R_NamesSymbol, names);

	double *lg = REAL(lower_sexp);
	double *ug = REAL(upper_sexp);

	for (ni = 0; ni < n_nodes; ni++)
	{
		double *max_upper = (double *)R_alloc(pow_nodes, sizeof(double));

		for (si = 0; si < pow_nodes; si++)
		{
			pos = si * n_nodes + ni;
			max_upper[si] = ifm[pos] ? up[pos] : R_NegInf;
		}

		for (si = 0; si < pow_nodes; si++)
			for (bi = 0; bi < n_nodes; bi++)
				if (si & (mask << bi))
				{
					pos2 = (si ^ (mask << bi));
					if (max_upper[pos2] > max_upper[si])
						max_upper[si] = max_upper[pos2];
				}

		for (si = 0; si < pow_nodes; si++)
		{
			pos = si * n_nodes + ni;
			int p_star = ps[pos];
			/* bnstruct_fbp leaves bps = 0 for candidate sets containing the node
			 * itself; those (i, U) entries are never read back (U = S \ {i} by
			 * construction), but indexing with p_star - 1 = -1 would read before
			 * the array. Note this is NOT the same as !ifm[pos]: a set excluded by
			 * max.fanin/layering still has a valid bps (the best feasible subset)
			 * and a meaningful Lg/Ug. */
			if (p_star == 0)
			{
				lg[pos] = R_NegInf;
				ug[pos] = R_NegInf;
				continue;
			}
			int p_star_pos = (p_star - 1) * n_nodes + ni;
			lg[pos] = low[p_star_pos];
			ug[pos] = max_upper[si];
		}
	}

	SET_VECTOR_ELT(res, 0, lower_sexp);
	SET_VECTOR_ELT(res, 1, upper_sexp);

	UNPROTECT(4);
	return res;
}

SEXP bnstruct_compute_lg_ug_tight(SEXP aflml, SEXP bps, SEXP lower, SEXP upper, SEXP imp_fam_mask)
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
	SEXP res, lower_sexp, upper_sexp, names;

	PROTECT(res = allocVector(VECSXP, 2));
	PROTECT(lower_sexp = allocMatrix(REALSXP, n_nodes, pow_nodes));
	PROTECT(upper_sexp = allocMatrix(REALSXP, n_nodes, pow_nodes));
	PROTECT(names = allocVector(STRSXP, 2));

	SET_STRING_ELT(names, 0, mkChar("lower"));
	SET_STRING_ELT(names, 1, mkChar("upper"));
	setAttrib(res, R_NamesSymbol, names);

	double *lg = REAL(lower_sexp);
	double *ug = REAL(upper_sexp);

	for (ni = 0; ni < n_nodes; ni++)
	{
		double *max_val = (double *)R_alloc(pow_nodes, sizeof(double));

		for (si = 0; si < pow_nodes; si++)
		{
			pos = si * n_nodes + ni;
			if (ifm[pos])
				max_val[si] = ml[pos] + up[pos];
			else
				max_val[si] = R_NegInf;
		}

		for (si = 0; si < pow_nodes; si++)
			for (bi = 0; bi < n_nodes; bi++)
				if (si & (mask << bi))
				{
					pos2 = (si ^ (mask << bi));
					if (max_val[pos2] > max_val[si])
						max_val[si] = max_val[pos2];
				}

		for (si = 0; si < pow_nodes; si++)
		{
			pos = si * n_nodes + ni;
			int p_star = ps[pos];
			/* see bnstruct_compute_lg_ug_cons: bps = 0 marks a candidate set that
			 * contains the node itself, for which p_star - 1 would index before the
			 * array. Such entries are never read back. */
			if (p_star == 0)
			{
				lg[pos] = R_NegInf;
				ug[pos] = R_NegInf;
				continue;
			}
			int p_star_pos = (p_star - 1) * n_nodes + ni;
			double s_star = ml[p_star_pos];
			lg[pos] = low[p_star_pos];
			ug[pos] = max_val[si] - s_star;
		}
	}

	SET_VECTOR_ELT(res, 0, lower_sexp);
	SET_VECTOR_ELT(res, 1, upper_sexp);

	UNPROTECT(4);
	return res;
}

SEXP bnstruct_compute_lf_uf(SEXP F, SEXP sinks, SEXP bps, SEXP aflml, SEXP Lg, SEXP Ug)
{
	R_len_t si, sink, nop;
	R_len_t n_nodes = nrows(aflml);
	R_len_t pow_nodes = ncols(aflml);
	R_len_t mask = 1;
	SEXP res, lower_sexp, upper_sexp, names;

	int *sinks_p = INTEGER(sinks);
	int *ps = INTEGER(bps);
	double *ml = REAL(aflml);
	double *f = REAL(F);
	double *lg = REAL(Lg);
	double *ug = REAL(Ug);

	PROTECT(res = allocVector(VECSXP, 2));
	PROTECT(lower_sexp = allocVector(REALSXP, pow_nodes));
	PROTECT(upper_sexp = allocVector(REALSXP, pow_nodes));
	PROTECT(names = allocVector(STRSXP, 2));

	SET_STRING_ELT(names, 0, mkChar("lower"));
	SET_STRING_ELT(names, 1, mkChar("upper"));
	setAttrib(res, R_NamesSymbol, names);

	double *lf = REAL(lower_sexp);
	double *uf = REAL(upper_sexp);

	for (si = 0; si < pow_nodes; si++)
	{
		lf[si] = R_NegInf;
		uf[si] = R_NegInf;
	}
	lf[0] = 0.0;
	uf[0] = 0.0;

	for (si = 1; si < pow_nodes; si++)
	{
		sink = sinks_p[si] - 1;
		if (sink < 0)
			continue;
		nop = si ^ (mask << sink);
		lf[si] = lf[nop] + lg[nop * n_nodes + sink];

		double max_uf = R_NegInf;
		R_len_t i;
		for (i = 0; i < n_nodes; i++)
			if (si & (mask << i))
			{
				R_len_t u = si ^ (mask << i);
				int p_star = ps[u * n_nodes + i];
				double A = f[u] + ml[(p_star - 1) * n_nodes + i];
				double UA = uf[u] + ug[u * n_nodes + i];
				double cand = A - f[si] + UA;
				if (cand > max_uf)
					max_uf = cand;
			}
		uf[si] = max_uf;
	}

	SET_VECTOR_ELT(res, 0, lower_sexp);
	SET_VECTOR_ELT(res, 1, upper_sexp);

	UNPROTECT(4);
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
	/* family scores on the distinct row patterns (see family_scores_bounds); identical
	 * to scoring every family with score_node_1 */
	SEXP all = PROTECT(family_scores_bounds(data, node_sizes, imp_fam_mask, *REAL(iss), *INTEGER(func), 1, 0));
	SEXP result = VECTOR_ELT(all, 0);
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

/* ------------------------------------------------------------------------------------
 * Coupled parent-set margins for every (node, candidate set): condition 3 of the CPDAG
 * certificate, computed once for all (i, U) instead of per reachable decision.
 *
 *   margin(i,U) = min over admissible P subset of U, P != P* = bps[i,U] - 1, of
 *                 s(P*) - s(P) + min_k [ lower_k(i,P*) - upper_k(i,P) ]
 *
 * The minima over P and over k commute, so for a common set K of child values
 *
 *   margin(i,U) = min_{k in K} [ s(P*) + lower_k(i,P*) - max_{P subset U, P != P*} (s(P) + upper_k(i,P)) ]
 *
 * and the maximum over the subsets of U excluding P* is a top-2 subset-maximum DP:
 * Theta(r p^2 2^p) for every (i, U). K must be the set of child values whose per-value
 * bounds are finite, and it must be the same for every admissible family of node i
 * (it is whenever each observed value occurs in some row, complete or not). A node
 * where that fails, or where K is empty, is handled by enumerating the submasks of
 * each U with the rule of sm_coupled_parent_witness. For the minimising competitor
 * the margin is recomputed with that rule, so values match the R witness.
 * The decoupled margin s(P*) - s(P) + lower(i,P*) - upper(i,P) (diagnostics) is
 * computed the same way.
 *
 * Returns list(margin, coupled, competitor, decoupled), each n_nodes x pow_nodes;
 * competitor is the 1-based column index of the minimising P; NA where the winning
 * family is infeasible or has no admissible competitor.
 * ------------------------------------------------------------------------------------ */
static void subset_top2(unsigned int n_nodes, unsigned int pow_nodes, unsigned int node,
						const double *v, const int *ifm, double *b1, int *a1, double *b2, int *a2)
{
	unsigned int U, b;
	for (U = 0; U < pow_nodes; U++)
	{
		if (U & (1U << node))
			continue;
		double x1 = R_NegInf, x2 = R_NegInf;
		int i1 = -1, i2 = -1;
		if (ifm[(R_xlen_t)U * n_nodes + node])
			update_best_pair(v[U], (int)U, &x1, &i1, &x2, &i2);
		for (b = 0; b < n_nodes; b++)
			if (U & (1U << b))
			{
				unsigned int V = U ^ (1U << b);
				update_best_pair(b1[V], a1[V], &x1, &i1, &x2, &i2);
				update_best_pair(b2[V], a2[V], &x1, &i1, &x2, &i2);
			}
		b1[U] = x1; a1[U] = i1; b2[U] = x2; a2[U] = i2;
	}
}

SEXP bnstruct_coupled_parent_margins(SEXP aflml_s, SEXP bps_s, SEXP ifm_s, SEXP node_sizes,
									 SEXP lower_k_s, SEXP upper_k_s, SEXP lower_s, SEXP upper_s)
{
	unsigned int n_nodes = nrows(aflml_s), pow_nodes = ncols(aflml_s);
	const double *s = REAL(aflml_s), *lk = REAL(lower_k_s), *uk = REAL(upper_k_s);
	const double *lo = REAL(lower_s), *up = REAL(upper_s);
	const int *bps = INTEGER(bps_s), *ifm = INTEGER(ifm_s), *ns = INTEGER(node_sizes);
	int max_r = INTEGER(getAttrib(lower_k_s, R_DimSymbol))[2];
	R_xlen_t npk = (R_xlen_t)n_nodes * pow_nodes;

	SEXP res, margin_s, coupled_s, comp_s, dec_s, names;
	PROTECT(margin_s = allocMatrix(REALSXP, n_nodes, pow_nodes));
	PROTECT(coupled_s = allocMatrix(REALSXP, n_nodes, pow_nodes));
	PROTECT(comp_s = allocMatrix(INTSXP, n_nodes, pow_nodes));
	PROTECT(dec_s = allocMatrix(REALSXP, n_nodes, pow_nodes));
	double *margin = REAL(margin_s), *coupled = REAL(coupled_s), *dec = REAL(dec_s);
	int *comp = INTEGER(comp_s);
	R_xlen_t t;
	for (t = 0; t < npk; t++)
	{
		margin[t] = NA_REAL; coupled[t] = NA_REAL; dec[t] = NA_REAL; comp[t] = NA_INTEGER;
	}

	double *v = (double *)R_alloc(pow_nodes, sizeof(double));
	double *b1 = (double *)R_alloc(pow_nodes, sizeof(double)), *b2 = (double *)R_alloc(pow_nodes, sizeof(double));
	int *a1 = (int *)R_alloc(pow_nodes, sizeof(int)), *a2 = (int *)R_alloc(pow_nodes, sizeof(int));
	double *best = (double *)R_alloc(pow_nodes, sizeof(double));
	int *best_arg = (int *)R_alloc(pow_nodes, sizeof(int));
	int valid[max_r];

	unsigned int i, U;
	int k;
#define LK(i, P, k) lk[(R_xlen_t)(i) + (R_xlen_t)(P) * n_nodes + (R_xlen_t)(k) * npk]
#define UK(i, P, k) uk[(R_xlen_t)(i) + (R_xlen_t)(P) * n_nodes + (R_xlen_t)(k) * npk]
	for (i = 0; i < n_nodes; i++)
	{
		int r = ns[i];
		/* common set of child values with finite per-value bounds */
		int common = 1, have_ref = 0, n_valid = 0;
		unsigned int P;
		for (P = 0; P < pow_nodes && common; P++)
		{
			if (!ifm[(R_xlen_t)P * n_nodes + i])
				continue;
			for (k = 0; k < r; k++)
			{
				int f = R_FINITE(LK(i, P, k)) && R_FINITE(UK(i, P, k));
				if (!have_ref)
					valid[k] = f;
				else if (valid[k] != f)
				{
					common = 0;
					break;
				}
			}
			have_ref = 1;
		}
		for (k = 0; k < r; k++)
			n_valid += valid[k];

		if (common && have_ref && n_valid > 0)
		{
			for (U = 0; U < pow_nodes; U++)
			{
				best[U] = R_PosInf;
				best_arg[U] = -1;
			}
			for (k = 0; k < r; k++)
			{
				if (!valid[k])
					continue;
				for (P = 0; P < pow_nodes; P++)
					v[P] = s[(R_xlen_t)P * n_nodes + i] + UK(i, P, k);
				subset_top2(n_nodes, pow_nodes, i, v, ifm, b1, a1, b2, a2);
				for (U = 0; U < pow_nodes; U++)
				{
					if (U & (1U << i))
						continue;
					int ps = bps[(R_xlen_t)U * n_nodes + i] - 1;
					double ss = s[(R_xlen_t)ps * n_nodes + i];
					if (!R_FINITE(ss))
						continue;
					double ex = (a1[U] == ps) ? b2[U] : b1[U];
					int exa = (a1[U] == ps) ? a2[U] : a1[U];
					if (exa < 0)
						continue;
					double cand = ss + LK(i, ps, k) - ex;
					if (cand < best[U])
					{
						best[U] = cand;
						best_arg[U] = exa;
					}
				}
			}
			for (U = 0; U < pow_nodes; U++)
			{
				if ((U & (1U << i)) || best_arg[U] < 0)
					continue;
				int ps = bps[(R_xlen_t)U * n_nodes + i] - 1, pc = best_arg[U];
				double cc = R_PosInf;
				for (k = 0; k < r; k++)
					if (valid[k] && LK(i, ps, k) - UK(i, pc, k) < cc)
						cc = LK(i, ps, k) - UK(i, pc, k);
				R_xlen_t pos = (R_xlen_t)U * n_nodes + i;
				coupled[pos] = cc;
				margin[pos] = s[(R_xlen_t)ps * n_nodes + i] - s[(R_xlen_t)pc * n_nodes + i] + cc;
				comp[pos] = pc + 1;
			}
		}
		else
		{
			/* explicit submask enumeration, rule of sm_coupled_parent_witness */
			for (U = 0; U < pow_nodes; U++)
			{
				if (U & (1U << i))
					continue;
				int ps = bps[(R_xlen_t)U * n_nodes + i] - 1;
				double ss = s[(R_xlen_t)ps * n_nodes + i];
				if (!R_FINITE(ss))
					continue;
				double worst = R_PosInf, cc_at = NA_REAL;
				int arg = -1;
				unsigned int Q = U;
				for (;;)
				{
					if (ifm[(R_xlen_t)Q * n_nodes + i] && (int)Q != ps)
					{
						double cc = R_PosInf;
						int any = 0;
						for (k = 0; k < r; k++)
							if (R_FINITE(LK(i, ps, k)) && R_FINITE(UK(i, Q, k)))
							{
								any = 1;
								if (LK(i, ps, k) - UK(i, Q, k) < cc)
									cc = LK(i, ps, k) - UK(i, Q, k);
							}
						if (!any)
							cc = lo[(R_xlen_t)ps * n_nodes + i] - up[(R_xlen_t)Q * n_nodes + i];
						double m = ss - s[(R_xlen_t)Q * n_nodes + i] + cc;
						if (m < worst)
						{
							worst = m;
							arg = (int)Q;
							cc_at = cc;
						}
					}
					if (Q == 0)
						break;
					Q = (Q - 1) & U;
				}
				if (arg >= 0 && R_FINITE(worst))
				{
					R_xlen_t pos = (R_xlen_t)U * n_nodes + i;
					margin[pos] = worst;
					coupled[pos] = cc_at;
					comp[pos] = arg + 1;
				}
			}
		}

		/* decoupled margin (diagnostics): same DP with family-wide bounds */
		{
			unsigned int P2;
			for (P2 = 0; P2 < pow_nodes; P2++)
				v[P2] = s[(R_xlen_t)P2 * n_nodes + i] + up[(R_xlen_t)P2 * n_nodes + i];
			subset_top2(n_nodes, pow_nodes, i, v, ifm, b1, a1, b2, a2);
			for (U = 0; U < pow_nodes; U++)
			{
				if (U & (1U << i))
					continue;
				int ps = bps[(R_xlen_t)U * n_nodes + i] - 1;
				double ss = s[(R_xlen_t)ps * n_nodes + i];
				if (!R_FINITE(ss))
					continue;
				int exa = (a1[U] == ps) ? a2[U] : a1[U];
				if (exa < 0)
					continue;
				dec[(R_xlen_t)U * n_nodes + i] = ss - s[(R_xlen_t)exa * n_nodes + i] +
												 lo[(R_xlen_t)ps * n_nodes + i] - up[(R_xlen_t)exa * n_nodes + i];
			}
		}
	}
#undef LK
#undef UK

	PROTECT(res = allocVector(VECSXP, 4));
	PROTECT(names = allocVector(STRSXP, 4));
	SET_VECTOR_ELT(res, 0, margin_s);
	SET_VECTOR_ELT(res, 1, coupled_s);
	SET_VECTOR_ELT(res, 2, comp_s);
	SET_VECTOR_ELT(res, 3, dec_s);
	SET_STRING_ELT(names, 0, mkChar("margin"));
	SET_STRING_ELT(names, 1, mkChar("coupled"));
	SET_STRING_ELT(names, 2, mkChar("competitor"));
	SET_STRING_ELT(names, 3, mkChar("decoupled"));
	setAttrib(res, R_NamesSymbol, names);
	UNPROTECT(6);
	return res;
}

/* ------------------------------------------------------------------------------------
 * Exact per-family coupling for the coupled sink condition of the CPDAG certificate:
 *
 *   c(i, R, P) = min over rows o of [ Delta(i, R, o) - Delta(i, P, o) ]
 *
 * where Delta is the one-row removal delta of the selected score (0 when the row is
 * incomplete in the family), for every node i, every admissible family P of i and
 * every reference parent set R = refs[i, t] (column t = one witness network).
 * c(i, R, R) = 0 exactly. Rows enter only through their K distinct patterns: each
 * family is counted once on the patterns and its delta evaluated at every pattern,
 * O(K) per family (subset_walk), then compared with the deltas of node i's distinct
 * references, O(D K) per family. Total O(N p log N + D K p 2^p); memory O(p 2^p D)
 * for the output plus O(p K + D K + largest count table).
 * Returns list(values, ref_index): values[i, P + 1, ref_index[i, t] + 1] is the
 * coupling of node i, family P and witness column t; inadmissible entries are 0.
 * ------------------------------------------------------------------------------------ */

/* removal delta of family (i, pa) at every pattern, counting on the patterns directly */
static void family_deltas_at_patterns(const row_patterns *rp, unsigned int n_nodes, const int *ns,
									  unsigned int i, unsigned int sset, double alpha, int scoring_func,
									  unsigned int *pa, double *out)
{
	unsigned int npa = get_bits(sset, pa, n_nodes), k;
	R_xlen_t u, K = rp->n, q = 1;
	int r_i = ns[i];
	for (k = 0; k < npa; k++)
		q *= ns[pa[k]];
	R_xlen_t *pc = (R_xlen_t *)R_alloc(K > 0 ? K : 1, sizeof(R_xlen_t));
	for (u = 0; u < K; u++)
	{
		R_xlen_t idx = 0, stride = 1;
		int xi = rp->x[u + (R_xlen_t)i * K];
		pc[u] = (xi == NA_INTEGER) ? -1 : 0;
		for (k = 0; k < npa && pc[u] >= 0; k++)
		{
			int e = rp->x[u + (R_xlen_t)pa[k] * K];
			if (e == NA_INTEGER)
				pc[u] = -1;
			else
			{
				idx += (R_xlen_t)(e - 1) * stride;
				stride *= ns[pa[k]];
			}
		}
		if (pc[u] >= 0)
			pc[u] = idx * r_i + (xi - 1); /* here pc holds the cell */
	}
	double *counts = (double *)calloc(q * r_i, sizeof(double));
	double *nij = (double *)calloc(q, sizeof(double));
	if (counts == NULL || nij == NULL)
		error("memory allocation failed in bnstruct_coupling_rows");
	double n_fam = 0.0;
	for (u = 0; u < K; u++)
		if (pc[u] >= 0)
		{
			counts[pc[u]] += rp->w[u];
			nij[pc[u] / r_i] += rp->w[u];
			n_fam += rp->w[u];
		}
	double a_ij = alpha / (double)q, a_ijk = alpha / (double)(q * r_i);
	double pen_relief = 0.0;
	if (scoring_func == 2 && n_fam >= 2.0)
		pen_relief = 0.5 * (double)q * (double)(r_i - 1) * log(n_fam / (n_fam - 1.0));
	for (u = 0; u < K; u++)
	{
		if (pc[u] < 0)
		{
			out[u] = 0.0;
			continue;
		}
		double n1 = nij[pc[u] / r_i], n2 = counts[pc[u]];
		out[u] = (scoring_func == 0) ? log(a_ij + n1 - 1.0) - log(a_ijk + n2 - 1.0)
									 : loo_f_cell(n1) - loo_f_cell(n2) + pen_relief;
	}
	free(counts);
	free(nij);
}

SEXP bnstruct_coupling_rows(SEXP data, SEXP node_sizes, SEXP imp_fam_mask, SEXP ess, SEXP func, SEXP refs_s)
{
	const unsigned int *d = (const unsigned int *)INTEGER(data);
	unsigned int n_ex = nrows(data), n_nodes = ncols(data);
	const int *ns = INTEGER(node_sizes), *ifm = INTEGER(imp_fam_mask), *refs = INTEGER(refs_s);
	unsigned int pow_nodes = ncols(imp_fam_mask), n_refs = ncols(refs_s);
	double alpha = *REAL(ess);
	int scoring_func = *INTEGER(func);
	unsigned int *pa = (unsigned int *)R_alloc(n_nodes, sizeof(unsigned int));

	row_patterns rp;
	build_row_patterns(d, n_ex, n_nodes, &rp);
	R_xlen_t K = rp.n;

	/* distinct reference sets per node: ref_index[i, t] = position of refs[i, t] among
	 * node i's distinct references (0-based); the value array holds one slice per
	 * distinct reference, so memory is n_nodes x pow_nodes x (max distinct references) */
	SEXP ref_index_s;
	PROTECT(ref_index_s = allocMatrix(INTSXP, n_nodes, n_refs));
	int *ref_index = INTEGER(ref_index_s);
	unsigned int *distinct = (unsigned int *)R_alloc((R_xlen_t)n_nodes * (n_refs > 0 ? n_refs : 1), sizeof(unsigned int));
	unsigned int *n_dist = (unsigned int *)R_alloc(n_nodes, sizeof(unsigned int));
	unsigned int i, m, max_dist = 1;
	R_xlen_t t, u;
	for (i = 0; i < n_nodes; i++)
	{
		n_dist[i] = 0;
		for (t = 0; t < (R_xlen_t)n_refs; t++)
		{
			unsigned int R = (unsigned int)refs[(R_xlen_t)i + t * n_nodes];
			for (m = 0; m < n_dist[i] && distinct[(R_xlen_t)i * n_refs + m] != R; m++)
				;
			if (m == n_dist[i])
				distinct[(R_xlen_t)i * n_refs + n_dist[i]++] = R;
			ref_index[(R_xlen_t)i + t * n_nodes] = (int)m;
		}
		if (n_dist[i] > max_dist)
			max_dist = n_dist[i];
	}

	SEXP vals_s, dim, res, names;
	PROTECT(vals_s = allocVector(REALSXP, (R_xlen_t)n_nodes * pow_nodes * max_dist));
	PROTECT(dim = allocVector(INTSXP, 3));
	INTEGER(dim)[0] = n_nodes;
	INTEGER(dim)[1] = pow_nodes;
	INTEGER(dim)[2] = max_dist;
	setAttrib(vals_s, R_DimSymbol, dim);
	double *out = REAL(vals_s);
	R_xlen_t npk = (R_xlen_t)n_nodes * pow_nodes;
	for (t = 0; t < npk * max_dist; t++)
		out[t] = 0.0;

	/* deltas of every node's distinct references at the patterns */
	double *dR = (double *)R_alloc((K > 0 ? K : 1) * (R_xlen_t)n_nodes * max_dist, sizeof(double));
	for (i = 0; i < n_nodes; i++)
		for (m = 0; m < n_dist[i]; m++)
			family_deltas_at_patterns(&rp, n_nodes, ns, i, distinct[(R_xlen_t)i * n_refs + m], alpha,
									  scoring_func, pa, dR + ((R_xlen_t)i * max_dist + m) * K);

	unsigned int maxsize;
	R_xlen_t maxcells, maxq;
	admissible_extent(ifm, ns, n_nodes, pow_nodes, &maxsize, &maxcells, &maxq);
	double *counts = (double *)R_alloc(maxcells, sizeof(double));
	double *nij = (double *)R_alloc(maxq, sizeof(double));
	double *dP = (double *)R_alloc(K > 0 ? K : 1, sizeof(double));
	double *dcell = (double *)R_alloc(maxcells, sizeof(double));
	R_xlen_t *touched = (R_xlen_t *)R_alloc(K > 0 ? K : 1, sizeof(R_xlen_t)), nt, tt;
	for (t = 0; t < maxcells; t++)
		counts[t] = 0.0;
	for (t = 0; t < maxq; t++)
		nij[t] = 0.0;

	subset_walk w;
	walk_init(&w, &rp, n_nodes, ns, maxsize);
	while (walk_next(&w))
	{
		unsigned int P = WALK_SET(&w);
		R_xlen_t q = WALK_Q(&w);
		const R_xlen_t *pc = WALK_PC(&w);
		for (i = 0; i < n_nodes; i++)
		{
			R_xlen_t pos = (R_xlen_t)P * n_nodes + i;
			if (!ifm[pos])
				continue;
			int r_i = ns[i];
			const int *xi = rp.x + (R_xlen_t)i * K;
			double n_fam = 0.0;
			nt = 0;
			for (u = 0; u < K; u++)
				if (xi[u] != NA_INTEGER && pc[u] >= 0)
				{
					R_xlen_t cell = pc[u] * r_i + (xi[u] - 1);
					if (counts[cell] == 0.0)
						touched[nt++] = cell;
					counts[cell] += rp.w[u];
					nij[pc[u]] += rp.w[u];
					n_fam += rp.w[u];
				}
			double a_ij = alpha / (double)q, a_ijk = alpha / (double)(q * r_i);
			double pen_relief = 0.0;
			if (scoring_func == 2 && n_fam >= 2.0)
				pen_relief = 0.5 * (double)q * (double)(r_i - 1) * log(n_fam / (n_fam - 1.0));
			/* one delta per occupied cell, then looked up at every pattern */
			for (tt = 0; tt < nt; tt++)
			{
				double n1 = nij[touched[tt] / r_i], n2 = counts[touched[tt]];
				dcell[touched[tt]] = (scoring_func == 0) ? log(a_ij + n1 - 1.0) - log(a_ijk + n2 - 1.0)
														 : loo_f_cell(n1) - loo_f_cell(n2) + pen_relief;
			}
			for (u = 0; u < K; u++)
				dP[u] = (xi[u] == NA_INTEGER || pc[u] < 0) ? 0.0 : dcell[pc[u] * r_i + (xi[u] - 1)];
			for (m = 0; m < n_dist[i]; m++)
			{
				double c = 0.0;
				if (distinct[(R_xlen_t)i * n_refs + m] != P && K > 0)
				{
					const double *dRm = dR + ((R_xlen_t)i * max_dist + m) * K;
					c = R_PosInf;
					for (u = 0; u < K; u++)
					{
						double x = dRm[u] - dP[u];
						if (x < c)
							c = x;
					}
				}
				out[pos + (R_xlen_t)m * npk] = c;
			}
			for (tt = 0; tt < nt; tt++)
			{
				counts[touched[tt]] = 0.0;
				nij[touched[tt] / r_i] = 0.0;
			}
		}
	}

	PROTECT(res = allocVector(VECSXP, 2));
	PROTECT(names = allocVector(STRSXP, 2));
	SET_VECTOR_ELT(res, 0, vals_s);
	SET_VECTOR_ELT(res, 1, ref_index_s);
	SET_STRING_ELT(names, 0, mkChar("values"));
	SET_STRING_ELT(names, 1, mkChar("ref_index"));
	setAttrib(res, R_NamesSymbol, names);
	UNPROTECT(5);
	return res;
}

/* ------------------------------------------------------------------------------------
 * Parent-set margins with the exact per-row coupling:
 *
 *   margin(i,U) = min over admissible P subset of U, P != R, of  s(R) - s(P) + c(i, R, P),
 *   R = P*(U) = bps[i,U] - 1,  c(i, R, P) = min over rows of Delta(i,R,o) - Delta(i,P,o),
 *
 * a lower bound on min_o [ s_{D\o}(R) - s_{D\o}(P) ] (so margin > 0 keeps R optimal among
 * the subsets of U under every single removal). For each node i and each distinct
 * reference R_m (values[i, , m], from bnstruct_coupling_rows), a top-2 subset-maximum DP
 * of s(P) - c(i, R_m, P) gives the margin at every U whose best parent set is R_m.
 * refs[i, t] / ref_index[i, t] list the references and their slices; only U with
 * bps[i,U] - 1 among node i's references are filled. Returns list(margin, competitor)
 * (competitor = 1-based column index), NA elsewhere.
 * ------------------------------------------------------------------------------------ */
SEXP bnstruct_parent_margins_rows(SEXP aflml_s, SEXP bps_s, SEXP ifm_s, SEXP values_s, SEXP refs_s, SEXP ref_index_s)
{
	unsigned int n_nodes = nrows(aflml_s), pow_nodes = ncols(aflml_s), n_refs = ncols(refs_s);
	const double *s = REAL(aflml_s), *vals = REAL(values_s);
	const int *bps = INTEGER(bps_s), *ifm = INTEGER(ifm_s), *refs = INTEGER(refs_s), *ridx = INTEGER(ref_index_s);
	R_xlen_t npk = (R_xlen_t)n_nodes * pow_nodes, t;

	SEXP res, margin_s, comp_s, names;
	PROTECT(margin_s = allocMatrix(REALSXP, n_nodes, pow_nodes));
	PROTECT(comp_s = allocMatrix(INTSXP, n_nodes, pow_nodes));
	double *margin = REAL(margin_s);
	int *comp = INTEGER(comp_s);
	for (t = 0; t < npk; t++)
	{
		margin[t] = NA_REAL;
		comp[t] = NA_INTEGER;
	}
	double *v = (double *)R_alloc(pow_nodes, sizeof(double));
	double *b1 = (double *)R_alloc(pow_nodes, sizeof(double)), *b2 = (double *)R_alloc(pow_nodes, sizeof(double));
	int *a1 = (int *)R_alloc(pow_nodes, sizeof(int)), *a2 = (int *)R_alloc(pow_nodes, sizeof(int));

	unsigned int i, P, U;
	for (i = 0; i < n_nodes; i++)
	{
		/* distinct slices used by node i */
		for (t = 0; t < (R_xlen_t)n_refs; t++)
		{
			int m = ridx[(R_xlen_t)i + t * n_nodes], R = refs[(R_xlen_t)i + t * n_nodes];
			R_xlen_t t2;
			int seen = 0;
			for (t2 = 0; t2 < t; t2++)
				if (ridx[(R_xlen_t)i + t2 * n_nodes] == m)
				{
					seen = 1;
					break;
				}
			if (seen)
				continue;
			double sR = s[(R_xlen_t)R * n_nodes + i];
			if (!R_FINITE(sR))
				continue;
			for (P = 0; P < pow_nodes; P++)
				v[P] = s[(R_xlen_t)P * n_nodes + i] - vals[(R_xlen_t)i + (R_xlen_t)P * n_nodes + (R_xlen_t)m * npk];
			subset_top2(n_nodes, pow_nodes, i, v, ifm, b1, a1, b2, a2);
			for (U = 0; U < pow_nodes; U++)
			{
				if ((U & (1U << i)) || bps[(R_xlen_t)U * n_nodes + i] - 1 != R)
					continue;
				int exa = (a1[U] == R) ? a2[U] : a1[U];
				if (exa < 0)
					continue;
				R_xlen_t pos = (R_xlen_t)U * n_nodes + i;
				margin[pos] = sR - s[(R_xlen_t)exa * n_nodes + i] +
							  vals[(R_xlen_t)i + (R_xlen_t)exa * n_nodes + (R_xlen_t)m * npk];
				comp[pos] = exa + 1;
			}
		}
	}
	PROTECT(res = allocVector(VECSXP, 2));
	PROTECT(names = allocVector(STRSXP, 2));
	SET_VECTOR_ELT(res, 0, margin_s);
	SET_VECTOR_ELT(res, 1, comp_s);
	SET_STRING_ELT(names, 0, mkChar("margin"));
	SET_STRING_ELT(names, 1, mkChar("competitor"));
	setAttrib(res, R_NamesSymbol, names);
	UNPROTECT(4);
	return res;
}
