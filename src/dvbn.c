/* dvbn.c — Algorithm 1 & 2 (DP discretization) with SEXP wrappers
 * Portable C (Linux/macOS/Windows), CRAN-friendly (uses lgammafn from Rmath.h)
 */

#include "dvbn.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <Rmath.h> /* for lgammafn */

#ifndef NA_INTEGER
#define NA_INTEGER (-2147483647 - 1)
#endif

/* ----------------------------- Utilities ----------------------------- */

static void *xmalloc(size_t n)
{
    void *p = malloc(n);
    if (!p)
        error("dvbn: OOM");
    return p;
}
static void *xcalloc(size_t a, size_t b)
{
    void *p = calloc(a, b);
    if (!p)
        error("dvbn: OOM");
    return p;
}

static inline int is_na_int(int v) { return v == NA_INTEGER; }

static inline double get_cont(const double *M, int n_cases, int col, int row)
{
    return M[(size_t)row + (size_t)col * (size_t)n_cases];
}
static inline int get_disc(const int *M, int n_cases, int col, int row)
{
    return M[(size_t)row + (size_t)col * (size_t)n_cases];
}

/* ---- portable argsort by values of X[,col] (column-major) ---- */
typedef struct
{
    double val;
    int idx;
} ValIdx;
static int cmp_ValIdx_by_val(const void *a, const void *b)
{
    const ValIdx *pa = (const ValIdx *)a, *pb = (const ValIdx *)b;
    if (pa->val < pb->val)
        return -1;
    if (pa->val > pb->val)
        return 1;
    if (pa->idx < pb->idx)
        return -1; /* tie-breaker for determinism */
    if (pa->idx > pb->idx)
        return 1;
    return 0;
}
static int *argsort_col(const double *X, int n_cases, int col)
{
    ValIdx *buf = (ValIdx *)xmalloc(sizeof(ValIdx) * (size_t)n_cases);
    for (int r = 0; r < n_cases; ++r)
    {
        buf[r].val = get_cont(X, n_cases, col, r);
        buf[r].idx = r;
    }
    qsort(buf, (size_t)n_cases, sizeof(ValIdx), cmp_ValIdx_by_val);
    int *order = (int *)xmalloc(sizeof(int) * (size_t)n_cases);
    for (int i = 0; i < n_cases; ++i)
        order[i] = buf[i].idx;
    free(buf);
    return order;
}

/* ---- Build tails[]: last position in sorted order for each run of equal values ---- */
static int *build_tails(const double *X, int n_cases, int col, const int *order, int *m_out)
{
    int *tails = (int *)xmalloc(sizeof(int) * (size_t)n_cases);
    int m = 0;
    for (int t = 0; t < n_cases; ++t)
    {
        int i = order[t];
        int same_next = (t + 1 < n_cases) &&
                        (get_cont(X, n_cases, col, i) == get_cont(X, n_cases, col, order[t + 1]));
        if (!same_next)
            tails[m++] = t;
    }
    *m_out = m;
    return tails;
}

static int max_card(const int *ns, int n_nodes)
{
    int mx = 1;
    for (int j = 0; j < n_nodes; ++j)
        if (ns[j] > mx)
            mx = ns[j];
    return mx;
}

/* split prior along sorted positions (last elem forced to 1) */
static double *compute_split_prior(const double *X, int n_cases, int col, const int *order, int lambda)
{
    double *split = (double *)xmalloc(sizeof(double) * (size_t)n_cases);
    double xmin = get_cont(X, n_cases, col, order[0]);
    double xmax = get_cont(X, n_cases, col, order[n_cases - 1]);
    double span = xmax - xmin;
    if (span <= 0.0)
        span = 1.0;
    for (int t = 0; t < n_cases - 1; ++t)
    {
        double a = get_cont(X, n_cases, col, order[t]);
        double b = get_cont(X, n_cases, col, order[t + 1]);
        double v = 1.0 - exp(-(double)lambda * (b - a) / span);
        if (v <= 0.0)
            v = DBL_MIN;
        if (v >= 1.0)
            v = 1.0 - 1e-15;
        split[t] = v;
    }
    split[n_cases - 1] = 1.0;
    return split;
}

/* ----------------- Helpers for category compression ----------------- */
static int cmp_int(const void *a, const void *b)
{
    int va = *(const int *)a, vb = *(const int *)b;
    return (va > vb) - (va < vb);
}
/* given raw int array arr[0..n-1], output:
   - uniq_vals[0..K-1] sorted unique values (allocated here)
   - K (# unique)
   - id[0..n-1] category id in 0..K-1 (allocated here)
*/
static void compress_to_ids(const int *arr, int n, int **uniq_vals, int *K, int **id)
{
    int *tmp = (int *)xmalloc(sizeof(int) * (size_t)n);
    memcpy(tmp, arr, sizeof(int) * (size_t)n);
    qsort(tmp, (size_t)n, sizeof(int), cmp_int);
    int k = 0;
    for (int i = 0; i < n; ++i)
        if (i == 0 || tmp[i] != tmp[i - 1])
            tmp[k++] = tmp[i];
    int *vals = (int *)xmalloc(sizeof(int) * (size_t)k);
    memcpy(vals, tmp, sizeof(int) * (size_t)k);
    free(tmp);

    int *ids = (int *)xmalloc(sizeof(int) * (size_t)n);
    for (int i = 0; i < n; ++i)
    {
        /* binary search in vals[0..k-1] */
        int v = arr[i], lo = 0, hi = k - 1, pos = 0;
        while (lo <= hi)
        {
            int mid = (lo + hi) >> 1;
            if (vals[mid] == v)
            {
                pos = mid;
                break;
            }
            if (vals[mid] < v)
                lo = mid + 1;
            else
                hi = mid - 1;
        }
        ids[i] = pos;
    }
    *uniq_vals = vals;
    *K = k;
    *id = ids;
}

/* --------- interval table for one categorical column over sorted order ----------
   tbl[i,j] = log(n_ij!) - sum_c log(n_ijc!)  (i..j inclusive), i<=j, else +inf
*/
static double *interval_table_single_cat(const int *col_sorted, int n)
{
    int *uniq = NULL, K = 0, *id = NULL;
    compress_to_ids(col_sorted, n, &uniq, &K, &id);
    // // Print unique values and their ids for debugging
    // for (int i = 0; i < K; ++i)
    //     printf("uniq[%d] = %d  ", i, uniq[i]);
    // // also print the n ids one by one separated by a space
    // printf("\nids: ");
    // for (int i = 0; i < n; ++i)
    //     printf("%d ", id[i]);
    // printf("\n");
    // End of debug print
    // Debug print first 10 values of joint_sorted
    // printf("joint_sorted before build: ");
    // for (int i = 0; i < 20 && i < n; ++i)
    //     printf("%d ", col_sorted[i]);
    // printf("\n");
    free(uniq);

    /* prefix counts per class */
    int **pref = (int **)xmalloc(sizeof(int *) * (size_t)K);
    for (int k = 0; k < K; ++k)
        pref[k] = (int *)xcalloc((size_t)n + 1, sizeof(int));
    for (int i = 0; i < n; ++i)
    {
        for (int k = 0; k < K; ++k)
            pref[k][i + 1] = pref[k][i];
        pref[id[i]][i + 1] += 1;
    }

    double *tbl = (double *)xmalloc(sizeof(double) * (size_t)n * (size_t)n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (j < i)
            {
                tbl[(size_t)i * (size_t)n + (size_t)j] = INFINITY;
                continue;
            }
            int tot = j - i + 1;
            double v = lgammafn((double)tot + 1.0);
            for (int k = 0; k < K; ++k)
            {
                int cnt = pref[k][j + 1] - pref[k][i];
                // if (i == 0 && j == 3)
                //     printf("cnt i: %d j: %d k: %d cnt: %d\n", i, j, k, cnt);
                v -= lgammafn((double)cnt + 1.0);
            }
            tbl[(size_t)i * (size_t)n + (size_t)j] = v;
        }
    }
    // // Debug print tbl in upper triangular form with nice formatting
    // int i = 0;
    // printf("tab_parents_contrib: \n");
    // for (int j = i; j < 20; ++j)
    // {
    //     if (j == i)
    //         printf("%10.4f ", tbl[(size_t)i * (size_t)n + (size_t)j]);
    //     else
    //         printf("%10.4f ", tbl[(size_t)i * (size_t)n + (size_t)j]);
    // }
    // printf("\n");

    // // End of debug print
    for (int k = 0; k < K; ++k)
        free(pref[k]);
    free(pref);
    free(id);
    return tbl;
}

/* spouse-child table over sorted order (child_sorted, sp_sorted both length n)
   For each spouse state s:
     stars&bars: ln C(tot + C - 1, C - 1) = lgamma(tot + C) - lgamma(C) - lgamma(tot+1)
     + multinomial: lgamma(tot+1) - sum_c lgamma(cnt_c+1)
   => sum over spouse states
*/
static double *interval_table_child_spouse(const int *child_sorted, const int *sp_sorted, int n)
{
    int *uniqC = NULL, C = 0, *cid = NULL;
    compress_to_ids(child_sorted, n, &uniqC, &C, &cid);
    free(uniqC);

    int *uniqS = NULL, S = 0, *sid = NULL;
    compress_to_ids(sp_sorted, n, &uniqS, &S, &sid);
    free(uniqS);

    /* prefix counts pref[S][C][n+1] */
    int ***pref = (int ***)xmalloc(sizeof(int **) * (size_t)S);
    for (int s = 0; s < S; ++s)
    {
        pref[s] = (int **)xmalloc(sizeof(int *) * (size_t)C);
        for (int c = 0; c < C; ++c)
            pref[s][c] = (int *)xcalloc((size_t)n + 1, sizeof(int));
    }
    for (int i = 0; i < n; ++i)
    {
        for (int s = 0; s < S; ++s)
            for (int c = 0; c < C; ++c)
                pref[s][c][i + 1] = pref[s][c][i];
        pref[sid[i]][cid[i]][i + 1] += 1;
    }

    double *tbl = (double *)xmalloc(sizeof(double) * (size_t)n * (size_t)n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (j < i)
            {
                tbl[(size_t)i * (size_t)n + (size_t)j] = INFINITY;
                continue;
            }
            double acc = 0.0;
            for (int s = 0; s < S; ++s)
            {
                int tot = 0;
                for (int c = 0; c < C; ++c)
                    tot += (pref[s][c][j + 1] - pref[s][c][i]);
                if (tot == 0)
                    continue;
                /* stars&bars */
                acc += lgammafn((double)tot + (double)C) - lgammafn((double)C) - lgammafn((double)tot + 1.0);
                /* multinomial over child categories given spouse */
                double v = lgammafn((double)tot + 1.0);
                for (int c = 0; c < C; ++c)
                {
                    int cnt = pref[s][c][j + 1] - pref[s][c][i];
                    v -= lgammafn((double)cnt + 1.0);
                }
                acc += v;
            }
            tbl[(size_t)i * (size_t)n + (size_t)j] = acc;
        }
    }

    for (int s = 0; s < S; ++s)
    {
        for (int c = 0; c < C; ++c)
            free(pref[s][c]);
        free(pref[s]);
    }
    free(pref);
    free(cid);
    free(sid);
    return tbl;
}

/* parents contribution: if approx=0, P_sorted is joint-coded parents; if approx=1, call per single parent */
static double *interval_table_parents(const int *P_sorted, int n, int approx, int n_pa)
{
    (void)approx;
    (void)n_pa;
    return interval_table_single_cat(P_sorted, n);
}

/* Build joint-coded parent column over sorted order (exact mode). */
static void build_joint_parents_sorted(const int *D, int n_cases, const int *order,
                                       const int *pa_idx, int n_pa, int *out_joint)
{
    for (int t = 0; t < n_cases; ++t)
    {
        int r = order[t];

        unsigned int h = 2166136261u;
        for (int p = 0; p < n_pa; ++p)
        {

            int v = get_disc(D, n_cases, pa_idx[p], r);
            h ^= (unsigned int)v;
            h *= 16777619u;
        }
        out_joint[t] = (int)(h & 0x7fffffff);
    }
}

/* Build Hfull (n x n) = parent term + sum child|spouses terms, over sorted order of target */
static double *build_Hfull(const double *X, const int *D,
                           int n_nodes, int n_cases, int target,
                           const int *order,
                           const int *ns,
                           const int *pa_idx, int n_pa,
                           const int *children, int n_ch,
                           SEXP parents_list, int approx_parents)
{
    (void)n_nodes;
    (void)ns;
    int n = n_cases;
    double *H = (double *)xcalloc((size_t)n * (size_t)n, sizeof(double));

    /* parents */
    if (n_pa > 0)
    {
        if (approx_parents)
        {
            for (int pi = 0; pi < n_pa; ++pi)
            {
                int pcol = pa_idx[pi];
                int *p_sorted = (int *)xmalloc(sizeof(int) * (size_t)n);
                for (int t = 0; t < n; ++t)
                    p_sorted[t] = get_disc(D, n_cases, pcol, order[t]);
                double *T = interval_table_single_cat(p_sorted, n);
                for (int i = 0; i < n; ++i)
                    for (int j = i; j < n; ++j)
                        H[(size_t)i * (size_t)n + (size_t)j] += T[(size_t)i * (size_t)n + (size_t)j];
                free(T);
                free(p_sorted);
            }
        }
        else
        {
            int *joint_sorted = (int *)xmalloc(sizeof(int) * (size_t)n);

            build_joint_parents_sorted(D, n_cases, order, pa_idx, n_pa, joint_sorted);
            double *T = interval_table_parents(joint_sorted, n, 0, n_pa);
            // printf("Hfull parents contribution: \n");
            for (int i = 0; i < n; ++i)
            {
                for (int j = i; j < n; ++j)
                {
                    H[(size_t)i * (size_t)n + (size_t)j] += T[(size_t)i * (size_t)n + (size_t)j];
                    // // Debug print H if i==0 and j<20
                    // if (i == 0 && j < 20)
                    // {
                    //     printf("%f, ", H[(size_t)i * (size_t)n + (size_t)j]);
                    // }
                }
            }
            free(T);
            free(joint_sorted);
        }
    }

    /* children + spouses(child) */
    for (int c = 0; c < n_ch; ++c)
    {
        int chcol = children[c];
        int *child_sorted = (int *)xmalloc(sizeof(int) * (size_t)n);
        for (int t = 0; t < n; ++t)
            child_sorted[t] = get_disc(D, n_cases, chcol, order[t]);

        /* spouses = parents(child) \ {target}; hash-combine them */
        SEXP par_ch = VECTOR_ELT(parents_list, chcol);
        int nsp_all = LENGTH(par_ch);
        const int *par_ch_idx = (nsp_all > 0) ? INTEGER(par_ch) : NULL;

        int *sp_sorted = (int *)xmalloc(sizeof(int) * (size_t)n);
        if (nsp_all == 0)
        {
            for (int t = 0; t < n; ++t)
                sp_sorted[t] = 0;
        }
        else
        {
            for (int t = 0; t < n; ++t)
            {
                int r = order[t];
                unsigned int h = 2166136261u;
                int any_na = 0, used = 0;
                for (int s = 0; s < nsp_all; ++s)
                {
                    int sidx = par_ch_idx[s];
                    if (sidx == target)
                        continue; /* exclude the target itself */
                    int v = get_disc(D, n_cases, sidx, r);

                    h ^= (unsigned int)v;
                    h *= 16777619u;
                    used = 1;
                }
                if (!used)
                {
                    sp_sorted[t] = 0;
                }
                else
                    sp_sorted[t] = (int)(h & 0x7fffffff);
            }
        }

        double *T = interval_table_child_spouse(child_sorted, sp_sorted, n);
        // Debug print T child-spouse contribution
        // printf("Hfull child-spouse contribution: \n");
        // int i = 0;

        // for (int j = i; j < 20; ++j)
        // {
        //     printf("%10.4f ", T[(size_t)i * (size_t)n + (size_t)j]);
        // }
        // printf("\n");

        // End of debug print
        for (int i = 0; i < n; ++i)
            for (int j = i; j < n; ++j)
                H[(size_t)i * (size_t)n + (size_t)j] += T[(size_t)i * (size_t)n + (size_t)j];
        free(T);
        free(child_sorted);
        free(sp_sorted);
    }

    //  Add -log(p_class_distribution_for_parent)
    int parent_size = 1;
    if (n_pa > 0)
    {
        for (int pi = 0; pi < n_pa; ++pi)
        {
            parent_size *= ns[pa_idx[pi]];
        }
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            H[(size_t)i * (size_t)n + (size_t)j] += lgammafn(j - i + (double)parent_size + 1.0) - lgammafn(j - i + 2.0) - lgammafn((double)parent_size);
        }
    }

    return H;
}

/* Map Hfull over sorted positions to Hblock over unique value blocks (m x m) using tails[] */
static double *map_to_blocks(const double *Hfull, int n, const int *tails, int m)
{
    double *Hb = (double *)xmalloc(sizeof(double) * (size_t)m * (size_t)m);
    for (int u = 0; u < m; ++u)
    {
        int tL = (u == 0) ? 0 : (tails[u - 1] + 1);
        for (int v = u; v < m; ++v)
        {
            int tR = tails[v];
            Hb[(size_t)u * (size_t)m + (size_t)v] = Hfull[(size_t)tL * (size_t)n + (size_t)tR];
        }
    }
    return Hb;
}
/* DP -> cuts (midpoints), excluding min/max, using full H matrix */
static void dp_and_edges_fullH(const double *Hfull, int n_cases,
                               const double *X, int col, const int *order,
                               const double *split, int lambda,
                               double **cuts_out, int *n_cuts_out)
{
    // Step 1: identify heads and tails
    int *heads = (int *)xmalloc(sizeof(int) * (size_t)n_cases);
    int *tails = (int *)xmalloc(sizeof(int) * (size_t)n_cases);
    int n_heads = 0, n_tails = 0;

    heads[n_heads++] = 0;
    for (int i = 1; i < n_cases; ++i)
    {
        double prev = get_cont(X, n_cases, col, order[i - 1]);
        double cur = get_cont(X, n_cases, col, order[i]);
        if (cur != prev)
        {
            heads[n_heads++] = i;
        }
    }
    // // debug print first 40 heads
    // if (col == 0)
    //     for (int i = 0; i < n_heads && i < 40; ++i)
    //     {
    //         printf("head[%d] = %d (value: %f)\n", i, heads[i], get_cont(X, n_cases, col, order[heads[i]]));
    //     }
    tails[n_tails++] = n_cases - 1;
    for (int i = n_cases - 2; i >= 0; --i)
    {
        double next = get_cont(X, n_cases, col, order[i + 1]);
        double cur = get_cont(X, n_cases, col, order[i]);
        if (cur != next)
        {
            tails[n_tails++] = i;
        }
    }
    // tails are collected backwards, reverse them
    for (int i = 0; i < n_tails / 2; ++i)
    {
        int tmp = tails[i];
        tails[i] = tails[n_tails - 1 - i];
        tails[n_tails - 1 - i] = tmp;
    }
    // // debug print first 40 tails
    // if (col == 0)
    //     for (int i = 0; i < n_tails && i < 40; ++i)
    //     {
    //         printf("tail[%d] = %d (value: %f)\n", i, tails[i], get_cont(X, n_cases, col, order[tails[i]]));
    //     }

    double xmin = get_cont(X, n_cases, col, order[0]);
    double xmax = get_cont(X, n_cases, col, order[n_cases - 1]);
    double span = (xmax > xmin) ? (xmax - xmin) : 1.0;
    // if (col == 0)
    // {
    //     // Debug print xmin, xmax, span
    //     printf("col: %d, xmin: %f, xmax: %f, span: %f\n", col, xmin, xmax, span);
    //     // ALSO print first 20 values of Hfull in the first row
    //     printf("Hfull first row: ");
    //     for (int i = 0; i < 20 && i < n_cases; ++i)
    //     {
    //         printf("%4.4f ", Hfull[0 * (size_t)n_cases + i]);
    //     }
    //     printf("\n");
    //     // Also fitst 20 values of split
    //     printf("split: ");
    //     for (int i = 0; i < 20 && i < n_cases; ++i)
    //     {
    //         printf("%4.4f ", split[i]);
    //     }
    // }

    // Step 2: DP arrays
    double *best = (double *)xmalloc(sizeof(double) * (size_t)n_tails);
    int **choice = (int **)xcalloc((size_t)n_tails, sizeof(int *));
    int *clen = (int *)xcalloc((size_t)n_tails, sizeof(int));

    for (int a = 0; a < n_tails; ++a)
    {
        int tail_a = tails[a];
        if (a == 0)
        {
            best[0] = -log(split[tail_a]) + Hfull[0 * (size_t)n_cases + tail_a];
            choice[0] = (int *)xmalloc(sizeof(int));
            choice[0][0] = tail_a;
            clen[0] = 1;
            // Debug print here
            // if (col == 0)
            // {
            //     printf("H: %f, split: %f, tail_a: %d, best[0]: %f, choice[0][0]: %d\n", Hfull[0 * (size_t)n_cases + tail_a], split[tail_a], tail_a, best[0], choice[0][0]);
            // }
        }
        else
        {
            double bv = INFINITY;
            int *bs = NULL;
            int bl = 0;
            for (int b = 0; b <= a; ++b)
            {
                double val;
                if (b == a)
                {
                    double frac = (get_cont(X, n_cases, col, order[tail_a]) - xmin) / span;
                    val = frac * (double)lambda - log(split[tail_a]) + Hfull[0 * (size_t)n_cases + tail_a];
                }
                else
                {
                    int head_next = heads[b + 1];
                    double frac = (get_cont(X, n_cases, col, order[tail_a]) -
                                   get_cont(X, n_cases, col, order[head_next])) /
                                  span;
                    val = best[b] + frac * (double)lambda - log(split[tail_a]) + Hfull[(size_t)head_next * (size_t)n_cases + tail_a];
                }
                if (val < bv)
                {
                    bv = val;
                    if (bs)
                        free(bs);
                    if (b == a)
                    {
                        bs = (int *)xmalloc(sizeof(int));
                        bs[0] = tail_a;
                        bl = 1;
                    }
                    else
                    {
                        bl = clen[b] + 1;
                        bs = (int *)xmalloc(sizeof(int) * (size_t)bl);
                        memcpy(bs, choice[b], sizeof(int) * (size_t)clen[b]);
                        bs[bl - 1] = tail_a;
                    }
                }
            }
            best[a] = bv;
            choice[a] = bs;
            clen[a] = bl;
        }
    }

    int K = clen[n_tails - 1]; /* number of bins */
    int ncuts = (K > 1) ? (K - 1) : 0;
    double *cuts = NULL;
    if (ncuts > 0)
    {
        cuts = (double *)xmalloc(sizeof(double) * (size_t)ncuts);
        for (int i = 0; i < ncuts; ++i)
        {
            int t = choice[n_tails - 1][i];
            double a_val = get_cont(X, n_cases, col, order[t]);
            double b_val = get_cont(X, n_cases, col, order[t + 1]);
            cuts[i] = 0.5 * (a_val + b_val);
        }
    }

    for (int i = 0; i < n_tails; ++i)
        if (choice[i])
            free(choice[i]);
    free(choice);
    free(clen);
    free(best);
    free(heads);
    free(tails);

    *cuts_out = cuts;
    *n_cuts_out = ncuts;
}

/* DP -> cuts (midpoints), excluding min/max */
static void dp_and_edges(const double *Hb, int m,
                         const double *X, int n_cases, int col, const int *order, const int *tails,
                         const double *split, int lambda,
                         double **cuts_out, int *n_cuts_out)
{
    int *heads = (int *)xmalloc(sizeof(int) * (size_t)m);
    int hcnt = 0;
    heads[hcnt++] = 0;
    for (int t = 1; t < n_cases && hcnt < m; ++t)
    {
        double prev = get_cont(X, n_cases, col, order[t - 1]);
        double cur = get_cont(X, n_cases, col, order[t]);
        if (cur != prev)
            heads[hcnt++] = t;
    }
    double xmin = get_cont(X, n_cases, col, order[0]);
    double xmax = get_cont(X, n_cases, col, order[n_cases - 1]);
    double span = xmax - xmin;
    if (span <= 0.0)
        span = 1.0;

    double *best = (double *)xmalloc(sizeof(double) * (size_t)m);
    int **choice = (int **)xcalloc((size_t)m, sizeof(int *));
    int *clen = (int *)xcalloc((size_t)m, sizeof(int));

    for (int a = 0; a < m; ++a)
    {
        int tail_a = tails[a];
        if (a == 0)
        {
            best[0] = -log(split[tail_a]) + Hb[0 * (size_t)m + 0];
            choice[0] = (int *)xmalloc(sizeof(int) * 1);
            choice[0][0] = tail_a;
            clen[0] = 1;
        }
        else
        {
            double bv = INFINITY;
            int *bs = NULL;
            int bl = 0;
            for (int b = 0; b <= a; ++b)
            {
                double val;
                if (b == a)
                {
                    double frac = (get_cont(X, n_cases, col, order[tail_a]) - xmin) / span;
                    val = frac * (double)lambda - log(split[tail_a]) + Hb[0 * (size_t)m + a];
                }
                else
                {
                    int head_next = heads[b + 1];
                    double frac = (get_cont(X, n_cases, col, order[tail_a]) - get_cont(X, n_cases, col, order[head_next])) / span;
                    val = best[b] + frac * (double)lambda - log(split[tail_a]) + Hb[(size_t)(b + 1) * (size_t)m + (size_t)a];
                }
                if (val < bv)
                {
                    bv = val;
                    if (bs)
                        free(bs);
                    if (b == a)
                    {
                        bs = (int *)xmalloc(sizeof(int) * 1);
                        bs[0] = tail_a;
                        bl = 1;
                    }
                    else
                    {
                        bl = clen[b] + 1;
                        bs = (int *)xmalloc(sizeof(int) * (size_t)bl);
                        memcpy(bs, choice[b], sizeof(int) * (size_t)clen[b]);
                        bs[bl - 1] = tail_a;
                    }
                }
            }
            best[a] = bv;
            choice[a] = bs;
            clen[a] = bl;
        }
    }

    int K = clen[m - 1]; /* bins */
    int ncuts = (K > 1) ? (K - 1) : 0;
    double *cuts = NULL;
    if (ncuts > 0)
    {
        cuts = (double *)xmalloc(sizeof(double) * (size_t)ncuts);
        for (int i = 0; i < ncuts; ++i)
        {
            int t = choice[m - 1][i];
            double a = get_cont(X, n_cases, col, order[t]);
            double b = get_cont(X, n_cases, col, order[t + 1]);
            cuts[i] = 0.5 * (a + b);
        }
    }

    for (int i = 0; i < m; ++i)
        if (choice[i])
            free(choice[i]);
    free(choice);
    free(clen);
    free(best);
    free(heads);

    *cuts_out = cuts;
    *n_cuts_out = ncuts;
}

/* Core Algorithm 1 for one variable; returns malloc'ed cuts array */
static void algo1_discretize_core(const double *X, const int *D,
                                  int n_nodes, int n_cases, const int *ns,
                                  int cont_col, int target,
                                  SEXP parents_list, SEXP children_list,
                                  int approx_parents,
                                  double **cuts_out, int *n_cuts_out)
{
    *cuts_out = NULL;
    *n_cuts_out = 0;

    int *order = argsort_col(X, n_cases, cont_col);
    // print the order for debugging
    // printf("order for column %d: ", cont_col);
    // for (int i = 0; i < 20; ++i)
    //     printf("%d ", order[i]);
    // printf("\n");
    // int m = 0;
    // int *tails = build_tails(X, n_cases, cont_col, order, &m);
    // if (m == 0)
    // {
    //     free(order);
    //     free(tails);
    //     return;
    // }

    SEXP pa_s = VECTOR_ELT(parents_list, target);
    int n_pa = LENGTH(pa_s);
    const int *pa_idx = (n_pa > 0) ? INTEGER(pa_s) : NULL;
    // Debug print of pa_idx
    // printf("parents_list of varcont[%d] and target[%d]: ", cont_col, target);
    // {
    //     for (int i = 0; i < n_pa; ++i)
    //         printf("%d ", pa_idx[i]);
    // }
    // printf("\n");
    // printf("Value of the parent of target at row 10: %d ", get_disc(D, n_cases, 3, 10));
    // End of debug print
    SEXP ch_s = VECTOR_ELT(children_list, target);
    int n_ch = LENGTH(ch_s);
    const int *children = (n_ch > 0) ? INTEGER(ch_s) : NULL;

    double *Hfull = build_Hfull(X, D, n_nodes, n_cases, target, order, ns,
                                pa_idx, n_pa, children, n_ch,
                                parents_list, approx_parents);
    // double *Hb = map_to_blocks(Hfull, n_cases, tails, m);

    int lambda = max_card(ns, n_nodes);
    double *split = compute_split_prior(X, n_cases, cont_col, order, lambda);

    // Debugging output
    // printf("Discretizing column %d (target=%d) with %d cases, %d parents, %d children\n",
    //        cont_col, target, n_cases, n_pa, n_ch);
    // printf("Hfull size: %d x %d\n", n_cases, n_cases);
    // printf("Hfull (first %d): ", 100);
    // for (int i = 0; i < n_cases * n_cases && i < 100; ++i)
    // {
    //     printf("%.8f", Hfull[i]);
    //     if (i < n_cases * n_cases - 1 && i < 99)
    //         printf(", ");
    // }
    // printf("\n");
    // printf("split: ");
    // for (int i = 0; i < n_cases; ++i)
    //     printf("%.8f ", split[i]);
    // printf("\n");
    // SEXP pa_s_dbg = VECTOR_ELT(parents_list, target);
    // int n_pa_dbg = LENGTH(pa_s_dbg);
    // printf("parents_list[%d]: ", target);
    // if (n_pa_dbg > 0)
    // {
    //     const int *pa_idx_dbg = INTEGER(pa_s_dbg);
    //     for (int i = 0; i < n_pa_dbg; ++i)
    //         printf("%d ", pa_idx_dbg[i]);
    // }
    // printf("\n");

    // SEXP ch_s_dbg = VECTOR_ELT(children_list, target);
    // int n_ch_dbg = LENGTH(ch_s_dbg);
    // printf("children_list[%d]: ", target);
    // if (n_ch_dbg > 0)
    // {
    //     const int *ch_idx_dbg = INTEGER(ch_s_dbg);
    //     for (int i = 0; i < n_ch_dbg; ++i)
    //         printf("%d ", ch_idx_dbg[i]);
    // }
    // printf("\n");

    // end of debugging output

    double *cuts = NULL;
    int ncuts = 0;
    // Try to change the dp_and_edges call to dp_and_edges_fullH
    // dp_and_edges(Hb, m, X, n_cases, cont_col, order, tails, split, lambda, &cuts, &ncuts);
    dp_and_edges_fullH(Hfull, n_cases, X, cont_col, order, split, lambda, &cuts, &ncuts);
    free(Hfull);
    // free(Hb);
    free(split);
    free(order);
    // free(tails);
    *cuts_out = cuts;
    *n_cuts_out = ncuts;
}

static int cuts_equal(const double *a, int na, const double *b, int nb)
{
    if (na != nb)
        return 0;
    for (int i = 0; i < na; ++i)
        if (fabs(a[i] - b[i]) > 1e-12)
            return 0;
    return 1;
}

/* ----------------------------- SEXP API ----------------------------- */

/* Algorithm 1 */
SEXP bnstruct_dvbn_discretize_one(SEXP data_cont, SEXP data_disc,
                                  SEXP n_nodes, SEXP n_cases,
                                  SEXP ns, SEXP cont_col, SEXP target,
                                  SEXP parents_list, SEXP children_list,
                                  SEXP approx_parents)
{
    if (!isReal(data_cont) || !isInteger(data_disc))
        error("data_cont must be REAL matrix; data_disc must be INTEGER matrix.");
    if (!isInteger(n_nodes) || !isInteger(n_cases))
        error("n_nodes/n_cases must be integer scalars.");
    if (!isInteger(ns))
        error("ns must be integer vector.");
    if (!isInteger(cont_col) || !isInteger(target))
        error("cont_col/target must be integer scalars.");
    if (!isNewList(parents_list) || !isNewList(children_list))
        error("parents_list and children_list must be lists.");
    if (!isInteger(approx_parents))
        error("approx_parents must be integer scalar.");

    const double *X = REAL(data_cont);
    const int *D = INTEGER(data_disc);
    int n_nodes_c = INTEGER(n_nodes)[0];
    int n_cases_c = INTEGER(n_cases)[0];
    const int *ns_c = INTEGER(ns);
    int cont_col_c = INTEGER(cont_col)[0];
    int target_c = INTEGER(target)[0];
    int approx_c = INTEGER(approx_parents)[0];

    double *cuts = NULL;
    int ncuts = 0;
    algo1_discretize_core(X, D, n_nodes_c, n_cases_c, ns_c,
                          cont_col_c, target_c,
                          parents_list, children_list,
                          approx_c,
                          &cuts, &ncuts);

    SEXP out = PROTECT(allocVector(REALSXP, ncuts));
    for (int i = 0; i < ncuts; ++i)
        REAL(out)
    [i] = cuts[i];
    if (cuts)
        free(cuts);
    UNPROTECT(1);
    return out;
}

/* Helper function to discretize a single column based on cut points */
static void discretize_column(const double *X, int *D, int n_cases, int col_x, int col_d,
                              const double *cuts, int ncuts)
{
    for (int i = 0; i < n_cases; ++i)
    {
        double val = get_cont(X, n_cases, col_x, i);
        int bin = 0;
        while (bin < ncuts && val > cuts[bin])
            bin++;
        D[i + (size_t)col_d * (size_t)n_cases] = bin;
    }
}
/* Algorithm 2 */
SEXP bnstruct_dvbn_discretize_all(SEXP data_cont, SEXP data_disc,
                                  SEXP n_nodes, SEXP n_cases,
                                  SEXP n_cont, SEXP cont_index,
                                  SEXP ns, SEXP parents_list, SEXP children_list,
                                  SEXP n_cycles, SEXP approx_parents)
{
    if (!isReal(data_cont) || !isInteger(data_disc))
        error("data_cont must be REAL matrix; data_disc must be INTEGER matrix.");
    if (!isInteger(n_nodes) || !isInteger(n_cases) || !isInteger(n_cont) || !isInteger(n_cycles))
        error("n_nodes, n_cases, n_cont, n_cycles must be integer scalars.");
    if (!isInteger(cont_index) || !isInteger(ns))
        error("cont_index/ns must be integer vectors.");
    if (!isNewList(parents_list) || !isNewList(children_list))
        error("parents_list and children_list must be lists.");
    if (!isInteger(approx_parents))
        error("approx_parents must be integer scalar.");

    const double *X = REAL(data_cont);
    int *D = INTEGER(data_disc);
    int n_nodes_c = INTEGER(n_nodes)[0];
    int n_cases_c = INTEGER(n_cases)[0];
    int n_cont_c = INTEGER(n_cont)[0];
    const int *cont_idx = INTEGER(cont_index);
    int *ns_c = INTEGER(ns);
    int max_cycles = INTEGER(n_cycles)[0];
    int approx_c = INTEGER(approx_parents)[0];

    double **cuts = (double **)xcalloc((size_t)n_cont_c, sizeof(double *));
    int *ncuts = (int *)xcalloc((size_t)n_cont_c, sizeof(int));

    // printf("ns_c: ");
    // for (int i = 0; i < n_nodes_c; ++i)
    //     printf("%d ", ns_c[i]);
    // printf("\n");
    int changed = 1;
    int iter = 0;
    while (changed && iter < max_cycles)
    {
        changed = 0;
        for (int j = 0; j < n_cont_c; ++j)
        {
            int col_d = cont_idx[j];
            int col_x = j; // position in X (data_cont)

            // Check if this variable has parents or children
            SEXP pa_s = VECTOR_ELT(parents_list, col_d);
            int n_pa = LENGTH(pa_s);
            SEXP ch_s = VECTOR_ELT(children_list, col_d);
            int n_ch = LENGTH(ch_s);

            double *new_cuts = NULL;
            int new_n = 0;

            if (n_pa == 0 && n_ch == 0)
            {
                // === Equal-width fallback (like Julia) ===
                int l_card = max_card(ns_c, n_nodes_c);
                double minv = R_PosInf, maxv = R_NegInf;
                for (int i = 0; i < n_cases_c; ++i)
                {
                    double v = get_cont(X, n_cases_c, col_x, i);
                    if (v < minv)
                        minv = v;
                    if (v > maxv)
                        maxv = v;
                }
                new_n = l_card - 1; // number of cuts
                if (new_n > 0)
                {
                    new_cuts = (double *)xcalloc((size_t)new_n, sizeof(double));
                    for (int k = 1; k < l_card; ++k)
                    {
                        new_cuts[k - 1] = minv + (maxv - minv) * ((double)k / (double)l_card);
                    }
                }
            }
            else
            {

                algo1_discretize_core(X, D, n_nodes_c, n_cases_c, ns_c,
                                      j, col_d,
                                      parents_list, children_list,
                                      approx_c,
                                      &new_cuts, &new_n);
                printf("iter %d col_d=%d ns=%d lambda=%d ncuts=%d\n",
                       iter, col_d, ns_c[col_d], max_card(ns_c, n_nodes_c), new_n);
            }
            if (!cuts_equal(cuts[j], ncuts[j], new_cuts, new_n))
            {
                if (cuts[j])
                    free(cuts[j]);
                cuts[j] = new_cuts;
                ncuts[j] = new_n;
                changed = 1;

                // Discretize the column in data_disc
                discretize_column(X, D, n_cases_c, j, col_d, new_cuts, new_n);
                // Update node size for that column
                ns_c[col_d] = new_n + 1; // +1 because cuts define bins, not categories
            }
            else
            {
                if (new_cuts)
                    free(new_cuts);
            }
        }
        iter++;
    }

    SEXP out = PROTECT(allocVector(VECSXP, n_cont_c));
    for (int j = 0; j < n_cont_c; ++j)
    {
        SEXP v = PROTECT(allocVector(REALSXP, ncuts[j]));
        for (int i = 0; i < ncuts[j]; ++i)
            REAL(v)
        [i] = cuts[j][i];
        SET_VECTOR_ELT(out, j, v);
        UNPROTECT(1);
        if (cuts[j])
            free(cuts[j]);
    }
    free(cuts);
    free(ncuts);

    UNPROTECT(1);
    return out;
}
