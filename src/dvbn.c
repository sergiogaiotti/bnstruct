/* dvbn.c — Algorithm 1 & 2 (DP discretization) with SEXP wrappers
 * Portable C (Linux/macOS/Windows), CRAN-friendly (uses lgammafn from Rmath.h)
 */

#include "dvbn.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <Rmath.h> /* for lgammafn */

#ifdef _OPENMP
#include <omp.h>
#endif

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
/* Returns a malloc'ed array of Markov blanket variable indices (excluding target), and sets *mb_size_out */
static int *markov_blanket(int n_nodes, SEXP parents_list, SEXP children_list, int target, int *mb_size_out)
{
    int *mb_flags = (int *)xcalloc((size_t)n_nodes, sizeof(int));
    // Add parents of target
    SEXP pa_s = VECTOR_ELT(parents_list, target);
    int n_pa = LENGTH(pa_s);
    if (n_pa > 0)
    {
        const int *pa_idx = INTEGER(pa_s);
        for (int i = 0; i < n_pa; ++i)
            if (pa_idx[i] != target)
                mb_flags[pa_idx[i]] = 1;
    }
    // Add children of target
    SEXP ch_s = VECTOR_ELT(children_list, target);
    int n_ch = LENGTH(ch_s);
    if (n_ch > 0)
    {
        const int *ch_idx = INTEGER(ch_s);
        for (int i = 0; i < n_ch; ++i)
            if (ch_idx[i] != target)
                mb_flags[ch_idx[i]] = 1;
        // Add co-parents of children (spouses)
        for (int i = 0; i < n_ch; ++i)
        {
            int child = ch_idx[i];
            SEXP child_pa_s = VECTOR_ELT(parents_list, child);
            int n_child_pa = LENGTH(child_pa_s);
            if (n_child_pa > 0)
            {
                const int *child_pa_idx = INTEGER(child_pa_s);
                for (int j = 0; j < n_child_pa; ++j)
                    if (child_pa_idx[j] != target)
                        mb_flags[child_pa_idx[j]] = 1;
            }
        }
    }
    // Count and collect indices
    int mb_size = 0;
    for (int i = 0; i < n_nodes; ++i)
        if (mb_flags[i] && i != target)
            mb_size++;
    int *mb = (int *)xmalloc(sizeof(int) * (size_t)mb_size);
    int k = 0;
    for (int i = 0; i < n_nodes; ++i)
        if (mb_flags[i] && i != target)
            mb[k++] = i;
    free(mb_flags);
    *mb_size_out = mb_size;
    return mb;
}

/* Returns the maximum cardinality among the Markov blanket of target */
static int max_card_mb(const int *ns, int n_nodes, SEXP parents_list, SEXP children_list, int target)
{
    int mb_size = 0;
    int *mb = markov_blanket(n_nodes, parents_list, children_list, target, &mb_size);
    int mx = 1;
    for (int i = 0; i < mb_size; ++i)
        if (ns[mb[i]] > mx)
            mx = ns[mb[i]];
    free(mb);
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

/* Optimized structure for sparse interval computation */
typedef struct
{
    int *uniq_vals;
    int K;
    int **prefix_counts; /* K x (n+1) */
    int n;
} CategoryData;

static CategoryData *create_category_data(const int *col_sorted, int n)
{
    CategoryData *cd = (CategoryData *)xmalloc(sizeof(CategoryData));

    int *uniq = NULL, K = 0, *id = NULL;
    compress_to_ids(col_sorted, n, &uniq, &K, &id);

    cd->uniq_vals = uniq;
    cd->K = K;
    cd->n = n;

    /* prefix counts per class */
    cd->prefix_counts = (int **)xmalloc(sizeof(int *) * (size_t)K);
    for (int k = 0; k < K; ++k)
        cd->prefix_counts[k] = (int *)xcalloc((size_t)n + 1, sizeof(int));

    for (int i = 0; i < n; ++i)
    {
        for (int k = 0; k < K; ++k)
            cd->prefix_counts[k][i + 1] = cd->prefix_counts[k][i];
        cd->prefix_counts[id[i]][i + 1] += 1;
    }

    free(id);
    return cd;
}

static void free_category_data(CategoryData *cd)
{
    if (!cd)
        return;
    if (cd->uniq_vals)
        free(cd->uniq_vals);
    if (cd->prefix_counts)
    {
        for (int k = 0; k < cd->K; ++k)
            if (cd->prefix_counts[k])
                free(cd->prefix_counts[k]);
        free(cd->prefix_counts);
    }
    free(cd);
}

/* Compute single interval score on-demand instead of full n×n table */
static double compute_interval_score_single_cat(CategoryData *cd, int i, int j)
{
    if (j < i)
        return INFINITY;

    int tot = j - i + 1;

    /* Skip very large intervals for efficiency - but make threshold much higher for small datasets */
    if (tot > 50000)
        return INFINITY;

    double v = lgammafn((double)tot + 1.0);

    for (int k = 0; k < cd->K; ++k)
    {
        int cnt = cd->prefix_counts[k][j + 1] - cd->prefix_counts[k][i];
        if (cnt > 0)
        { /* Only compute lgamma for non-zero counts */
            v -= lgammafn((double)cnt + 1.0);
        }
    }

    return v;
}

/* Legacy function for compatibility - now uses sparse computation */
static double *interval_table_single_cat(const int *col_sorted, int n)
{
    CategoryData *cd = create_category_data(col_sorted, n);

    double *tbl = (double *)xmalloc(sizeof(double) * (size_t)n * (size_t)n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            tbl[(size_t)i * (size_t)n + (size_t)j] = compute_interval_score_single_cat(cd, i, j);
        }
    }

    free_category_data(cd);
    return tbl;
}

/* Optimized child-spouse data structure */
typedef struct
{
    int C, S;             /* number of child and spouse categories */
    int ***prefix_counts; /* S x C x (n+1) */
    int n;
} ChildSpouseData;

static ChildSpouseData *create_child_spouse_data(const int *child_sorted, const int *sp_sorted, int n)
{
    ChildSpouseData *csd = (ChildSpouseData *)xmalloc(sizeof(ChildSpouseData));

    int *uniqC = NULL, *cid = NULL;
    compress_to_ids(child_sorted, n, &uniqC, &csd->C, &cid);
    free(uniqC);

    int *uniqS = NULL, *sid = NULL;
    compress_to_ids(sp_sorted, n, &uniqS, &csd->S, &sid);
    free(uniqS);

    csd->n = n;

    /* prefix counts pref[S][C][n+1] */
    csd->prefix_counts = (int ***)xmalloc(sizeof(int **) * (size_t)csd->S);
    for (int s = 0; s < csd->S; ++s)
    {
        csd->prefix_counts[s] = (int **)xmalloc(sizeof(int *) * (size_t)csd->C);
        for (int c = 0; c < csd->C; ++c)
            csd->prefix_counts[s][c] = (int *)xcalloc((size_t)n + 1, sizeof(int));
    }

    for (int i = 0; i < n; ++i)
    {
        for (int s = 0; s < csd->S; ++s)
            for (int c = 0; c < csd->C; ++c)
                csd->prefix_counts[s][c][i + 1] = csd->prefix_counts[s][c][i];
        csd->prefix_counts[sid[i]][cid[i]][i + 1] += 1;
    }

    free(cid);
    free(sid);
    return csd;
}

static void free_child_spouse_data(ChildSpouseData *csd)
{
    if (!csd)
        return;
    if (csd->prefix_counts)
    {
        for (int s = 0; s < csd->S; ++s)
        {
            if (csd->prefix_counts[s])
            {
                for (int c = 0; c < csd->C; ++c)
                    if (csd->prefix_counts[s][c])
                        free(csd->prefix_counts[s][c]);
                free(csd->prefix_counts[s]);
            }
        }
        free(csd->prefix_counts);
    }
    free(csd);
}

static double compute_interval_score_child_spouse(ChildSpouseData *csd, int i, int j)
{
    if (j < i)
        return INFINITY;

    /* Skip very large intervals for efficiency - but make threshold much higher for small datasets */
    if (j - i + 1 > 50000)
        return INFINITY;

    double acc = 0.0;
    for (int s = 0; s < csd->S; ++s)
    {
        int tot = 0;
        for (int c = 0; c < csd->C; ++c)
            tot += (csd->prefix_counts[s][c][j + 1] - csd->prefix_counts[s][c][i]);
        if (tot == 0)
            continue;

        /* stars&bars */
        acc += lgammafn((double)tot + (double)csd->C) - lgammafn((double)csd->C) - lgammafn((double)tot + 1.0);
        /* multinomial over child categories given spouse */
        double v = lgammafn((double)tot + 1.0);
        for (int c = 0; c < csd->C; ++c)
        {
            int cnt = csd->prefix_counts[s][c][j + 1] - csd->prefix_counts[s][c][i];
            if (cnt > 0)
            { /* Only compute lgamma for non-zero counts */
                v -= lgammafn((double)cnt + 1.0);
            }
        }
        acc += v;
    }
    return acc;
}

/* Legacy function for compatibility */
static double *interval_table_child_spouse(const int *child_sorted, const int *sp_sorted, int n)
{
    ChildSpouseData *csd = create_child_spouse_data(child_sorted, sp_sorted, n);

    double *tbl = (double *)xmalloc(sizeof(double) * (size_t)n * (size_t)n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            tbl[(size_t)i * (size_t)n + (size_t)j] = compute_interval_score_child_spouse(csd, i, j);
        }
    }

    free_child_spouse_data(csd);
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

/* Simple hash table for memoization */
typedef struct CacheEntry
{
    int i, j;
    double score;
    struct CacheEntry *next;
} CacheEntry;

typedef struct
{
    CacheEntry **buckets;
    int size;
    int hits, misses;
} ScoreCache;

/* Optimized sparse scoring structure that avoids O(n²) memory */
typedef struct
{
    CategoryData **parent_data;      /* for approximation mode */
    CategoryData *joint_parent_data; /* for exact mode */
    ChildSpouseData **child_spouse_data;
    int n_pa, n_ch;
    int target, n_cases;
    const int *order;
    const int *ns;
    const int *pa_idx;
    int parent_size;
    int approx_parents;
    ScoreCache *cache; /* memoization cache */
} SparseScorer;

/* Cache management functions */
static ScoreCache *create_cache(int cache_size)
{
    ScoreCache *cache = (ScoreCache *)xmalloc(sizeof(ScoreCache));
    cache->size = cache_size;
    cache->buckets = (CacheEntry **)xcalloc(cache_size, sizeof(CacheEntry *));
    cache->hits = 0;
    cache->misses = 0;
    return cache;
}

static void free_cache(ScoreCache *cache)
{
    if (!cache)
        return;
    for (int i = 0; i < cache->size; i++)
    {
        CacheEntry *entry = cache->buckets[i];
        while (entry)
        {
            CacheEntry *next = entry->next;
            free(entry);
            entry = next;
        }
    }
    free(cache->buckets);
    free(cache);
}

static inline int hash_pair(int i, int j, int size)
{
    return ((unsigned int)(i * 1009 + j * 1013)) % size;
}

static double cache_get(ScoreCache *cache, int i, int j, int *found)
{
    int hash = hash_pair(i, j, cache->size);
    CacheEntry *entry = cache->buckets[hash];
    while (entry)
    {
        if (entry->i == i && entry->j == j)
        {
            cache->hits++;
            *found = 1;
            return entry->score;
        }
        entry = entry->next;
    }
    cache->misses++;
    *found = 0;
    return 0.0;
}

static void cache_put(ScoreCache *cache, int i, int j, double score)
{
    int hash = hash_pair(i, j, cache->size);
    CacheEntry *entry = (CacheEntry *)xmalloc(sizeof(CacheEntry));
    entry->i = i;
    entry->j = j;
    entry->score = score;
    entry->next = cache->buckets[hash];
    cache->buckets[hash] = entry;
}

static SparseScorer *create_sparse_scorer(const double *X, const int *D,
                                          int n_nodes, int n_cases, int target,
                                          const int *order, const int *ns,
                                          const int *pa_idx, int n_pa,
                                          const int *children, int n_ch,
                                          SEXP parents_list, int approx_parents)
{
    (void)n_nodes;
    SparseScorer *ss = (SparseScorer *)xmalloc(sizeof(SparseScorer));
    ss->n_pa = n_pa;
    ss->n_ch = n_ch;
    ss->target = target;
    ss->n_cases = n_cases;
    ss->order = order;
    ss->ns = ns;
    ss->pa_idx = pa_idx;
    ss->approx_parents = approx_parents;

    /* Initialize parent data */
    ss->parent_data = NULL;
    ss->joint_parent_data = NULL;

    if (n_pa > 0)
    {
        if (approx_parents)
        {
            ss->parent_data = (CategoryData **)xmalloc(sizeof(CategoryData *) * (size_t)n_pa);
            for (int pi = 0; pi < n_pa; ++pi)
            {
                int pcol = pa_idx[pi];
                int *p_sorted = (int *)xmalloc(sizeof(int) * (size_t)n_cases);
                for (int t = 0; t < n_cases; ++t)
                    p_sorted[t] = get_disc(D, n_cases, pcol, order[t]);
                ss->parent_data[pi] = create_category_data(p_sorted, n_cases);
                free(p_sorted);
            }
        }
        else
        {
            int *joint_sorted = (int *)xmalloc(sizeof(int) * (size_t)n_cases);
            build_joint_parents_sorted(D, n_cases, order, pa_idx, n_pa, joint_sorted);
            ss->joint_parent_data = create_category_data(joint_sorted, n_cases);
            free(joint_sorted);
        }
    }

    /* Initialize child-spouse data */
    ss->child_spouse_data = NULL;
    if (n_ch > 0)
    {
        ss->child_spouse_data = (ChildSpouseData **)xmalloc(sizeof(ChildSpouseData *) * (size_t)n_ch);
        for (int c = 0; c < n_ch; ++c)
        {
            int chcol = children[c];
            int *child_sorted = (int *)xmalloc(sizeof(int) * (size_t)n_cases);
            for (int t = 0; t < n_cases; ++t)
                child_sorted[t] = get_disc(D, n_cases, chcol, order[t]);

            /* spouses = parents(child) \ {target} */
            SEXP par_ch = VECTOR_ELT(parents_list, chcol);
            int nsp_all = LENGTH(par_ch);
            const int *par_ch_idx = (nsp_all > 0) ? INTEGER(par_ch) : NULL;

            int *sp_sorted = (int *)xmalloc(sizeof(int) * (size_t)n_cases);
            if (nsp_all == 0)
            {
                for (int t = 0; t < n_cases; ++t)
                    sp_sorted[t] = 0;
            }
            else
            {
                for (int t = 0; t < n_cases; ++t)
                {
                    int r = order[t];
                    unsigned int h = 2166136261u;
                    int used = 0;
                    for (int s = 0; s < nsp_all; ++s)
                    {
                        int sidx = par_ch_idx[s];
                        if (sidx == target)
                            continue;
                        int v = get_disc(D, n_cases, sidx, r);
                        h ^= (unsigned int)v;
                        h *= 16777619u;
                        used = 1;
                    }
                    sp_sorted[t] = used ? (int)(h & 0x7fffffff) : 0;
                }
            }

            ss->child_spouse_data[c] = create_child_spouse_data(child_sorted, sp_sorted, n_cases);
            free(child_sorted);
            free(sp_sorted);
        }
    }

    /* Compute parent size */
    ss->parent_size = 1;
    if (n_pa > 0)
    {
        for (int pi = 0; pi < n_pa; ++pi)
        {
            ss->parent_size *= ns[pa_idx[pi]];
        }
    }

    /* Initialize memoization cache - size based on expected interval queries */
    int cache_size = (n_cases < 1000) ? 1009 : (n_cases < 10000) ? 10007
                                                                 : 50021;
    ss->cache = create_cache(cache_size);

    return ss;
}

static void free_sparse_scorer(SparseScorer *ss)
{
    if (!ss)
        return;

    if (ss->parent_data)
    {
        for (int pi = 0; pi < ss->n_pa; ++pi)
        {
            free_category_data(ss->parent_data[pi]);
        }
        free(ss->parent_data);
    }

    if (ss->joint_parent_data)
    {
        free_category_data(ss->joint_parent_data);
    }

    if (ss->child_spouse_data)
    {
        for (int c = 0; c < ss->n_ch; ++c)
        {
            free_child_spouse_data(ss->child_spouse_data[c]);
        }
        free(ss->child_spouse_data);
    }

    /* Print cache statistics for debugging */
    if (ss->cache)
    {
        int total = ss->cache->hits + ss->cache->misses;
        if (total > 0)
        {
            double hit_rate = (double)ss->cache->hits / total * 100.0;
            printf("Cache stats: %d hits, %d misses (%.1f%% hit rate)\n",
                   ss->cache->hits, ss->cache->misses, hit_rate);
        }
        free_cache(ss->cache);
    }

    free(ss);
}

/* Compute interval score on-demand without storing full matrix */
static double compute_sparse_score(SparseScorer *ss, int i, int j)
{
    if (j < i)
        return INFINITY;

    /* Check cache first */
    int found = 0;
    double cached_score = cache_get(ss->cache, i, j, &found);
    if (found)
    {
        return cached_score;
    }

    double score = 0.0;

    /* Parents contribution */
    if (ss->n_pa > 0)
    {
        if (ss->approx_parents)
        {
            for (int pi = 0; pi < ss->n_pa; ++pi)
            {
                score += compute_interval_score_single_cat(ss->parent_data[pi], i, j);
            }
        }
        else
        {
            score += compute_interval_score_single_cat(ss->joint_parent_data, i, j);
        }
    }

    /* Children + spouses contribution */
    for (int c = 0; c < ss->n_ch; ++c)
    {
        score += compute_interval_score_child_spouse(ss->child_spouse_data[c], i, j);
    }

    /* Prior term */
    score += lgammafn(j - i + (double)ss->parent_size + 1.0) - lgammafn(j - i + 2.0) - lgammafn((double)ss->parent_size);

    /* Cache the result for future use */
    cache_put(ss->cache, i, j, score);

    return score;
}

/* Legacy function for compatibility - now redirects to sparse implementation */
static double *build_Hfull(const double *X, const int *D,
                           int n_nodes, int n_cases, int target,
                           const int *order,
                           const int *ns,
                           const int *pa_idx, int n_pa,
                           const int *children, int n_ch,
                           SEXP parents_list, int approx_parents)
{
    double *H = (double *)xcalloc((size_t)n_cases * (size_t)n_cases, sizeof(double));
    SparseScorer *ss = create_sparse_scorer(X, D, n_nodes, n_cases, target, order, ns,
                                            pa_idx, n_pa, children, n_ch,
                                            parents_list, approx_parents);

    for (int i = 0; i < n_cases; ++i)
    {
        for (int j = i; j < n_cases; ++j)
        {
            H[(size_t)i * (size_t)n_cases + (size_t)j] = compute_sparse_score(ss, i, j);
        }
    }

    free_sparse_scorer(ss);
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

/* Optimized DP that works directly with sparse scorer with early pruning */
static void dp_and_edges_sparse(SparseScorer *ss, int m,
                                const double *X, int cont_col, const int *order, const int *tails,
                                const double *split, int lambda,
                                double **cuts_out, int *n_cuts_out)
{
    int n_cases = ss->n_cases;

    int *heads = (int *)xmalloc(sizeof(int) * (size_t)m);
    int hcnt = 0;
    heads[hcnt++] = 0;
    for (int t = 1; t < n_cases && hcnt < m; ++t)
    {
        double prev = get_cont(X, n_cases, cont_col, order[t - 1]);
        double cur = get_cont(X, n_cases, cont_col, order[t]);
        if (cur != prev)
            heads[hcnt++] = t;
    }

    double xmin = get_cont(X, n_cases, cont_col, order[0]);
    double xmax = get_cont(X, n_cases, cont_col, order[n_cases - 1]);
    double span = (xmax > xmin) ? (xmax - xmin) : 1.0;

    double *best = (double *)xmalloc(sizeof(double) * (size_t)m);
    int **choice = (int **)xcalloc((size_t)m, sizeof(int *));
    int *clen = (int *)xcalloc((size_t)m, sizeof(int));

    for (int a = 0; a < m; ++a)
    {
        int tail_a = tails[a];
        if (a == 0)
        {
            double score = compute_sparse_score(ss, 0, tail_a);
            best[0] = -log(split[tail_a]) + score;
            choice[0] = (int *)xmalloc(sizeof(int));
            choice[0][0] = tail_a;
            clen[0] = 1;
        }
        else
        {
            double bv = INFINITY;
            int *bs = NULL;
            int bl = 0;

            /* Parallel computation of candidate values for each breakpoint */
            int num_candidates = a + 1;
            double *candidate_vals = (double *)xmalloc(sizeof(double) * (size_t)num_candidates);
            double **candidate_choices = (double **)xmalloc(sizeof(double *) * (size_t)num_candidates);
            int *candidate_lens = (int *)xmalloc(sizeof(int) * (size_t)num_candidates);

            /* Initialize arrays */
            for (int b = 0; b <= a; ++b)
            {
                candidate_vals[b] = INFINITY;
                candidate_choices[b] = NULL;
                candidate_lens[b] = 0;
            }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for (int b = 0; b <= a; ++b)
            {
                double val;
                if (b == a)
                {
                    double frac = (get_cont(X, n_cases, cont_col, order[tail_a]) - xmin) / span;
                    double score = compute_sparse_score(ss, 0, tail_a);
                    val = frac * (double)lambda - log(split[tail_a]) + score;

                    /* Store single-element choice */
                    candidate_choices[b] = (double *)xmalloc(sizeof(double));
                    candidate_choices[b][0] = (double)tail_a;
                    candidate_lens[b] = 1;
                }
                else
                {
                    int head_next = heads[b + 1];
                    double frac = (get_cont(X, n_cases, cont_col, order[tail_a]) -
                                   get_cont(X, n_cases, cont_col, order[head_next])) /
                                  span;

                    double score = compute_sparse_score(ss, head_next, tail_a);
                    val = best[b] + frac * (double)lambda - log(split[tail_a]) + score;

                    /* Store extended choice */
                    candidate_lens[b] = clen[b] + 1;
                    candidate_choices[b] = (double *)xmalloc(sizeof(double) * (size_t)candidate_lens[b]);
                    for (int k = 0; k < clen[b]; k++)
                    {
                        candidate_choices[b][k] = (double)choice[b][k];
                    }
                    candidate_choices[b][candidate_lens[b] - 1] = (double)tail_a;
                }
                candidate_vals[b] = val;
            }

            /* Sequential reduction to find minimum */
            int best_b = -1;
            for (int b = 0; b <= a; ++b)
            {
                if (candidate_vals[b] < bv)
                {
                    bv = candidate_vals[b];
                    best_b = b;
                }
            }

            /* Copy best choice */
            if (best_b >= 0)
            {
                bl = candidate_lens[best_b];
                bs = (int *)xmalloc(sizeof(int) * (size_t)bl);
                for (int k = 0; k < bl; k++)
                {
                    bs[k] = (int)candidate_choices[best_b][k];
                }
            }

            /* Cleanup candidate arrays */
            for (int b = 0; b <= a; ++b)
            {
                if (candidate_choices[b])
                {
                    free(candidate_choices[b]);
                }
            }
            free(candidate_vals);
            free(candidate_choices);
            free(candidate_lens);

            best[a] = bv;
            choice[a] = bs;
            clen[a] = bl;
        }
    }

    int K = clen[m - 1];
    int ncuts = (K > 1) ? (K - 1) : 0;
    double *cuts = NULL;
    if (ncuts > 0)
    {
        cuts = (double *)xmalloc(sizeof(double) * (size_t)ncuts);
        for (int i = 0; i < ncuts; ++i)
        {
            int t = choice[m - 1][i];
            double a_val = get_cont(X, n_cases, cont_col, order[t]);
            double b_val = get_cont(X, n_cases, cont_col, order[t + 1]);
            cuts[i] = 0.5 * (a_val + b_val);
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
                                  int approx_parents, int lambda,
                                  double **cuts_out, int *n_cuts_out)
{
    *cuts_out = NULL;
    *n_cuts_out = 0;

    int *order = argsort_col(X, n_cases, cont_col);
    int m = 0;
    int *tails = build_tails(X, n_cases, cont_col, order, &m);
    if (m == 0)
    {
        free(order);
        free(tails);
        return;
    }

    SEXP pa_s = VECTOR_ELT(parents_list, target);
    int n_pa = LENGTH(pa_s);
    const int *pa_idx = (n_pa > 0) ? INTEGER(pa_s) : NULL;

    SEXP ch_s = VECTOR_ELT(children_list, target);
    int n_ch = LENGTH(ch_s);
    const int *children = (n_ch > 0) ? INTEGER(ch_s) : NULL;

    /* Use optimized sparse scorer instead of full matrices */
    SparseScorer *ss = create_sparse_scorer(X, D, n_nodes, n_cases, target, order, ns,
                                            pa_idx, n_pa, children, n_ch,
                                            parents_list, approx_parents);

    if (lambda <= 0)
        lambda = max_card(ns, n_nodes);
    double *split = compute_split_prior(X, n_cases, cont_col, order, lambda);

    double *cuts = NULL;
    int ncuts = 0;
    /* Use optimized sparse DP algorithm */
    dp_and_edges_sparse(ss, m, X, cont_col, order, tails, split, lambda, &cuts, &ncuts);

    free_sparse_scorer(ss);
    free(split);
    free(order);
    free(tails);
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
                                  SEXP approx_parents, SEXP l_card)
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
    if (!isInteger(l_card))
        error("l_card must be integer scalar.");
    const double *X = REAL(data_cont);
    const int *D = INTEGER(data_disc);
    int n_nodes_c = INTEGER(n_nodes)[0];
    int n_cases_c = INTEGER(n_cases)[0];
    const int *ns_c = INTEGER(ns);
    int cont_col_c = INTEGER(cont_col)[0];
    int target_c = INTEGER(target)[0];
    int approx_c = INTEGER(approx_parents)[0];
    int lambda = INTEGER(l_card)[0];

    double *cuts = NULL;
    int ncuts = 0;
    algo1_discretize_core(X, D, n_nodes_c, n_cases_c, ns_c,
                          cont_col_c, target_c,
                          parents_list, children_list,
                          approx_c, lambda,
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

/* Optimized discretize_column with early termination for identical values */
static void discretize_column(const double *X, int *D, int n_cases, int col_x, int col_d,
                              const double *cuts, int ncuts)
{
    if (ncuts == 0)
    {
        /* No cuts - all values go to bin 0 */
        for (int i = 0; i < n_cases; ++i)
        {
            D[i + (size_t)col_d * (size_t)n_cases] = 0;
        }
        return;
    }

    for (int i = 0; i < n_cases; ++i)
    {
        double val = get_cont(X, n_cases, col_x, i);

        /* Binary search for efficiency with many cuts */
        int bin = 0;
        if (ncuts > 8)
        { /* Use binary search for many cuts */
            int low = 0, high = ncuts - 1;
            while (low <= high)
            {
                int mid = (low + high) / 2;
                if (val <= cuts[mid])
                {
                    high = mid - 1;
                }
                else
                {
                    low = mid + 1;
                    bin = low;
                }
            }
        }
        else
        { /* Linear search for few cuts */
            while (bin < ncuts && val > cuts[bin])
                bin++;
        }

        D[i + (size_t)col_d * (size_t)n_cases] = bin;
    }
}
/* Algorithm 2 */
SEXP bnstruct_dvbn_discretize_all(SEXP data_cont, SEXP data_disc,
                                  SEXP n_nodes, SEXP n_cases,
                                  SEXP n_cont, SEXP cont_index,
                                  SEXP ns, SEXP parents_list, SEXP children_list,
                                  SEXP n_cycles, SEXP approx_parents, SEXP l_card)
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
    if (!isInteger(l_card))
        error("l_card must be integer scalar.");

    const double *X = REAL(data_cont);
    int *D = INTEGER(data_disc);
    int n_nodes_c = INTEGER(n_nodes)[0];
    int n_cases_c = INTEGER(n_cases)[0];
    int n_cont_c = INTEGER(n_cont)[0];
    const int *cont_idx = INTEGER(cont_index);
    int *ns_c = INTEGER(ns);
    int max_cycles = INTEGER(n_cycles)[0];
    int approx_c = INTEGER(approx_parents)[0];
    int lambda = INTEGER(l_card)[0];

    double **cuts = (double **)xcalloc((size_t)n_cont_c, sizeof(double *));
    int *ncuts = (int *)xcalloc((size_t)n_cont_c, sizeof(int));

    int changed = 1;
    int iter = 0;
    int consecutive_no_change = 0; /* Early convergence detection */

    while (changed && iter < max_cycles)
    {
        changed = 0;
        int vars_changed_this_iter = 0;

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
                                      approx_c, lambda,
                                      &new_cuts, &new_n);
                // Debug print for iteration tracking
                printf("iter %d col_d=%d ns=%d lambda=%d ncuts=%d\n",
                       iter, col_d, ns_c[col_d], lambda, new_n);
            }

            if (!cuts_equal(cuts[j], ncuts[j], new_cuts, new_n))
            {
                if (cuts[j])
                    free(cuts[j]);
                cuts[j] = new_cuts;
                ncuts[j] = new_n;
                changed = 1;
                vars_changed_this_iter++;

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

        /* Early convergence detection */
        if (vars_changed_this_iter == 0)
        {
            consecutive_no_change++;
            if (consecutive_no_change >= 2)
            {
                break; /* Converged early */
            }
        }
        else
        {
            consecutive_no_change = 0;
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
