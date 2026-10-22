#include <random>
#include <chrono>

#include <Rcpp.h>

// OpenMP parallel regions (guarded: builds single-threaded without OpenMP).
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

 // Use R random number generator to control seeds
 // Otherwise use the following:
 // static std::random_device rd;
 // initialize Mersennes' twister using rd to generate the seed
 // static std::mt19937 gen{rd()};
 // std::uniform_real_distribution<double> dist(0, 1);

 // calls an arbitrary R function `f` on a numeric vector `x` and returns its result.
 NumericVector callFunction(NumericVector x, Function f) {
   NumericVector res = f(x);
   return res;
 }

 // Inverse-CDF categorical draw: index in [0, n) of the first cumulative
 // weight exceeding `u * total`, with `u` a pre-drawn Uniform(0,1). Zero-weight
 // entries are never selected (unless every weight is zero, in which case the
 // last index is returned as a numerical fallback). Consuming pre-drawn
 // uniforms instead of calling R's sampler keeps the kernels thread-safe and
 // `set.seed()`-deterministic.
 static int sampleCdf(const double *w, int n, double u)
 {
   double total = 0.0;
   for (int i = 0; i < n; i++) total += w[i];
   double threshold = u * total;
   double acc = 0.0;
   for (int i = 0; i < n; i++)
   {
     acc += w[i];
     if (threshold < acc) return i;
   }
   return n - 1;
 }

 //' indexPatterns
 //'
 //' @param U IntegerVector with factor values corresponding to the patterns of Partially Identifying Variables (PIVs) observed among records in the concerned source
 //' @param nvals integer for the total number of possible patterns (among all sources)
 //'
 //' @return List: for each pattern in key, a vector of indices of the records having the pattern in the concerned source
 //' @export
 // [[Rcpp::export]]
 List indexPatterns(IntegerVector U, int nvals)
 {
   // Create one empty "bucket" (integer vector) for each pattern
   List out(nvals);
   for (int i = 0; i < nvals; i++)
   {
     IntegerVector tmp;
     out[i] = tmp;
   }
   // Which value does observation i have
   // Assign it to the appropriate "bucket" in the list
   for (int i = 0; i < U.length(); i++)
   {
     IntegerVector tmpX = out[U[i]-1];
     tmpX.push_back(i+1);
     out[U[i]-1] = tmpX;
   }
   return out;
 }

 //' pairPatterns
 //'
 //' @param A List with for each pattern in value, count of the records having the pattern in the concerned source
 //' @param B List with for each pattern in value, count of the records having the pattern in the concerned source
 //' @param nvals integer for the total number of possible patterns (among all sources)
 //'
 //' @return IntegerMatrix: nrow=nbr of potential links, ncol = 2, indicates the indices (of records from A, records from B) of records with matching patterns to be considered in the likelihood as potential links. It represents the cartesian product of both lists passed as parameters to represent all possible linked pairs of records.
 //' @export
 // [[Rcpp::export]]
 IntegerMatrix pairPatterns(List A, List B, int nvals)
 {
   // Figure out how many candidate pairs there will be in total
   int ntotal = 0;
   for (int k = 0; k < nvals; k++)
   {
     // Check a bucket to see whether and how many combos are possible
     IntegerVector tmpA = A[k];
     IntegerVector tmpB = B[k];
     ntotal += tmpA.length() * tmpB.length();
   }
   IntegerMatrix tmpC(ntotal,2);
   int counter = 0;
   // For each pattern (index in list A or B), give all possible pairs of indices carrying that pattern
   for (int k = 0; k < nvals; k++)
   {
     IntegerVector tmpA = A[k];
     IntegerVector tmpB = B[k];
     for (int i = 0; i < tmpA.length(); i++)
     {
       for (int j = 0; j < tmpB.length(); j++)
       {
         tmpC(counter,0) = tmpA[i];
         tmpC(counter,1) = tmpB[j];
         counter += 1;
       }
     }
   }
   return tmpC;
 }

 //' buildCandidatePairs
 //'
 //' Compute the candidate-pair matrix: all pairs of records from A and B whose
 //' stable-PIV true-value patterns are identical, i.e. the cartesian product of
 //' the A-records and B-records within every pattern group. This is a
 //' single-call, allocation-light replacement for the R-side string pipeline
 //' `pasteIntoPattern -> unique -> factor -> indexPatterns -> pairPatterns`
 //' (per-record ostringstream building plus factor() over millions of strings
 //' per Gibbs sweep). Patterns are compared as integer tuples via a
 //' lexicographic sort, so no hashing collisions or string churn.
 //'
 //' @param patsA IntegerMatrix (nA x m) of true stable-PIV values for file A (m may be 0: every record then shares the single empty pattern)
 //' @param patsB IntegerMatrix (nB x m) of true stable-PIV values for file B
 //'
 //' @return IntegerMatrix: 2 columns with 1-based indices into A and B; groups
 //' ordered lexicographically by pattern, records in ascending index order
 //' within a group (deterministic output order).
 //' @export
 // [[Rcpp::export]]
 IntegerMatrix buildCandidatePairs(IntegerMatrix patsA, IntegerMatrix patsB)
 {
   int nA = patsA.nrow(), nB = patsB.nrow(), m = patsA.ncol();
   if (patsB.ncol() != m)
     Rcpp::stop("`patsA` and `patsB` must have the same number of columns.");
   // Row-major copies for cache-friendly lexicographic comparison
   std::vector<int> flatA((size_t)nA * m), flatB((size_t)nB * m);
   for (int i = 0; i < nA; i++)
     for (int c = 0; c < m; c++)
       flatA[(size_t)i * m + c] = patsA(i, c);
   for (int j = 0; j < nB; j++)
     for (int c = 0; c < m; c++)
       flatB[(size_t)j * m + c] = patsB(j, c);
   // Sort record indices by (pattern tuple, original index): the index
   // tie-break keeps the output deterministic across std::sort implementations.
   std::vector<int> ordA(nA), ordB(nB);
   for (int i = 0; i < nA; i++) ordA[i] = i;
   for (int j = 0; j < nB; j++) ordB[j] = j;
   auto cmpA = [&](int a, int b) {
     for (int c = 0; c < m; c++)
     {
       int va = flatA[(size_t)a * m + c], vb = flatA[(size_t)b * m + c];
       if (va != vb) return va < vb;
     }
     return a < b;
   };
   auto cmpB = [&](int a, int b) {
     for (int c = 0; c < m; c++)
     {
       int va = flatB[(size_t)a * m + c], vb = flatB[(size_t)b * m + c];
       if (va != vb) return va < vb;
     }
     return a < b;
   };
   std::sort(ordA.begin(), ordA.end(), cmpA);
   std::sort(ordB.begin(), ordB.end(), cmpB);
   auto rowEqA = [&](int a, int b) {
     for (int c = 0; c < m; c++)
       if (flatA[(size_t)a * m + c] != flatA[(size_t)b * m + c]) return false;
     return true;
   };
   auto rowEqB = [&](int a, int b) {
     for (int c = 0; c < m; c++)
       if (flatB[(size_t)a * m + c] != flatB[(size_t)b * m + c]) return false;
     return true;
   };
   auto rowABcmp = [&](int a, int b) -> int {  // A row `a` vs B row `b`
     for (int c = 0; c < m; c++)
     {
       int va = flatA[(size_t)a * m + c], vb = flatB[(size_t)b * m + c];
       if (va != vb) return va < vb ? -1 : 1;
     }
     return 0;
   };
   // First pass: walk both sorted orders jointly, count the pairs per group
   long long ntotal = 0;
   int pa = 0, pb = 0;
   while (pa < nA && pb < nB)
   {
     int c = rowABcmp(ordA[pa], ordB[pb]);
     if (c < 0) { pa++; continue; }
     if (c > 0) { pb++; continue; }
     int ea = pa + 1, eb = pb + 1;
     while (ea < nA && rowEqA(ordA[pa], ordA[ea])) ea++;
     while (eb < nB && rowEqB(ordB[pb], ordB[eb])) eb++;
     ntotal += (long long)(ea - pa) * (long long)(eb - pb);
     pa = ea; pb = eb;
   }
   if (ntotal > 2147483647LL)
     Rcpp::stop("Too many candidate pairs (", ntotal,
                "): the stable-PIV patterns leave one giant group; more discriminative PIVs are needed.");
   // Second pass: emit the cartesian product of each group (A-major)
   IntegerMatrix out((int)ntotal, 2);
   int w = 0;
   pa = 0; pb = 0;
   while (pa < nA && pb < nB)
   {
     int c = rowABcmp(ordA[pa], ordB[pb]);
     if (c < 0) { pa++; continue; }
     if (c > 0) { pb++; continue; }
     int ea = pa + 1, eb = pb + 1;
     while (ea < nA && rowEqA(ordA[pa], ordA[ea])) ea++;
     while (eb < nB && rowEqB(ordB[pb], ordB[eb])) eb++;
     for (int ii = pa; ii < ea; ii++)
       for (int jj = pb; jj < eb; jj++)
       {
         out(w, 0) = ordA[ii] + 1;
         out(w, 1) = ordB[jj] + 1;
         w++;
       }
     pa = ea; pb = eb;
   }
   return out;
 }

 //' pasteIntoPattern
 //'
 //' @param A IntegerMatrix with values to form patterns
 //'
 //' @return CharacterVector: "pasteIntoPattern(A)" is a faster version in Rcpp of "do.call(paste, c(as.data.frame(A), list(sep="_")))"
 //' @export
 // [[Rcpp::export]]
 CharacterVector pasteIntoPattern(IntegerMatrix A)
 {
   int i = 0, j = 0;
   int sz = A.nrow();
   CharacterVector res(sz);
   // Guard against datasets with no stable PIVs: a 0-column matrix has no
   // elements to read, and A(i,0) below would be an out-of-bounds access.
   // All records then share a single pattern (blocking degenerates to the
   // full Cartesian product of candidate pairs).
   if (A.ncol() == 0)
   {
     for (i = 0; i < sz; i++) res[i] = "";
     return res;
   }
   for (std::ostringstream oss; i < sz; i++, oss.str(""))
   {
     oss << A(i,0);
     for (j = 1; j < A.ncol(); j++)
     {
       oss << "_" << A(i,j);
     }
     res[i] = oss.str();
   }
   return res;
 }

 //' pairLogLik
 //'
 //' Compute the complete-data log-likelihood pieces for the D-sweep in a single
 //' pass over the records and candidate pairs (the C++ implementation of
 //' simulateD's likelihood block):
 //' - LLA[i]: log-likelihood contribution of A-record i if it is NOT linked,
 //' - LLB[j]: idem for B-record j,
 //' - LLLq[q]: contribution of candidate pair q if it IS linked, aligned with
 //'   the rows of `select`.
 //'
 //' Per PIV k each contribution is
 //'   log(eta_k[h]) + log P(registered g | true h)
 //' with P(g|h) = pMissing if g == 0, pAgree if g == h, pTypo otherwise, plus,
 //' for an unstable PIV, the survival term
 //'   log( pSameH ) if hA == hB, else log( (1 - pSameH)/(nval - 1) ),
 //' where pSameH = S(cov_k, alpha_k, times) is evaluated once per PIV, so
 //' user-defined survival models keep working.
 //'
 //' @param select IntegerMatrix (nPairs x 2, 1-based record indices) of candidate pairs
 //' @param truepivsA,truepivsB IntegerMatrix (n x nPIVs) of latent true values
 //' @param encodedA,encodedB List: the encoded data frames (columns looked up by name)
 //' @param PIVs CharacterVector of PIV column names (k order)
 //' @param PIVs_stable LogicalVector: TRUE where the PIV does not change over time
 //' @param n_values IntegerVector of the number of unique values per PIV
 //' @param eta,phi,alpha Lists of per-PIV numeric vectors (phi = [agreeA, agreeB, missA, missB])
 //' @param covA_names,covB_names Lists of per-PIV CharacterVectors with hazard-covariate column names (may be empty)
 //' @param S R function model_dynamics$S(cov, alpha, times), called once per unstable PIV
 //' @param threads Number of OpenMP threads for the per-pair loops (default 1 =
 //'   sequential). Ignored — the code runs sequentially — if the package was
 //'   built without OpenMP. Results are identical for any thread count: every
 //'   pair writes its own slot, no reduction crosses threads.
 //' @return List with elements LLA (nA), LLB (nB) and LLLq (nPairs)
 //' @export
 // [[Rcpp::export]]
 List pairLogLik(const IntegerMatrix & select,
                 const IntegerMatrix & truepivsA,
                 const IntegerMatrix & truepivsB,
                 List encodedA,
                 List encodedB,
                 CharacterVector PIVs,
                 LogicalVector PIVs_stable,
                 IntegerVector n_values,
                 List eta,
                 List phi,
                 List alpha,
                 List covA_names,
                 List covB_names,
                 Function S,
                 int threads = 1)
 {
   int nA = truepivsA.nrow(), nB = truepivsB.nrow();
   int nPairs = select.nrow();
   int nPIVs = n_values.length();
   if (threads < 1) threads = 1;
#ifndef _OPENMP
   (void)threads;  // built without OpenMP: silently sequential
#endif
   NumericVector LLA(nA, 0.0);
   NumericVector LLB(nB, 0.0);
   NumericVector LLLq(nPairs, 0.0);
   // Raw pointers for the hot loops (Rcpp matrices are column-major, like R)
   const int *sel = select.begin();     // sel[q] = A index, sel[nPairs + q] = B index
   const int *tpA = truepivsA.begin();  // tpA[i + k*nA]
   const int *tpB = truepivsB.begin();  // tpB[j + k*nB]
   double *lla = LLA.begin(), *llb = LLB.begin(), *lllq = LLLq.begin();
   for (int k = 0; k < nPIVs; k++)
   {
     NumericVector eta_k = eta[k];
     NumericVector phi_k = phi[k];
     int nval = n_values[k];
     double pMissingA = phi_k[2], pMissingB = phi_k[3];
     double pTypoA = (1-pMissingA) * (1-phi_k[0]) / (nval-1);
     double pAgreeA = (1-pMissingA) * phi_k[0];
     double pTypoB = (1-pMissingB) * (1-phi_k[1]) / (nval-1);
     double pAgreeB = (1-pMissingB) * phi_k[1];
     // log(eta) and the log error probabilities, once per PIV
     std::vector<double> logeta(nval);
     for (int h = 0; h < nval; h++) logeta[h] = std::log(eta_k[h]);
     double lMissA = std::log(pMissingA), lTypoA = std::log(pTypoA), lAgreeA = std::log(pAgreeA);
     double lMissB = std::log(pMissingB), lTypoB = std::log(pTypoB), lAgreeB = std::log(pAgreeB);
     // Registered-value columns for this PIV
     std::string nmA = Rcpp::as<std::string>(PIVs[k]), nmB = nmA;
     NumericVector colA = Rcpp::as<NumericVector>(encodedA[nmA]);
     NumericVector colB = Rcpp::as<NumericVector>(encodedB[nmB]);
     const double *gA = colA.begin(), *gB = colB.begin();
     // ---- Unlinked-record contributions (LLA, LLB) ----
     for (int i = 0; i < nA; i++)
     {
       int gi = (int)gA[i], hi = tpA[i + (size_t)k * nA];
       double lc = (gi == 0) ? lMissA : ((gi == hi) ? lAgreeA : lTypoA);
       lla[i] = lla[i] + logeta[hi - 1] + lc;
     }
     for (int j = 0; j < nB; j++)
     {
       int gj = (int)gB[j], hj = tpB[j + (size_t)k * nB];
       double lc = (gj == 0) ? lMissB : ((gj == hj) ? lAgreeB : lTypoB);
       llb[j] = llb[j] + logeta[hj - 1] + lc;
     }
     // ---- Candidate-pair contributions (LLLq) ----
     if (PIVs_stable[k])
     {
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
       for (int q = 0; q < nPairs; q++)
       {
         int ia = sel[q] - 1, ib = sel[nPairs + q] - 1;
         int ha = tpA[ia + (size_t)k * nA], hb = tpB[ib + (size_t)k * nB];
         int ga = (int)gA[ia], gb = (int)gB[ib];
         double la = (ga == 0) ? lMissA : ((ga == ha) ? lAgreeA : lTypoA);
         double lb = (gb == 0) ? lMissB : ((gb == hb) ? lAgreeB : lTypoB);
         lllq[q] = lllq[q] + logeta[ha - 1] + la + lb;
       }
     }else
     {
       // Unstable PIV: assemble the covariate matrix and times, evaluate the
       // survival function once, then fold in the extra term per pair
       CharacterVector cnA = covA_names[k], cnB = covB_names[k];
       int nca = cnA.length(), ncb = cnB.length();
       NumericMatrix cov(nPairs, 1 + nca + ncb);
       double *covp = cov.begin();
       for (int q = 0; q < nPairs; q++) covp[q] = 1.0;  // column 0: intercept
       for (int c = 0; c < nca; c++)
       {
         NumericVector col = Rcpp::as<NumericVector>(encodedA[std::string(cnA[c])]);
         const double *cp = col.begin();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
         for (int q = 0; q < nPairs; q++)
           covp[(size_t)(1 + c) * nPairs + q] = cp[sel[q] - 1];
       }
       for (int c = 0; c < ncb; c++)
       {
         NumericVector col = Rcpp::as<NumericVector>(encodedB[std::string(cnB[c])]);
         const double *cp = col.begin();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
         for (int q = 0; q < nPairs; q++)
           covp[(size_t)(1 + nca + c) * nPairs + q] = cp[sel[nPairs + q] - 1];
       }
       NumericVector dateA = Rcpp::as<NumericVector>(encodedA["date"]);
       NumericVector dateB = Rcpp::as<NumericVector>(encodedB["date"]);
       const double *dA = dateA.begin(), *dB = dateB.begin();
       NumericVector times(nPairs);
       double *tsp = times.begin();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
       for (int q = 0; q < nPairs; q++)
         tsp[q] = std::abs(dB[sel[nPairs + q] - 1] - dA[sel[q] - 1]);
       NumericVector pSameH = S(cov, alpha[k], times);
       const double *ps = pSameH.begin();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
       for (int q = 0; q < nPairs; q++)
       {
         int ia = sel[q] - 1, ib = sel[nPairs + q] - 1;
         int ha = tpA[ia + (size_t)k * nA], hb = tpB[ib + (size_t)k * nB];
         int ga = (int)gA[ia], gb = (int)gB[ib];
         double la = (ga == 0) ? lMissA : ((ga == ha) ? lAgreeA : lTypoA);
         double lb = (gb == 0) ? lMissB : ((gb == hb) ? lAgreeB : lTypoB);
         double lh = (ha == hb) ? std::log(ps[q]) : std::log((1 - ps[q]) / (nval - 1));
         lllq[q] = lllq[q] + logeta[ha - 1] + la + lb + lh;
       }
     }
   }
   List ret;
   ret["LLA"] = LLA;
   ret["LLB"] = LLB;
   ret["LLLq"] = LLLq;
   return ret;
 }

 //' sampleD
 //'
 //' Performs one Gibbs-sampling sweep over all candidate record pairs in `S`,
 //' proposing for each pair to either:
 //'  - add a link, by comparing the log-likelihood of the "linked" vs "not linked" state, or
 //'  - remove an existing link (if the pair is currently linked)
 //' Each proposal is accepted stochastically.
 //'
 //' @param S IntegerMatrix where each row correspond to the indices (from source A and source B) of records for which the true values matches (representing the potential links)
 //' @param LLA NumericVector gives the likelihood contribution of each non linked record from A
 //' @param LLB NumericVector gives the likelihood contribution of each non linked record from B
 //' @param LLL NumericVector gives the likelihood contribution of each potential linked records (from select)
 //'   — the vector `pairLogLik()` returns as `LLLq`
 //' @param gamma Numeric, the parameter gamma (proportion of linked records)
 //'   — a single scalar, as everywhere else in the model
 //' @param nlinkrec integer for the current number of linked records
 //' @param links IntegerMatrix (0-based, k x 2) of the currently linked (A, B)
 //'   record pairs: pass the `links` element returned by the previous sweep,
 //'   or a 0x2 matrix when nothing is linked yet. A sweep carries no state
 //'   between calls — everything it needs is in the arguments.
 //' @param threads Number of OpenMP threads for the per-pair delta precompute
 //'   pass (default 1 = sequential; ignored if built without OpenMP; results
 //'   are identical for any value).
 //' @return
 //' List:
 //' - new set of links (0-based, ordered by ascending A index)
 //' - linksQ: 1-based row indices in S of the finally linked pairs (same set as
 //'   `links`, unordered) — lets R look up their LLL contributions without a
 //'   sparse-matrix lookup
 //' - new sumRowD
 //' - new sumColD
 //' - new value of the complete log likelihood
 //' - new number of linked records
 //' @export
 // [[Rcpp::export]]
 List sampleD(const IntegerMatrix & S,
              const NumericVector & LLA,
              const NumericVector & LLB,
              const NumericVector & LLL,
              double gamma,
              int nlinkrec,
              const IntegerMatrix & links,
              int threads = 1)
 {
   int nB = LLB.length();
   int nA = LLA.length();
   if (threads < 1) threads = 1;
#ifndef _OPENMP
   (void)threads;  // built without OpenMP: silently sequential
#endif
   // One-to-one link state as flat partner arrays: linkOfA[i] is the 0-based
   // B partner of A-record i, or -1 when unlinked; linkOfB[j] symmetrically.
   // Array reads keep the scan's linked/unlinked tests branch-cheap and
   // updates are O(1) stores. Built from the `links` argument (the previous
   // sweep's output) and handed back to R at the end.
   std::vector<int> linkOfA(nA, -1), linkOfB(nB, -1);
   for (int r = 0; r < links.nrow(); r++)
   {
     int ia = links(r, 0), ib = links(r, 1);
     if (ia < 0 || ia >= nA || ib < 0 || ib >= nB)
       Rcpp::stop("sampleD: `links` references record indices outside 0..n-1.");
     linkOfA[ia] = ib;
     linkOfB[ib] = ia;
   }
   // The acceptance probability of every flip reduces to
   //   pswitch = exp(delta) / (1 + exp(delta)),
   // with delta the pair's own log-likelihood change: the running log-likelihood
   // cancels in the ratio, so it is not needed as an input. The reported value
   // is reconstructed from scratch as
   //   (matching + unlinked-record terms, accumulated up front)
   // + sum of LLL over finally-linked pairs (accumulated during the scan)
   // + sum of the deltas of accepted flips.
   double loglik = 0.0;
   // gamma is the model's scalar proportion of linked records; its two logs
   // are loop-invariant over the scan, so they are computed once.
   double logg = std::log(gamma), log1mg = std::log(1 - gamma);
   // log(nB - nlinkrec) (and its +1 variant for the removal branch) only
   // change when a flip is accepted, so they are maintained as scalars too.
   double logRem = std::log(nB - nlinkrec);
   double logRemP = std::log(nB - nlinkrec + 1);
   for (int i = 0; i < nA; i++)
     loglik += linkOfA[i] >= 0 ? logg : log1mg;
   for (int r = 0; r < nlinkrec; r++)
     loglik -= log(nB - r);
   // Sum the unlinked-record contributions over the overlap of the record and
   // likelihood vectors (equal lengths in the normal pipeline; the overlap
   // guard only tolerates degenerate direct calls, e.g. scalar LLA/LLB).
   int nAterm = LLA.length() < nA ? LLA.length() : nA;
   for (int i = 0; i < nAterm; i++) if (linkOfA[i] < 0) loglik += LLA(i);
   int nBterm = LLB.length() < nB ? LLB.length() : nB;
   for (int j = 0; j < nBterm; j++) if (linkOfB[j] < 0) loglik += LLB(j);
   double llLinked = 0.0;
   double llDelta = 0.0;
   std::vector<int> linkedQ;
   int nPairs = S.nrow();
   // Two-pass design: pass 1 (parallel) computes, for every candidate pair q,
   // the pair-specific part of the log-likelihood change of linking it,
   //   deltaBase[q] = -LLB(j) - LLA(i) - log(1-gamma) + log(gamma),
   // a pure gather with no coupling between pairs. Pass 2 (the sequential
   // scan below) reads it back, adds the bipartite-matching count term and
   // evaluates the acceptance ratio, so the per-pair work left on the
   // critical path is one add, one exp, one division and one comparison.
   NumericVector deltaBase(nPairs);
   double *db = deltaBase.begin();
   const int *sel = S.begin();  // column-major: sel[q] = A index, sel[nPairs + q] = B index
   const double *llap = LLA.begin(), *llbp = LLB.begin();
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
   for (int q = 0; q < nPairs; q++)
     db[q] = - llbp[sel[nPairs + q] - 1] - llap[sel[q] - 1] - log1mg + logg;
   // One Uniform(0,1) deviate per candidate pair, pre-drawn from R's RNG in a
   // single block and consumed by pair index (actionable pairs only). All
   // sampling kernels draw from such pre-drawn blocks: draws stay
   // set.seed()-controlled while the kernels themselves never call the RNG.
   NumericVector unif = Rcpp::runif(nPairs);
   const double *up = unif.begin();
   // Iterate over every candidate pair
   for (int q = 0; q < nPairs; q++)
   {
     int i = sel[q] - 1;
     int j = sel[nPairs + q] - 1;
     // If non-linked -> possibly becomes linked
     if((linkOfA[i] < 0) && (linkOfB[j] < 0))
     {
       // Non-LLL part of the log-likelihood change if this pair were linked.
       // The full delta (used for the acceptance ratio) also carries + LLL(q);
       // llLinked accounts for LLL separately, so the accumulated delta must
       // not include it (it would double-count the linked contribution).
       double delta = db[q] - logRem;
       double pswitch = 1.0 / (1.0 + exp(-(delta + LLL(q))));
       // Random number smaller than pswitch -> generates binomial value
       // Accept the new link with probability pswitch
       bool link = up[q]  < pswitch;
       if(link)
       {
         linkOfA[i] = j;
         linkOfB[j] = i;
         nlinkrec = nlinkrec + 1;
         llDelta += delta;
         llLinked += LLL(q);
         linkedQ.push_back(q + 1);
         logRem = std::log(nB - nlinkrec);
         logRemP = std::log(nB - nlinkrec + 1);
       }
     }else if(linkOfA[i] == j)
     {
       // If linked -> possibly becomes non-linked (analog). The full delta
       // carries - LLL(q); the pair's LLL was never added to llLinked, so the
       // accumulated delta excludes it.
       double delta = LLB(j) + LLA(i) + log1mg - logg + logRemP;
       double pswitch = 1.0 / (1.0 + exp(-(delta - LLL(q))));
       // Accept breaking the link with probability pswitch
       bool nolink = up[q]  < pswitch;
       if(nolink)
       {
         linkOfA[i] = -1;
         linkOfB[j] = -1;
         nlinkrec = nlinkrec - 1;
         llDelta += delta;
         logRem = std::log(nB - nlinkrec);
         logRemP = std::log(nB - nlinkrec + 1);
       }else
       {
         // Pair stays linked: its LLL contribution remains in the likelihood
         llLinked += LLL(q);
         linkedQ.push_back(q + 1);
       }
     }
   }
   loglik += llLinked + llDelta;
   // Materialise the updated links: one linear pass over linkOfA, ordered by
   // ascending A index.
   int nlinks = 0;
   for (int i = 0; i < nA; i++) if (linkOfA[i] >= 0) nlinks++;
   IntegerMatrix linksOut(nlinks, 2);
   int w = 0;
   for (int i = 0; i < nA; i++)
   {
     if (linkOfA[i] >= 0)
     {
       linksOut(w, 0) = i;
       linksOut(w++, 1) = linkOfA[i];
     }
   }
   // Linked-record flags rebuilt from the flat state (logical vectors, one
   // entry per record in A / B).
   LogicalVector rowD(nA);
   for (int i = 0; i < nA; i++) rowD(i) = linkOfA[i] >= 0;
   LogicalVector colD(nB);
   for (int j = 0; j < nB; j++) colD(j) = linkOfB[j] >= 0;
   // Return to R
   List ret;
   ret["links"] = linksOut;
   ret["linksQ"] = wrap(linkedQ);
   ret["sumRowD"] = rowD;
   ret["sumColD"] = colD;
   ret["loglik"] = loglik;
   ret["nlinkrec"] = nlinkrec;
   return ret;
 }

 //' sampleNL
 //'
 //' Samples the latent "true" value of a single PIV for non-linked records

 //' @param G IntegerVector of registered values for a certain Partially Identifying Variable (PIV) for non linked records
 //' @param eta NumericVector parameter for the distribution of the PIV concerned
 //' @param phi NumericVector parameter for the registration errors for the PIV concerned
 //' @param unif NumericVector of one Uniform(0,1) deviate per record, pre-drawn from R's RNG
 //'
 //' @return IntegerVector: of latent true values underlying G
 //' @export
 // [[Rcpp::export]]
 IntegerVector sampleNL(IntegerVector G, NumericVector eta, NumericVector phi, NumericVector unif)
 {
   IntegerVector H(G.length());
   // Number of possible values
   int nval = eta.length();
   // Possible registration errors
   double pMissing = phi[1];
   double pTypo = (1-pMissing) * (1-phi[0]) / (eta.length()-1);
   double pAgree = (1-pMissing) * phi[0];
   // Reusable weight buffer: for each record, P(True=h) for h = 1..nval
   std::vector<double> w(nval);
   // Iterate over all elements
   for(int i = 0; i < G.length(); i++)
   {
     if(G(i) == 0)
     {
       // Registered value is missing: sample from the PIV distribution
       for (int h = 0; h < nval; h++) w[h] = eta[h];
     }else
     {
       // Registered value is observed: P(True=h) proportional to eta[h] * P(G=h|True=h)
       for (int h = 0; h < nval; h++) w[h] = eta[h] * pTypo;
       w[G(i)-1] = eta[G(i)-1] * pAgree;
     }
     H(i) = sampleCdf(w.data(), nval, unif[i]) + 1;
   }
   return H;
 }

 //' sampleL
 //'
 //' Samples the latent joint "true" values of a single PIV for linked records
 //'
 //' @param GA IntegerVector of registered values for a certain Partially Identifying Variable (PIV) for linked records from A
 //' @param GB IntegerVector of registered values for a certain PIV for linked records from B
 //' @param survivalpSameH NumericVector of probabilities that the concerned PIV values coincide between file A and file B
 //' @param nval integer for the number of unique values in the PIV concerned
 //' @param phikA NumericVector parameter for the registration errors in A for the PIV concerned
 //' @param phikB NumericVector parameter for the registration errors in B for the PIV concerned
 //' @param eta NumericVector parameter for the distribution of the true value in A: the joint conditional is
 //'        proportional to eta[a] * pjoint(a,b) * P(GA=a|gA) * P(GB=b|gB) (paper Eq. 4/12)
 //' @param unif NumericVector of two Uniform(0,1) deviates per pair, pre-drawn from R's RNG
 //'
 //' @return IntegerMatrix: two columns, the latent joint true values (H^A, H^B) underlying GA and GB
 //' @export
 // [[Rcpp::export]]
 IntegerMatrix sampleL(IntegerVector GA, IntegerVector GB, NumericVector survivalpSameH,
                       int nval, NumericVector phikA, NumericVector phikB, NumericVector eta,
                       NumericVector unif)
 {
   IntegerMatrix H(GA.length(), 2);
   // Possible actions for file A
   double pMissingA = phikA[1];
   double pTypoA = (1-pMissingA) * (1-phikA[0]) / (nval-1);
   double pAgreeA = (1-pMissingA) * phikA[0];
   // Possible actions for file B
   double pMissingB = phikB[1];
   double pTypoB = (1-pMissingB) * (1-phikB[0]) / (nval-1);
   double pAgreeB = (1-pMissingB) * phikB[0];
   // Reused scratch buffers (no per-pair allocations)
   std::vector<double> wa(nval), wb(nval);
   // Iterate over all matches
   for(int i = 0; i < GA.length(); i++)
   {
     // Prob that both TRUE values are the same
     double pSameH = survivalpSameH[i];
     double pDiffH = (1 - pSameH) / (nval - 1);
     int gA = GA(i), gB = GB(i);
     // ---- Stage 1: draw H^A = a from its exact marginal.
     // Joint conditional: P(a,b) ∝ eta[a] * pjoint(a,b) * P(GA=a|gA) * P(GB=b|gB)
     // with pjoint(a,b) = pSameH if a==b else pDiffH, hence
     // P(a) ∝ eta[a] * P(GA=a|gA) * [pSameH*P(GB=a|gB) + pDiffH*(1-P(GB=a|gB))].
     for (int a = 1; a <= nval; a++)
     {
       double pGAa = (gA == 0) ? pMissingA : ((gA == a) ? pAgreeA : pTypoA);
       double pGBa = (gB == 0) ? pMissingB : ((gB == a) ? pAgreeB : pTypoB);
       wa[a-1] = eta[a-1] * pGAa * (pSameH * pGBa + pDiffH * (1 - pGBa));
     }
     int a = sampleCdf(wa.data(), nval, unif[2*i]) + 1;
     // ---- Stage 2: draw H^B = b | a. P(b|a) ∝ pjoint(a,b) * P(GB=b|gB).
     double pGBa = (gB == 0) ? pMissingB : ((gB == a) ? pAgreeB : pTypoB);
     for (int b = 1; b <= nval; b++)
     {
       double pGBb = (gB == 0) ? pMissingB : ((gB == b) ? pAgreeB : pTypoB);
       wb[b-1] = (b == a) ? (pSameH * pGBa) : (pDiffH * pGBb);
     }
     H(i,0) = a;
     H(i,1) = sampleCdf(wb.data(), nval, unif[2*i+1]) + 1;
   }
   return H;
 }

 //' ExpandGrid
 //'
 //' Compute the cartesian product of vectors
 //'
 //' @param vector1 first IntegerVector of values to compute the cartesian product
 //' @param vector2 second IntegerVector of values to compute the cartesian product
 //'
 //' @return IntegerMatrix: of 2 columns with the cartesian product of vec1 and vec2
 //' @export
 // [[Rcpp::export]]
 IntegerMatrix ExpandGrid(IntegerVector vec1, IntegerVector vec2) {
   int n = vec1.size();
   int m = vec2.size();
   IntegerMatrix result(n * m, 2);
   for (int i = 0; i < n; i++) {
     for (int j = 0; j < m; j++) {
       result(i * m + j, 0) = vec1[i];
       result(i * m + j, 1) = vec2[j];
     }
   }
   return result;
 }

 //' generateSequence
 //'
 //' @param n integer superior to 1
 //'
 //' @return IntegerVector: with n values from 1 to n
 //' @export
 // [[Rcpp::export]]
 IntegerVector generateSequence(int n) {
   IntegerVector result(n);
   for (int i = 0; i < n; i++) {
     result[i] = i + 1;
   }
   return result;
 }

 //' sampleH
 //'
 //' @param nA IntegerVector of dimensions of registered values of the Partially Identifying Variables (PIVs) in A
 //' @param nB IntegerVector of dimensions of registered values of the PIVs in B
 //' @param links IntegerMatrix of 2 columns with the indices of the linked records
 //' @param survivalpSameH NumericMatrix with for each PIV the probability that true values coincide (if stable: filled with 1)
 //' @param pivs_stable LogicalVector indicating for each PIV whether it is stable of not (if not we expect survivalpSameH for that same element to not be filled with 1 but with lower values)
 //' @param pivsA List ith registered data from A
 //' @param pivsB List with registered data from B
 //' @param nvalues IntegerVector with number of unique values of each PIV
 //' @param nonlinkedA LogicalVector indicating for all records in A whether they are linked or not
 //' @param nonlinkedB LogicalVector indicating for all records in B whether they are linked or not
 //' @param eta List parameters of the PIVs distributions
 //' @param phi List parameters of the PIVs registration errors
 //'
 //' @return List:
 //' - truePIVsA, true values underlying data in A
 //' - truePIVsB, true values underlying data in B
 //' @export
 // [[Rcpp::export]]
 List sampleH(IntegerVector nA, IntegerVector nB, IntegerMatrix links, NumericMatrix survivalpSameH, LogicalVector pivs_stable, List pivsA, List pivsB, IntegerVector nvalues, LogicalVector nonlinkedA, LogicalVector nonlinkedB, List eta, List phi)
 {
   // Output matrices of true PIV values, same shape as the registered data
   IntegerMatrix truepivsA(nA[0],nA[1]);
   IntegerMatrix truepivsB(nB[0],nB[1]);
   int nphi = 2;
   // Pre-draw all Uniform(0,1) deviates for this sweep from R's RNG (set.seed()
   // controlled, sequential). The sampling kernels consume their own slots, so
   // they never touch the RNG themselves — thread-safe by construction.
   int nNL_A = 0, nNL_B = 0;
   for (int i = 0; i < nonlinkedA.length(); i++) if (nonlinkedA[i]) nNL_A++;
   for (int i = 0; i < nonlinkedB.length(); i++) if (nonlinkedB[i]) nNL_B++;
   int nPIVs = nvalues.length();
   // Two deviates per linked pair for sampleL (marginal draw + conditional draw)
   int nlinks = links.nrow();
   int nLunif = (nlinks > 0) ? nPIVs * 2 * nlinks : 0;
   NumericVector unif = Rcpp::runif(nPIVs * (nNL_A + nNL_B) + nLunif);
   int off = 0;
   for (int k = 0; k < nvalues.length(); k++) {
     IntegerVector truepivsA_k = truepivsA(_,k);
     IntegerVector truepivsB_k = truepivsB(_,k);
     NumericVector eta_k = eta[k];
     NumericVector phi_k = phi[k];
     // Unpack the registration-error parameters for this PIV:
     // phi_k = [agree_A, agree_B, missing_A, missing_B]
     NumericVector phi_k_A(nphi);
     phi_k_A[0] = phi_k[0];
     phi_k_A[1] = phi_k[2];
     NumericVector phi_k_B(nphi);
     phi_k_B[0] = phi_k[1];
     phi_k_B[1] = phi_k[3];
     IntegerVector pivsA_k = pivsA[k];
     IntegerVector pivsB_k = pivsB[k];
     IntegerVector pivsA_k_L = pivsA_k[links(_,0)];
     IntegerVector pivsB_k_L = pivsB_k[links(_,1)];
     IntegerVector pivsA_k_NL = pivsA_k[nonlinkedA];
     IntegerVector pivsB_k_NL = pivsB_k[nonlinkedB];
     // MArginally sample true values for non-linked pairs
     IntegerVector truepivsA_k_NL = sampleNL(pivsA_k_NL, eta_k, phi_k_A,
                                             nNL_A > 0 ? unif[Range(off, off + nNL_A - 1)] : NumericVector(0));
     off += nNL_A;
     IntegerVector truepivsB_k_NL = sampleNL(pivsB_k_NL, eta_k, phi_k_B,
                                             nNL_B > 0 ? unif[Range(off, off + nNL_B - 1)] : NumericVector(0));
     off += nNL_B;
     truepivsA_k[nonlinkedA] = truepivsA_k_NL;
     truepivsB_k[nonlinkedB] = truepivsB_k_NL;
     if (links.nrow()>0)
     {
       int nlinks = links.nrow();
       NumericVector survivalpSameH_k = survivalpSameH(_,k);
       // Jointly sample true values for linked pairs. eta stays indexed by the
       // true value of A (paper Eq. 4: marginal of H^A is eta; H^B | H^A
       // follows the survival mechanism).
       IntegerMatrix out = sampleL(pivsA_k_L, pivsB_k_L, survivalpSameH_k, nvalues[k], phi_k_A, phi_k_B, eta_k,
                                   unif[Range(off, off + 2*nlinks - 1)]);
       off += 2*nlinks;
       IntegerVector outA = out(_,0);
       IntegerVector outB = out(_,1);
       truepivsA_k[links(_,0)] = outA;
       truepivsB_k[links(_,1)] = outB;
     }
     truepivsA(_,k) = truepivsA_k;
     truepivsB(_,k) = truepivsB_k;
   }
   List ret;
   ret["truepivsA"] = truepivsA;
   ret["truepivsB"] = truepivsB;
   return ret;
 }

 //' survivalS_cpp
 //'
 //' Fast C++ evaluation of the survival function S(t | X, alpha) for the
 //' built-in survival models: the exact formulas of the R closures built by
 //' [survival_model()] (see its documentation), computed in one fused pass
 //' over the rows (dot product, then the family's exp/pow terms), with the
 //' row loop parallelisable over OpenMP threads. Reached through the
 //' internal R dispatcher `survival_S()`, which custom ("custom"-type)
 //' models bypass entirely — they keep calling the user's R function.
 //'
 //' Differences from the R closures are at BLAS-rounding level (~1e-16:
 //' R's `%*%` is dgemv with its own summation order; R's `sum()` accumulates
 //' in extended precision), far inside every tolerance the package uses.
 //' Results are bitwise identical for any `threads` value (each row is
 //' computed by exactly one thread).
 //'
 //' @param cov NumericMatrix (n x (1 + n_cov)): intercept column first, then
 //'   the covariates — the same X the R S functions receive. The piecewise
 //'   family ignores the intercept column, exactly like its R closure.
 //' @param alpha NumericVector of the model's parameters (layout per family,
 //'   as documented in [survival_model()])
 //' @param times NumericVector of the time gaps (one per row)
 //' @param type One of "exponential", "weibull", "gompertz", "piecewise"
 //' @param cuts NumericVector of the piecewise cut points (empty for the
 //'   other families)
 //' @param threads Number of OpenMP threads (default 1 = sequential; ignored
 //'   if the package was built without OpenMP)
 //' @return NumericVector of survival probabilities, one per row
 //' @export
 // [[Rcpp::export]]
 NumericVector survivalS_cpp(const NumericMatrix & cov,
                             const NumericVector & alpha,
                             const NumericVector & times,
                             const std::string & type,
                             const NumericVector & cuts,
                             int threads = 1)
 {
   int n = times.length();
   if (cov.nrow() != n)
     Rcpp::stop("survivalS_cpp: cov has %d rows but there are %d times.", cov.nrow(), n);
   int nCov = cov.ncol();
   int na = alpha.length();
   if (threads < 1) threads = 1;
#ifndef _OPENMP
   (void)threads;  // built without OpenMP: silently sequential
#endif
   // Column-major (like R): cov(i, c) = cp[(size_t)c * n + i]. Each thread
   // streams a contiguous chunk of each column, so every element is read
   // exactly once regardless of the thread count.
   const double *cp = cov.begin();
   const double *tp = times.begin();
   const double *a = alpha.begin();
   NumericVector S(n);
   double *sp = S.begin();

   if (type == "exponential")
   {
     // R: exp(-ph(X, alpha) * times), ph = exp(X %*% alpha)
     if (na != nCov)
       Rcpp::stop("survivalS_cpp: exponential expects %d alphas (one per cov column), got %d.", nCov, na);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
     for (int i = 0; i < n; i++)
     {
       double d = 0.0;
       for (int c = 0; c < nCov; c++) d += cp[(size_t)c * n + i] * a[c];
       sp[i] = std::exp(-std::exp(d) * tp[i]);
     }
   }else if (type == "weibull")
   {
     // R: exp(-(ph(X, alpha[-1]) * times)^exp(alpha[1]))
     if (na != nCov + 1)
       Rcpp::stop("survivalS_cpp: weibull expects %d alphas (log-shape + one per cov column), got %d.", nCov + 1, na);
     double k = std::exp(a[0]);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
     for (int i = 0; i < n; i++)
     {
       double d = 0.0;
       for (int c = 0; c < nCov; c++) d += cp[(size_t)c * n + i] * a[c + 1];
       sp[i] = std::exp(-std::pow(std::exp(d) * tp[i], k));
     }
   }else if (type == "gompertz")
   {
     // R: g <- alpha[1]; if (abs(g) < 1e-8) exp(-ph(X, alpha[-1]) * times)
     //    else exp(-ph(X, alpha[-1]) * (exp(g * times) - 1) / g)
     if (na != nCov + 1)
       Rcpp::stop("survivalS_cpp: gompertz expects %d alphas (g + one per cov column), got %d.", nCov + 1, na);
     double g = a[0];
     bool small_g = std::abs(g) < 1e-8;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
     for (int i = 0; i < n; i++)
     {
       double d = 0.0;
       for (int c = 0; c < nCov; c++) d += cp[(size_t)c * n + i] * a[c + 1];
       double lam = std::exp(d);
       // Same expression shape as the R closure: (exp(g t) - 1) / g, not
       // expm1, so the values stay at rounding-level distance from R even
       // for small g.
       sp[i] = small_g ? std::exp(-lam * tp[i])
                       : std::exp(-lam * (std::exp(g * tp[i]) - 1.0) / g);
     }
   }else if (type == "piecewise")
   {
     // R: h <- exp(alpha[seq_len(n_int)])
     //    H0 <- sum(h * pmax(0, pmin(t, edges[-1]) - edges[-length(edges)]))
     //    lam <- if (length(beta)) ph(X[, -1], beta) else 1
     //    exp(-lam * H0), edges = c(0, cuts, Inf)
     int n_int = cuts.length() + 1;
     if (na != n_int + (nCov - 1))
       Rcpp::stop("survivalS_cpp: piecewise expects %d alphas (%d interval log-hazards + %d covariate coefficients), got %d.",
                  n_int + nCov - 1, n_int, nCov - 1, na);
     for (int k = 1; k < cuts.length(); k++)
       if (cuts[k] <= cuts[k - 1])
         Rcpp::stop("survivalS_cpp: piecewise `cuts` must be increasing.");
     std::vector<double> edges(n_int + 1);
     edges[0] = 0.0;
     edges[n_int] = INFINITY;
     for (int k = 1; k < n_int; k++) edges[k] = cuts[k - 1];
     std::vector<double> h(n_int);
     for (int k = 0; k < n_int; k++) h[k] = std::exp(a[k]);
     const double *beta = a + n_int;  // coefficients for cov columns 1..nCov-1
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static)
#endif
     for (int i = 0; i < n; i++)
     {
       double t = tp[i];
       double H0 = 0.0;
       for (int k = 0; k < n_int; k++)
       {
         double spent = std::min(t, edges[k + 1]) - edges[k];
         if (spent > 0.0) H0 += h[k] * spent;
       }
       double lam = 1.0;
       if (nCov > 1)
       {
         double d = 0.0;
         for (int c = 1; c < nCov; c++) d += cp[(size_t)c * n + i] * beta[c - 1];
         lam = std::exp(d);
       }
       sp[i] = std::exp(-lam * H0);
     }
   }else
   {
     Rcpp::stop("survivalS_cpp: unknown type '%s' (built-in families only; custom models keep their R function).", type.c_str());
   }
   return S;
 }
