#include <random>
#include <chrono>

#include <Rcpp.h>

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

 // Draws a single random number from Uniform(0,1). Separated to control from R.
 double drawUniform01() {
   double randomnumber = R::runif(0,1);
   return randomnumber;
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

 std::map<int, std::set<int>> _DeltaMap;

 //' initDeltaMap
 //'
 //' @return void: Initialise the cpp map _DeltaMap representing the sparse linkage matrix Delta.
 //' @export
 // [[Rcpp::export]]
 void initDeltaMap()
 {
   _DeltaMap = {};
 }

 //' Deltafind()
 //'
 //' @return IntegerMatrix: find the indices of the elements in the cpp map _DeltaMap representing the sparse linkage matrix Delta.
 //' @export
 // [[Rcpp::export]]
 IntegerMatrix Deltafind()
 {
   // Count total number of links
   int length = 0;
   for(const auto & rowColSet: _DeltaMap)
   {
     length += rowColSet.second.size();
   }
   IntegerMatrix found(length, 2);
   int index = 0;
   // Flatten the map
   for(const auto & rowColSet: _DeltaMap)
   {

     for (int c : rowColSet.second)
     {
       found(index,0)   = rowColSet.first;
       found(index++,1) = c;
     }
   }
   if(index != length)
   {
     Rcpp::stop("Something went wrong creating Δ.");
   }
   return found;
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
 //' @param gamma NumericVector repeats the value of the parameter gamma (proportion of linked records) number of potential linked records (nrow of S) times
 //' @param loglik double for the value of the current complete log likelihood of the model
 //' @param nlinkrec integer for the current number of linked records
 //' @param sumRowD A LogicalVector vector indicating, for each row of the linkage matrix, i.e. for each record in the smallest file A, whether the record has a link in B or not.
 //' @param sumColD A LogicalVector vector indicating, for each column of the linkage matrix, i.e. for each record in the largest file B, whether the record has a link in A or not.
 //'
 //' @return
 //' List:
 //' - new set of links
 //' - new sumRowD
 //' - new sumColD
 //' - new value of the complete log likelihood
 //' - new number fo linked records
 //' @export
 // [[Rcpp::export]]
 List sampleD(const IntegerMatrix & S,
              const NumericVector & LLA,
              const NumericVector & LLB,
              const NumericVector & LLL,
              const NumericVector & gamma,
              double loglik,
              int nlinkrec,
              LogicalVector & sumRowD,
              LogicalVector & sumColD)
 {
   // Work on copies: when a caller passes a length-1 logical (e.g. the result
   // of a scalar comparison), R hands us its interned TRUE/FALSE constant and
   // writing through the reference corrupts R's global constant (seen as
   // "internal FALSE value has been modified" and flipped comparison results
   // session-wide, e.g. `1 > 2` evaluating to TRUE).
   LogicalVector rowD = Rcpp::clone(sumRowD);
   LogicalVector colD = Rcpp::clone(sumColD);
   // Iterate over every candidate pair
   for (int q = 0; q < S.nrow(); q++)
   {
     int i = S(q,0)-1;
     int j = S(q,1)-1;
     // If non-linked -> possibly becomes linked
     if((rowD(i)==false) && (colD(j)==false))
     {
       // Log-likelihood if this pair were linked:
       // remove the two "unlinked" contributions, add the "linked" contribution,
       // and adjust for the bipartite-matching
       double loglikNew = loglik
       // Comparisons
       - LLB(j) - LLA(i)
       + LLL(q)
       // Bipartite matching
       - log(1-gamma(i)) + log(gamma(i))
       - log(LLB.length() - nlinkrec);
       double sumlogdensity = log(1 + exp(loglik-loglikNew)) + loglikNew;
       double pswitch = exp(loglikNew - sumlogdensity);
       // Random number smaller than pswitch -> generates binomial value
       // Accept the new link with probability pswitch
       bool link = drawUniform01()  < pswitch;
       if(link)
       {
         loglik = loglikNew;
         _DeltaMap[i].insert(j);
         rowD(i) = true;
         colD(j) = true;
         nlinkrec = nlinkrec + 1;
       }
     }else if(_DeltaMap.count(i) && _DeltaMap.at(i).count(j))
     {
       // If linked -> possibly becomes non-linked (analog)
       double loglikNew = loglik
       // Comparisons
       + LLB(j) + LLA(i)
       - LLL(q)
       // Bipartite matching
       + log(1-gamma(i)) - log(gamma(i))
       + log(LLB.length() - nlinkrec+1);
       double sumlogdensity = log(1 + exp(loglik-loglikNew)) + loglikNew;
       double pswitch = exp(loglikNew - sumlogdensity);
       // Accept breaking the link with probability pswitch
       bool nolink = drawUniform01()  < pswitch;
       if(nolink)
       {
         loglik = loglikNew;
         if(_DeltaMap.count(i) && _DeltaMap[i].count(j))
           _DeltaMap[i].erase(j);
         rowD(i) = false;
         colD(j) = false;
         nlinkrec = nlinkrec - 1;
       }
     }
   }
   // Convert the updated sparse linkage map into the (idxA, idxB) matrix form
   IntegerMatrix links=Deltafind();
   // Return to R
   List ret;
   ret["links"] = links;
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
