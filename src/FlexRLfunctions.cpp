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

//' @export
// [[Rcpp::export]]
NumericVector callFunction(NumericVector x, Function f) {
  NumericVector res = f(x);
  return res;
}

double ControlRandomNumberGen() {
  double randomnumber = R::runif(0,1);
  return randomnumber;
}

double ControlRandomSampleGen(IntegerVector choiceset, NumericVector probavec) {
  double generatedrandom = Rcpp::sample(choiceset, 1, false, probavec)[0];
  return generatedrandom;
}

//' @export
 // [[Rcpp::export]]
List F2(IntegerVector U, int nvals)
{
  List out(nvals);
  for (int i = 0; i < nvals; i++)
  {
    IntegerVector tmp;
    out[i] = tmp;
  }
  // Which value does observation i have. Assign it to the appropriate "bucket" in the list
  for (int i = 0; i < U.length(); i++)
  {
    IntegerVector tmpX = out[U[i]-1];
    tmpX.push_back(i+1);
    out[U[i]-1] = tmpX;
  }
  return out;
}

// Cartesian product of both lists. Which combos are possible based on the true values?
//' @export
 // [[Rcpp::export]]
IntegerMatrix F33(List A, List B, int nvals)
{
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

//' @export
// [[Rcpp::export]]
CharacterVector sspaste2(IntegerMatrix A)
{
  int i = 0, j = 0;
  int sz = A.nrow();
  CharacterVector res(sz);
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

List F11(const List & F, const int & nvals)
{
  List out(nvals);
  for (int k = 0; k < nvals; k++)
  {
    List tmp = F[k];
    IntegerVector tmpA = tmp[0];
    IntegerVector tmpB = tmp[1];
    IntegerMatrix tmpC(tmpA.length() * tmpB.length(),2);
    int counter = 0;
    for (int i = 0; i < tmpA.length(); i++)
    {
      for (int j = 0; j < tmpB.length(); j++)
      {
        tmpC(counter,0) = tmpA[i];
        tmpC(counter,1) = tmpB[j];
        counter += 1;
      }
    }
    out[k] = tmpC;
  }
  return out;
}

List F1(const IntegerVector & HA, const IntegerVector & HB, const int & nvals)
{
  List out(nvals);
  for (int i = 0; i < nvals; i++)
  { IntegerVector tmpA;
    IntegerVector tmpB;
    List both(2);
    both[0] = tmpA;
    both[1] = tmpB;
    out[i] = both;
  }
  for (int i = 0; i < HA.length(); i++)
  {
    List tmp = out[HA[i]-1];
    IntegerVector tmpA = tmp[0];
    tmpA.push_back(i+1);
    tmp[0] = tmpA;
    out[HA[i]-1] = tmp;
  }
  for (int j = 0; j < HB.length(); j++)
  {
    List tmp = out[HB[j]-1];
    IntegerVector tmpB = tmp[1];
    tmpB.push_back(j+1);
    tmp[1] = tmpB;
    out[HB[j]-1] = tmp;
  }
  return out;
}

std::map<int, std::set<int>> _ΔMap;

//' @export
 // [[Rcpp::export]]
void initΔMap()
{
  _ΔMap = {};
}

//' @export
 // [[Rcpp::export]]
IntegerMatrix Δfind()
{
  int length = 0;
  for(const auto & rowColSet: _ΔMap)
  {
    length += rowColSet.second.size();
  }
  IntegerMatrix found(length, 2);
  int index = 0;
  for(const auto & rowColSet: _ΔMap)
  {

    for (int c : rowColSet.second)
    {
      found(index,0)   = rowColSet.first;
      found(index++,1) = c;
    }
  }
  if(index != length)
  {
    Rf_error("Something went wrong creating Δ.");
  }
  return found;
}

//' @export
// [[Rcpp::export]]
List sampleD(const Rcpp::IntegerMatrix & S,
             const Rcpp::NumericVector & LLA,
             const Rcpp::NumericVector & LLB,
             const Rcpp::NumericVector & LLL,
             const Rcpp::NumericVector & gamma,
             double loglik,
             int nlinkrec,
             LogicalVector & sumRowD,
             LogicalVector & sumColD)
{
  for (int q = 0; q < S.nrow(); q++)
  {
    int i = S(q,0)-1;
    int j = S(q,1)-1;
    // If non match and possible match -> check if match
    if((sumRowD(i)==false) && (sumColD(j)==false))
    {

      double loglikNew = loglik
      // Comparison vectors
      - LLB(j) - LLA(i)
      + LLL(q)
      // Bipartite matching
      - log(1-gamma(i)) + log(gamma(i))
      - log(LLB.length() - nlinkrec);
      double sumlogdensity = log(1 + exp(loglik-loglikNew)) + loglikNew;
      double pswitch = exp(loglikNew - sumlogdensity);
      // Random number smaller than prob -> generate binomial value
      bool link = ControlRandomNumberGen()  < pswitch;
      if(link)
      {
        loglik = loglikNew;
        _ΔMap[i].insert(j);
        sumRowD(i) = true;
        sumColD(j) = true;
        nlinkrec = nlinkrec + 1;
      }
    }else if(_ΔMap.count(i) && _ΔMap.at(i).count(j))
    {
      // If match -> check if nonmatch
      double loglikNew = loglik
      // Comparison vectors
      + LLB(j) + LLA(i)
      - LLL(q)
      // Bipartite matching
      + log(1-gamma(i)) - log(gamma(i))
      + log(LLB.length() - nlinkrec+1);
      double sumlogdensity = log(1 + exp(loglik-loglikNew)) + loglikNew;
      double pswitch = exp(loglikNew - sumlogdensity);
      bool nolink = ControlRandomNumberGen()  < pswitch;
      if(nolink)
      {
        loglik = loglikNew;
        if(_ΔMap.count(i) && _ΔMap[i].count(j))
          _ΔMap[i].erase(j);
        sumRowD(i) = false;
        sumColD(j) = false;
        nlinkrec = nlinkrec - 1;
      }
    }
  }
  IntegerMatrix links=Δfind();
  // Return to R
  List ret;
  ret["links"] = links;
  ret["sumRowD"] = sumRowD;
  ret["sumColD"] = sumColD;
  ret["loglik"] = loglik;
  ret["nlinkrec"] = nlinkrec;
  return ret;
}

//' @export
// [[Rcpp::export]]
IntegerVector sampleNL(IntegerVector G, NumericVector eta, NumericVector phi)
{
  IntegerVector H(G.length());
  // Number of possible values
  int nval = eta.length();
  // Create the possible values to sample from
  IntegerVector choice_set = seq_len(nval);
  // Possible registration errors
  double pMissing = phi[1];
  double pTypo = (1-pMissing) * (1-phi[0]) / (eta.length()-1);
  double pAgree = (1-pMissing) * phi[0];
  // Iterate over all elements
  for(int i = 0; i < G.length(); i++)
  {
    // Create a vector indicating P(Registered=X|True)
    // First value is for the missings
    // What happens if missing:
    if(G(i) == 0)
    {
      // All equally likely
      H(i) = ControlRandomSampleGen(choice_set, eta);
    }else
    {
      NumericVector help1(nval, pTypo);
      help1(G(i)-1) = pAgree;
      // Joint probability to have the registered and true value
      NumericVector prob = eta * help1;
      H(i) = ControlRandomSampleGen(choice_set, prob);
    }
  }
  return H;
}

//' @export
// [[Rcpp::export]]
IntegerVector sampleL(IntegerVector GA, IntegerVector GB, NumericVector survivalpSameH,
                      IntegerMatrix choice_set, IntegerVector choice_equal,
                      int nval, NumericVector phikA, NumericVector phikB, NumericVector eta)
{
  IntegerVector H(GA.length());
  int size_choice_set = choice_set.nrow();
  IntegerVector choice_index = seq_len(size_choice_set);
  // Possible actions for file A
  double pMissingA = phikA[1];
  double pTypoA = (1-pMissingA) * (1-phikA[0]) / (nval-1);
  double pAgreeA = (1-pMissingA) * phikA[0];
  // Possible actions for file B
  double pMissingB = phikB[1];
  double pTypoB = (1-pMissingB) * (1-phikB[0]) / (nval-1);
  double pAgreeB = (1-pMissingB) * phikB[0];
  // Iterate over all matches
  for(int i = 0; i < GA.length(); i++)
  {
    // Prob that both TRUE values are the same
    double pSameH = survivalpSameH[i];
    // Dummy Vectors
    // Define P(Hb|Ha)
    NumericVector probH(size_choice_set, pSameH);
    // Define P(Ga|Ha) and P(Gb|Hb)
    NumericVector helpA(size_choice_set, pTypoA);
    NumericVector helpB(size_choice_set, pTypoB);
    for(int j = 0; j < size_choice_set; j++)
    {
      // Prob to observe HB=b|HA=a where a!=b
      if(choice_equal(j)==0)
      {
        probH(j) = (1-pSameH)/(nval-1); //
      }
      if(GA(i) == choice_set(j,0))
        helpA(j) = pAgreeA;
      if(GB(i) == choice_set(j,1))
        helpB(j) = pAgreeB;
    }
    // There are four options possible
    // Both missing
    if(GA(i)==0 && GB(i)==0)
    {
      NumericVector prob = eta * probH;
      H(i) = ControlRandomSampleGen(choice_index, prob);
    }else if(GA(i)>0 && GB(i)==0)
    {
      // Joint probability to have the registered and true value
      NumericVector prob = eta * probH * helpA;
      H(i) = ControlRandomSampleGen(choice_index, prob);
    }else if(GA(i)==0 && GB(i)>0)
    {
      // Joint probability to have the registered and true value
      NumericVector prob = eta * probH * helpB;
      H(i) = ControlRandomSampleGen(choice_index, prob);
    }else if(GA(i)>0 && GB(i)>0)
    {
      // None missing
      // Create vectors indicating P(Registered=X|True)
      // Joint probability to have the registered and true value
      NumericVector prob = eta * probH * helpA * helpB;
      H(i) = ControlRandomSampleGen(choice_index, prob);
    }
  }
  return H;
}

//' @export
// [[Rcpp::export]]
IntegerMatrix cartesianProduct(IntegerVector vec1, IntegerVector vec2) {
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

//' @export
// [[Rcpp::export]]
IntegerMatrix ExpandGrid(IntegerVector vector1, IntegerVector vector2) {
  return cartesianProduct(vector1, vector2);
}

//' @export
// [[Rcpp::export]]
IntegerVector generateSequence(int n) {
  IntegerVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = i + 1;
  }
  return result;
}

//' @export
// [[Rcpp::export]]
List sampleH(IntegerVector nA, IntegerVector nB, IntegerMatrix links, NumericMatrix survivalpSameH, LogicalVector pivs_stable, List pivsA, List pivsB, IntegerVector nvalues, LogicalVector nonlinkedA, LogicalVector nonlinkedB, List eta, List phi)
{
  IntegerMatrix truepivsA(nA[0],nA[1]);
  IntegerMatrix truepivsB(nB[0],nB[1]);
  int nphi = 2;
  for (int k = 0; k < nvalues.length(); k++) {
    IntegerVector truepivsA_k = truepivsA(_,k);
    IntegerVector truepivsB_k = truepivsB(_,k);
    NumericVector eta_k = eta[k];
    NumericVector phi_k = phi[k];
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
    IntegerVector truepivsA_k_NL = sampleNL(pivsA_k_NL, eta_k, phi_k_A);
    IntegerVector truepivsB_k_NL = sampleNL(pivsB_k_NL, eta_k, phi_k_B);
    truepivsA_k[nonlinkedA] = truepivsA_k_NL;
    truepivsB_k[nonlinkedB] = truepivsB_k_NL;
    IntegerMatrix choice_set;
    IntegerVector choice_equal;
    NumericVector eta_choice;
    if (links.nrow()>0)
    {
      IntegerVector values = generateSequence(nvalues[k]);
      choice_set = ExpandGrid(values, values);
      choice_equal = choice_set(_,1) == choice_set(_,0);
      eta_choice = eta_k[choice_set(_,1) - 1];
      NumericVector survivalpSameH_k = survivalpSameH(_,k);
      IntegerVector out = sampleL(pivsA_k_L, pivsB_k_L, survivalpSameH_k, choice_set, choice_equal, nvalues[k], phi_k_A, phi_k_B, eta_choice);
      IntegerVector choice_set_A = choice_set(_,0);
      IntegerVector choice_set_B = choice_set(_,1);
      truepivsA_k[links(_,0)] = choice_set_A[out - 1];
      truepivsB_k[links(_,1)] = choice_set_B[out - 1];
    }
    truepivsA(_,k) = truepivsA_k;
    truepivsB(_,k) = truepivsB_k;
  }
  List ret;
  ret["truepivsA"] = truepivsA;
  ret["truepivsB"] = truepivsB;
  return ret;
}
