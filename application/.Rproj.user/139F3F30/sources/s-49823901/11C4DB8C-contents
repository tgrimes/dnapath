// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec d_edgesC(mat nw1, mat nw2, double lp) {
  int p = nw1.n_cols;
  int N = p * (p - 1) / 2;
  vec diff(N);
  
  int k = 0;
  for(int j = 0; j < p - 1; j++) {
    for(int i = j + 1; i < p; i++) {
      diff[k] = pow(abs(nw1(j, i) - nw2(j, i)), lp);
      k += 1;
    }
  }
  
  return diff;
}

// [[Rcpp::export]]
double d_pathwayC(mat nw1, mat nw2, double lp) {
  int p = nw1.n_cols;
  int N = p * (p - 1) / 2;
  double diff = accu(d_edgesC(nw1, nw2, lp)) / (double) N;
  
  if(lp > 1) {
    diff = pow(diff, 1 / lp);
  }
  
  return diff;
}

// [[Rcpp::export]]
vec d_genesC(mat nw1, mat nw2, double lp) {
  int p = nw1.n_cols;
  vec diff(p);
  
  for(int i = 0; i < p; i++) {
    diff[i] = 0;
  }
  
  for(int i = 0; i < p; i++) {
    for(int j = 0; j < p; j++) {
      if(i != j) {
        diff[i] += pow(abs(nw1(i, j) - nw2(i, j)), lp);
      }
    }
    diff[i] = diff[i] / (p - 1);
    if(lp < 1) {
      diff[i] = pow(diff[i], 1 / lp);
    }
  }
  return diff;
}

// [[Rcpp::export]]
uvec setdiffC(uvec x, const uvec y){
  
  for(int i = 0; i < y.n_elem; i++) {
    uword q1 = conv_to<uword>::from(find(x == y[i]));
    x.shed_row(q1);
  }
  return x;
}

// [[Rcpp::export]]
mat scaleC(mat x) {
  for(int i = 0; i < x.n_cols; i++) {
    x.col(i) = (x.col(i) - mean(x.col(i))) / stddev(x.col(i));
  }
  return x;
}

// [[Rcpp::export]]
mat corC(mat x) {
  return cor(x);
}

// [[Rcpp::export]]
mat centerC(mat x) {
  for(int i = 0; i < x.n_cols; i++) {
    x.col(i) = x.col(i) - mean(x.col(i));
  }
  return x;
}

double sqr_errC(vec x, vec y) {
  vec vec_prod = x % y;
  return accu(pow(vec_prod - mean(vec_prod), 2));
}

// [[Rcpp::export]]
mat cor_shrinkC(mat x) {
  int p = x.n_cols;
  double n = (double)x.n_rows;
  
  // Calculate r and var(r), using centered and scaled x.
  x = scaleC(x);
  mat r = cor(x);
  double numer_1 = 0;
  double denom_1 = 0;
  vec vec_prod(p);
  for(int j = 0; j < p - 1; j++) {
    for(int i = j + 1; i < p; i++) {
      vec_prod = x.col(i) % x.col(j);
      numer_1 += accu(pow(vec_prod - mean(vec_prod), 2));
      denom_1 += pow(r(i, j), 2);
    }
  }
  numer_1 = numer_1 * n * pow(n - 1, -3);
  
  double lambda_1 = numer_1 / denom_1;
  if(denom_1 == 0) lambda_1 = 0;
  if(lambda_1 > 1) lambda_1 = 1;
  
  mat r_star = (1 - lambda_1) * r;
  r_star.diag().ones();
  
  return r_star;
}


// [[Rcpp::export]]
mat cov2cor(mat x) {
  vec v = x.diag();
  mat v_mat = sqrt(v * v.t());
  x = x / v_mat;
  
  return x;
}

// [[Rcpp::export]]
mat pcor_shrinkC(mat x) {
  mat r = cor_shrinkC(x);
  r = -inv(r);
  r.diag() = -r.diag();
  vec v = r.diag();
  mat v_mat = sqrt(v * v.t());
  r = r / v_mat;
  
  return r;
}

// [[Rcpp::export]]
List dna_pcorC(mat counts, umat permutations, double lp) {
  int n = counts.n_rows;
  int n_genes = counts.n_cols;
  int n_edges = n_genes * (n_genes - 1) / 2;
  int n_perm = permutations.n_cols;
  
  uvec all_samples(n);
  for(int i = 0; i < n; i++) {
    all_samples(i) = i + 1;
  }
  
  uvec index1 = permutations.col(0);
  uvec index2 = setdiffC(all_samples, index1);
  
  mat scores1 = pcor_shrinkC(counts.rows(index1 - 1));
  mat scores2 = pcor_shrinkC(counts.rows(index2 - 1));
  
  vec d_gene_scores = d_genesC(scores1, scores2, lp);
  vec d_edge_scores = d_edgesC(scores1, scores2, lp);
  
  vec p_value_genes;
  p_value_genes.ones(n_genes);
  vec p_value_edges;
  p_value_edges.ones(n_edges);
  
  for(int j = 1; j < n_perm; j++) {
    index1 = permutations.col(j);
    index2 = setdiffC(all_samples, index1);
    
    scores1 = pcor_shrinkC(counts.rows(index1 - 1));
    scores2 = pcor_shrinkC(counts.rows(index2 - 1));
    
    p_value_genes = p_value_genes + (d_genesC(scores1, scores2, lp) >= d_gene_scores);
    p_value_edges = p_value_edges + (d_edgesC(scores1, scores2, lp) >= d_edge_scores);
  }
  
  p_value_genes = p_value_genes / (double)n_perm;
  p_value_edges = p_value_edges / (double)n_perm;
  
  return List::create(Named("p_value_genes") = p_value_genes,
                      Named("p_value_edges") = p_value_edges,
                      Named("d_genes") = d_gene_scores,
                      Named("d_edges") = d_edge_scores,
                      Named("lp") = lp);
}

// [[Rcpp::export]]
List dna_corC(mat counts, umat permutations, double lp) {
  int n = counts.n_rows;
  int n_genes = counts.n_cols;
  int n_edges = n_genes * (n_genes - 1) / 2;
  int n_perm = permutations.n_cols;
  
  uvec all_samples(n);
  for(int i = 0; i < n; i++) {
    all_samples(i) = i + 1;
  }
  
  uvec index1 = permutations.col(0);
  uvec index2 = setdiffC(all_samples, index1);
  
  mat scores1 = cor(counts.rows(index1 - 1));
  mat scores2 = cor(counts.rows(index2 - 1));
  
  vec d_gene_scores = d_genesC(scores1, scores2, lp);
  vec d_edge_scores = d_edgesC(scores1, scores2, lp);
  
  vec p_value_genes;
  p_value_genes.ones(n_genes);
  vec p_value_edges;
  p_value_edges.ones(n_edges);
  
  for(int j = 1; j < n_perm; j++) {
    index1 = permutations.col(j);
    index2 = setdiffC(all_samples, index1);
    
    scores1 = cor(counts.rows(index1 - 1));
    scores2 = cor(counts.rows(index2 - 1));
    
    p_value_genes = p_value_genes + (d_genesC(scores1, scores2, lp) >= d_gene_scores);
    p_value_edges = p_value_edges + (d_edgesC(scores1, scores2, lp) >= d_edge_scores);
  }
  
  p_value_genes = p_value_genes / (double)n_perm;
  p_value_edges = p_value_edges / (double)n_perm;
  
  return List::create(Named("p_value_genes") = p_value_genes,
                      Named("p_value_edges") = p_value_edges,
                      Named("d_genes") = d_gene_scores,
                      Named("d_edges") = d_edge_scores,
                      Named("lp") = lp);
}