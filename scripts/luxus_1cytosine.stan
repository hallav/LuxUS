data {
  //Third version of LuxUS, version for one cytosine. Final version: removed unnecessary lines.

  //Dimensions
  int<lower=1> n_replicates; // number of replicates for each cytosine, name is kept compact
  int<lower=1> n_predictors; // number of predictors

  //Experimental parameters (check dimensions! should be a matrix/longer vector as there are multiple cytosines)
  real<lower=0,upper=1> bsEff[n_replicates];
  real<lower=0,upper=1> bsBEff[n_replicates];
  real<lower=0,upper=1> seqErr[n_replicates];

  //The total number of reads and number of reads that are C
  int<lower=0> bsC[n_replicates];
  int<lower=0> bsTot[n_replicates];

  //Design matrix X consists of block matrices, each representing the X matrix for individual cytosine. Rest of the matrix is zeros.
  matrix[n_replicates, n_predictors] X;

  //Indicator vector for block effect: replicate
  //int Z_R[n_replicates];

  //Fixed variance parameter for the coefficients vector B
  real<lower=0> sigmaB2;

  // prior for sigmaE2
  real<lower=0> alpha;
  real<lower=0> beta;

  //prior for sigmaR2
  real<lower=0> alphaR;
  real<lower=0> betaR;

}
transformed data {

}
parameters {
  vector[n_predictors] B;
  vector[n_replicates] ranint_rep;
  real<lower=0> sigmaE2;
  vector[n_replicates] Y;
  real<lower=0> sigmaR2;
  
}

transformed parameters {

  vector<lower=0,upper=1>[n_replicates] theta;
  vector[n_replicates] Y_hat;
  

  for (i in 1:n_replicates){
    Y_hat[i] = dot_product(X[i,:], B) + ranint_rep[i];
  }


  //The softmax function must be changed into logit function (or something)
  for (i in 1:n_replicates){
    //same as sigmoid function
    theta[i] = inv_logit(Y[i]);
  }
  

}

model {

  sigmaE2 ~ gamma(alpha,beta);
  sigmaR2 ~ gamma(alphaR,betaR);

  B ~ multi_normal(rep_vector(0,n_predictors),diag_matrix(rep_vector(sigmaB2,n_predictors)));

  ranint_rep ~ multi_normal(rep_vector(0,n_replicates),diag_matrix(rep_vector(sigmaR2,n_replicates)));


  Y ~ multi_normal(Y_hat, diag_matrix(rep_vector(sigmaE2, n_replicates)));

  for (i in 1:n_replicates){
    bsC[i] ~ binomial(bsTot[i],
            theta[i] * ((1.0 - seqErr[i]) * (1.0 - bsEff[i]) + seqErr[i] * bsEff[i]) +
            (1-theta[i]) * ((1.0 - bsBEff[i]) * (1.0 - seqErr[i]) + seqErr[i] * bsBEff[i]));
  }
}
