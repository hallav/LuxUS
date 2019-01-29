data {

  int<lower=1> n_cytosines;
  int<lower=1> n_replicates;
  int<lower=1> n_predictors;

  real<lower=0,upper=1> bsEff[n_replicates * n_cytosines];
  real<lower=0,upper=1> bsBEff[n_replicates * n_cytosines];
  real<lower=0,upper=1> seqErr[n_replicates * n_cytosines];

  int<lower=0> bsC[n_replicates * n_cytosines];
  int<lower=0> bsTot[n_replicates * n_cytosines];

  matrix[n_cytosines * n_replicates, n_predictors] X;

  int Z_R[n_replicates * n_cytosines];
  int Z_C[n_replicates * n_cytosines];

  real<lower=0> sigmaB2;

  real<lower=0> alpha;
  real<lower=0> beta;

  real<lower=0> alphaR;
  real<lower=0> betaR;

  real<lower=0> alphaC;
  real<lower=0> betaC; 

  real<lower=0> alphal;
  real<lower=0> betal;

  real<lower=0> coordinates[n_cytosines];

}
parameters {

  vector[n_predictors] B;
  vector[n_cytosines] ranint_cyt;
  vector[n_replicates] ranint_rep;
  real<lower=0> sigmaE2;
  vector[n_cytosines * n_replicates] Y;
  real<lower=0> l;
  real<lower=0> sigmaC2;
  real<lower=0> sigmaR2;
  
}

transformed parameters {

  vector<lower=0,upper=1>[n_cytosines * n_replicates] theta;
  cov_matrix[n_cytosines] SigmaC;
  vector[n_cytosines * n_replicates] Y_hat;
  
  for (i in 1:(n_cytosines * n_replicates)){
    theta[i] = inv_logit(Y[i]);
  }

  for (i2 in 1:(n_cytosines)){
    SigmaC[i2,i2] = 1;
  }

  for (i in 1:n_cytosines){
    for (j in 1:n_cytosines){
      if (i != j) {
        SigmaC[i,j] = exp(-(fabs(coordinates[j]-coordinates[i])) / pow(l,2));
      }
    }
  }

  for (i in 1:n_cytosines * n_replicates){
    Y_hat[i] = dot_product(X[i,:], B) + ranint_cyt[Z_C[i]] + ranint_rep[Z_R[i]];
  }
  
}

model {

  l ~ gamma(alphal,betal);

  sigmaE2 ~ gamma(alpha,beta);
  sigmaR2 ~ gamma(alphaR,betaR);
  sigmaC2 ~ gamma(alphaC,betaC);

  B ~ multi_normal(rep_vector(0,n_predictors),diag_matrix(rep_vector(sigmaB2,n_predictors)));

  ranint_cyt ~ multi_normal(rep_vector(0,n_cytosines), sigmaC2 * SigmaC);
  ranint_rep ~ multi_normal(rep_vector(0,n_replicates),diag_matrix(rep_vector(sigmaR2,n_replicates)));

  Y ~ multi_normal(Y_hat, diag_matrix(rep_vector(sigmaE2, n_cytosines * n_replicates)));

  for (i in 1:n_cytosines * n_replicates){
    bsC[i] ~ binomial(bsTot[i],
            theta[i] * ((1.0 - seqErr[i]) * (1.0 - bsEff[i]) + seqErr[i] * bsEff[i]) +
            (1-theta[i]) * ((1.0 - bsBEff[i]) * (1.0 - seqErr[i]) + seqErr[i] * bsBEff[i]));
  }
}
