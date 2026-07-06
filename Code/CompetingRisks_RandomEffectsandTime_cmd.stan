data {
  int<lower=1> T;                   // Number of time steps
  int<lower=1> N;                   // Number of data rows
  int<lower=1> P;                   // Number of predictors 
  array[N] int<lower=1> Time;             // Time interval (e.g., 1, 2, 3,...)
  int<lower=1> K;                   // Number of categories
  array[N] int<lower=1, upper = K> Event; // Response: 1=none, 2=metamorphosis, 3=death
  matrix[N,P] X;                    // Predictors metamorphosis
  int<lower=1> G;                   // Number of tanks
  array[N] int<lower=1, upper=G> tank;    // Tank index for each observation
  int<lower=1> E;                   // Number of Enclosures
  array[N] int<lower=1, upper=E> enclosure;// Enclosure index for each observation
}

parameters {
  vector[T] alpha;                 // Time-specific intercepts for metamorphosis
  vector[T] theta;                 // Time-specific intercepts for death
  real<lower = 0> sigma_alpha;       // Enclosure-level standard deviations
  real<lower = 0> sigma_theta;       // Enclosure-level standard deviation
   
  vector[P] beta;                  // Coefficients for predictors for metamorphosis
  vector[P] gamma;                 // Coefficients for predictors for death
  
  vector[G] tau;                   // Standardized tank-level effects on metamorphosis
  vector[G] phi;                   // Standardized tank-level effects on death
  real<lower = 0> sigma_tau;       // Tank-level standard deviations
  real<lower = 0> sigma_phi;       // Tank-level standard deviations
  
  vector[E] rho;                   // Standardized Enclosure-level effects on metamorphosis
  vector[E] psi;                   // Standardized Enclosure-level effects on death
  real<lower = 0> sigma_rho;       // Enclosure-level standard deviations
  real<lower = 0> sigma_psi;       // Enclosure-level standard deviations
}

transformed parameters {
  vector[N] lambda;                // Metamorphosis hazard rate
  vector[N] eta;                   // Death hazard rate
  matrix[N,K] pi;
  for (i in 1:N) {
    int t = Time[i];
    int g = tank[i];
    int en = enclosure[i];
    lambda[i] = exp(alpha[t] + X[i]* beta + tau[g] + rho[en]);
    eta[i] = exp(theta[t] + X[i]* gamma + phi[g] + psi[en]);
    pi[i,1] = exp(-(lambda[i]+eta[i]));
    pi[i,2] = lambda[i]/(lambda[i]+ eta[i])*(1-exp(-(lambda[i]+eta[i])));
    pi[i,3] = eta[i]/(lambda[i]+ eta[i])*(1-exp(-(lambda[i]+eta[i])));
  }
}

model {
  // Priors (modify these as appropriate)
  alpha ~ normal(0, sigma_alpha);
  theta ~ normal(0, sigma_theta);
  sigma_alpha ~ exponential(1);
  sigma_theta ~ exponential(1);  
    
  beta  ~ normal(0, 1);
  gamma ~ normal(0, 1);
  
  tau ~ normal(0, sigma_tau);
  phi ~ normal(0, sigma_phi);
  sigma_tau ~ exponential(1);
  sigma_phi ~ exponential(1);
  
  rho ~ normal(0, sigma_rho);
  psi ~ normal(0, sigma_psi);
  sigma_rho ~ exponential(1);
  sigma_psi ~ exponential(1);
  
  // Likelihood
  for (i in 1:N){
    Event[i] ~ categorical(pi[i,]');
  }
}

generated quantities {
  vector[N] log_lik;                // For LOO cross-validation
  array[N] int M_rep; 
  for (i in 1:N) {
    log_lik[i] = categorical_lpmf(Event[i] | pi[i,]');
    M_rep[i] = categorical_rng(pi[i,]'); 
  }
}
