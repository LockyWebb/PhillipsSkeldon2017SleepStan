//
// This Stan program is the  Phillips sleep model from 
// Skeldon 2017
// Based on code provided by Andrew Phillips  in October 2020
// using cmdstan

//vector phillipssleep07(real Q_max, real theta, real sigma, real nu_ma_Q_ao, real nu_vm, real nu_mv, 
//real nu_vc, real nu_vh, real chi, real mu, real tau_m, real tau_v, real c_0, real omega, real alpha)



// 
// function dY = odefun(t, Y, Q_max, theta, sigma, nu_ma_Q_ao, nu_vm, nu_mv, nu_vc, nu_vh, chi, mu, tau_m, tau_v, c_0, omega, alpha)
// %function dY = odefun(t, Y)
// 	
// 	%dY(1) = (nu_vm*(Q_max/(1 + exp(-(Y(2)-theta)/sigma))) + nu_vc*(c_0 + cos(2*pi/24*t + alpha)) + nu_vh*Y(3) - Y(1))/(tau_v/3600);
// 	dY(1) = (nu_vm*Sigm(Y(2),Q_max,theta,sigma) + nu_vc*C_drv(t,c_0, omega,alpha) + nu_vh*Y(3) - Y(1))/(tau_v/3600);
// 	
// 	%dY(2) = (nu_ma_Q_ao + nu_mv*(Q_max/(1 + exp(-(Y(1)-theta)/sigma))) - Y(2))/(tau_m/3600);
// 	dY(2) = (nu_ma_Q_ao + nu_mv*Sigm(Y(1),Q_max,theta,sigma) - Y(2))/(tau_m/3600);
// 
// 	%dY(3) = (mu*(Q_max/(1 + exp(-(Y(2)-theta)/sigma))) - Y(3))/chi;
// 	dY(3) = (mu*Sigm(Y(2),Q_max,theta,sigma) - Y(3))/chi;
// 	
// 	dY = dY(:);
// end

functions {
  
  real Sigm(real x,
                real Q_max,
                real theta,
                real sigma) {
    real Q = (Q_max/(1 + exp(-(x-theta)/sigma)));
    return Q;
  }
  
  real C_drv(real t,
                 real c_0, 
                 real omega,
                 real alpha) {
    real C = c_0 + cos(omega*t + alpha);
    return C;
  }
  
  vector phillipssleep07(real t, 
                         vector Y, 
                         real Q_max, 
                         real theta, 
                         real sigma, 
                         real nu_ma_Q_ao, 
                         real nu_vm, 
                         real nu_mv,
                         real nu_vc, 
                         real nu_vh, 
                         real chi, 
                         real mu, 
                         real tau_m, 
                         real tau_v, 
                         real c_0, 
                         real omega, 
                         real alpha) {
    vector[3] dydt;
    
    dydt[1] = (nu_vm*Sigm(Y[2],Q_max,theta,sigma) + nu_vc*C_drv(t,c_0, omega,alpha) + nu_vh*Y[3] - Y[1])/(tau_v/3600); // Vv
	
    dydt[2] = (nu_ma_Q_ao + nu_mv*Sigm(Y[1],Q_max,theta,sigma) - Y[2])/(tau_m/3600); // Vm

  	dydt[3] = (mu*Sigm(Y[2],Q_max,theta,sigma) - Y[3])/chi; // H
  	
  	return dydt;
    
  }
  
  real get_state(vector Y) {
    real state;
    state = 1.0*(Y[2] > Y[1]);
    //state = to_array_1d(state);
    return state;
  }

  
  
}

data {
  int<lower = 1>    T_n;       // number of points
  real<lower = 0>  ts[T_n];   // times values
  vector[3]         y0;      // inital value
  //real              Y[T_n];    // function values
  //real            
  //real             state_obs[T_n];
  real             prop_obs;
}

transformed data{
  real t0 = 0;      // inital time 
  real Q_max = 100;
  real theta = 10; 
  real sigma = 3; 
  real nu_ma_Q_ao = 1.3;
  real nu_vm = -2.1;
  real nu_mv = -1.8;
  real nu_vc = -2.9; 
  real nu_vh = 1.0; 
  //real chi = 45; 
  real mu = 4.4; 
  real tau_m = 10; 
  real tau_v = 10; 
  real c_0 = 4.5; 
  real omega = 0.2617994; // 2*pi/24
  real alpha = 0;
}

parameters {
 real chi;
}

transformed parameters{
  vector[3] Y[T_n] = ode_rk45(phillipssleep07, y0, t0, ts,
                            Q_max, theta, sigma, nu_ma_Q_ao, nu_vm, 
                            nu_mv, nu_vc, nu_vh, chi, mu,  tau_m, 
                            tau_v, c_0, omega, alpha);
                            
  
  real state[T_n];
  for (n in 1:T_n){
    if (Y[n,2]>Y[n,1])
      state[n] = 1;
    else 
      state[n] = 0;
  }
    //state[n] = 1.0*(Y[n,2]>Y[n,1]); // 1/TRUE is wake, 0/FALSE is sleep
    
  real prop_sleep = sum(state)/T_n;
    
}

model {
  //chi ~ normal(45,1);
  chi ~ cauchy(40,2);
  
  //for (n in 1:T_n){
  //  state_obs[n] ~ normal(state[n], 0.01);
  //}
  
  prop_obs ~ normal(prop_sleep, 0.01);
   
  

}

generated quantities {
//real state[T_n];
//  for (n in 1:T_n){
//    if (Y[n,2]>Y[n,1])
//      state[n] = 1;
//    else 
//      state[n] = 0;
//  }
}