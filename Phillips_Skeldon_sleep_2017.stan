//
// This Stan program is the  Phillips sleep model from 
// Skeldon 2017
// Based on code provided by Andrew Phillips  in October 2020
// using cmdstan

//vector phillipssleep07(real Q_max, real theta, real sigma, real nu_ma_Q_ao, real nu_vm, real nu_mv, 
//real nu_vc, real nu_vh, real chi, real mu, real tau_m, real tau_v, real c_0, real omega, real alpha)




// 	
// dy = [p{1}*(-y(1) + (wake_func(t)==0).*(p{2}*sigmoid(y(2)) + p{3}) + (wake_func_LW(t)==1).*(Dv>2.4).*(wake_effort(Dv)) + (wake_func_LW(t)==1).*(Dv<=2.4).*(0.5*(1+tanh(-1-y(1))).*1 + (p{2}*sigmoid(y(2)) + p{3})));...
//     p{1}*(-y(2) + p{4}*sigmoid(y(1)) + p{5}*y(3) + p{6}*C + p{22}*((C-0.5).^2) + p{7});...
//     p{8}*(-y(3) + p{9}*sigmoid(y(1)))
//     p{15}*(alpha*(1-y(4))-p{16}*y(4));
//     p{17}*(p{18}*(y(5)-p{19}*y(5).^3)-y(6)*(p{20}+p{21}*B));
//     p{17}*(y(5)+B)];

// p{1} = sph/10; % 1/tau_v=1/tau_m
// p{2} = -1.8/sph; % numv
// p{3} = 1.3; % A
// p{4} = -2.1/sph; % nuvm
// p{5} = 1.0; % nuvh
// p{6} = -2;%-3.37; %-1.1; %-3.37; % nuvc
// p{7} = -10.2; %-3;%-10.2; % D_0
// p{8} = 1/45; %1/20; %1/45; % 1/chi
// p{9} = 4.20/sph; %5/sph; %4.20/sph; % mu

// % light processing parameters
// p{10} = 0.16; % alpha0
// p{11} = 9500; % I0
// p{12} = 0.6; % p
// p{13} = 0.4; % b
// p{14} = 19.9; % G
// p{15} = 60; % lambda
// p{16} = 0.013; % beta
// 
// % circadian pacemaker parameters
// tauc = 24.2; % Define intrinsic circadian period (h)
// p{17} = pi/12; %1/kappa
// p{18} = 0.23; %gamma
// p{19} = 4/3; % coefficient for x^3 term
// p{20} = (24/(0.99669*tauc))^2; % tau term
// p{21} = 0.55; % k
// 
// % circadian squared parameter
// p{22} = -5.5; %nu_c2


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
  
  real light_func(real t) {
    real dshift = 0; // Amount by which to shift the solar curve (hours)
    real evelight = 40; // Set level (lux) of evening light (home lighting during wake, outside solar day) 

    real dawn=(8+dshift)*3600;     //Clock time when light is approx (day+evening)/2
    real dusk=(17+dshift)*3600;    //Clock time when light is approx (day+evening)/2
    real daylight=500; // Max daylight level (lux)
    real steepness=1/6000; // Steepness of solar curve
    real a=evelight;
    real b=(daylight-evelight)/2;
    // Combine two sigmoid curves to approximate the solar curve
    real time = mod(t,24);
    //lightf=@(a,b,dawn,dusk,time)a+b*(tanh(steepness*(time*3600-dawn))-tanh(steepness*(time*3600-dusk)));
    real L = a+b*(tanh(steepness*(time*3600-dawn))-tanh(steepness*(time*3600-dusk)));
    //L = lightf(a,b,dawn,dusk,mod(t,24));
  }
  
  vector phillipssleep17(real t, 
                         vector Y, 
                         real Q_max, 
                         real theta, 
                         real sigma, 
                         real nu_ma_Q_ao, // also called A
                         real nu_vm, 
                         real nu_mv,
                         real nu_vc, 
                         real nu_vh, 
                         real D_0,
                         real chi, 
                         real mu, 
                         real tau_m, 
                         real tau_v, 
                         //real c_0, 
                         //real omega, 
                         //real alpha  // alpha is a different parameter in this model
                         
                         real alpha0,
                         real I0,
                         real p,
                         real b,
                         real G,
                         real lambda,
                         real beta,
                         
                         real kappa,
                         real gamma,
                         real x3coef,
                         real taut
                         real k
                         
                         real nu_c2
                         
                         ) {
    real C_drv = 0.5 * (1 + 0.8 * Y[6] - 0.47 * Y[5]);
    //real C_drv = C_drv(t,c_0, omega,alpha); // real c_0 = 4.5; real omega = 0.2617994; // 2*pi/24 real alpha = 0;
    
    real Drv = nu_vh * Y[3] + nu_vc * C_drv + nu_c2 * ((C_drv-0.5).^2) + D_0;
    
    real I = light_func(t) .* (Y[1]>Y[2]);
    
    real alpha = alpha0 * (I/I0).^p; 
    
    real B = (1 - b*Y[5]).*(1 - b*Y[6])*G*alpha*(1 - Y[4]);
                           
                           
    vector[6] dydt;
    
    // 2007 has Y[1] being Vv, and this one has Y[1] as Vm
    
    dydt[1] = (nu_ma_Q_ao + nu_mv*Sigm(Y[2],Q_max,theta,sigma) - Y[1])/(tau_m/3600); // Vm
    
    dydt[2] = (nu_vm*Sigm(Y[1],Q_max,theta,sigma) + Drv - Y[2])/(tau_v/3600); // Vv
	
  	dydt[3] = (mu*Sigm(Y[1],Q_max,theta,sigma) - Y[3])/chi; // H
  	
  	dydt[4] = lambda * (alpha * (1 - Y[4]) - beta *  Y[4]); // 
  	
  	dydt[5] = kappainv * (gamma * (Y[5] - (x3coef * (Y[5].^3)) ) - Y[6] * (taut + k * B) ); //
  	
  	dydt[6] = kappainv * (Y[5] + B)
  	
  	return dydt;
    
  }
  
  real get_state(vector Y) {
    real state;
    state = 1.0*(Y[1] > Y[2]);
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
  real D_0 = -10.2;
  //real chi = 45; 
  real mu = 4.4; 
  real tau_m = 10; 
  real tau_v = 10; 
  //real c_0 = 4.5; 
  //real omega = 0.2617994; // 2*pi/24
  //real alpha = 0;
  
  real alpha0 = 0.16;
  real I0  =  9500;
  real p = 0.6;
  real b = 0.4;
  real G = 19.9;
  real lambda = 60;
  real beta = 0.013,
                         
  real kappa = pi/12;
  real gamma = 0.23;
  real x3coef = 4/3;
  real taut = (24/(0.99669*24.2))^2;
  real k = 0.55;
                         
  real nu_c2 = -3.5; // make 0 to remove squared term
                         
  
  
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