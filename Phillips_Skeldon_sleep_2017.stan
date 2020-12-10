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
    real time = fmod(t,24);
    //lightf=@(a,b,dawn,dusk,time)a+b*(tanh(steepness*(time*3600-dawn))-tanh(steepness*(time*3600-dusk)));
    real L = a+b*(tanh(steepness*(time*3600-dawn))-tanh(steepness*(time*3600-dusk)));
    //L = lightf(a,b,dawn,dusk,mod(t,24));
    return L;
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
                         
                         real kappainv,
                         real gamma,
                         real x3coef,
                         real taut,
                         real k,
                         
                         real nu_c2
                         
                         ) {
    real C = 0.5 * (1 + 0.8 * Y[6] - 0.47 * Y[5]);
    //real C = C_drv(t,c_0, omega,alpha); // real c_0 = 4.5; real omega = 0.2617994; // 2*pi/24 real alpha = 0;
    
    real Drv = nu_vh * Y[3] + nu_vc * C + nu_c2 * ((C-0.5).^2) + D_0;
    
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
  	
  	dydt[6] = kappainv * (Y[5] + B);
  	
  	return dydt;
    
  }
  
  real get_state(vector Y) {
    real state;
    state = 1.0*(Y[1] > Y[2]);
    //state = to_array_1d(state);
    return state;
  }

  int num_matches(real[] x, int a) {
    int n = 0;
    for (i in 1:size(x))
      if (x[i] == a)
        n += 1;
    return n;
  }
  
  int[] which_equal(int[] x, int a, int n) {
    int match_positions[n];
    int pos = 1;
    for (i in 1:size(x)) {
      if (x[i] == a) {
        match_positions[pos] = x[i];
        pos += 1;
      }
    }
    return match_positions;
  }
  
  
}

data {
  int<lower = 1>    T_n;       // number of points
  real<lower = 0>  ts[T_n];   // times values
  vector[6]         y0;      // inital value
  //real              Y[T_n];    // function values
  //real            
  //real             state_obs[T_n];
//  real             prop_obs;
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
  real beta = 0.013;
                         
  real kappainv = pi()/12;
  real gamma = 0.23;
  real x3coef = 4/3;
  real taut = (24/(0.99669*24.2))^2;
  real k = 0.55;
                         
  real nu_c2 = -3.5; // make 0 to remove squared term
                         
  
  int<lower = 1> max_bout = 10;
}

parameters {
 real chi;
}

transformed parameters{
  vector[6] Y[T_n] = ode_rk45(phillipssleep17, y0, t0, ts,
                            Q_max, theta, sigma, nu_ma_Q_ao, nu_vm, 
                            nu_mv, nu_vc, nu_vh, D_0, chi, mu,  tau_m,
                            tau_v, alpha0, I0, p, b, G, lambda, beta,
                            kappainv, gamma, x3coef, taut, k, nu_c2);
                            
  
  real state[T_n];
  
  //real sleep_ontemp;
  //real sleep_offtemp;
  
  for (n in 1:T_n){
    if (Y[n,1]>Y[n,2])
      state[n] = 1;
    else 
      state[n] = 0;
    
  //  if (n >= 2)
  //   statedif[n-1] = state[n] - state[n-1];
     //sleep_ontemp[n-1] = (statedif[n-1] == -1);
     //sleep_offtemp[n-1] = (statedif[n-1] == 1);
  }
    //state[n] = 1.0*(Y[n,2]>Y[n,1]); // 1/TRUE is wake, 0/FALSE is sleep
    
  real prop_sleep = sum(state)/T_n;
  
  // real sleep_ontemp = ts[statedif == -1];
  // real sleep_offtemp = ts[statedif == 1];
  // 
  // //if (sleep_ontemp[T_n-1] > sleep_offtemp[T_n-1])
  // //  sleep_on = sleep_ontemp[1:(Y_n-2)];
  //   
  // //if (sleep_offtemp[1] < sleep_ontemp[1])
  // //  sleep_off = sleep_offtemp[2:(T_n-1)];
  //   
  // real midsleep = 0.5 * (sleep_ontemp + sleep_offtemp);
  // real midsleepclk = fmod(midsleep, 24);
  // 
  // real duration = sleep_off-sleep_on;  
    
}

model {
  chi ~ normal(45,1);
  //chi ~ cauchy(40,2);
  
  //for (n in 1:T_n){
  //  state_obs[n] ~ normal(state[n], 0.01);
  //}
  
//  prop_obs ~ normal(prop_sleep, 0.01);
   
  

}

generated quantities {
  
  
  // int statedif[T_n-1];
  // int idx_on[T_n-1];
  // int idx_off[T_n-1];
  // for (n in 2:T_n){
  //   if ((state[n] - state[n-1]) == 0){
  //     statedif[n-1] = 0;
  //     idx_on[T_n-1] = 0;
  //     idx_off[T_n-1] = 0;
  //   }
  //   else if ((state[n] - state[n-1]) == -1){
  //     statedif[n-1] = -1;
  //     idx_on[n-1] = 1;
  //     idx_off[T_n-1] = 0;
  //   }
  //   else if ((state[n] - state[n-1]) == 1){
  //     statedif[n-1] = 1;
  //     idx_on[n-1] = 0;
  //     idx_off[T_n-1] = 1;
  //   }
  //   
  //  
  // }
  
//   //int num_on = sum(idx_on);
//   // can't make an integer array based on number of events
//   // make it a laceholder of 8 (four per  day id run for two days)
//   int idx_on_int[4];
//   int nnon = 1;
//   for(n_on in 1:(T_n-1)){
//     if (idx_on[n_on] == 1){
//       idx_on_int[nnon] = n_on;
//       nnon += 1;
//     }
//   }
// 
//   //int num_off = sum(idx_off);
//   int idx_off_int[4];
//   int nnoff = 1;
//   for(n_off in 1:(T_n-1)){
//     if (idx_off[n_off] == 1){
//       idx_off_int[nnoff] = n_off;
//       nnoff += 1;
//     }
//   }
// 
// 
// 
// //  for(nn in 1:(T_n - 1)){
// //    idx_on[nn] = statedif[nn] == -1;
// //    idx_off[nn] = statedif[nn] == 1;
// //  }
// 
// //  int num_sleep_on = num_matches(statedif,-1);
// 
// //vector[num_sleep_on] idx_on = which_equal(statedif, -1);
// 
// 
//   real sleep_ontemp[8] = ts[idx_on_int];
//   real sleep_offtemp[8] = ts[idx_off_int];
// 
//   real sleep_on[8];
//   real sleep_off[8];
// 
//   if (sleep_ontemp[8] > sleep_offtemp[8])
//     sleep_on = sleep_ontemp[1:(7)];
// 
//   if (sleep_offtemp[1] < sleep_ontemp[1])
//     sleep_off = sleep_offtemp[2:(8)];
// 
// //  real midsleep = 0.5 * (sleep_ontemp + sleep_offtemp);
// //  real midsleepclk = fmod(midsleep, 24);
// 
// //  real duration = sleep_off-sleep_on;

// try a more smiple way

  int nn_on = 1;
  int nn_off = 1;
  
  real sleep_on[max_bout];
  real sleep_off[max_bout];
  
  for (n in 1:(T_n-1)){
    if ((state[n] - state[n-1]) == -1){
      sleep_on[nn_on] = ts[n];
      nn_on += 1;
    }
    else if ((state[n] - state[n-1]) == 1){
      sleep_off[nn_off] = ts[n];
      nn_off += 1;
    }
  }
  
  real sleep_on_use[max_bout-1];
  real sleep_off_use[max_bout-1];
  
  sleep_on_use = sleep_on[1:(max_bout-1)];
  
  if(sleep_off[1] < sleep_on[1]){
    sleep_off_use = sleep_off[2:];
  }else{
    sleep_off_use = sleep_off[1:(max_bout-1)];
  }
  
  
  real midsleep[(max_bout-1)];
  real sleep_on_clk[(max_bout-1)];
  real sleep_off_clk[(max_bout-1)];
  real midsleep_clk[(max_bout-1)];
  
  for (ii in 1:(max_bout-1)){
    if (!is_nan(sleep_on_use[ii])){
      midsleep[ii] = 0.5 * (sleep_on_use[ii] - sleep_off_use[ii]);
      sleep_on_clk[ii] = fmod(sleep_on_use[ii],24);
      sleep_off_clk[ii] = fmod(sleep_off_use[ii],24);
      midsleep_clk[ii] = fmod(midsleep[ii],24);
    }
    
  }





}