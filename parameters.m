% Parameters Settting

% PSF estimation parameters
maxROIsize = [1024 1024];
MSlevels = 4;

% parameters
PAR.verbose = 2; 


% data term weight gamma
PAR.gamma = 1e2;
PAR.Lp = 0.3;

% PSFs estimation
PAR.beta_h = 1e4*PAR.gamma;
PAR.alpha_h = 1e1*PAR.gamma;
PAR.centering_threshold = 20/255; 

% image estimation
PAR.beta_u = 1e0*PAR.gamma;
PAR.alpha_u = 1e-2*PAR.gamma; 

% non-blind image estimation (final E-step)

PAR.gamma_nonblind = 2e1*PAR.gamma;
PAR.beta_u_nonblind = 1e-2*PAR.gamma_nonblind;
PAR.Lp_nonblind = 1;

% number of iterations
PAR.maxiter_u = 10;
PAR.maxiter_h = 10;
PAR.maxiter = 10;
PAR.ccreltol = 1e-3;