function U = LagrangianApproach(G,H,PAR)
% fast nonblind deconv using augmented Lagrangian approach 


% number of iterations
maxiter = PAR.maxiter_u;

% alpha
alpha = PAR.alpha_u;

% rel. tolerace used as convergence criterion
ccreltol = PAR.ccreltol;

% use gamma if gamma_nonblind is not defined, otherwise use gamma_nonblind
if isfield(PAR,'gamma_nonblind')
    gamma = PAR.gamma_nonblind;
else
    gamma = PAR.gamma;
end


if isfield(PAR,'beta_u_nonblind')
    beta = PAR.beta_u_nonblind;
else
    beta = PAR.beta_u;
end


if isfield(PAR,'Lp_nonblind')
    Lp = PAR.Lp_nonblind;
else
    Lp = PAR.Lp;
end

% size of image U
usize = size(G);
usize(3) = size(G,3);

% vrange ... range of intensity values in each color channel
vrange = zeros(usize(3),2);
for c=1:usize(3)
    vrange(c,:) = [min(reshape(G(:,:,c),[],1)), max(reshape(G(:,:,c),[],1))];
end

% If we work with FFT, we have to move H center into the origin
hshift = zeros(size(H));
hshift(floor(size(H,1)/2)+1, floor(size(H,2)/2)+1) = 1;

% FU ... FFT of u
FU = 0;

% FDx, FDx ... FFT of x and y derivative operators
FDx = repmat(fft2([1 -1],usize(1),usize(2)),[1 1 usize(3)]);
FDy = repmat(fft2([1; -1],usize(1),usize(2)),[1 1 usize(3)]);

% FH ... FFT of PSF
FH = repmat(conj(fft2(hshift,usize(1),usize(2))).*fft2(H,usize(1),usize(2)), [1 1 usize(3)]); 
FHTH = conj(FH).*FH; % FT of (H^T)H (RGB)

% FGs ... FFT of H^T*g
% Note that we use edgetaper to reduce border effect
eG = edgetaper(G,H);
FGu = fft2(eG);
FGs = conj(FH).*FGu;

DTD = conj(FDx).*FDx + conj(FDy).*FDy;

if PAR.verbose > 1

else

end

% extra variables for Bregman iterations
Bx = zeros(usize);
By = zeros(usize);
Vx = zeros(usize);
Vy = zeros(usize);

% main iteration loop, do everything in the FT domain
for i = 1:maxiter
    if PAR.verbose
        disp(['nonblind deconv step ',num2str(i)]);
    end
    
    FUp = FU;
    b = FGs + beta/gamma*(conj(FDx).*fft2(Vx+Bx) + conj(FDy).*fft2(Vy+By));
	FU = b./(FHTH + beta/gamma*DTD);
    
    % Prepare my Lp prior
    Pr = PriorNorm(Lp,alpha,beta);
    xD = real(ifft2(FDx.*FU));
    yD = real(ifft2(FDy.*FU));
    xDm = xD - Bx;
    yDm = yD - By;
    nDm = repmat(sqrt(sum(xDm.^2,3) + sum(yDm.^2,3)),[1 1 usize(3)]);
    Vy = Pr.fh(yDm,nDm);
    Vx = Pr.fh(xDm,nDm);
    
    % update Bregman variables
    Bx = Bx + Vx - xD;
    By = By + Vy - yD;
   
    if PAR.verbose>1
     
        U = real(ifft2(FU));
        % impose constraints on U 
        U = uConstr(U,vrange);
        E = sqrt(Vy.^2+Vx.^2);
    end
    
  
    relcon = sqrt(sum(abs(FUp(:)-FU(:)).^2))/sqrt(sum(abs(FU(:)).^2));
    if PAR.verbose
        disp(['relcon:',num2str([relcon])]);
    end
    if relcon < ccreltol
        break;
    end
end

if PAR.verbose<2
    U = real(ifft2(FU));
    U = uConstr(U,vrange);
end
end

function newU = uConstr(U,vrange)


newU = U;
for c = 1:size(U,3)
	m = false(size(U));
	m(:,:,c) = U(:,:,c)<vrange(c,1);
	newU(m) = vrange(c,1);
	m(:,:,c) = U(:,:,c)>vrange(c,2);
	newU(m) = vrange(c,2);
end
end