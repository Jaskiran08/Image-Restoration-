function [H, U, Report] = ERSteps(G, iH, PAR, Hstar)
% PSFestim


Report = [];

gamma = PAR.gamma;
Lp = PAR.Lp;
ccreltol = PAR.ccreltol;

% size of H
hsize = [size(iH,1) size(iH,2)];
gsize = [size(G,1), size(G,2)];
usize = gsize;

if PAR.verbose > 1

end

% if true PSF Hstar is provided -> calculate MSE
if exist('Hstar','var') && ~isempty(Hstar)
    doMSE = 1;
    Report.hstep.mse =  zeros(1,PAR.maxiter+1);
else
    doMSE = 0;
end

U = zeros(usize);
H = iH;

% Initialization of variables for min_U step, which do not change
% FU ... FFT of u
FU = fft2(U);
% FDx, FDx ... FFT of x and y derivative operators
FDx = fft2([1 -1],usize(1),usize(2));
FDy = fft2([1; -1],usize(1),usize(2));
DTD = conj(FDx).*FDx + conj(FDy).*FDy;

% FUx, FUy ... FFT of x (y) derivatives of U
FUx = zeros(usize);
FUy = zeros(usize);

% auxiliary variables for image gradient and blurs
% initialize to zeros
Vx = zeros(usize);
Vy = zeros(usize);
Vh = zeros([usize]);
% extra variables for Bregman iterations
Bx = zeros(usize);
By = zeros(usize);
Bh = zeros([usize]);

if doMSE
	Report.hstep.mse(1) = calculateMSE(H,Hstar);
end

eG = edgetaper(G,ones(hsize)/prod(hsize));

FeGu = fft2(eG);
FeGx = FDx.*FeGu;
FeGy = FDy.*FeGu;

% main loop which alternates between E-estimation and R-estimation
for mI = 1:PAR.maxiter
    % E-estimation
    Estep;
   
	% R-estimation
    Rstep;
	
	% reporting
	if doMSE
        Report.hstep.mse(mI+1) = calculateMSE(H,Hstar);
	end
	
    % increasing gamma helps
    gamma = gamma*1.5;
end

% PSF centering
if(PAR.centering_threshold > 0)
	H = centerPSF(H, PAR.centering_threshold);
end

% =========================
% E STEP
% =========================
function Estep
	FHS = fft2(H,usize(1),usize(2)); 
	FHTH = conj(FHS).*FHS; 
	FGs = sum(conj(FHS).*FeGu,3); 

	beta = PAR.beta_u;
	alpha = PAR.alpha_u;

	% main iteration loop, do everything in the FT domain
	for i = 1:PAR.maxiter_u
		FUp = FU;
		b = FGs + beta/gamma*(conj(FDx).*fft2(Vx+Bx) + conj(FDy).*fft2(Vy+By));
		FU = b./( FHTH + beta/gamma*DTD);

		% Prepare my Lp prior
		Pr = PriorNorm(Lp,alpha,beta);
		FUx = FDx.*FU;
		FUy = FDy.*FU;
		xD = real(ifft2(FUx));
		yD = real(ifft2(FUy));
		xDm = xD - Bx;
		yDm = yD - By;
		nDm = sqrt(xDm.^2 + yDm.^2);
		Vy = Pr.fh(yDm,nDm);
		Vx = Pr.fh(xDm,nDm);
		% update Bregman variables
		Bx = Bx + Vx - xD;
		By = By + Vy - yD;

		E = sqrt(Vy.^2+Vx.^2);

		% Calculate relative convergence criterion
		relcon = sqrt(sum(abs(FUp(:)-FU(:)).^2))/sqrt(sum(abs(FU(:)).^2));

		if relcon < ccreltol
			break;
		end
	end

	if PAR.verbose
		disp(['Steps E: ',num2str(i),' relcon:',num2str([relcon])]);
	end
	U = real(ifft2(FU));

end % end of Estep

% =======================
% R STEP
% =======================
function Rstep
	FUD = FeGx.*conj(FUx) + FeGy.*conj(FUy); 
	FUTU = conj(FUx).*FUx + conj(FUy).*FUy; 
	FH = fft2(H,usize(1),usize(2));
	
	beta = PAR.beta_h;
	alpha = PAR.alpha_h;
	
	% main loop
	for i = 1:PAR.maxiter_h
		FHp = FH; 
		b = beta/gamma*fft2(Vh+Bh) + FUD;
		FH = b./(FUTU + beta/gamma);
		
		% Calculate relative convergence criterion
		relcon = sqrt(sum(abs(FHp(:)-FH(:)).^2))/sqrt(sum(abs(FH(:)).^2));
		
		Pr = PriorNorm(1,alpha,beta);
		hI = real(ifft2(FH));
		hIm = hI - Bh;
		nIm = abs(hIm);
		Vh = Pr.fh(hIm,nIm); 
		Vh(Vh<0) = 0;
		
		Vh(hsize(1)+1:end,:,:) = 0; Vh(1:hsize(1),hsize(2)+1:end,:) = 0;
		% update Bregman variables
		Bh = Bh + Vh - hI;

		H = hI(1:hsize(1),1:hsize(2),:); % new H estimate

		E = abs(Vh);

		% convergence test
		if relcon < ccreltol
			break;
		end
	end

	if PAR.verbose
		disp([' Step R',num2str(i),' relcon:',num2str([relcon])]);
	end
end % end of Rstep
end

function r = calculateMSE(h,hs)
    hsize = size(hs);
    i = size(h)-hsize+1;
	h = h/sum(h(:))*sum(hs(:));
	R = im2col(h,[size(hs,1) size(hs,2)],'sliding');
    
    s = sqrt(sum((R-repmat(reshape(hs,prod(hsize(1:2)),1),1,prod(i(1:2)))).^2,1));
    r = s(ceil(prod(i(1:2))/2));
end

function [H] = centerPSF(H, thresh)

hsize = [size(H,1) size(H,2)];

for i=1:size(H,3)
	h = mat2gray(H(:,:,i)); 
	
	% nonzero mask
	m = h >= thresh;
	m2 = bwmorph(m, 'clean');
	if(any(m2(:))) m = m2; end 
	
	% determine mask support
	sum1 = sum(m, 1);
	sum2 = sum(m, 2);
	L = [find(sum2, 1, 'first') find(sum1, 1, 'first')];
	R = [find(sum2, 1, 'last') find(sum1, 1, 'last')];
	topleft = fix((L+R+1-hsize)/2); 
	
	% indexing (=shifting)
	h = h(max(topleft(1),1):min(topleft(1)+hsize(1)-1,end), max(topleft(2),1):min(topleft(2)+hsize(2)-1,end)); % get the 'existing' data, then pad borders with zeros
	h = padarray(h, max(topleft-[1,1],0), 0, 'post'); % pad with zeros to end up with the same size
	h = padarray(h, max([1,1]-topleft,0), 0, 'pre');
	
	% normalize sum to 1
	H(:,:,i) = h/sum(h(:));
end
end

