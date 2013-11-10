%	10/31/2008: modified line 130


%This file is for computing the self-organizing maps of the V1. We are 
%	interested in the relative relationship between five maps:retinotopic map, 
%	orientation, ocular dominance, and spatial freuqency maps. 
%	We use the Kohonen self-organizing map algorithm for the map formation.
%	Reference: K.Obermayer et al. Phys Rev A, 45, 7568 (1992).

%rand('state',sum(100*clock));
rand('state', 22);

%---WEIGHTS for features
xWeight = 1;
yWeight = 1;
qrWeight = 1.3;
qsWeight = 1.3;
zWeight = 1.3;
fWeight = 1.3;

%---flags-----------------------------------------

iContinue = 0;		% if 1, continue on the previous simulation.
iPlot = 0;		% if 1, plot on the way. otherwise no plots.
iFixQ = 0;		% if 1, orientation tuning strength is fixed and not sampled.
iSaveIntermediateFiles=1;	% if 1, 
iInterFile = 1

if iContinue == 1 %%%NOTE: SOME HARD CODE in else block line 103-110, pay attention to runFlag, ifixQ,etc
	fileContinue = 'runData/mapInter-11.mat';
	iInterFile = 11;
end

if iContinue ~= 1

iStart = 1;

%---parameters------------------------------------

N = 513;		% size of the cortical grid.
D = 6;			% dimenstion of the features (x,y,z,f,phi), (x,y) retinotopic, phi orientation, z oclular dominance, f spatial frequency.
ratioXY = 4;		% ratio of magnification in x and y directions
iVPosition =1;		% the positions of the plots.
d1 = N*ratioXY;		% range of the retinotopic map in x direction.
d2 = N;			% range of the retinotopic map in y direction.
qPat  = 40;		% the range of the orientation tuning stregth.
zPat  = 60;		% [0 zPat] is the range of the ocular dominance.
fPat  = 60;		% [0 fPat] is the range of the frequency tuning.

alpha  = 0.02;		% learning rate.
sigma1 = 5;		% the Gaussian parameter of the neighbourhood function in direction 1. 
sigma2 = 5;		% the Gaussian parameter of the neighbourhood function in direction 2.
maxIter = 6000000;	% maximum number of iterations.


runFlag = 0;		% 0, default runs.
			% 1, fix Q, and change sigma1 and sigma2 from large to small
			% 2, round visual field to oval cortical field.
iSampleFixQ = 0;

if runFlag == 1
	sigmaStart = 25;
	stepScale = 100000;	% sigma is scaled down according to sigma1,2 = sigma1,2/2^exponent at each stepScale.
	iSampleFixQ = 1;
	sigma1 = sigmaStart;
	sigma2 = sigmaStart;
	exponent = log(sigmaStart)/log(2)/(maxIter/stepScale-1);
end

%---analyze the stability criterion-----------------------

Tth1 = 1/2*sqrt(exp(1))*d1/N*min(sigma1,sigma2);
Tth2 = 1/2*sqrt(exp(1))*d2/N*min(sigma1,sigma2);
T3 = qPat/2;
T4 = qPat/2;
T5 = zPat/sqrt(12);
T6 = fPat/sqrt(12);
disp(['Tth1 = ',num2str(Tth1)]);
disp(['Tth2 = ',num2str(Tth2)]);
disp(['T3   = ',num2str(T3)]);
disp(['T4   = ',num2str(T4)]);
disp(['T5   = ',num2str(T5)]);
disp(['T6   = ',num2str(T6)]);

%---variables-------------------------------------

M = zeros(N,N,D);	% the map.
V = zeros(N,N,D);	% the intermediate for finding the vector.

[X Y] = meshgrid(1:N,1:N);
disp('Running the simulation ...');

%---now do the simulation-------------------------

% initialize.

M(:,:,1) = X/N*d1;
M(:,:,2) = Y/N*d2;
dang=2*pi/N;
ang = [0:dang:(2*pi-dang)]; ang=ang' * ones(1,N) + ones(N,1) * ang;
M(:,:,3) = qPat * sin(ang);
M(:,:,4) = qPat * cos(ang);

Range = [zPat fPat];
for i=5:D %%% initialize values to half of the range
	M(:,:,i) = Range(i-4)/2;
end

else
	load(fileContinue);
	disp(['Continuing simulation on ',fileContinue]);
	iStart = iter;
    runFlag = 0;
    maxIter = 6000000;
    iSampleFixQ = 0;
    
    

end % if iContinue ~= 1


plotIter = 500;		% plotting steps.
saveIter = 100 * plotIter;

if runFlag == 1	
	filename=['runData/mapFixQScaleSigma-N-',num2str(N),'-rxy-',num2str(ratioXY),'-d1-',num2str(d1),'-qPat-',num2str(qPat),'-zPat-',num2str(zPat),'-fPat-',num2str(fPat),'-sigmaStart-',num2str(sigmaStart),'-stepScale-',num2str(stepScale),'-maxIter',num2str(maxIter),'-aplha-',num2str(alpha),'.mat'];
elseif iFixQ == 1
	filename=['runData/mapFixQ-N-',num2str(N),'-rxy-',num2str(ratioXY),'-d1-',num2str(d1),'-qPat-',num2str(qPat),'-zPat-',num2str(zPat),'-fPat-',num2str(fPat),'-sigma1-',num2str(sigma1),'-sigma2-',num2str(sigma2),'-aplha-',num2str(alpha),'.mat'];
else
	filename=['runData/map-N-',num2str(N),'-rxy-',num2str(ratioXY),'-d1-',num2str(d1),'-qPat-',num2str(qPat),'-zPat-',num2str(zPat),'-fPat-',num2str(fPat),'-sigma1-',num2str(sigma1),'-sigma2-',num2str(sigma2),'-aplha-',num2str(alpha),'.mat'];
end
disp(['ratioXY = ',num2str(ratioXY),' d1=',num2str(d1),' qPat=',num2str(qPat),' zPat=',num2str(zPat),' fPat=',num2str(fPat)]);

%---orientation tuning colormap------------------

ONES = ones(N,N);		% useful matrix.
DIST = zeros(N,N);
%%CHANGE SIGMA - HARDCODE
%%sigma1 = 8; %sigma for ODSF
%%igma2 = 8; %sigma for ODSF

NF1 = 3*sigma1; % 
NF2 = 3*sigma2;
H = zeros(2*NF1+1,2*NF2+1);		% the neiborghood function. (ORIGINAL)
for i=-NF1:NF1
for j=-NF2:NF2
	H(i+NF1+1,j+NF2+1) = exp(-i^2/sigma1^2-j^2/sigma2^2); 
end
end


% now self organize
if iPlot == 1
ww = 231;
hh = 207;
x0 = 7;
y0 = 897 - (iVPosition-1)*hh*1.4;
for i=1:D
	hand(i)=figure;
	set(gcf,'doublebuffer','on','position',[x0   y0   ww   hh]);
	x0 = x0 + ww;
end
end

M0 = M;
Range = [zPat fPat];
diffs = [];
nums = 0;

for iter=iStart:maxIter

	if runFlag == 2		% round visual field mapping to oval cortex
		sampleAng = rand * 2 * pi;
		sampleRad = N * rand;
		V(:,:,1) = ONES * sampleRad * cos(sampleAng);
		V(:,:,2) = ONES * sampleRad * sin(sampleAng);
	else
		V(:,:,1) = ONES * rand * d1;
		V(:,:,2) = ONES * rand * d2;
	end

	if iFixQ == 1 | iSampleFixQ == 1
		V(:,:,3) = ONES * qPat * cos(2*pi*rand);
		V(:,:,4) = ONES * qPat * sin(2*pi*rand);
	else
		V(:,:,3) = ONES * qPat * rand * cos(2*pi*rand);
		V(:,:,4) = ONES * qPat * rand * sin(2*pi*rand);
	end
		
	for i=5:D  %%%pick a random value within the range for Z, F.
		V(:,:,i) = Range(i-4) * rand;
	end

	DIST = zeros(N,N);
	if iFixQ == 1 
		for i=[1:2 5:D]	
			DIST = DIST + (M(:,:,i)-V(:,:,i)).^2;
		end
		AnglesSampled = atan2(V(:,:,4),V(:,:,3));
		AnglesMap = atan2(M(:,:,4),M(:,:,3));
		DIST = DIST + qPat^2 *(atan(tan(AnglesSampled-AnglesMap))).^2;
    else 
        %% change the weight (scale down) for the orientation, ocular
        %% dominance, and spatial frequency features (Q, Z, F). The find
        %% the distance between the stimulus, and receptive field vectors
        

        V(:,:,1) = xWeight * V(:,:,1);
        V(:,:,2) = yWeight * V(:,:,2);
        V(:,:,3) = qrWeight * V(:,:,3);
        V(:,:,4) = qsWeight * V(:,:,4);
        V(:,:,5) = zWeight * V(:,:,5);
        V(:,:,6) = fWeight * V(:,:,6);
     
        
		%for i=1:D	%%%FIND the euclidean distance between the vectors
			%DIST = DIST + (M(:,:,i)-V(:,:,i)).^2; %%%old version: take the
			%norm distance
            %DIST(:) = reshape(M(:,:,D), N*N, D)*V(1,1,D); %%%new version 1: take dot product            
		%end

        DIST = sum((M-V).^2, 3); %%%new version 2: take norm distance, add it to current DIST, on dimension 3       
	end	
	[mm mi] = min(DIST); %%%old version: take the MIN of all norm distances
    %[mm mi] = max(DIST); %%%new version 1: take the MAX of all dot products
	%[mm mi] = max(DIST); %%%new version 2: take the MIN of all dot products
    [mm i2] = min(mm);
    
	i1 = mi(i2);

	if i1-NF1 < 1
		idm = 1;
		idh = NF1+2-i1;
	else
		idm = i1-NF1;
		idh = 1;
	end
	if i1+NF1 > N
		ium = N;
		iuh = NF1+1+N-i1;
	else
		ium = i1+NF1;
		iuh = 2*NF1+1;
	end

	if i2-NF2 < 1
		jdm = 1;
		jdh = NF2+2-i2;
	else
		jdm = i2-NF2;
		jdh = 1;
	end
	if i2+NF2 > N
		jum = N;
		juh = NF2+1+N-i2;
	else
		jum = i2+NF2;
		juh = 2*NF2+1;
	end
	
	if iFixQ == 1 
		for i=[1:2 5:D]
			M(idm:ium,jdm:jum,i) = M(idm:ium,jdm:jum,i) + alpha * H(idh:iuh,jdh:juh) .* (V(idm:ium,jdm:jum,i)-M(idm:ium,jdm:jum,i));
		end
		AnglesChange = 0.01 * alpha * H(idh:iuh,jdh:juh) .* atan(tan(AnglesSampled(idm:ium,jdm:jum)-AnglesMap(idm:ium,jdm:jum)));
		AnglesMap(idm:ium,jdm:jum) = AnglesMap(idm:ium,jdm:jum) + AnglesChange;
		M(idm:ium,jdm:jum,3) = qPat * cos(AnglesMap(idm:ium,jdm:jum));
		M(idm:ium,jdm:jum,4) = qPat * sin(AnglesMap(idm:ium,jdm:jum));		
    else
        %%MODIFIED: take into account sigmaXY = 3 instead of 5.
        %%(neighborhood function is Hmod)
		for i=1:D
			M(idm:ium,jdm:jum,i) = M(idm:ium,jdm:jum,i) + alpha * H(idh:iuh,jdh:juh) .* (V(idm:ium,jdm:jum,i)-M(idm:ium,jdm:jum,i));
        end
    end
    
%     %%Revisions Scott: for "if iFixQ == 1" loop
%     iidx = round(idm:ium);
%     jidx = round(jdm:jum);
%  
%     h = H_Mat(round(idh:iuh),round(jdh:juh), :);
%     if iFixQ == 1
%         didx = [1:2 5:D];
%         M(iidx,jidx,didx) = M(iidx,jidx,didx) + alpha * h(:,:,didx) .* (V(iidx, jidx,didx)-M(idx, jidx,didx));
%        
%         AnglesChange = 0.01 * alpha * h(:,:,1) .* atan(tan(AnglesSampled(idx,jidx)-AnglesMap(iidx,jidx)));
%         AnglesMap(iidx,jidx) = AnglesMap(iidx, jidx) + AnglesChange;
%         M(iidx,jidx,3) = qPat * cos(AnglesMap(idx, jidx));
%         M(iidx, jidx,4) = qPat * sin(AnglesMap(iidx, jidx));       
%     else
%         M(iidx, jidx, :) = M(iidx, jidx, :) + alpha*h.*(V(iidx, jidx, :)-M(iidx, jidx, :));
%     end
%     %%/Revisions Scott
    
		
	if mod(iter,plotIter) == 0 & iPlot == 1
		nums = nums+1;
	
		disp(['Fraction complete = ',num2str(iter/maxIter)]);
		colormap('default');
		figure(hand(1));
		subplot(2,1,1);imagesc(M(:,:,1));  
		subplot(2,1,2);imagesc(M(:,:,2));

		qMap = sqrt(M(:,:,3).^2+M(:,:,4).^2);
		phiMap = atan2(M(:,:,4),M(:,:,3))/2+pi/2; 
		figure(hand(2));imagesc(qMap);  
		figure(hand(3));imagesc(phiMap); colormap('hsv');

		for i=5:D
			figure(hand(i-1));
			imagesc(M(:,:,i)); colormap('default');
		end
		
		for i=1:D
			diffs(i,nums) = sum(sum(abs(M(:,:,i)-M0(:,:,i))))/(N*N);
		end
		figure(hand(D));
		for i=1:D
			subplot(D,1,i); plot(diffs(i,:));
		end
		M0 = M;		
		pause(0.5);
		if nums==200
			nums = 0;
		end
	end
	if mod(iter,saveIter) == 0
		disp(['At iter = ',num2str(iter),' saving the data ...']);
		save(filename,'N','d1','d2','qPat','zPat','fPat','ratioXY','sigma1','sigma2','alpha','M','iter','D','diffs');
		if iSaveIntermediateFiles == 1
			iInterFile=iInterFile+1
			fn=['runData/mapInter-',num2str(iInterFile),'.mat'];
			save(fn,'N','d1','d2','qPat','zPat','fPat','ratioXY','sigma1','sigma2','alpha','M','iter','D','diffs');
		end
	end	
	
	if runFlag == 1 & mod(iter,stepScale) == 0		% aneal sigma1 and sigma2		
		sigma1 = max(1,round(sigma1/2^exponent));
		sigma2 = sigma1;
		
		NF1 = 3*sigma1;
		NF2 = 3*sigma2;
		H = zeros(2*NF1+1,2*NF2+1);		% the neiborghood function.
		for i=-NF1:NF1
		for j=-NF2:NF2
			H(i+NF1+1,j+NF2+1) = exp(-i^2/sigma1^2-j^2/sigma2^2);
		end
		end
	
	end
end

%---save the result-----------------------------------------------------

save(filename,'N','d1','d2','qPat','zPat','fPat','ratioXY','sigma1','sigma2','alpha','M','iter','D','diffs');
 
%---plot the space projection of the topographic map--------------

if iPlot == 1
figure;
hold on;
iskip = 5;
for i=1:iskip:N
for j=1:iskip:N-iskip
	plot([M(i,j,1) M(i,j+iskip,1)],[M(i,j,2) M(i,j+iskip,2)]);
end
end
for j=1:iskip:N
for i=1:iskip:N-1
	plot([M(i,j,1) M(i+iskip,j,1)],[M(i,j,2) M(i+iskip,j,2)]);
end
end
hold off;
axis equal;
end
