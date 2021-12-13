% [HIST, LDPIM] = LDP(IMAGE, SUBSZ, ORDER)  
% Calculates the LDP images and histogram of IMAGE. The concatenation  
% of the LDP histograms of each subblock (of size SUBSZ) in the image 
% is returned in HIST, and the LDP code images are returned in LDPIM. 
% The LDP order is specified by ORDER. 
% This code is based on B Zhang, Y Gao, S Zhao & J Liu, "Local Derivative Pattern Versus Local Binary Pattern: Face Recognition With High-Order Local Pattern Descriptor," Image Processing, IEEE Transactions on, vol.19, no.2, pp.533,544, Feb. 2010 
function [hout,ldpim]=LDP(im, subsz, order) 
numBinsPerBlock = 32;       % number of bins per direction per block 
 
assert(order >=2 && order <= 10); 
assert(all(floor(size(im) ./ subsz(:)') == (size(im) ./ subsz(:)'))); 
 
im = double(im); 
der0 = im; 
der45 = im; 
der90 = im; 
der135 = im; 
 
% calc derivative image of appropriate order 
for ii=2:order 
    der0 = [der0(:, 1:end-1)-der0(:, 2:end) zeros(size(der0,1),1)];    % subtract right 
    der45 = [zeros(1,size(der45,2)); der45(2:end, 1:end-1)-der45(1:end-1, 2:end) zeros(size(der45,1)-1,1)];    % subtract above-right 
    der90 = [zeros(1,size(der90,2)); der90(2:end, :)-der90(1:end-1, :)];    % subtract above 
    der135 = [zeros(1,size(der135,2)); zeros(size(der135,1)-1,1) der135(2:end, 2:end)-der135(1:end-1, 1:end-1)];    % subtract above-left 
end 
 
% calc LDP code using sign(dot product), combine adjacent codes to get 32 codes/bins 
div = 256 / numBinsPerBlock; 
codes0 = floor(doDotProdCode(der0) / div); 
codes45 = floor(doDotProdCode(der45) / div); 
codes90 = floor(doDotProdCode(der90) / div); 
codes135 = floor(doDotProdCode(der135) / div); 
 
% make histograms for each subimage (block) within the code images 
bins = 0:numBinsPerBlock-1; 
hout = zeros(numBinsPerBlock * 4, prod(size(im) ./ subsz)); 
hidx = 1; 
for ii=0:size(im,1)/subsz(1)-1 
    bii = ii*subsz(1) + (1:subsz(1)); 
    for jj=0:size(im,2)/subsz(2)-1 
        bjj = jj*subsz(2) + (1:subsz(2)); 
        block0 = codes0(bii, bjj); 
        block45 = codes45(bii, bjj); 
        block90 = codes90(bii, bjj); 
        block135 = codes135(bii, bjj); 
         
        hout(:, hidx) = [intHist(block0(:), bins) intHist(block45(:), bins) intHist(block90(:), bins) intHist(block135(:), bins)]'; 
%         hout(:, hidx) = [hist(block0(:), bins) hist(block45(:), bins) hist(block90(:), bins) hist(block135(:), bins)]'; 
         
        hidx = hidx + 1; 
    end 
end 
 
hout = hout(:); 
 
if nargout > 1 
	ldpim.b0 = block0; 
	ldpim.b45 = block45; 
	ldpim.b90 = block90; 
	ldpim.b135 = block135; 
end 
 
% figure(1)		% show the derivative images 
% subplot(221), imshow(der0,[]); 
% subplot(222), imshow(der45,[]); 
% subplot(223), imshow(der90,[]); 
% subplot(224), imshow(der135,[]); 
%  
% figure(2)		% show 
% subplot(221), imshow(codes0,[]); 
% subplot(222), imshow(codes45,[]); 
% subplot(223), imshow(codes90,[]); 
% subplot(224), imshow(codes135,[]); 
 
end 
 
function yim=doDotProdCode(im) 
[m,n] = size(im); 
shim = zeros(m+2,n+2); 
shim(2:end-1, 2:end-1) = im;    % zero-pad im 
 
dp1 = (shim(1:end-2,1:end-2) .* im) <= 0;   % Z1 . Z0 
dp2 = (shim(1:end-2,2:end-1) .* im) <= 0;   % Z2 . Z0 
dp3 = (shim(1:end-2,3:end) .* im) <= 0;     % Z3 . Z0 
dp4 = (shim(2:end-1,3:end) .* im) <= 0;     % Z4 . Z0 
dp5 = (shim(3:end,3:end) .* im) <= 0;       % Z5 . Z0 
dp6 = (shim(3:end,2:end-1) .* im) <= 0;     % Z6 . Z0 
dp7 = (shim(3:end,1:end-2) .* im) <= 0;     % Z7 . Z0 
dp8 = (shim(2:end-1,1:end-2) .* im) <= 0;   % Z8 . Z0 
 
yim = 1*dp1 + 2*dp2 + 4*dp3 + 8*dp4 + 16*dp5 + 32*dp6 + 64*dp7 + 128*dp8; 
end 
 
% note that all entries must belong to a value in bins, or it won't be 
% counted 
function ho=intHist(v, bins) 
ho = zeros(1,length(bins)); 
for ii=1:length(bins) 
    ho(ii) = sum(v == bins(ii)); 
end 
end