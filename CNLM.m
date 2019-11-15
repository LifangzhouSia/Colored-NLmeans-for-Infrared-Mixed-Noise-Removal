function [dI]  = CNLM(Insy, r, d, nl)

if(~exist('nl', 'var')),
    [rnd_str, fpn_str] = NoiseLevelEstimation(Insy);
    nl = [rnd_str, fpn_str];
    fprintf('\tEstimated Random Noise Level is %s\n', num2str(nl(1)))
    fprintf('\tEstimated Stripe Noise Level is %s\n', num2str(nl(2)))
end

[dI]  = CNLMFilter(Insy, r, d, nl);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Noise Level Estimation %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rnd_str, fpn_str] = NoiseLevelEstimation(Insy)

dx = diff(Insy, 1, 2);
dy = diff(Insy, 1, 1);

sx = median(abs(dx(:) - median(dx(:))))/0.6745;
sy = median(abs(dy(:) - median(dy(:))))/0.6745;

rnd_str  = sy/1.414;
fpn_str = sqrt(sx^2 - sy^2)/1.414;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Noise Removal Filter %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dI] = CNLMFilter(Insy, r, d, nl)

Insy = double(Insy);
[hei, wid] = size(Insy);
rnd_str = nl(1); fpn_str = nl(2);
vn = 1/rnd_str^2;
vs = fpn_str^2/(rnd_str^2 + (2*r+1)*fpn_str^2);
gamma = 10;

L = padarray(Insy, [d, d], 'symmetric', 'both');
N = boxfilter(ones(hei, wid), r);

cL = zeros(hei, wid);
wL = zeros(hei, wid);
wmax = zeros(hei, wid);
for k = -d:d
    for g = -d:d
        
        % Check for middle patch
        if k==0 && g==0
            continue;
        end
        
        % Check for vertical patch
        if g == 0
            F = 0;
        else
            F = 0.81;
%             F = 0;
        end
        
        % Construct residual image
        tL = L(1+d+k:end-d+k,1+d+g:end-d+g);
        dL = Insy - tL;
        
        % Construct weight matrix
        meanDL = boxfilter(dL, r)./N;
        xN     = boxfilter(dL.*dL, r)./N;
        varDL  = xN - meanDL.*meanDL;
        Phi    = F.*exp(- varDL./(gamma*rnd_str^2 + gamma*fpn_str^2));
        xS     = stripeDistance(dL, r);
        thisW = exp(- vn.*(xN - Phi.*xS));
        cL = cL + thisW.*tL;
        wL = wL + thisW;
        
        wmax = max(wmax, thisW); wmax(wmax == 0) = 0.000001; 
    end
end

dI = (cL + wmax.*Insy) ./ (wL + wmax);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Some Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I2 = stripeDistance(I, r)

[hei, wid] = size(I);

I2 = verticalBoxfilter(I, r);
I2 = (verticalBoxfilter((I2.*I2).', r)).';

end

function mI = verticalBoxfilter(I, r)

[hei, wid] = size(I);
mI = zeros(hei, wid);
Nv = D1boxfilter(ones(hei, 1), r);
for k = 1:wid
    mI(:, k) = D1boxfilter(I(:, k), r)./Nv;
end

end

function imDst = D1boxfilter(I, r)

imLine = I(:);
lth = length(imLine);
imDst = zeros(size(imLine));

imCum = cumsum(imLine, 1);
imDst(1:r+1) = imCum(1+r:2*r+1);
imDst(r+2:lth-r) = imCum(2*r+2:lth) - imCum(1:lth-2*r-1);
imDst(lth-r+1:lth) = repmat(imCum(lth), [r, 1]) - imCum(lth-2*r:lth-r-1);

end

function imDst = boxfilter(imSrc, r)
%   BOXFILTER   O(1) time box filtering using cumulative sum
%
%   - Definition imDst(x, y)=sum(sum(imSrc(x-r:x+r,y-r:y+r)));
%   - Running time independent of r; 
%   - Equivalent to the function: colfilt(imSrc, [2*r+1, 2*r+1], 'sliding', @sum);
%   - But much faster.

[hei, wid] = size(imSrc);
imDst = zeros(size(imSrc));

%cumulative sum over Y axis
imCum = cumsum(imSrc, 1);
%difference over Y axis
imDst(1:r+1, :) = imCum(1+r:2*r+1, :);
imDst(r+2:hei-r, :) = imCum(2*r+2:hei, :) - imCum(1:hei-2*r-1, :);
imDst(hei-r+1:hei, :) = repmat(imCum(hei, :), [r, 1]) - imCum(hei-2*r:hei-r-1, :);

%cumulative sum over X axis
imCum = cumsum(imDst, 2);
%difference over Y axis
imDst(:, 1:r+1) = imCum(:, 1+r:2*r+1);
imDst(:, r+2:wid-r) = imCum(:, 2*r+2:wid) - imCum(:, 1:wid-2*r-1);
imDst(:, wid-r+1:wid) = repmat(imCum(:, wid), [1, r]) - imCum(:, wid-2*r:wid-r-1);
end




