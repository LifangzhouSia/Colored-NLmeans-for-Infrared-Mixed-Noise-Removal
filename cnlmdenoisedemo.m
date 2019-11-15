function cnlmdenoisedemo
%CNLMDENOISEDEMO CNLM denoising method.
%  CNLMDENOISEDEMO reads an image, adds random stripe and white noise and denoises it
%  using CNLM denoising. 
%
%  To run the demo, type CNLMDENOISEDEMO from the Matlab prompt.

%  Fangzhou Li
%
%  November 2019

disp(' ');
disp('  **********  CNLM Denoising Demo  **********');
disp(' ');
disp('  This demo reads an image, adds random Stripe and White mixed Gaussian noise.');
disp('  The mixed noise will be removed using an improved NLM method (CNLM).');
disp('  The denoised image will be shown.');
disp(' ');

%% prompt user for image %%

im = readImage('cnlmdenoisedemo');

%% generate noisy image %%

sigma_white = 5;
sigma_stripe = 5;

disp(' ');
disp('Generating mixed noisy image...');

n = randn(size(im)) * sigma_white;
n = n + repmat(randn(1, size(im, 2)), size(im, 1), 1) .* sigma_stripe;
imnoise = im + n;

% denoise!
disp('Performing CNLM denoising...');
[dI] = CNLM(imnoise, 2, 5);

% show results %

figure; imshow(newlp(im));
title('Original image');

figure; imshow(newlp(imnoise)); 
title('Noisy image')
% title(sprintf('Noisy image, PSNR = %.2fdB', 20*log10(params.maxval * sqrt(numel(im)) / norm(im(:)-imnoise(:))) ));

figure; imshow(newlp(dI));
% title(sprintf('Denoised image, PSNR: %.2fdB', 20*log10(params.maxval * sqrt(numel(im)) / norm(im(:)-imout(:))) ));
title('Denoised image')

figure; imshow(newlp(n))
title('Extracted mixed noise')
