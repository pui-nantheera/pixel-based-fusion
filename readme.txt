**************************************************************************
*                  Pixel-based image fusion tools                        *
*                                                                        *
*             Copyright (C)2012 University of Bristol                    *
**************************************************************************

These tools provide MATLAB functions for image fusion (IFT) and quality 
metrics for image fusion assessment. These functions have been tested with
MATLAB R2012a(7.14).

Three major groups of the IFTs and quality metrics are listed below:


List of pixel-based image fusion functions
------------------------------------------

I. Image Fusion in Spatial Domain
I.1) Average intensity:       				fuse_avg.m
I.2) PCA method:              				fuse_pca.m
I.3) Min or Max intensity:    				selc.m
I.4) Weighted average with spatial frequency		fuse_spafrq.m

--------------------------------------------------------------------------
II. Image Fusion in Pyramid Transform Domain
II.1) Laplacian pyramid       				fuse_lap.m
II.2) Filter-Substract-Decimate Laplacian Pyramid    	fuse_fsd.m   
II.3) Ratio Pyramid           				fuse_rat.m
II.4) Contrast Pyramid        				fuse_con.m
II.5) Gradient Pyramid        				fuse_gra.m
II.6) Morphological Pyramid				fuse_mod.m

    In pyramid transform domain, 
	select base coefficients using 		selb.m
	select other scale coefficients using 	selc.m

--------------------------------------------------------------------------
III. Image Fusion in Wavelet Transform Domain
III.1) Discrete Wavelet Transform				fuse_dwb.m
III.2) Lifting Scheme on Quincunx Grids			fuse_lisq.m
III.3) Shift-Invariant Discrete Wavelet Transform	fuse_sih.m
III.4a) Dual-tree complex wavelet transfrom  (simple)   fuse_cwt.m
III.4b) Dual-tree complex wavelet transfrom  (model1)   fuse_model.m
	Note: Only two images
III.4c) Dual-tree complex wavelet transfrom  (model2)   fuse_dtcwt.m
	Note: Denoising is an option.

    In wavelet transform domain, 
	select low-pass coefficients using 	selb.m
	select high-pass coefficients using 	selcComplex.m (for simple)
	select high-pass coefficients using 	selcModel.m   (for model1)
	select high-pass coefficients using 
			       functions inside fuse_dtcwt.m  (for model2)

==========================================================================

List of quality metrics for objective assessment
------------------------------------------------
I.  Mutual information			mif.m
II. Petrovic and Xydeas Metric		petmetric.m
III.Piella's Quality Index		imqmet.m
IV. Cvejic's Quality Index		Cvejic_metric.m

==========================================================================

Usage and Examples
------------------
All image fusion methods are included in fuse.m and the usage of fuse.m can 
be found in demo.m.

Input images can be stored in Cell array if they are colour or greyscale images
or 3 dimension matrix if all of them are grayscale images.
Examples: 	> origImgs{1} = im2double(imread('imageColour1.png'));
		> origImgs{2} = im2double(imread('image2.png'));
Or		> origGImgs(:,:,1) = im2double(imread('imageGrey1.png'));
		> origGImgs(:,:,2) = im2double(imread('imageGrey2.png'));

If any of them is colour image, the fuse.m will convert to greyscale using 
HSV transformation. The intensity channel is used in fusion process, whilst 
the others are fused by averaging amongst the available channels. The fused 
results are converted back to RGB at the end of the process. 
Note that the rgb2ycbcr can be also employed.

The first output is the greyscale fused image. The second output is the colour 
fused image. If the colours are not available, this output will be same as the first output. The third output is the computational time of the fusion process.

The fusion method is identified with the name (string) or number. This can be 
found by using help.
Examples:	> help fuse
		> [fusedGImg, fusedCImg, processTime] = fuse(OrigImgs,'SIDWT','high',1, 'low',3);
		> figure; imshow(fusedGImg);

For more flexible usage, please refer to each image fusion function listed above.
Usage of each function can be found by type help xxx.m
Examples:	> help fuse_dtcwt
		> fusedGImg = fuse_dtcwt(origGImgs); % use all default parameters
If you want to use some default parameters, you can use the empty matrix ([]) to
define the value.
Examples:	> model = 'CAU'   % using Cauchy distribution
		> variate = 2;    % Bivariate
		> gain = [2 1.4]; % sharpening
		> fusedGImg = fuse_dtcwt(origGImgs, [], model, [], [], [], variate, [], gain);

Usage of quality metrics can be found by using help with that function.
The first input is the 3D matrix of the grayscale input images. The second input is the fused image.
Examples: 	> help petmetric
		> score = petmetric(OrigGImgs, fusedGImg);

==========================================================================
More questions or some advises are welcome, 
please contact N.Anantrasirichai@bristol.ac.uk
