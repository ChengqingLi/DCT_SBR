# Directory

- Objectives and motivations

  - What is DCT transform?

  - Sign bit recovery.

  - Sign bit recovery model.

- Installation

  - Matlab installation

  - Cplex installation

- Main function

  - cplex_signbit

  - cplex_signbitHierarchyDCACMIP

# Objectives and motivations

## What is DCT transform ?

Focus on a grey scale image with width w and height h. A w\*h matrix can be obtained, in which each pixel value varies from 0 to 255. Calculate the contributions of each of our cosine waves to this image and determine the bits of cosine needed to add together to create an image that looks exactly like this. To start with, centre all of these values to 0, which are currently centred at 128, for cosine wave goes from 1 to -1 instead of 1 to 0. Simply take away 128 off every value, and the shifted block centred to 0 will be obtained. Then, calculate DCT2 coefficients, and a new matrix will be obtained, the number in which varies from -1024 to 1024, representing the weight or the amount of each of our cosine waves. Later we use quantization table to get the quantized matrix and omit the elements whose value is close to 0. Finally, reverse the procedure and a compressed imaged is created.

##  Sign bit recovery

It is possible to recovery the missing sign bits, since there are relations between adjacent pixel values. The difference between any two neighboring pixels is a Laplacian variate with zero mean and a small variance, which is a well-known feature of most natural images. Our proposed model is to minimize the objective function defined as the summation of the absolute values between adjacent pixels, which is a liner problem and can be solved in linear time.

## Sign bit recovery model

### Objective function:

Minimize the summation of h variables (summation of all absolute values between adjacent pixels.

### Variables:

In order to build the system of 2-D DCT sign bits recovery model, it is necessary to introduce several variables, namely
- x-pixel values(all nonnegative),

- y-DCT coefficients (without DC prediction),

- QZ-quantization number (quantized z * quantization table),

- z-DCT coefficients (with DC prediction),

- h-absolute values between adjacent pixels (vertical and horizontal).

- Constraints: (rhs and lhs for sparse matrix A)

### Constraints on prediction mode:
There are four prediction mode to be chosen:

0. No prediction, just return x directly

1. Predicted from previously coded block in the same block-row

2. Predicted from the previously coded blocks in the same image scanned using the raster order, as defined in JPEG standard

3. Predicted from one or two previously coded blocks in the same image scanned using the raster order, as defined in MPEG-1/2 standards

For each prediction mode, the constraints differ a little bit.

Constraints on QZ

0. DCT coefficients with DC prediction is centred with quantization number varying from minus or plus half of the number on the quantization table, which is to say QZ(i)-quan_table(i)/2<=z(i)<=QZ(i)+quan_table(i)/2

Constraints on x

0. X is ruled by the DCT transform formula, which is to say x(k,l) =sum_u\sum_v A(u,v,k,l)*z(u,v)

Constraints on h

0. H is set to be equal to the absolute value between adjacent pixels, which is to say x(i)-x(j) <= h(x(i),x(j)) and x(j)-x(i) <= h(x(i),x(j)).

Constraints on y

0. For known y, set them to be equal

1. For unknown y, we do a ‘relaxation’, we relax y to be in-between its minus absolute value and plus absolute value.

Lower and upper bounds: (up and lb for variables)

- Each x itself is between its minimum and maximum, which is to say x min<= x <= x max.

- The distance between two values is between 0 and maximum of x, which is to say 0<=h<=x max.

- y and z's lower bound and upper bound are determined implicitly via the DCT constraints, therefore, no need to add upper bound and lower bound.

Remark:

- Error tolerance, a value smaller than 0.5 is sufficient to make sure rounding will work properly for x-variables (pixel values).

- eps = 0.499;


# Installation

Software used:

·    Matlab_R 2015b

o  Matlab is used as an interface only.

·    Cplex_studio1263

o  CPLEX for MATLAB is an extension to IBM ILOG CPLEX Optimizers that allows a user to define optimization problems and solve them with in MATLAB. Thus a student or practitioner who is using MATLAB can easily solve optimization problems within that framework. ---From IBM ILOG CPLEX Optimization Studio Getting Started with CPLEX for MATLAB.

Installation:

·    Installing Matlab

1.   Download MATLAB for Mac

2.   Double-click the downloaded file to extract the files and open the Matlab folder.

3.   Double-click the InstallForMacOSX file to begin installation.

4.   On the Select Installation Method screen, choose Log in with a MathWorks Account and click next.

·    Installing CPLEX

1.   Download cplex_studio1263.osx.bin

2.   Open ‘terminal’ and enter ‘chomd a+x cplex_studio1263.osx.bin’

3.   Execute the .bin file by entering Cplex_studio1263.osx.bin

4.   The whole procedure may cost approximate 4-6 minutes.

Path setting:

o  Using addpath function to add cplex in Matlab and save function to save it.

1.   Add path /Users/myAccount/Applications/IBM/ILOG/CPLEX_Studio1263/cplex/Matlab/x86-64_osx

2.   Save path


# Main function

## cplex_signbit:

This function invokes IBM ILOG CPLEX’s MATLAB Cplex class to solve the problem of recovering unknown sign bits of DCT coefficients of a DCT-transformed image with a linear programming model.

Example:

[imgs, Ts, imgs_dct2,PSNRs,SSIMs,BERs] = cplex_signbit…

(img_dct2, mask, DCpre_mode, DCTabs_min, relaxX, relaxZ, recovery_method, bDClevelshifted, H, W, Xmax)

Explanation:

The very things that you need to input are the image and the mask and other parameters all have default setting. You can directly call the image itself or enter the file name. As for the mask, it depends on the number of unknown coefficients. For a file name as the input, the program will convert the image to its block wise 8*8 2D DCT matrix automatically.

For mac os x:

Image_dct2 can be replaced by ‘images/image1.pgm’.

mask can be replaced by ‘masks/mask_01.pgm’.

For windows:

Image_dct2 can be replaced by ‘images\image1.pgm’.

mask can be replaced by ‘masks\mask_01.pgm’.

Input arguments:

Img_dct2: the DCT-transformed image with unknown coefficients (can also take an image or a file name as the input which will be read/converted to its block wise 8x8 2D DCT matrix)

mask: a logical mask array showing which coefficients are unknown (default: sign bits of all DC coefficients are unknown) (can take a binary image or a filename as well)

Other input arguments:

DCpred_mode: mode of DC coefficient prediction (default: 0)

0.   No prediction, just return x directly

1.   Predicted from previously coded block in the same block-row

2.   Predicted from the previously coded blocks in the same image scanned using the raster order, as defined in JPEG standard

3.   Predicted from one or two previously coded blocks in the same image scanned using the raster order, as defined in MPEG-1/2 standards


DCTabs_min: threshold for removing DC coefficients with unknown sign bits but smaller amplitude (if amplitude <= DCTabs_min, the coeffient will be set to 0 and not be recovered; default: 0)

NB: When DCTabs_min is too large it may create too strong so that a feasible solution does not exist!

relaxX: true: x-variables are relaxed by xError for JPEG images. xError is calculated from the quantization errors of DCT coefficients. false: x-variables are not relaxed. relaxX is ignored if the input image is not JPEG. (default: true)

relaxZ: true: z-variables are relaxed by quantization errors of JPEG images. false: z-variables are not relaxed. relaxZ is ignored if the input image is not JPEG. (default: true)

recovery_method: the method for recovering the unknown sign bits

0. all the methods (default)

1. when the recovered DCT coefficient is 0, the result DCT coefficient is 0 as well

2. when the recovered DCT coefficient is 0, the result DCT coefficient's sign is set to 1

3. when the recovered DCT coefficient is 0, the result DCT coefficient's sign is set to 0

4. when the recovered DCT coefficient is 0, the result DCT coefficient's sign is randomly assigned

5. all signs are randomly guessed following random bounding, using a uniform distribution

6. all signs are randomly guessed following random bounding, using a distribution learned from some test images


bDClevelshifted: true = level shift DC coefficients by -round(xmax/2) as defined in JPEG standard (default) false = do not level shift DC coefficients (in this case DC coefficients' sign bits are always 1 so effectively known in all cases)

H, W: height and width of the 2-D DCT block (default: H = 8, W = H)

Xmax: maximal pixel value (default: 255) minimum pixel value Xmin is assumed to be always 0.

Output arguments:

imgs: a cell array holding a number of images recovered using different methods

- imgs.x0: the original image as input

- imgs.x: the image recovered directly by CPLEX

- imgs.x1: the image recovered from recovery method 1

- imgs.x2: the image recovered from recovery method 2

- imgs.x3: the image recovered from recovery method 3

- imgs.x4: the image recovered from recovery method 4

- imgs.x5: the image recovered from recovery method 5

- imgs.x6: the image recovered from recovery method 6

- imgs.x00: the image recovered from a naive method (all unknown signs are set to 0)

- imgs.x01: the image recovered from a naive method (all unknown signs are set to 1)

- imgs.x02: the image recovered from a naive method (all unknown signs are randomly assigned to 0 or 1)

Ts: the time consumed by different part of the function

- Ts.seconds_prep: time consumed by preprocessing input arguments in seconds

- Ts.seconds_x0123: time consumed by creating initial condition(s)

- Ts.seconds_model: time consumed by creating the CPLEX optimization model

- Ts.seconds_cplex: time consumed by solving the CPLEX optimization model in seconds

- Ts.ticks_cplex: deterministic time consumed by solving the CPLEX optimization model in ticks

- Ts.seconds_solutions: time consumed by creating the final solutions in seconds

- Ts.seconds_all: time consumed by the whole function in seconds



imgs_dct2: DCT coefficients of some images (with differentially encoded DC coefficients if DC prediction is involved)

- imgs_dct2.y0: DCT coefficients of the input (original) image

- imgs_dct2.y: DCT coefficients of img.x (the one recovered directly by CPLEX)


The above are useful for producing the figure on accuracy sign bit estimation using the optimization model.

PSNRs, SSIMs: indicators of visual quality of recovered images as measured by PSNR and SSIM

BERs: Bit Error Rates of different recovery methods. BER = Number of wrong bits / Total number of bits


## cplex_signbitHierarchyDCACMIP:

This function solves the problem of recovering unknown sign bits of DCT coefficients of a DCT-transformed image ‘image_dct2’ region by region (each by a mixed integer programming model).

Example

[imgs, Ts, imgs_dct2, SB, PSNRs, SSIMs, BERs] =cplex_signbitHierarchyDCACMIP(img_dct2, mask, regionSize_H, regionSize_W, DCpred_mode, DCdepend, relaxXBound, relaxZBound, DCpassByRegion, DCTabs_min, timeLimit, relatGap, DCpass_method, recovery_method, relaxXBound_LP, relaxZBound_LP, bDClevelshifted, H, W, Xmax, fileLabel)

Input arguments:

img_dct2: the DCT-transformed image with unknown coefficients (can also take an image or a file name as the input which will be read/converted to its blockwise 8x8 2D DCT matrix)

mask: a logical mask array showing which coefficients are unknown

(default: sign bits of all DC coefficients are unknown)

(can take a binary image or a filename as well)

regionSize_H, regionSize_W: height or width of each region, if the image size cannot not be divided by the region height or width, region in the last region-column will be of width=mod(w,regionSize_W), region in the last region-row will be of height=mod(h,regionSize_H).

(Default: regionSize_H=64, regionSize_W=regionSize_H). Note if DCpred_mode\==2 && (DCdepend\==1 || DCdepend\==2), regionSize_W will be forced to be equal to image width.

DCpred_mode: mode of DC coefficient prediction (default: 0)

DCdepend:

The DC dependency mode used in 1st pass to deal with the missing information from the DC coefficients of previous regions

DCdepend=0: all DC coefficients that have dependency on previous regions will be relaxed to [-inf inf], the dependency within the region is maintained though.

DCdepend=1: as the regions are solved row by row, the region extended to be a rectangle from the upper left corner to the lower right corner of the current unsolved region, the sign bits from previous solved regions are marked as known and use the previously solved values as their value, i.e. the DC prediction will be maintained and all the pixel differences in this extended region will be concerned.

DCdepend=2: the current solved region will use the DC coefficients from previously solved regions as an offset to the current DC coefficent, but only the pixel differences in the current unsolved region is concerned. If the input is JPEG file, as the quantization error accumulates along the DC prediction path, the accumulated quantization error is added to the relaxation of z-varibles when relaxZBound is set to true relaxXBound: when set to true, the x-variables are relaxed to [-inf inf] in 1st pass when DCdepend=1,2 (x-variables bounds are NOT relaxed when DCdepend=0 regardless of the value of relaxXBound. Note x-variables are relaxed to [-inf inf] in 2nd pass in all cases.

relaxZBound: when set to true and when the input is a JPEG file, the bounds of z-variables will be relaxed by the quantization error

DCpassByRegion: mode of LP 2nd pass.

- 0: try all methods (default)

- 1: the brightness are aligned in a region-based manner

- 2: all blocks participate in the LP 2nd pass individually

DCTabs_min: threshold for removing DC coefficients with unknown sign bits but smaller amplitude (if amplitude <= DCTabs_min, the coefficient will be set to 0 and not be recovered; default: 0)

NB: When DCTabs_min is too large it may create too strong constraints so that a feasible solution does not exist!

timeLimit: time limit for CPLEX MIP solver. If time limit=[], the time limit is determined by CPLEX default (default: []). Note the LP 2nd pass is also limited by this parameter, as it's passed to CPLEX as a MIP model with no value-changeable binaries.

relatGap: the CPLEX relative MIP gap tolerance, if relatGap is [],relative MIP gap tolerance is determined by CPLEX default (default: []). Note the LP 2nd pass is not limited by this parameter

DCpass_method: method of 2nd pass

- 0: try all (default)

- 1: MIP

- 2: LP

recovery_method: the method in each MIP pass for regions or MIP/LP 2nd pass for recovering the unknown sign bits with non-zero absolute values < DCTabs_min. Note the results will be in the mat file saved from each pass, but only the one generated directly as CPLEX output is used in this function

- 0: try all

- 1: set the SB to 1

- 2: set the SB to 0

- 3: set the SB randomly

- other number: try none, leave those coefficients as 0

The results of setting those coefficients to zeros is generated as CPLEX output (default)

relaxXBound_LP: relaxX used in the pure LP method (default: true)

relaxZBound_LP: relaxZ used in the pure LP method (default: true)

bDClevelshifted:

- true = level shift DC coefficients by -round(xmax/2) as defined in JPEG standard (default)

- false = do not level shift DC coefficients (in this case DC coefficients' sign bits are always 1 so effectively known in all cases)

H, W: height and width of the 2-D DCT block (default: H = 8, W = H)

Xmax: maximal pixel value (default: 255)

The minimum pixel value Xmin is assumed to be always 0.

fileLabel: the (imaginary) file path for img_dct2 to be used for naming the mat data file generated. Used only when the img_dct2 is not a filename



Output arguments:

imgs: a cell array holding a number of images recovered using different methods

- imgs.x0: the original image as input

- imgs.x: 1st pass

- imgs.x1: 2nd MIP pass

- imgs.x2: 2nd LP pass by blocks

- imgs.x3: 2nd LP pass by regions

- imgs.x00: the image recovered from a naive method (all unknown signs are set to 0)

- imgs.x01: the image recovered from a naive method (all unknown signs are set to 1)

- imgs.x02: the image recovered from a naive method (all unknown signs are randomly assigned to 0 or 1)

- imgs.{'x000' 'x001' 'x002' 'x003' 'x004', 'x005', 'x006'}: pure LP methods, please refer to cplex_signbit for details

Ts: the time consumed by different part of the function

   Ts.seconds_prep: time consumed by preprocessing input arguments in seconds (summed up in all passes)

- Ts.seconds_x0123: time consumed by creating initial condition(s)

- Ts.seconds_model: time consumed by creating the Gurobi optimization model (summed up in all passes)

- Ts.seconds_cplex: time consumed by solving the Gurobi optimization model in seconds (sumed up in all passes)

- Ts.ticks_cplex: deterministic time consumed by solving the CPLEX optimization model in ticks (sumed up in all passes)

- Ts.seconds_solutions: time consumed by creating the final solutions in seconds (summed up in all passes)

- Ts.seconds_all: time consumed by the whole function in seconds

- Ts.seconds_1stPass: time consumed in 1st pass

- Ts.seconds_cplex_1stPass, Ts.ticks_cplex_1stPass: a matrix containing cplex seconds, ticks for each region

- Ts.ticks_cplex_2ndPassMIP, Ts.seconds_cplex_2ndPassMIP: time consumed by CPLEX in 2nd MIP pass

- Ts.seconds_2ndPassMIP: total time in 2nd MIP pass

- Ts.seconds_cplex_2ndPassLP, Ts.seconds_ticks_2ndPassLP: time consumed by CPLEX in 2nd LP pass by blocks

- Ts.seconds_cplex_2ndPassLP_region, Ts.seconds_ticks_2ndPassLP_region: time consumed by CPLEX in 2nd LP pass by regions

- Ts.seconds_2ndPassLP: total time in 2nd LP pass

- Ts.seconds_cplex_LP, Ts.ticks_cplex_LP: time consumed by CPLEX in pure LP methods

imgs_dct2: DCT coefficients of some images (with differentially encoded DC coefficients if DC prediction is involved)

- imgs_dct2.y0: DCT coefficients of the input (original) image

- imgs_dct2.y: 1st pass

- imgs_dct2.y1: MIP 2nd pass

- imgs_dct2.y2: LP 2nd pass by blocks

The above are useful for producing the figure on accuracy of sign bit estimation using the optimization model.

SB: the recovered sign bit of DCT coefficient

  sign bit = 1 if DCT coefficient>0

  sign bit = 0 if DCT coefficient<=0

- SB.x: 1st pass

- SB.x1: MIP 2nd pass

- SB.x2: LP 2nd pass by blocks(SB is a floating point number with 0<=SB<=1 in this case)

PSNRs, SSIMs: indicators of visual quality of recovered images as measured by PSNR and SSIM BERs:

Bit Error Rates of different recovery methods

    BER = Number of wrong bits / Total number of bits

Note BERs are calculated with SB even for LP 2nd passes
