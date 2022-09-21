# Recovery of Sign Bits of DCT Coefficients

Developed by Ruiyuan Lin, Jun Jiang, Sheng Liu, and Shujun Li (in order of contribution).

## Installation

1. Install MATLAB

2. Install CPLEX >= 12.6 and its MATLAB interface

3. Add path '/MIP_private' and '/jpeg_toolbox' in MATLAB:

   ```matlab
   addpath('./MIP_private');
   addpath('./jpeg_toolbox');
   ```

## Main function

- `cplex_signbit`: This function invokes IBM ILOG CPLEX’s MATLAB Cplex class to solve the problem of recovering unknown sign bits of DCT coefficients of a DCT-transformed image with a **linear programming** model. 
  Example:

  ```matlab
  [imgs, Ts, imgs_dct2, PSNRs, SSIMs, BERs] = cplex_signbit(img_dct2, mask, DCpre_mode, DCTabs_min, relaxX, relaxZ, recovery_method, bDClevelshifted, H, W, Xmax)
  ```

- `cplex_signbitHierarchyDCACMIP`: This function solves the problem of recovering unknown sign bits of DCT coefficients of a DCT-transformed image ‘image_dct2’ region by region (each by a **mixed integer programming** model).

  Example:

  ```matlab
  [imgs, Ts, imgs_dct2, SB, PSNRs, SSIMs, BERs] = cplex_signbitHierarchyDCACMIP(img_dct2, mask, regionSize_H, regionSize_W, DCpred_mode, DCdepend, relaxXBound, relaxZBound, DCpassByRegion, DCTabs_min, timeLimit, relatGap, DCpass_method, recovery_method, relaxXBound_LP, relaxZBound_LP, bDClevelshifted, H, W, Xmax, fileLabel)
  ```

## Batch scripts

Running the following scripts with some modification to reproduce the experimental result:

- `batch_LP.m`
- `batch_LP_JPEG.m`
- `batch_MIP.m`

For more detailed description, see [document.md](/document.md).

## License

This project is released under the [Apache 2.0 license](LICENSE).

## Citation
<To be added>
