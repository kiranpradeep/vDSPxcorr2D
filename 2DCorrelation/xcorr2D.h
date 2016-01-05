// correlation(img, kernel) = inv( fft(img) * conjugate(fft(kernel)) )
#ifndef correlate_h
#define correlate_h

#include <stdio.h>

#include <math.h>
#include "Accelerate/Accelerate.h"

//set up fft weights array
FFTSetup setup_fft_weights( int max_dimension );

// correlation as O(n log n) operation.
void xcorr2D( FFTSetup fft_weights, const int img_rows, const int img_cols, const float img[img_rows][img_cols],  const int kernel_rows, const int kernel_cols, const float kernel[kernel_rows][kernel_cols], float corr[img_rows-kernel_rows+1][img_cols-kernel_cols+1] );

// correlation as O(N^2) operation.
void convolve2D_slow( const vImage_Buffer* img, const vImage_Buffer* corr, const float* kernel, const int kernel_rows, const int kernel_cols );

//tear down
void teardown_fft_weights( FFTSetup fft_setup );


#endif /* correlate_h */
