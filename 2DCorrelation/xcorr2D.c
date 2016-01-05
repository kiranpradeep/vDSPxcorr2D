#include "xcorr2D.h"

FFTSetup setup_fft_weights( int max_dimension )
{
    return vDSP_create_fftsetup(max_dimension, FFT_RADIX2);
}

void teardown_fft_weights( FFTSetup fft_setup )
{
    vDSP_destroy_fftsetup(fft_setup);
}

void xcorr2D( FFTSetup fft_weights, const int img_rows, const int img_cols, const float img[img_rows][img_cols],  const int kernel_rows, const int kernel_cols, const float kernel[kernel_rows][kernel_cols], float corr[img_rows-kernel_rows+1][img_cols-kernel_cols+1] )
{
    //dft size - use these rows/cols if correlation is needed rows/cols near border
    //const int dftRows = img_rows + kernel_rows - 1 ;
    //const int dftCols = img_cols + kernel_cols - 1 ;
    
    //expand to power of 2
    //const int rows = 1 << (int)ceilf(log2f(dftRows));
    //const int cols = 1 << (int)ceilf(log2f(dftCols));
    const int rows = 1 << (int)ceilf(log2f(img_rows));
    const int cols = 1 << (int)ceilf(log2f(img_cols));
    
    //log to base 2 of rows and cols
    const int log2_cols = log2( cols ) ;
    const int log2_rows = log2( rows ) ;
    
    //total number of elements
    const int n = rows * cols;
    
    //half of total nuber of elements
    const int nOver2 = n / 2;
    
    //FFT row and column strides
    const int row_stride = 1;
    const int column_stride = 0;
    
    //correlation = inv( fft(img) * conjugate(fft(kernel)) )
    const int conjugate = -1;
    
    //compensations in scaling for forward and revert fft
    const float scaleForwardFft =  0.5f;
    const float scaleInverseFft = 1.0f / n;
    
    //allocate mem for kernel(n), image(n) and result of first two columns(2 * rows) and zero initialize
    float* const fft_mem = ( float* ) calloc ( ( 2 * n +  2 * rows ), sizeof( float ) );
    
    //pointers to store kernel in split odd even format
    float* const kernel_real = fft_mem;
    float* const kernel_imag = kernel_real + nOver2;
    
    //pointers to store image in split odd even format
    float* const img_real = fft_mem + n;
    float* const img_imag = img_real + nOver2;
    
    //pointers to store first columns in split odd even format
    float* const temp_real = fft_mem + 2 * n;
    float* const temp_imag = temp_real + rows;
    
    const COMPLEX_SPLIT kernel_fft = { kernel_real, kernel_imag };
    const COMPLEX_SPLIT img_fft = { img_real, img_imag };
    const COMPLEX_SPLIT temp_fft = { temp_real, temp_imag };
    
    //store kernel in odd even split format.
    for( int i=0; i<kernel_rows; ++i) {
        //process each row of kernel
        const float* const row = kernel[i];
        const COMPLEX_SPLIT splitKernel = { &kernel_real[i * cols/2], &kernel_imag[i * cols/2]};
        vDSP_ctoz(( const DSPComplex *)row, 2, &splitKernel, 1, kernel_cols/2);
        //if num of columns is odd, store the last element into real part.
        if( kernel_cols % 2 == 1) {
            kernel_real[i*cols/2+(kernel_cols/2)] = kernel[i][kernel_cols-1];
        }
    }
    
    //store image in odd even split format.
    for( int i=0; i<img_rows; ++i) {
        //process each row of image
        const float* const row = img[i];
        
        const COMPLEX_SPLIT splitImg = { &img_real[i * cols/2], &img_imag[i * cols/2]};
        vDSP_ctoz((const DSPComplex *)row, 2, &splitImg, 1, img_cols/2);
        //if num of columns is odd, store the last element into real part.
        if( img_cols % 2 == 1) {
            img_real[i*cols/2+(img_cols/2)] = img[i][img_cols-1];
        }
    }
    
    //fft kernel
    vDSP_fft2d_zrip(fft_weights, &kernel_fft, row_stride, column_stride, log2_rows, log2_cols, FFT_FORWARD);
    
    //fft image
    vDSP_fft2d_zrip(fft_weights, &img_fft, row_stride, column_stride, log2_rows, log2_cols, FFT_FORWARD);
    
    //correct vDSP forward scaling of kernel - divide by 2
    vDSP_vsmul(kernel_fft.realp, 1, &scaleForwardFft, kernel_fft.realp, 1, nOver2);
    vDSP_vsmul(kernel_fft.imagp, 1, &scaleForwardFft, kernel_fft.imagp, 1, nOver2);
    
    //correct vDSP forward scaling of image - divide by 2
    vDSP_vsmul(img_fft.realp, 1, &scaleForwardFft, img_fft.realp, 1, nOver2);
    vDSP_vsmul(img_fft.imagp, 1, &scaleForwardFft, img_fft.imagp, 1, nOver2);
    
    //We have now completed taking FFT of image and kernel and now need to multilply. But that will produce garbage in the first column but complete all the other columns correctly. So overwrite the bad first column of the output with correct results: Do the four real multiplications for the four real elements and use a couple of vDSP_zvmul calls with strides to handle the weirdly-arranged elements in the first column of the reals and the first column of the imaginaries. [ Thanks super helpful Scitech list].
    
    const int one_row  = cols/2;
    const int two_row_stride = cols;
    
    // to store - first column of reals
    const COMPLEX_SPLIT temp_real_col= { temp_real, temp_real + 1  };
    const COMPLEX_SPLIT img_real_col = { img_real, img_real + one_row } ;
    const COMPLEX_SPLIT kernel_real_col= { kernel_real, kernel_real  + one_row  };
    
    //to store - first column of imaginary
    const COMPLEX_SPLIT temp_imag_col= { temp_imag, temp_imag + 1  };
    const COMPLEX_SPLIT img_imag_col = { img_imag, img_imag + one_row } ;
    const COMPLEX_SPLIT kernel_imag_col= { kernel_imag, kernel_imag + one_row  };
    
    // Since we are not using a seperate result array, temporarily save first column to be copied later. Store first column of real/image as below.
    //store - first column of reals
    vDSP_zvmul(&img_real_col, two_row_stride, &kernel_real_col, two_row_stride, &temp_real_col, 2, rows/2, conjugate);
    //store - first column of imaginary
    vDSP_zvmul(&img_imag_col, two_row_stride, &kernel_imag_col, two_row_stride, &temp_imag_col, 2, rows/2, conjugate);
    
    // multiply four reals
    temp_fft.realp[0] = img_fft.realp[0] * kernel_fft.realp[0];
    temp_fft.imagp[0] = img_fft.imagp[0] * kernel_fft.imagp[0];
    temp_fft.realp[1] = img_fft.realp[cols/2] * kernel_fft.realp[cols/2];
    temp_fft.imagp[1] = img_fft.imagp[cols/2] * kernel_fft.imagp[cols/2];
    
    //multiply entire 2D array
    vDSP_zvmul(&img_fft, 1, &kernel_fft, 1, &img_fft, 1, nOver2, conjugate);
    
    //copy first two cols stored in temp
    vDSP_zvmov(&temp_fft, 1, &img_fft, cols/2, rows);
    
    //inverse fft
    vDSP_fft2d_zrip(fft_weights, &img_fft, row_stride, column_stride, log2_rows, log2_cols, kFFTDirection_Inverse);
    
    //correct vDSP inverse scaling - divide by n
    vDSP_vsmul(img_fft.realp, 1, &scaleInverseFft, img_fft.realp, 1, nOver2);
    vDSP_vsmul(img_fft.imagp, 1, &scaleInverseFft, img_fft.imagp, 1, nOver2);
    
    //vertical flip on result
    const int rowsOver2 = rows/2;
    const int colsOver2 = cols/2;
    for(int idx=1; idx<rowsOver2; ++idx) {
        const int j = idx * colsOver2;
        const int k = (rows - idx) * colsOver2;
        vDSP_vswap( &img_fft.realp[k], 1, &img_fft.realp[j], 1, colsOver2);
        vDSP_vswap( &img_fft.imagp[k], 1, &img_fft.imagp[j], 1, colsOver2);
    }
    
    //horizontal flip on result except first real col
    for(int idx=0; idx<colsOver2/2; ++idx) {
        const int k = colsOver2 - idx - 1;
        vDSP_vswap( &img_fft.imagp[idx], colsOver2, &img_fft.imagp[k], colsOver2, rows);
        if( idx > 0 ) {
            vDSP_vswap( &img_fft.realp[idx], colsOver2, &img_fft.realp[k+1], colsOver2, rows);
        }
    }
    
    const int corr_rows = img_rows - kernel_rows + 1;
    const int corr_cols = img_cols - kernel_cols + 1;
    
    //convert correlation result in odd even split format to normal 2d array representation
    for( int i=0; i<corr_rows; ++i) {
        float* const row = corr[i];
        const COMPLEX_SPLIT splitResult = { &img_real[i * cols/2], &img_imag[i * cols/2] };
        vDSP_ztoc(&splitResult, 1, (COMPLEX *)row, 2, corr_cols/2);
        if( corr_cols % 2 == 1) {
            corr[i][corr_cols-1] = img_real[ i * cols/2 + corr_cols/2];
        }
    }
    
    free( fft_mem );
}

void convolve2D_slow( const vImage_Buffer* img, const vImage_Buffer* corr, const float* kernel, const int kernel_rows, const int kernel_cols )
{
    vImageConvolve_PlanarF( img, corr, NULL, 0, 0, kernel, kernel_rows, kernel_cols, 0, kvImageCopyInPlace);
}
