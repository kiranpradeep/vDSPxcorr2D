// sample for 2D cross orrelation

#include <stdio.h>
#include <sys/time.h>
#include "xcorr2D.h"

#define TEST_PERFORMANCE
//#define SAMPLE_RUN

void demo_run_with_results();
void test_performance();
uint64_t get_time_usec();
void print2D( int rows, int cols, float corr[rows][cols] );

int main()
{

#ifdef TEST_PERFORMANCE
    //test performance by comparing timing with vImage_convole method O(n^2)
    test_performance();
#else
    //sample usage with results
    demo_run_with_results();
#endif
    
    return 0;
}

void demo_run_with_results()
{
    //max dimension for setting fft weights
    //2^11 = 2048 and so adjust MAX_DIMENESION as your needs
    const int MAX_DIMENSION = 11;
    
    //dimensions of image
    const int IMG_ROWS = 5;
    const int IMG_COLS = 5;
    
    //dimensions of kernel
    const int KERNEL_ROWS = 3;
    const int KERNEL_COLS = 3;
    
    //dimensions of correlation result
    const int CORR_ROWS = IMG_ROWS - KERNEL_ROWS + 1;
    const int CORR_COLS = IMG_COLS - KERNEL_COLS + 1;
    float corr[CORR_ROWS][CORR_COLS];
    
    const float img[IMG_ROWS][IMG_COLS]= { {17,24,1,8,15}, {23,5,7,14,16}, {4,6,13,20,22}, {10,12,19,21,3}, {11,18,25,2,9} };
    const float kernel[KERNEL_ROWS][KERNEL_COLS]= { {8,1,6}, {3,5,7}, {4,9,2}  };
    
    //const float img[IMG_ROWS][IMG_COLS]= { {1,1,1,1}, {1,1,1,1}, {1,1,1,1}, {1,1,1,1} };
    //const float kernel[KERNEL_ROWS][KERNEL_COLS]= { {1,2,3} } ;
    
    //setup weights
    FFTSetup fft_weights = setup_fft_weights( MAX_DIMENSION );
    
    //do correlation
    xcorr2D( fft_weights, IMG_ROWS, IMG_COLS, img, KERNEL_ROWS, KERNEL_COLS, kernel, corr );
    
    //print result
    print2D(CORR_ROWS, CORR_COLS, corr);
    
    //cleanup
    teardown_fft_weights( fft_weights );
}

void test_performance()
{
    const int MAX_TESTS=7;
    //dimensions at which performances are compared
    float sizeImg[MAX_TESTS] =    {20, 32,  128, 512,  1024, 2048, 2048 };
    float sizeKernel[MAX_TESTS] = {8-1,10-1,32-1,128-1,256-1,256-1,1024-1};
    
    for( int i=0; i<MAX_TESTS; i++ ) {
        const int img_rows = sizeImg[i];
        const int img_cols = sizeImg[i];
        
        const int kernel_rows = sizeKernel[i];
        const int kernel_cols = sizeKernel[i];
        
        const int corr_rows = img_rows - kernel_rows + 1;
        const int corr_cols = img_cols - kernel_cols + 1;
        
        //alloc mem for image, kernel and result - actually 2D arrays
        float* const img = (float*)malloc(img_rows * img_cols * sizeof(float));
        float* const kernel= (float*)malloc(kernel_rows * kernel_cols * sizeof(float));
        float* const corr= (float*)calloc(corr_rows * corr_cols, sizeof(float));
        
        //random init image array
        for(int i=0; i<img_rows * img_cols; i++) {
            img[i] = rand() % 1000;
        }
        
        //random init kernel array
        for(int i=0; i<kernel_rows * kernel_cols; i++) {
            kernel[i] = rand() % 1000;
        }

        //setup weights
        const int max = (int)ceilf(log2f(img_rows));
        FFTSetup fft_weights = setup_fft_weights( max );
        
        //do correlation - O(n log n)
        const uint64_t t1 = get_time_usec();
        xcorr2D( fft_weights, img_rows, img_cols, img, kernel_rows, kernel_cols, kernel, corr );
        const uint64_t t2 = get_time_usec();
        
        //cleanup
        teardown_fft_weights( fft_weights );
        
        //do correlation - O(n^2) - with vImage_convole method
        float* corr_slow = (float*)calloc(img_rows * img_cols, sizeof(float));
        const vImage_Buffer imgBuf  = { img, img_rows, img_cols, img_cols * sizeof(float) };
        const vImage_Buffer resBuf  = { corr_slow, img_rows, img_cols, img_cols * sizeof(float) };
        
        const uint64_t t3 = get_time_usec();
        convolve2D_slow( &imgBuf, &resBuf, kernel, kernel_rows, kernel_cols );
        const uint64_t t4 = get_time_usec();
        
        printf("%5d - %5d - %9.6f- %9.6f\n", img_rows, kernel_rows, (t2-t1)/1e6, (t4-t3)/1e6 );
        
        //free image, kernel and result
        free(img);
        free(kernel);
        free(corr);
        free(corr_slow);
    }

}

void print2D( int rows, int cols, float result[rows][cols] )
{
    for(int i=0; i<rows; ++i) {
        for(int j=0; j<cols; ++j) {
            printf("%9.1f ", result[i][j]);
        }
        printf("\n");
    }
    printf("-------------\n");
}

uint64_t get_time_usec()
{
    struct timeval time;
    if(gettimeofday(&time, NULL) < 0) {
        return 0;
    }
    return (time.tv_sec * 1000 * 1000) + time.tv_usec;
}

