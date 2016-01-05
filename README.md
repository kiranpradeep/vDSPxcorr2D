##vDSPxcorr2D##

A small sample to do FFT based 2D cross-correlation on OSX/iOS using Accelerate framework.  Starting point of this sample is [vDSPxcorr](https://github.com/kunalkandekar/vDSPxcorr/) project.

Methods  in Accelerate framework like `vDSP_imgfir` or `vImage_convolve` does 2D convolutions. But both methods currently are O(n<sup>2</sup>) algorithms. This sample uses FFT based calculations and hence is O( n log<sub>2</sub>n ) 

>crosscorr(img, kernel) = invFFT( FFT(img) * CONJUGATE(FFT(kernel)) )


And so for computing cross correlation from FFT based calculations, this sample internally uses Accelerate framework provided methods like `vDSP_fft2d_zip`, `vDSP_zmmul` etc. Code would have been cleaner if we used complex-complex FFTs, or out of place operations or different arrays for input and output to FFT. But each of those would double the memory and there may be earlier than necessary overflows from L2 cache. vDSP routines work best when data fits in L1/L2 caches. So, a real to real in place FFT with same arrays for input and output are used. This effects readability as we have to deal with 2D array special packing.

##Timing Comparison##

These timings( code inside ) are with out repetitions. But the time difference is so huge( at higher dimensions ) that reps are unnecessary. Tests were taken on a 2010 iMac 3.06GHz Intel Core i3. Dimensions of image and kernel given. Image and kernel are considered 2D square matrices for tests. Time in seconds.

    +-------+-------+------------+----------------+
    |  IMG  | KERNEL| vDSPxcorr2D|vImage_convolve |
    +-------+-------+------------+----------------+
    |    20 |     7 |   0.000099 |       0.000028 |
    |    32 |     9 |   0.000032 |       0.000134 |
    |   128 |    31 |   0.000570 |       0.001395 |
    |   512 |   127 |   0.015879 |       0.220793 |
    |  1024 |   255 |   0.069091 |       2.677988 |
    |  2048 |   255 |   0.506878 |      16.996459 |
    |  2048 |  1023 |   0.515079 |      76.712189 |
    +-------+-------+------------+----------------+


##TODO##
* The sample is only for single precision data. Yet to do for double precision.
* Kernel FFTs could be cached for future operations - improves performance on further runs
* No validation for wrong input. E.g. kernel size higher than image size etc
* Not well tested code. Github issues are welcome.

##Thanks##
To the really helpful Vector and Numerics Group who provided replies[ in Apple lists ] to every query and most of this based on their inputs. But, any wrong information here would be for me only.
