//=================================================================================
//  slicpython.c
//=================================================================================


#include <stdbool.h>
#include "slic.c"
    

void SLICmain(  double* img, const int width, const int height,
                const int nchannels, const int numSuperpixels, const double compactness,
                const bool doRGBtoLAB, int* klabels, int* numlabels)
{
    int sz = width*height;
    double** channels = (double**)malloc(sizeof(double*)*nchannels);
    for(int c = 0; c < nchannels; c++)
    {
        channels[c] = img + c*sz;
    }
    //---------------------------
    // Perform color conversion
    //---------------------------
    if(doRGBtoLAB && nchannels==3)
    {
        rgbtolab(channels[0],channels[1],channels[2],sz,channels[0],channels[1],channels[2]);
    }
    //---------------------------
    // Find seeds
    //---------------------------
    int numseeds = 0;
    int step = sqrt((double)(sz)/(double)(numSuperpixels))+0.5;

    int* seedIndices = (int*)malloc( sizeof(int)* (int)(numSuperpixels * 1.1 + 10) ); // in case the number of seeds is greater than requested
    getSeeds(numSuperpixels,width,height,seedIndices,&numseeds);

    double*  kseedsx    = (double* )malloc( sizeof(double)      * numseeds ) ;
    double*  kseedsy    = (double* )malloc( sizeof(double)      * numseeds ) ;
    double** kseedsc    = (double**)malloc( sizeof(double*)     * nchannels) ;
    for(int c = 0; c < nchannels; c++)
    {
        kseedsc[c] = (double*)malloc( sizeof(double) * numseeds );
    }
    for(int k = 0; k < numseeds; k++)
    {
        kseedsx[k] = seedIndices[k]%width;
        kseedsy[k] = seedIndices[k]/width;
        for(int c = 0; c < nchannels; c++)
        {
            kseedsc[c][k] = channels[c][seedIndices[k]];
        }
    }
    //---------------------------
    // Create superpixels
    //---------------------------
    PerformSuperpixelSLIC(channels,kseedsx,kseedsy,kseedsc,width,height,nchannels,numseeds,klabels,step,compactness);
    //---------------------------
    // Enforce connectivity
    //---------------------------
    int finalNumberOfLabels = 0;
    int* clabels = (int*)malloc( sizeof(int) * sz );//corrected labels after enforcing connectivity

    EnforceSuperpixelConnectivity(klabels,width,height,numSuperpixels,clabels,&finalNumberOfLabels);
    //---------------------------
    // Copy connected labels back
    //---------------------------
    for(int i = 0; i < sz; i++)
    {
        klabels[i] = clabels[i];
    }
    *numlabels = finalNumberOfLabels;

    for(int c = 0; c < nchannels; c++)
    {
        free(kseedsc[c]);
    }
    free(kseedsc);
    free(kseedsx);
    free(kseedsy);
    free(clabels);
    free(seedIndices);
    free(channels);
}

