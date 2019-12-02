//=================================================================================
//  slicmex.c
//=================================================================================

#include<mex.h>
#include "slic.c"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int width;
    int height;
    int nchannels;
    int sz;
    int i, ii;
    int x, y;
    int* rin; int* gin; int* bin;
    int* klabels;
    int* clabels;
    double* lvec; double* avec; double* bvec;
    double** channels;
    int step;
    int* seedIndices;
    int numseeds;
    double* kseedsx;double* kseedsy;
    double** kseedsc;
    int k;
    const mwSize* dims;//int* dims;
    mwSize numdims;
    int* outputNumSuperpixels;
    int* outlabels;
    int finalNumberOfLabels;
    unsigned char* imgbytes;
    int numelements;
    int numSuperpixels = 200;//default value
    double compactness = 10;//default value
    bool doRGBtoLAB = false;

    if (nrhs < 1) {
        mexErrMsgTxt("At least one argument is required.\n") ;
    } else if(nrhs > 4) {
        mexErrMsgTxt("Too many input arguments.\n");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("SLIC:nlhs","Two outputs required, a labels and the number of labels, i.e superpixels.");
    }
    //---------------------------
    numelements   = mxGetNumberOfElements(prhs[0]); // total number of pixels, for instance
    numdims = mxGetNumberOfDimensions(prhs[0]); // numdims is 2 for a grayscale image; the value is always 2 or greater.
    dims  = mxGetDimensions(prhs[0]) ;
    imgbytes  = (unsigned char*)mxGetData(prhs[0]) ;//mxGetData returns a void pointer, so cast it
    height = dims[0];//Note: first dimension provided is height and second is width
    width = dims[1];
    nchannels = 1;
    if(numdims > 2)
    {
        nchannels = dims[2];
    }

    sz = width*height;
    //---------------------------
    numSuperpixels  = mxGetScalar(prhs[1]);
    compactness     = mxGetScalar(prhs[2]);
    doRGBtoLAB      = mxGetScalar(prhs[3]);
    
    //---------------------------
    // Allocate memory
    //---------------------------
    channels = mxMalloc( sizeof(double*)*nchannels);
    for(int c = 0; c < nchannels; c++)
    {
        channels[c] = mxMalloc( sizeof(double) * sz );
    }

    klabels = mxMalloc( sizeof(int)         * sz );//original k-means labels
    clabels = mxMalloc( sizeof(int)         * sz );//corrected labels after enforcing connectivity
    seedIndices = mxMalloc( sizeof(int)     * (int)(numSuperpixels * 1.1 + 10) );
    
    
    for(x = 0, ii = 0; x < width; x++)//reading data from column-major MATLAB matrics to row-major C matrices (i.e perform transpose)
    {
        for(y = 0; y < height; y++,ii++)
        {
            for(int c = 0; c < nchannels; c++)
            {
                i = y*width+x;
                channels[c][i] = imgbytes[ii+sz*c];
            }
        }
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
    step = sqrt((double)(sz)/(double)(numSuperpixels))+0.5;
    getSeeds(numSuperpixels,width,height,seedIndices,&numseeds);
    
    kseedsx    = mxMalloc( sizeof(double)      * numseeds ) ;
    kseedsy    = mxMalloc( sizeof(double)      * numseeds ) ;
    kseedsc    = mxMalloc( sizeof(double*)     * nchannels) ;
    for(int c = 0; c < nchannels; c++)
    {
        kseedsc[c] = mxMalloc( sizeof(double) * numseeds );
    }
    for(k = 0; k < numseeds; k++)
    {
        kseedsx[k] = seedIndices[k]%width;
        kseedsy[k] = seedIndices[k]/width;
        for(int c = 0; c < nchannels; c++)
        {
            kseedsc[c][k] = channels[c][seedIndices[k]];
        }
    }
    //---------------------------
    // Compute superpixels
    //---------------------------
    PerformSuperpixelSLIC(channels,kseedsx,kseedsy,kseedsc,width,height,nchannels,numseeds,klabels,step,compactness);
    //---------------------------
    // Enforce connectivity
    //---------------------------
    EnforceSuperpixelConnectivity(klabels,width,height,numSuperpixels,clabels,&finalNumberOfLabels);
    //---------------------------
    // Assign output labels
    //---------------------------
    plhs[0] = mxCreateNumericMatrix(height,width,mxINT32_CLASS,mxREAL);
    outlabels = mxGetData(plhs[0]);
    for(x = 0, ii = 0; x < width; x++)//copying data from row-major C matrix to column-major MATLAB matrix (i.e. perform transpose)
    {
        for(y = 0; y < height; y++)
        {
            i = y*width+x;
            outlabels[ii] = clabels[i];
            ii++;
        }
    }
    //---------------------------
    // Assign number of labels/seeds
    //---------------------------
    plhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    outputNumSuperpixels = (int*)mxGetData(plhs[1]);//gives a void*, cast it to int*
    *outputNumSuperpixels = finalNumberOfLabels;
    //---------------------------
    // Deallocate memory
    //---------------------------
    for(int c = 0; c < nchannels; c++)
    {
        mxFree(channels[c]);
        mxFree(kseedsc[c]);
    }
    mxFree(channels);
    mxFree(kseedsc);
    mxFree(kseedsx);
    mxFree(kseedsy);
    mxFree(klabels);
    mxFree(clabels);
    mxFree(seedIndices);
}

