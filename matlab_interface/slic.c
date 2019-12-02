//=================================================================================
//  slic.c
//=================================================================================

//------------------------------------------------------------------------
// The following should make all the code to work for both C and C++
// cases, i.e., when the compiler is gcc or g++.
//------------------------------------------------------------------------
#ifdef __cplusplus
extern "C"{
#endif

#include <stdlib.h>
#include <math.h>
#include <float.h>

//=================================================================================
/// rgbtolab
//=================================================================================
//void rgbtolab(int* rin, int* gin, int* bin, int sz, double* lvec, double* avec, double* bvec)
void rgbtolab(double* rin, double* gin, double* bin, int sz, double* lvec, double* avec, double* bvec)
{
    double sR, sG, sB;
    double R,G,B;
    double X,Y,Z;
    double r, g, b;
    const double epsilon = 0.008856;	//actual CIE standard
    const double kappa   = 903.3;		//actual CIE standard
    
    const double Xr = 0.950456;	//reference white
    const double Yr = 1.0;		//reference white
    const double Zr = 1.088754;	//reference white
    double xr,yr,zr;
    double fx, fy, fz;
    double lval,aval,bval;
    
    for(int i = 0; i < sz; i++)
    {
        sR = rin[i]; sG = gin[i]; sB = bin[i];
        R = sR/255.0;
        G = sG/255.0;
        B = sB/255.0;
        
        if(R <= 0.04045)	r = R/12.92;
        else				r = pow((R+0.055)/1.055,2.4);
        if(G <= 0.04045)	g = G/12.92;
        else				g = pow((G+0.055)/1.055,2.4);
        if(B <= 0.04045)	b = B/12.92;
        else				b = pow((B+0.055)/1.055,2.4);
        
        X = r*0.4124564 + g*0.3575761 + b*0.1804375;
        Y = r*0.2126729 + g*0.7151522 + b*0.0721750;
        Z = r*0.0193339 + g*0.1191920 + b*0.9503041;
        
        //------------------------
        // XYZ to LAB conversion
        //------------------------
        xr = X/Xr;
        yr = Y/Yr;
        zr = Z/Zr;
        
        if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
        else				fx = (kappa*xr + 16.0)/116.0;
        if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
        else				fy = (kappa*yr + 16.0)/116.0;
        if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
        else				fz = (kappa*zr + 16.0)/116.0;
        
        lval = 116.0*fy-16.0;
        aval = 500.0*(fx-fy);
        bval = 200.0*(fy-fz);
        
        lvec[i] = lval; avec[i] = aval; bvec[i] = bval;
    }
}

//=================================================================================
/// getSeeds
///
/// The function computes the grid step and uses it to find seed indices. The
/// indices are used to extract the feature values (like [l,a,b,,x,y]).
//=================================================================================
void getSeeds(int numk, int width, int height, int* seedIndices, int* numseeds)
{
    const int sz = width*height;
    double gridstep = sqrt((double)(sz)/(double)(numk)) + 0.5;
    if(1)
    {
        int minerr = abs( (int)(0.5 + width/gridstep)*(int)(0.5 + height/gridstep) - numk);
        double minstep = gridstep-1.0; double maxstep = gridstep+1.0;
        for(double x = minstep; x <= maxstep; x += 0.1)
        {
            int err = abs( (int)(0.5 + width/x)*(int)(0.5 + height/x) - numk);
            if(err < minerr)
            {
                minerr = err;
                gridstep = x;
            }
        }
    }

    double halfstep = gridstep/2.0;

    int n = 0;
    for(double y = halfstep; y <= height; y += gridstep)
    {
        int yval = y+0.5;
        if(yval < height)
        {
            for(double x = halfstep; x <= width; x += gridstep)
            {
                int xval = x+0.5;
                if(xval < width)
                {
                    seedIndices[n] = yval*width + xval;
                    n++;
                }
            }
        }
    }

    *numseeds = n;
}

//=================================================================================
/// PerformSuperpixelSLIC
//=================================================================================
void PerformSuperpixelSLIC(double** channels, double* kseedsx, double* kseedsy, double** kseedsc, int width, int height, int nchannels, int numseeds, int* klabels, int STEP, double compactness)
{
    int x1, y1, x2, y2;
	double l, a, b;
	double distc;
	double distxy;
    double dist;
    double diff;
    int itr;
    int n;
    int x,y;
    int i;
    int ind;
    //int r,c;
    int c;
    int k;
    int sz = width*height;
	const int numk = numseeds;
	int offset = STEP;
    
    double* clustersize = (double*)malloc(sizeof(double)*numk);
    double* inv         = (double*)malloc(sizeof(double)*numk);
    double* sigmax      = (double*)malloc(sizeof(double)*numk);
    double* sigmay      = (double*)malloc(sizeof(double)*numk);
    double** sigmac     = (double**)malloc(sizeof(double*)*nchannels);
    for(c = 0; c < nchannels; c++)
    {
        sigmac[c] = (double*)malloc(sizeof(double)*numk);
    }
	double invwt = 1.0/((STEP/compactness)*(STEP/compactness));
    double* distvec     = (double*)malloc(sizeof(double)*sz);
    
	for( itr = 0; itr < 10; itr++ )
	{
		for(i = 0; i < sz; i++){distvec[i] = DBL_MAX;}
     
		for( n = 0; n < numk; n++ )
		{
            x1 = kseedsx[n]-offset; if(x1 < 0) x1 = 0;
            y1 = kseedsy[n]-offset; if(y1 < 0) y1 = 0;
            x2 = kseedsx[n]+offset; if(x2 > width)  x2 = width;
            y2 = kseedsy[n]+offset; if(y2 > height) y2 = height;
            
			for( y = y1; y < y2; y++ )
			{
				for( x = x1; x < x2; x++ )
				{
					i = y*width + x;
                    
                    distxy = (x - kseedsx[n])*(x - kseedsx[n]) + (y - kseedsy[n])*(y - kseedsy[n]);
					distc =	0;
                    for(c = 0; c < nchannels; c++)
                    {
                        diff = channels[c][i]-kseedsc[c][n];
                        distc += diff*diff;
                    }
					//--------------------------
					dist = distc + distxy*invwt;
                    //--------------------------
					if(dist < distvec[i])
					{
						distvec[i] = dist;
						klabels[i]  = n;
					}
				}
			}
		}
        
		//-----------------------------------------------------------------
		// Reset the centroid in order to compute the new ones
		//-----------------------------------------------------------------
        {for(k = 0; k < numk; k++)
        {
            for(c = 0; c < nchannels; c++)
            {
                sigmac[c][k] = 0;
            }
            sigmax[k] = 0;
            sigmay[k] = 0;
            clustersize[k] = 0;
        }}
        
        //-----------------------------------------------------------------
        // Recalculate the centroid and store in the seed values
        //-----------------------------------------------------------------
		
        {for( y = 0, ind = 0; y < height; y++ )
        {
            for( x = 0; x < width; x++, ind++ )
            {
                if(klabels[ind] >= 0)
                {
                    for(c = 0; c < nchannels; c++)
                    {
                        sigmac[c][klabels[ind]] += channels[c][ind];
                    }
                    sigmax[klabels[ind]] += x;
                    sigmay[klabels[ind]] += y;
                    clustersize[klabels[ind]] += 1.0;
                }
            }
        }}
        
		{for( k = 0; k < numk; k++ )
		{
			if( clustersize[k] <= 0 ) clustersize[k] = 1;
			inv[k] = 1.0/clustersize[k];//computing inverse now to multiply, than divide later
		}}
		
		{for( k = 0; k < numk; k++ )
		{
            for(c = 0; c < nchannels; c++)
            {
                kseedsc[c][k] = sigmac[c][k]*inv[k];
            }
			kseedsx[k] = sigmax[k]*inv[k];
			kseedsy[k] = sigmay[k]*inv[k];
		}}
        
	}
    for(c = 0; c < nchannels; c++)
    {
        free(sigmac[c]);
    }
    free(sigmac);
    free(sigmax);
    free(sigmay);
    free(clustersize);
    free(inv);
    free(distvec);

}

//=================================================================================
/// EnforceSuperpixelConnectivity
///
/// This function uses the queue based flood-fill algorithm to detect segments
/// smaller than a threshold and absorb them directly into a previous superpixel.
/// This step is needed to avoid "orphan" clusters, which are pixels belonging to
/// a certain superpixel, but lying outside of it.
//=================================================================================
void EnforceSuperpixelConnectivity(int* labels, int width, int height, int numSuperpixels,int* nlabels, int* finalNumberOfLabels)
{
    int i,j,k;
    int n,c,count;
    int x,y;
    int ind;
    int oindex, adjlabel;
    int label;
    const int dx4[4] = {-1,  0,  1,  0};
	const int dy4[4] = { 0, -1,  0,  1};
    const int sz = width*height;
    const int SUPSZ = sz/numSuperpixels;
    int* xvec = (int*)malloc(sizeof(int)*SUPSZ*10);
	int* yvec = (int*)malloc(sizeof(int)*SUPSZ*10);

	for( i = 0; i < sz; i++ ) nlabels[i] = -1;
    oindex = 0;
    adjlabel = -1;//adjacent label
    label = 0;
	for( j = 0; j < height; j++ )
	{
		for( k = 0; k < width; k++ )
		{
			if( 0 > nlabels[oindex] )
			{
				nlabels[oindex] = label;
				//--------------------
				// Start a new segment
				//--------------------
				xvec[0] = k;
				yvec[0] = j;
				//-------------------------------------------------------
				// Quickly find an adjacent label for use later if needed
				//-------------------------------------------------------
				for( n = 0; n < 4; n++ )
				{
					int x = xvec[0] + dx4[n];
					int y = yvec[0] + dy4[n];
					// if( ((x >= 0) && (x < width)) && ((y >= 0) && (y < height)) )
                    if(!(x < 0 || x >= width || y < 0 || y >= height))
					{
						int nindex = y*width + x;
						if(nlabels[nindex] >= 0)
                        {
                            adjlabel = nlabels[nindex];
                        }
					}
				}
                
				count = 1;
				for( c = 0; c < count; c++ )
				{
					for( n = 0; n < 4; n++ )
					{
						x = xvec[c] + dx4[n];
						y = yvec[c] + dy4[n];
                        
						// if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                        if(!(x < 0 || x >= width || y < 0 || y >= height))
						{
							int nindex = y*width + x;
                            
							if( (0 > nlabels[nindex]) && (labels[oindex] == labels[nindex]) )
							{
								xvec[count] = x;
								yvec[count] = y;
								nlabels[nindex] = label;
								count++;
							}
						}
                        
					}
				}
				//-------------------------------------------------------
				// If segment size is less then a limit, assign an
				// adjacent label found before, and decrement label count.
				//-------------------------------------------------------
				if(count <= (SUPSZ >> 2))
				{
					for( c = 0; c < count; c++ )
					{
                        ind = yvec[c]*width+xvec[c];
						nlabels[ind] = adjlabel;
					}
					label--;
				}
				label++;
			}
			oindex++;
		}
	}
	*finalNumberOfLabels = label;
    
	free(xvec);
	free(yvec);
}

#ifdef __cplusplus
}
#endif