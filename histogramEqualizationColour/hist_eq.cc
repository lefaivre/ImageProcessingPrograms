/*
	Adam Lefaivre - 001145679
	Assignment 4 - CPSC 5990 (D.I.P.)
	Programming Question 
	Histogram Equalization - Color 
	Dr. Howard Cheng
*/

#include <vector>
#include <pam.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <limits>

using namespace std;

// Here we define a mapping between the intensity itself,
// the number of pixels in the image with that intensity,
// and the intensity that we map each input intensity to.
struct mapping
{
	int r;
	int n;
	int s;
};

// globals
const int MAXINTENSITY = 256;
const double PI = 3.14159265;
vector<mapping> mappingVals;
double epsilon = numeric_limits<double>::epsilon();

// Function stubs.
tuple ** read_image(char *filename, pam &inpam);
void write_image(char *filename, const pam &inpam, tuple **array);
vector<vector<vector<double> > > RGB_HSI(const pam &inpam, tuple **array);
void HSI_RGB(const pam &inpam, tuple **array,
 vector<vector<vector<double> > > & HSI);
void fillMappingValues(const pam &inpam, tuple **array,
 vector<vector<vector<double> > > & HSI);
void snapRGB(double & R, double & G, double & B);
void mapValuesToNewImage(const pam &inpam, tuple **array,
 vector<vector<vector<double> > > & HSI);
tuple ** hist_eq(const pam &inpam, tuple **array);

int main(int argc, char *argv[])
{
  /* structures for input image */
  pam inpam;

  /* a dynamic two-dimensional array to store the pixels...note that
     pnm uses a tuple (for color images with multiple planes) for
     each pixel.  For PGM files it will only be one plane. */
  tuple **array;
  tuple **outArray;

  /* initializes the library */
  pm_init(argv[0], 0);

  /* read the image */
  array = read_image(argv[1], inpam);
  outArray = hist_eq(inpam, array);

  /* write the output */
  write_image(argv[2], inpam, outArray);

  /* clean up */
  pnm_freepamarray(array, &inpam);
  pnm_freepamarray(outArray, &inpam);

  return 0;
}

// This is the driver function where we do both the conversion 
// to and from RGB and HSI as well as the histogram 
// equalization on intensities.
tuple ** hist_eq(const pam &inpam, tuple **array)
{
	tuple **outArray;
	outArray = pnm_allocpamarray(&inpam);

	vector<vector<vector<double> > > HSI = RGB_HSI(inpam, array);
	fillMappingValues(inpam, array, HSI);
	mapValuesToNewImage(inpam, array, HSI);
	HSI_RGB(inpam, outArray, HSI);

	return outArray;
}

// This applies the formulas from pg 410-411 and returns a 3d vector
// that contains the hue, saturation and intensity values.
vector<vector<vector<double> > > RGB_HSI(const pam &inpam, tuple **array)
{
	double R, G, B = 0.0;
	double H, S, I = 0.0;
	double minRGB = 0.0;
	
	// The vector to return, with the same dimensions as the original
	vector<vector<vector<double> > > HSI(inpam.height,
	 vector<vector<double> >(inpam.width, vector<double>(3)));
 
	for(int row = 0; row < inpam.height; row++)
	{
		for(int col = 0; col < inpam.width; col++)
		{
			R = array[row][col][0];
			G = array[row][col][1];
			B = array[row][col][2];

			// When all RGB components are equal, S = 0, furthermore
			// this prevents division by a very small number, when
			// calculating the theta component
			if ((fabs(R - G) < epsilon) && 
				(fabs(R - B) < epsilon) && 
				(fabs(B - G) < epsilon))
			{	
				H = 0.0;
			}
			else
			{
				if(B <= G)
				{
					// H = theta
					H = acos(
						((0.5)*((R-G)+(R-B)))/
						sqrt((pow((R-G), 2)) + ((R-B)*(G-B))));
				}
				else
				{
					// H = 360 - theta
					H = (360.0*(PI/180)) - acos(
						((0.5)*((R-G)+(R-B)))/
						sqrt((pow((R-G), 2)) + ((R-B)*(G-B))));
				}
			}
			
			// Finish calculating the S and I components
			minRGB = min(min(R,G), B);
			S = 1 - ((3.0 * minRGB)/(R+B+G));
			I = ((1.0/3.0)*(R+G+B));
			
			HSI[row][col][0] = H;
			HSI[row][col][1] = S;
			HSI[row][col][2] = I;
		}
	}
	return HSI;
}

void fillMappingValues(const pam &inpam, tuple **array,
 vector<vector<vector<double> > > & HSI)
{

	// fill r_i
	// -----------------------------------------------------------
	// Just store each intensity value and generate the 
	// corresponding mapping elements for each intensity
	// Note that it is imperative to initialize n and s
	// or else the image will be filled with noise intensities
	for(int k = 0; k < MAXINTENSITY; k++)
	{
		mapping temp;
		temp.r = k;
		temp.n = 0; 
		temp.s = 0;  
		mappingVals.push_back(temp);
	}

	// fill n_k
	// ----------------------------------------------------
	// Check out what the intensity value is at each pixel.
	// Increment the count for that intensity value.
	double oldIntensity = 0.0;
	for(int row = 0; row < inpam.height; row++)
	{
		for(int col = 0; col < inpam.width; col++)
		{
			oldIntensity = HSI[row][col][2];
			mappingVals[static_cast<int>(round(oldIntensity))].n ++;
		}
	}

	// fill s_k
	// -----------------------------------------------------------
	// We now calculate the summation of n for the s at k value.
	// With equation 3.3-8 we notice that we can just store 
	// each sum of n for the intensity k, and then for k+1 we just 
	// look at the previous sum, and add to it the next n 
	// value (the value at that intensity).  After storing
	// these values we just loop through them and multiply
	// each sum by (L-1)/MN to get the final s result for 
	// each intensity.

	double constant = 255.0/(inpam.height * inpam.width);
	mappingVals[0].s = mappingVals[0].n;

	for(int k = 1; k < MAXINTENSITY; k++)
	{
		mappingVals[k].s = static_cast<int>(round((mappingVals[k - 1].s +
		 mappingVals[k].n)));
	}

	for(int k = 0; k < MAXINTENSITY; k++)
	{
		mappingVals[k].s = static_cast<int>(round(mappingVals[k].s * constant));
	}
}

void mapValuesToNewImage(const pam & inpam, tuple ** inArray,
 vector<vector<vector<double> > > & HSI)
{
	// So now that we have all of the mappings from r to s,
	// We loop through the image and get the intensity at that
	// pixel location, then we use that intensity to find the
	// mapping to the output s intensity.  We assign that s
	// intensity value to the pixel. 
	int oldIntensity = 0;
	for(int row = 0; row < inpam.height; row++)
	{
		for(int col = 0; col < inpam.width; col++)
		{
			oldIntensity = HSI[row][col][2];
			HSI[row][col][2] = mappingVals[oldIntensity].s;
		}
	}
}

// This applies the formulas from pg 411-412 and returns 
// a new array with the new RGB values
void HSI_RGB(const pam &inpam, tuple ** outArray,
 vector<vector<vector<double> > > & HSI)
{

	double H, S, I = 0.0;
	double R, G, B = 0.0;

	for (int row = 0; row < inpam.height; row++)
	{
		for(int col = 0; col < inpam.width; col++)
		{
			H = HSI[row][col][0];
			S = HSI[row][col][1];
			I = HSI[row][col][2];
			
			// If S is 0 then the conversion back from HSI means that all 
			// of the components from RGB have the same value.
			if(fabs(S) <= epsilon)
			{
				R = G = B = I;
			}
			else
			{
				// RG sector --> 0 <= H < 120
				if ((0.0 <= H) && (H < (120.0 * (PI/180))))
				{
					R = I * (1.0 + ((S * cos(H))/(cos((60.0 * (PI/180)) - H))));
					B = I * (1.0 - S);
					G = (3.0 * I) - (R + B);
				}
				// GB sector --> 120 <= H < 240
				else if (((120.0 * (PI/180)) <= H) && (H < (240.0 * (PI/180))))
				{
					H = H - (120.0 * (PI/180));
					G = I * (1.0 + ((S* cos(H))/(cos((60.0 * (PI/180)) - H))));
					R = I * (1.0 - S);
					B = (3.0 * I) - (R + G);
				}
				// BR sector --> 240 <= H <= 360
				else 
				{
					H = H - (240.0 * (PI/180));
					B = I * (1.0 + ((S * cos(H))/(cos((60.0 * (PI/180)) - H))));
					G = I * (1.0 - S);
					R = (3.0 * I) - (G + B);
				}
			}

			// Snapping is also in case of very small denominators
			// (i.e. in the case cos(60 - H) in each of the sectors.
			snapRGB(R, G, B);

			outArray[row][col][0] = static_cast<int>(round(R));
			outArray[row][col][1] = static_cast<int>(round(G));
			outArray[row][col][2] = static_cast<int>(round(B));
		}
	}
}

void snapRGB(double & R, double & G, double & B)
{
	if(R < 0)
		R = 0.0;
	else if(R > 255.0)
		R = 255.0;	
	if(G < 0.0)
		G = 0.0;
	else if(G > 255.0)
		G = 255.0;
	if(B < 0.0)
		B = 0.0;
	else if(B > 255.0)
		B = 255.0;
}

/* reads the image into the netpbm structure */
tuple **read_image(char *filename, pam &inpam)
{
  FILE *f;
  tuple **A;

  if ((f = pm_openr(filename)) == NULL) {
    cerr << "Cannot open file \"" << filename << "\" for reading." << endl;
    exit(1);
  }

  if ((A = pnm_readpam(f, &inpam, PAM_STRUCT_SIZE(tuple_type))) == NULL) {
    cerr << "Cannot read image \"" << filename << "\"." << endl;
    exit(1);
  }
  
  pm_close(f);
  return A;
}

/* writes the image to the given file */
void write_image(char *filename, const pam &inpam, tuple **array)
{
  FILE *f;
  pam outpam = inpam;
  
  if ((f = pm_openw(filename)) == NULL) {
    cerr << "Cannot open file \"" << filename << "\" for writing." << endl;
    exit(1);
  }

  /* NOTE: if you change other attributes such as height and width, you
     should change it here too. */
  outpam.file = f;
  
  pnm_writepam(&outpam, array);

  pm_close(f);
}