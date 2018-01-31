/*
	Adam Lefaivre - 001145679
	Assignment 2 - CPSC 5990 (D.I.P.)
	Programming Question 1 
	Histogram Equalization
	Dr. Howard Cheng
*/

#include <vector>
#include <pam.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <math.h>

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

const int MAXINTENSITY = 256;
std::vector<mapping> mappingVals;

tuple **read_image(char *filename, pam &inpam);
void write_image(char *filename, const pam &inpam, tuple **array);
void fillMappingVals(tuple **array, const pam &inpam);
tuple ** mapValuesToNewImage(tuple **array, const pam &inpam);

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

  //-------------------------------------------------
  // Here is where the histogram equalization occurs.
  //-------------------------------------------------
  fillMappingVals(array, inpam);
  outArray = mapValuesToNewImage(array, inpam);

  /* write the output */
  write_image(argv[2], inpam, outArray);

  /* clean up */
  pnm_freepamarray(array, &inpam);
  pnm_freepamarray(outArray, &inpam);

  return 0;
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

void fillMappingVals(tuple **array, const pam &inpam)
{

  // fill r_i
  // -----------------------------------------------------------
  // Just store each intensity value and generate the 
  // corresponding mapping elements for each intensity
  for(int k = 0; k < MAXINTENSITY; k++)
  {
  	mapping temp;
  	temp.r = k;
  	mappingVals.push_back(temp);
  }

  // fill n_k
  // -----------------------------------------------------------
  // Loop through the image and check out what the intensity
  // value is at each pixel.  Increment the count for that
  // intensity value.
  for(int row = 0; row < inpam.height; row++)
  {
  	for(int col = 0; col < inpam.width; col++)
  	{
  		//get the n_k values
  		mappingVals[*(array[row][col])].n ++;
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
  	mappingVals[k].s = static_cast<int>(round((mappingVals[k - 1].s + mappingVals[k].n)));
  }

  for(int k = 0; k < MAXINTENSITY; k++)
  {
  	mappingVals[k].s = static_cast<int>(round(mappingVals[k].s * constant));
  }
}

tuple ** mapValuesToNewImage(tuple ** inArray, const pam & inpam)
{
	// So now that we have all of the mappings from r to s,
	// We loop through the image and get the intensity at that
	// pixel location, then we use that intensity to find the
	// mapping to the output s intensity.  We assign that s
	// intensity value to the pixel. 

	tuple ** outArray = pnm_allocpamarray(&inpam);
	int oldIntensity = 0;
	for(int row = 0; row < inpam.height; row++)
	{
		for(int col = 0; col < inpam.width; col++)
		{
			oldIntensity = *(inArray[row][col]);
			*(outArray[row][col]) = mappingVals[oldIntensity].s;
		}
	}

	return outArray;
}

