/*
	Adam Lefaivre - 001145679
	Assignment 3 - CPSC 5990 (D.I.P.)
	Motion Blurring
	Dr. Howard Cheng
*/

#include <pam.h>
#include <vector>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include "fft.h"

using namespace std;

const double PI = 3.14159265359;

// Function stubs
void write_image(char *filename, const pam &inpam, tuple **array);
tuple **read_image(char *filename, pam &inpam);

// functions for the fourier transform pipeline (in order)
void applyFTRecipe(const pam &originalPam, pam & paddedPam, 
	tuple ** originalArray, double a, double b, double T);
tuple ** padImage(const pam &inpam, pam & outpam, tuple **array);
vector<Array> shift(const pam &paddedPam, tuple** inArr);
vector<Array> FT(vector<Array> inArr, const pam & paddedPam);
void applyMotionBlur(vector<Array> & ftResult, const pam & paddedPam,
 double a, double b, double T);
void IFT(tuple** inArr, const pam &inpam, const pam & paddedPam, vector<Array> & F_by_H);
void shift(const pam &paddedPam, vector<Array> freqArray);
void unpadImageAndGetResult(const pam &originalPam, tuple ** paddedArray, 
	tuple ** originalArray, const vector<Array> & ftResult);

// A function to return the nearest power of 2 (we need this for the FFT)
int nearestPowerOfTwo(int input);

int main(int argc, char *argv[])
{
	// The height is given by the number of rows in the
	// image and the width is given by the number of columns.

	// pam: structures for input/output images
	// tuple: a dynamic two-dimensional array to store the pixels...note that
	// pnm uses a tuple (for color images with multiple planes) for
	// each pixel.  For PGM files it will only be one plane.

  	// structures for input image
 	 pam originalPam;

  	// a dynamic two-dimensional array to store the pixels...note that
  	// pnm uses a tuple (for color images with multiple planes) for
    // each pixel.  For PGM files it will only be one plane. 

  	tuple **originalArray;

 	// initializes the library 
  	pm_init(argv[0], 0);

  	// read the image 
  	originalArray = read_image(argv[1], originalPam);

  	// Get arguments for motion blur equation.
  	double a, b, T = 0.0;
  	a = atof(argv[3]);
  	b = atof(argv[4]);
  	T = atof(argv[5]);

  	// Create a new pam data structure to be used throughout the program
  	// for padding throughout.
  	pam paddedPam = originalPam;

  	// This applies the steps needed for frequency domain filtering.
	applyFTRecipe(originalPam, paddedPam, originalArray, a, b, T);
	
	// write the output
	write_image(argv[2], originalPam, originalArray);
	
	// clean up
	pnm_freepamarray(originalArray, &originalPam);
	
	return 0;
}

// We can now swap out the motion blur function for any other filter
// that we want to apply in the frequency domain! =)
void applyFTRecipe(const pam &originalPam, pam & paddedPam, 
	tuple ** originalArray, double a, double b, double T)
{
	tuple ** paddedArray;
	paddedArray = padImage(originalPam, paddedPam, originalArray);
	vector<Array> shiftedVector = shift(paddedPam, paddedArray);
	vector<Array> ftResult = FT(shiftedVector, paddedPam);
	applyMotionBlur(ftResult, paddedPam, a, b, T);
	IFT(paddedArray, originalPam, paddedPam, ftResult);
	shift(paddedPam, ftResult);
	unpadImageAndGetResult(originalPam, paddedArray, originalArray, ftResult);

	for(int row = 0; row < paddedPam.height; row++)
	{
		delete [] ftResult[row];
	}
}

// Get the real portion of the IFT and snap the intensities to between
// 255 and 0.  Note that we are only looping up to the original arrays
// width and height and not the padded width and height.
void unpadImageAndGetResult(const pam &originalPam, tuple ** paddedArray,
 tuple ** originalArray, const vector<Array> & ftResult)
{
	int resultIntensity = 0;
	for(int row = 0; row < originalPam.height; row++)
	{
		for(int col = 0; col < originalPam.width; col++)	
		{
			// We take the real portion here and omit the imaginary portion
			resultIntensity = 
			static_cast<int>(round(ftResult[row][col].real()));

			if(resultIntensity > 255)
			{
				resultIntensity = 255;
			}
			else if (resultIntensity < 0)
			{
				resultIntensity = 0;	
			}

			originalArray[row][col][0] = resultIntensity;
		}
	}
}

// This will apply the motion blur equation, if pi*(ua+vb) == 0
// then we just have sin(0)/0 = 1.  The N/2 and M/2 are 
// offsets, so that the filter works with respect to the 
// centered FT.
void applyMotionBlur(vector<Array> & ftResult, const pam & paddedPam, 
	double a, double b, double T)
{	
	double H_Real = 0.0;
	int shifted_u = 0.0;
	int shifted_v = 0.0;

	for(int u = 0; u < paddedPam.height; u++)
	{
		for(int v = 0; v < paddedPam.width; v++)
		{
			shifted_u = (u-(paddedPam.height/2));
			shifted_v = (v-(paddedPam.width/2));

			Coeff H_Imaginary = 
			exp(Coeff(0, -1*PI*((shifted_u*a) + (shifted_v*b))));

			if(static_cast<int>(
				round(PI * ((shifted_u*a) + (shifted_v*b)))) == 0)
			{
				H_Real = T;
			}
			else
			{
				H_Real = (T/(PI * ((shifted_u*a) + (shifted_v*b)))) * 
				(sin(PI*((shifted_u*a) + (shifted_v*b))));
			}

			Coeff H_Product(H_Real, 0);
			H_Product *= H_Imaginary;	
			ftResult[u][v] *= H_Product;		
		}
	}
}

// This shift function applies -1^x+y to center the transform
vector<Array> shift(const pam &paddedPam, tuple** inArr)
{
	vector<Array> toReturn;
	for(int row = 0; row < paddedPam.height; row++)
	{
		Array currRow = new Coeff [paddedPam.width];
		for(int col = 0; col < paddedPam.width; col++)
		{
			Coeff valAfterShift(inArr[row][col][0], 0);
			valAfterShift *= pow(-1, row+col);
			currRow[col] = valAfterShift;
		}
		toReturn.push_back(currRow);
	}
	return toReturn;
}

// This shift function applies -1^x+y to obtain the result 
// after applying the filter
void shift(const pam &paddedPam, vector<Array> freqArray)
{
	for(int row = 0; row < paddedPam.height; row++)
	{
		for(int col = 0; col < paddedPam.width; col++)
		{
			Coeff shifter(pow(-1, row+col), 0);
			freqArray[row][col] *= shifter;
		}
	}
}

// We use the separability property of the 2-D DFT here.
// Loop through rows and create a new row with
// the right number of elements (i.e. populating the
// row with the corresponding column values).
// We then call the FFT function and store the result
// in an output vector.  We then repeat this process
// for columns and then copy the FFT results of that column
// into the resultant output vector.
vector<Array> FT(vector<Array> inArr, const pam & paddedPam) 
{
	vector<Array> ftResult;
	Coeff rowBuffer[paddedPam.width];  
	for(int row = 0; row < paddedPam.height; row++)
	{
		// Populate new row before sending it to FFT
		Array currRow = new Coeff [paddedPam.width];
		for(int col = 0; col < paddedPam.width; col++)
		{
			currRow[col] = inArr[row][col];
		}

		FFT(currRow, paddedPam.width, rowBuffer);
		ftResult.push_back(currRow);
	}
	
	Coeff colBuffer[paddedPam.height];  
	for(int col = 0; col < paddedPam.width; col++)
	{
		// Populate column before sending it to FFT
		Coeff currCol[paddedPam.height];
		for(int row = 0; row < paddedPam.height; row++)
		{
			currCol[row] = ftResult[row][col];
		}

		FFT(currCol, paddedPam.height, colBuffer);
		for(int row = 0; row < paddedPam.height; row++)
		{
			ftResult[row][col] = currCol[row];
		}
	}
	return ftResult;
}

// The IFT has much of the same logic as for the FT.
void IFT(tuple** inArr, const pam &inpam, 
	const pam & paddedPam, vector<Array> & F_by_H) 
{
	Coeff rowBuffer[paddedPam.width];
	for(int row = 0; row < paddedPam.height; row++)
	{
		// Send rows to IFFT
		inverseFFT(F_by_H[row], paddedPam.width, rowBuffer);
	}
	
	Coeff colBuffer[paddedPam.height];  
	for(int col = 0; col < paddedPam.width; col++)
	{
		// Populate column before sending it to FFT
		Coeff currCol[paddedPam.height];
		for(int row = 0; row < paddedPam.height; row++)
		{
			currCol[row] = F_by_H[row][col];
		}
		
		inverseFFT(currCol, paddedPam.height, colBuffer);
		for(int row = 0; row < paddedPam.height; row++)
		{
			F_by_H[row][col] = currCol[row];
		}
	}
}

// A helper function for the FFT algorithm
int nearestPowerOfTwo(int input)
{
	return pow(2, ceil(log(input)/log(2)));
}

// We use the above function to pad excess zeros up to the nearest power of two 
// that is also at least twice the size of the image.  This way both the FFT 
// algorithm and the proper image padding in the frequency domain
// will be satisfied.
tuple ** padImage(const pam & originalPam, pam & paddedPam, tuple ** oldArray)
{
	int paddedWidth, paddedHeight = 0;

	paddedWidth = nearestPowerOfTwo(originalPam.width * 2);
	paddedHeight = nearestPowerOfTwo(originalPam.height * 2);

	paddedPam.width = paddedWidth;
	paddedPam.height = paddedHeight;

	tuple ** paddedArray;
	paddedArray = pnm_allocpamarray(&paddedPam);

	for (int row = 0; row < originalPam.height; row++)
	{
		for (int col = 0; col < originalPam.width; col++)
		{
			paddedArray[row][col][0] = oldArray[row][col][0];
		}
	}

	for (int row = originalPam.height; row < paddedHeight; row++)
	{
		for (int col = originalPam.width; col < paddedWidth; col++)
		{
			paddedArray[row][col][0] = 0;
		}
	}

	return paddedArray;
}

// The usual...
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

void write_image(char *filename, const pam &inpam, tuple **array)
{
	FILE *f;
	pam outpam = inpam;

	if ((f = pm_openw(filename)) == NULL) {
		cerr << "Cannot open file \"" << filename << "\" for writing." << endl;
		exit(1);
	}

	// NOTE: if you change other attributes such as height and width, you
	// should change it here too.
	outpam.file = f;

	pnm_writepam(&outpam, array);

	pm_close(f);
}