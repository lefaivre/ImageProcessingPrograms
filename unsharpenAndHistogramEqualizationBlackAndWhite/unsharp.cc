/*
	Adam Lefaivre - 001145679
	Assignment 2 - CPSC 5990 (D.I.P.)
	Programming Question 2 
	Unsharp Masking
	Dr. Howard Cheng
*/

#include <pam.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

// Function stubs
tuple **read_image(char *filename, pam &inpam);
void write_image(char *filename, const pam &inpam, tuple **array);
tuple ** unsharpMasking(const pam & inpam, tuple ** originalArray); 
double getAveragedIntensity(int col, int row, const pam & inpam, tuple ** originalArray);

int main(int argc, char *argv[])
{
	// The height is given by the number of rows in the
	// image and the width is given by the number of columns.

	// pam: structures for input/output images
	// tuple: a dynamic two-dimensional array to store the pixels...note that
	// pnm uses a tuple (for color images with multiple planes) for
	// each pixel.  For PGM files it will only be one plane.

  	// structures for input image
 	 pam inpam;

  	// a dynamic two-dimensional array to store the pixels...note that
  	// pnm uses a tuple (for color images with multiple planes) for
    // each pixel.  For PGM files it will only be one plane. 

  	tuple **oldArray;
  	tuple **newArray;

 	// initializes the library 
  	pm_init(argv[0], 0);

  	// read the image 
  	oldArray = read_image(argv[1], inpam);

  	newArray = unsharpMasking(inpam, oldArray);

	// write the output
	write_image(argv[2], inpam, newArray);

	// clean up
	pnm_freepamarray(oldArray, &inpam);
	pnm_freepamarray(newArray, &inpam);

	return 0;
}

tuple ** unsharpMasking(const pam & inpam, tuple ** originalArray)
{
	tuple ** outArray;
	outArray = pnm_allocpamarray(&inpam);

	double averagedIntensity = 0.0;
	int result = 0;
	
	for(int col = 0; col < inpam.width; col++)
	{
		for(int row = 0; row < inpam.height; row++)
		{

			averagedIntensity = getAveragedIntensity(col, row, inpam, originalArray);

			// Note: Since we are doing unsharp masking and not highboost filtering
			// we can just double the original intensity and subtract the average

			// Had there been a different constant multiplying the mask, then it 
			// should be accounted for, but we don't need to worry about it.

			result = static_cast<int>(round((2 * originalArray[row][col][0]) - averagedIntensity));

			// Now since it is possible that we get intensities that are larger than
			// the highest intensity and lower than 0, we just set these values, which 
			// are out of bounds, to 255 or 0, depending on the value.
			// If we don't do this, we get some grainy looking distortion.
			
			if (result < 0)
			{
				result = 0;
			}
			else if (result > 255)
			{
				result = 255;
			}

			// Here it is imperative to assign the intensity to a new image array
			// and not the original one!

			outArray[row][col][0] = result;
		}
	}

	return outArray;
}

double getAveragedIntensity(int col, int row, const pam & inpam, tuple ** originalArray)
{
	double averagedIntensity = 0.0;

	// We start by checking all of the boundary conditions!
	// Checks 1 through 4 are corner points 
	// (WE NEED TO DO THESE FIRST FOR THE LOGIC TO WORK)
	// Checks 5 through 8 are edges.

	if ((col == 0) && (row == 0)) //origin
	{
		averagedIntensity = originalArray[row][col][0] + 
		originalArray[row + 1][col][0] + 
		originalArray[row][col + 1][0] + 
		originalArray[row + 1][col + 1][0];

		averagedIntensity = averagedIntensity/4.0;
	}
	else if((col == inpam.width - 1) && (row == 0)) //top right
	{
		averagedIntensity = originalArray[row][col][0] + 
		originalArray[row][col - 1][0] + 
		originalArray[row + 1][col][0] + 
		originalArray[row + 1][col - 1][0];

		averagedIntensity = averagedIntensity/4.0;
	}
	else if((col == inpam.width - 1) && (row == inpam.height - 1)) //bottom right 
	{
		averagedIntensity = originalArray[row][col][0] + 
		originalArray[row][col - 1][0] + 
		originalArray[row - 1][col][0] + 
		originalArray[row - 1][col - 1][0];

		averagedIntensity = averagedIntensity/4.0;
	}
	else if((col == 0) && (row == inpam.height - 1)) //bottom left 
	{
		averagedIntensity = originalArray[row][col][0] + 
		originalArray[row][col + 1][0] + 
		originalArray[row - 1][col][0] + 
		originalArray[row - 1][col + 1][0];

		averagedIntensity = averagedIntensity/4.0;
	}
	else if(row == 0) //TOP EDGE
	{
		averagedIntensity = originalArray[row][col][0] +
		originalArray[row][col -1][0] + 
		originalArray[row][col + 1][0] + 
		originalArray[row + 1][col][0] + 
		originalArray[row + 1][col - 1][0] + 
		originalArray[row + 1][col + 1][0];

		averagedIntensity = averagedIntensity/6.0;
	}
	else if(col == inpam.width - 1) //RIGHT EDGE
	{
		averagedIntensity = originalArray[row][col][0] + 
		originalArray[row - 1][col][0] + 
		originalArray[row + 1][col][0] + 
		originalArray[row][col - 1][0] + 
		originalArray[row + 1][col - 1][0] +
		originalArray[row - 1][col - 1][0];

		averagedIntensity = averagedIntensity/6.0;
	}
	else if(row == inpam.height - 1) //BOTTOM EDGE
	{
		averagedIntensity = originalArray[row][col][0] + 
		originalArray[row][col - 1][0] + 
		originalArray[row][col + 1][0] + 
		originalArray[row - 1][col - 1][0] + 
		originalArray[row - 1][col][0] +
		originalArray[row -1][col + 1][0];

		averagedIntensity = averagedIntensity/6.0;
	}
	else if(col == 0) //LEFT EDGE
	{
		averagedIntensity = originalArray[row][col][0] + 
		originalArray[row - 1][col][0] + 
		originalArray[row + 1][col][0] + 
		originalArray[row][col + 1][0] + 
		originalArray[row - 1][col + 1][0] +
		originalArray[row + 1][col + 1][0];

		averagedIntensity = averagedIntensity/6.0;
	}
	else // NOW WE DO THE GENERAL CASE!
	{
		averagedIntensity = originalArray[row][col][0] + 
		originalArray[row - 1][col][0] + 
		originalArray[row + 1][col][0] + 
		originalArray[row][col + 1][0] + 
		originalArray[row][col - 1][0] +
		originalArray[row - 1][ col - 1][0] + 
		originalArray[row - 1][col + 1][0] + 
		originalArray[row + 1][col + 1][0] + 
		originalArray[row + 1][col - 1][0];

		averagedIntensity = averagedIntensity/9.0;
	}

	return averagedIntensity;
}


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