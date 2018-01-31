/*
	Adam Lefaivre - 001145679
	Assignment 1 - CPSC 5990 (D.I.P.)
	Dr. Howard Cheng
	Resources: 
		http://supercomputingblog.com/graphics/coding-bilinear-interpolation/
		Accessed May, 2017
*/

#include <pam.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <limits>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

const double PI  = 3.141592653;
const double DEFAULTINTENSITY = 100.0; // Set this to 0.0 for black 

using namespace std;

// Function stubs
tuple **read_image(char *filename, pam &inpam);
void write_image(char *filename, const pam &inpam, tuple **array);
double CW_RowEq(double v, double w, double theta);
double CW_ColEq(double v, double w, double theta);
double CCW_RowEq(double v, double w, double theta);
double CCW_ColEq(double v, double w, double theta);
std::vector<double> findFrameSizeForNewImage(const pam &inpam, double theta);
std::vector<double> getHeightPoints(const pam &inpam, double theta);
std::vector<double> getWidthPoints(const pam &inpam, double theta);
void inverseMapping(const pam & inpam, pam &outpam,
 tuple **oldArray, tuple ** rotatedArray, double theta);
void bilerp(const pam &inpam, tuple **oldArray, 
	tuple ** rotatedArray, double oldRow, double oldCol, int row, int col);

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

 	// initializes the library 
  	pm_init(argv[0], 0);

  	// read the image 
  	oldArray = read_image(argv[1], inpam);

	// Convert to radians
	double theta = atof(argv[3]) * (PI/180.0);

	// Assign the total row size to height, and the total column size to width.
	std::vector<double> widthAndHeight = 
	findFrameSizeForNewImage(inpam, theta);

	// Allocate a new pam and array object for our output!
	pam outpam = inpam;
	tuple ** rotatedArray;
	outpam.width = widthAndHeight[0];
	outpam.height = widthAndHeight[1];
	rotatedArray = pnm_allocpamarray(&outpam);

	// Apply inverse mapping using bilinear interpolation.
	inverseMapping(inpam, outpam, oldArray, rotatedArray, theta);

	// write the output
	write_image(argv[2], outpam, rotatedArray);

	// clean up
	pnm_freepamarray(oldArray, &inpam);
	pnm_freepamarray(rotatedArray, &outpam);

	return 0;
}

// After say a 45 degree rotation of a rectangular image,
// there will need to be a new height and width that allows
// for the image to fit in the "frame".  So we need to calculate
// the difference between the ceil(max) - floor(min) rows and columns.
// Floor and ceiling are used because the rotation equation gives us 
// double values and we should allow for this edge condition row/col 
// to be accessed if need be. 

vector<double> findFrameSizeForNewImage(const pam &inpam, double theta)
{
	vector<double> rowVector = getHeightPoints(inpam, theta);
	vector<double> colVector = getWidthPoints(inpam, theta);
	vector<double> outVector;

	double maxRow = *(max_element(rowVector.begin(), rowVector.end()));
	double minRow = *(min_element(rowVector.begin(), rowVector.end()));
	double maxCol = *(max_element(colVector.begin(), colVector.end()));
	double minCol = *(min_element(colVector.begin(), colVector.end()));

	outVector.push_back(ceil(maxCol) - floor(minCol));
	outVector.push_back(ceil(maxRow) - floor(minRow));

	return outVector;
}

// Utility functions to help find the new size for the rotated image
// and for further shifting of points in other functions.

vector<double> getHeightPoints(const pam &inpam, double theta)
{
	std::vector<double> outVect;
	outVect.push_back(CCW_RowEq(0, 0, theta));
	outVect.push_back(CCW_RowEq(0, inpam.width - 1, theta));
	outVect.push_back(CCW_RowEq(inpam.height - 1, 0, theta));
	outVect.push_back(CCW_RowEq(inpam.height - 1, inpam.width - 1, theta));

	return outVect;
}

vector<double> getWidthPoints(const pam &inpam, double theta)
{
	std::vector<double> outVect;
	outVect.push_back(CCW_ColEq(0, 0, theta));
	outVect.push_back(CCW_ColEq(0, inpam.width - 1, theta));
	outVect.push_back(CCW_ColEq(inpam.height - 1, 0, theta));
	outVect.push_back(CCW_ColEq(inpam.height - 1, inpam.width - 1, theta));

	return outVect;
}

// The function to look back at the original image using CW rotation
// equations and execute interpolation

void inverseMapping(const pam & inpam, pam &outpam, 
	tuple **oldArray, tuple ** rotatedArray, double theta)
{
	// Find the old images min width and height points for proper iterpolation.

	vector<double> widthPoints = getWidthPoints(inpam, theta);
	vector<double> heightPoints = getHeightPoints(inpam, theta);
	
	double shiftWidth = floor(*(std::min_element(widthPoints.begin(),
	 widthPoints.end())));

	double shiftHeight = floor(*(std::min_element(heightPoints.begin(),
	 heightPoints.end())));
	
	double oldRow, oldCol, actualRow, actualCol = 0.0;

	// For each pixel in the new frame, look back at the original image.
	// and do bilinear interpolation 

	for(int row = 0; row < outpam.height; row++)
	{
		for(int col = 0; col < outpam.width; col++)
		{
			// We are rotating about the origin, so we need to
			// add the min row and col values to match the location
			// of the intensity values in the new image.  Looping 
			// through (0,0) to (max height, max width) of the new image
			// is not correct given the coordinate system.

			actualRow = shiftHeight + row;
			actualCol = shiftWidth + col;

			// Go backwards to the original image and look at the old row/column
			// value and then perform the interpolation on that.

			oldRow = CW_RowEq(actualRow, actualCol, theta);
			oldCol = CW_ColEq(actualRow, actualCol, theta);

			// Now call the actual bilinear interpolation function
			// Assigning a new intensity value based on the four nearest
			// pixel coordinates in the original image.

			bilerp(inpam, oldArray, rotatedArray, oldRow, oldCol, row, col);
		}
	}
}

// The function to execute the bilinear interpolation
void bilerp(const pam &inpam, tuple **oldArray, tuple ** rotatedArray,
 double oldRow, double oldCol, int row, int col)
{
	// The Q values represent the 4 pixels nearest the original pixel.
	// The result values are averages of the Q values.
	// P is the NEW intensity value that is returned.

	double Q11, Q12, Q21, Q22, P, result1, result2 = 0.0;
	int R1, R2, C1, C2 = 0;

	// Ex:
	// oldRow,oldCol = 180.6,60.4
	// R1,C1 = 180,60
	// R2,C2 = 181,61
	// We need to cast to int since
	// these are actually pixel coordinates
	
	R1 = static_cast<int>(floor(oldRow));
	R2 = static_cast<int>(ceil(oldRow));
	C1 = static_cast<int>(floor(oldCol));
	C2 = static_cast<int>(ceil(oldCol));

	// If the R,C (row, column) values are outside of the old image after,
	// applying the rotation equations, then we set them to a default
	// intensity, otherwise we get the pixel's neighbouring 4 values.
	// 
	// Imagine we rotate 45 degrees on a rectangular image
	// and do the min row,col shifting in the inverse 
	// mapping function, then we'll have an image where:
	//		new image width > old image width
	//		new image height < old image height
	// So some value in the new image will need to be set to default
	// after rotation.

	if (((R1 < 0) || (!(R1 < inpam.height))) || 
		((C1 < 0) || (!(C1 < inpam.width))))
		Q11 = DEFAULTINTENSITY;
	else
		Q11 = *(oldArray[R1][C1]);
	
	if (((R1 < 0) || (!(R1 < inpam.height))) ||
		((C2 < 0) || (!(C2 < inpam.width))))
		Q12 = DEFAULTINTENSITY;
	else
		Q12 = *(oldArray[R1][C2]);

	if (((R2 < 0) || (!(R2 < inpam.height))) ||
		((C1 < 0) || (!(C1 < inpam.width))))
		Q21 = DEFAULTINTENSITY; 
	else
		Q21 = *(oldArray[R2][C1]);

	if (((R2 < 0) || (!(R2 < inpam.height))) ||
		((C2 < 0) || (!(C2 < inpam.width))))
		Q22 = DEFAULTINTENSITY;
	else
		Q22 = *(oldArray[R2][C2]);

	// Now we will do the actual bilinear interpolation part!
	// We should check for division by 0 too!

	if ((R2 - R1) == 0)
	{
		result1 = 0.0;
		result2 = 0.0;
	}
	else
	{	
		// Calculate weightings
		result1 = ((R2 - oldRow)/(R2 - R1))*Q11 
			+ ((oldRow - R1)/(R2 - R1))*Q21;
		result2 = ((R2 - oldRow)/(R2 - R1))*Q12 
			+ ((oldRow - R1)/(R2 - R1))*Q22;
	}

	if ((C2 - C1) == 0)
	{
		P = 0.0;
	}
	else
	{
		P = ((C2 - oldCol)/(C2 - C1))*result1
		+ ((oldCol - C1)/(C2 - C1))*result2;
	}

	// Assign the new array at r,c the interpolated intensity!
	*(rotatedArray[row][col]) = P;
}

// Rotation Equations
// We are rotating about the origin here.  
// Note: (0,0) stays as (0,0).

double CCW_RowEq(double v, double w, double theta)
{
	return (v * cos(theta)) - (w * sin(theta));
}

double CCW_ColEq(double v, double w, double theta)
{
	return (v * sin(theta)) + (w * cos(theta));
}

double CW_RowEq(double v, double w, double theta)
{
	return (v * cos(theta)) + (w * sin(theta));
}

double CW_ColEq(double v, double w, double theta)
{
	return (v * (-(sin(theta)))) + (w * cos(theta));
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