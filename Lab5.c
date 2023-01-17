/*
	ECE 4310
	Lab 3
	Roderick "Rance" White
	
	This program 
	1) Reads a PPM image, msf image, and ground truth files.
	2) For a range of T, loop through: 
		a) Loop through the ground truth letter locations.
			i)		Check a 9 x 15 pixel area centered at the ground truth location. If
						any pixel in the msf image is greater than the threshold, consider
						the letter “detected”. If none of the pixels in the 9 x 15 area are
						greater than the threshold, consider the letter “not detected”.
			ii) 	If the letter is "not detected" continue to the next letter.
			iii)	Create a 9 x 15 pixel image that is a copy of the area centered at 
						the ground truth location (center of letter) from the original image.
			iv) 	Threshold this image at 128 to create a binary image.
			v) 		Thin the threshold image down to single-pixel wide components.
			vi) 	Check all remaining pixels to determine if they are branchpoints or 
						endpoints.
			vii) 	If there are not exactly 1 branchpoint and 1 endpoint, do not 
						further consider this letter (it becomes "not detected").
		b) Count up the number of FP (letters detected that are not 'e') and TP (number of 
			 letters detected that are 'e').
		c) Output the total FP and TP for each T. 


	The program also demonstrates how to time a piece of code.
	To compile, must link using -lrt  (man clock_gettime() function).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define SQR(x) ((x)*(x))

/*
typedef struct ActiveContour {
	int RowCord;			// Row coordinate of the contour
	int ColCord;			// Column coordinate of the contour
	int IntEnergy1;		// Square distance between itself and previous point
	int IntEnergy2;		// Square distance between itself and the following point
	int ExtEnergy;
} ActCont;
*/


/* A distance calculation function. This is mostly for readability */
float Distance_Function(float x1, float x2, float y1, float y2)
{
	float DistanceValue;
	DistanceValue = sqrt(SQR(x2-x1)+SQR(y2-y1));
	return DistanceValue;
}


/* This function takes in an array of floats and returns an array of unsigned chars */
void Float_to_UnsignedChar(float *imageFl, unsigned char *imageUnChar, int ROWS, int COLS)
{
	int i;
	for(i=0; i<(ROWS*COLS); i++)
		imageUnChar[i] = imageFl[i];
}

/* This function serves to read and open the image, given it's name. 
 * It will return the values within the file, the number of rows, the number of columns, and
 * the number of bytes.
 */
unsigned char *Image_Read(char *FileName, char *header, int *r, int *c, int *b)
{
	FILE						*fpt;
	unsigned char		*image;
	int							ROWS,COLS,BYTES;
	
	/* read image */
	if ((fpt=fopen(FileName,"rb")) == NULL)
	{
		printf("Unable to open %s for reading\n", FileName);
		exit(0);
	}

	/* read image header (simple 8-bit greyscale PPM only) */
	fscanf(fpt,"%s %d %d %d",header,&COLS,&ROWS,&BYTES);	
	if (strcmp(header,"P5") != 0  ||  BYTES != 255)
	{
		printf("Not a greyscale 8-bit PPM image\n");
		exit(0);
	}
	
	/* allocate dynamic memory for image */
	image=(unsigned char *)calloc(ROWS*COLS,sizeof(unsigned char));
	header[0]=fgetc(fpt);	/* read white-space character that separates header */
	fread(image,1,COLS*ROWS,fpt);
	fclose(fpt);
	
	*r = ROWS;		//Return number of rows for the image
	*c = COLS;		//Return number of columns for the image
	*b = BYTES;
	
	return image;
}


/* This function serves to write the image to the appropriate file for more concise code */
void Image_Write(unsigned char *image, char *FileName, int ROWS, int COLS)
{
	FILE						*fpt;
	fpt=fopen(FileName,"w");
	fprintf(fpt,"P5 %d %d 255\n",COLS,ROWS);
	fwrite(image,COLS*ROWS,1,fpt);
	fclose(fpt);
}


/* This function serves to put the 7x7 + shaped contour points onto the image at the 
 * appropriate coordinates */
void Draw_Contour_Points(unsigned char *image, int COLS, int ContRow[], 
	int ContCol[], int ContNum)
{
	int i, j;
	/* Loop through the list of coordinates and replace the values with black + */
	for(i=0; i<ContNum; i++)
	{
		//Replace the pixels around the coordinate to make the + shape		
		for(j=-3; j<=3; j++) image[(ContRow[i]+j)*COLS+ContCol[i]]=0;		//Row line
		for(j=-3; j<=3; j++) image[ContRow[i]*COLS+(ContCol[i]+j)]=0;		//Col line
	}
}




/* Finds the minimum and maximum values within the image */
void Find_Max_and_Min(float *image, int ROWS, int COLS, float *minimumV, float *maximumV)
{
	int i;
	float min, max;
	
	//Set initial values
	min = image[0];
	max = image[0];
	
	//Compare all the values in the image to find the maximum and minimum values.
	for(i=0; i<(ROWS*COLS); i++)
	{
		//If the value is less than the current minimum, it is the new minimum
		if(min>image[i])
			min = image[i];
		
		//If the value is more than the current maximum, it is the new maximum
		if(max<image[i])
			max = image[i];
	}
	
	//Set the values at the end to return them both;
	*minimumV = min;
	*maximumV = max;
}


/* The following takes in an image, its size, current min and max values, and desired min and max
 * values and normalizes the image according to these values.
 */
float *Normalize_Image(float *image, int ROWS, int COLS, float Min, float Max, 
	int NewMin, int NewMax)
{
	int i;
	float *imageNorm;
	
	/* allocate dynamic memory for image */
	imageNorm=(float *)calloc(ROWS*COLS,sizeof(float));
	
	//Loop through all the values and normalize them using it and the new and old min/max
	for(i=0; i<(ROWS*COLS); i++)
		imageNorm[i] = (image[i]-Min)*(NewMax-NewMin)/(Max-Min)+NewMin;
	return imageNorm;	
}


/* Finds the Sobel Convolution of Two Images */
float *Sobel_Convolve(unsigned char *image, int ROWS, int COLS)
{
	int r, c, dr, dc, sumx, sumy;
	int f1[] = {-1,0,1,-2,0,2,-1,0,1}, f2[] = {-1,-2,-1,0,0,0,1,2,1};
	float *SobelConvl;
	
	/* allocate dynamic memory for resulting MSF */
	SobelConvl=(float *)calloc(ROWS*COLS,sizeof(float));
	
	//Loop through every row,column coordinate and convolve for each
	for(r=1; r<ROWS-1; r++)
	{
		for(c=1; c<=COLS-1; c++)
		{
			//Convolve for each value
			sumx=0;
			sumy=0;
			for(dr=-1; dr<=1; dr++)
			{
				for(dc=-1; dc<=1; dc++)
				{
					//Adding the proper convolution for each coordinate
					sumx+= image[(r+dr)*COLS+(c+dc)] * f1[(dr+1)*3+(dc+1)];
					sumy+= image[(r+dr)*COLS+(c+dc)] * f2[(dr+1)*3+(dc+1)];
				}
			}
			SobelConvl[r*COLS+c]=Distance_Function(0, sumx, 0, sumy);		//Value is the distance
		}
	}
	return SobelConvl;
}


/*
Creates the sobel image of the input image
The convolution will have a window of 3x3 due to this being the gradient size
*/
float *Sobel_Edge_Gradient_Make(unsigned char *image, int ROWS, int COLS)
{
	float SobelMin, SobelMax;
	float *SobelImageInit, *SobelImageNorm;
	
	//Allocate space for the initial sobel and the normalized sobel
	SobelImageInit = (float *)calloc(ROWS*COLS, sizeof(float));
	SobelImageNorm = (float *)calloc(ROWS*COLS, sizeof(float));
	
	SobelImageInit = Sobel_Convolve(image, ROWS, COLS);									//Find the initial convolution
	
	//Normalize the Sobel Image
	Find_Max_and_Min(SobelImageInit, ROWS, COLS, &SobelMin, &SobelMax);	//Needed for normalizing
	SobelImageNorm = Normalize_Image(SobelImageInit, ROWS, COLS, SobelMin, SobelMax, 0, 255);
	
	free(SobelImageInit);

	return SobelImageNorm;
}


/* Function to find the first internal energy 
 * This is done by squaring the distance between points
 * Internal energy 1 is the square of the distance between current coord and the next
*/
float Internal_Energy_1(float x1, float x2, float y1, float y2)
{
	float InternalEnergy1;
	InternalEnergy1 = SQR(Distance_Function(x1, x2, y1, y2));
	return InternalEnergy1;
}

/* Function to find the second internal energy 
 * This is done by taking the distance between points and subtracting the average distance
 * Internal energy 2 is the distance minus the average distance overall
*/
float Internal_Energy_2(float InternalEnergy1, float AverageDistance)
{
	float InternalEnergy2;
	InternalEnergy2 = SQR(sqrt(InternalEnergy1) - AverageDistance);
	return InternalEnergy2;
}

/* Function to find the external energy 
 * This is done by finding the square of the inverse value of the sobel input
*/
float External_Energy(float SobInput, float Max)
{
	float ExternalEnergy;
	ExternalEnergy = SQR(Max - SobInput);
	return ExternalEnergy;
}


void Active_Contour(float *SobIm, int ROWS, int COLS, int ContRow[], int ContCol[], int ContNum)
{
	float Min, Max, MinEn, MaxEn;
	float *InEn1, *InEn2, *ExEn;
	
	int i, j, r, c, W, EnIndex, MinEnergySumLocation;
	int Window = 11, IterCur = 0, IterNum = 7;
	float AvgDist, EnergySum, MinEnergySum;
	
	W = Window/2;			//Half the size of the window.
	
	/* Allocate space for the arrays of energies*/
	InEn1 = (float *)calloc(SQR(Window), sizeof(float));
	InEn2 = (float *)calloc(SQR(Window), sizeof(float));
	ExEn = (float *)calloc(SQR(Window), sizeof(float));
	
	/* Find the energies */
	printf("Beginning active contour algorithm.\nWindow size: %d\nRunning %d times\n", 
		Window, IterNum);
	
	Find_Max_and_Min(SobIm, ROWS, COLS, &Min, &Max);		//Need max for the external energy calculations
	
	/* Run for the number of total iterations */
	for(IterCur=0; IterCur<IterNum; IterCur++){
		AvgDist = 0;
		/* Add together the distances of all points */
		//For loop doesn't perform the last distance to loop data
		for(i=0; i<ContNum-1; i++)
			//Distance between current point and the next
			AvgDist += Distance_Function(ContRow[i], ContRow[i+1], ContCol[i], ContCol[i+1]);
			
		//Final calculations for average distance
		AvgDist += Distance_Function(ContRow[i], ContRow[0], ContCol[i], ContCol[0]);
		AvgDist = AvgDist/ContNum;							//Divide by number of points
		
		/* Loop to find the energies for each coordinate */
		for(i=0; i<ContNum; i++)
		{
			EnIndex=0;
			/* Loop through the window around the coordinate for find all the energies in the window */
			for(r=(ContRow[i]-W); r<=(ContRow[i]+W); r++)
			{
				for(c=(ContCol[i]-W); c<=(ContCol[i]+W); c++)
				{
					//Energy calculations
					InEn1[EnIndex] = Internal_Energy_1(r, ContRow[(i+1)%ContNum], c, ContCol[(i+1)%ContNum]);
					InEn2[EnIndex] = Internal_Energy_2(InEn1[EnIndex], AvgDist);
					ExEn[EnIndex] = External_Energy(SobIm[r*COLS+c], Max);
					EnIndex++;
				}
			}
			
			/* Normalize all the energies to have a minimum of 0 and a maximum of 1 */
			//Normalize Internal Energy 1
			Find_Max_and_Min(InEn1, Window, Window, &MinEn, &MaxEn);
			InEn1 = Normalize_Image(InEn1, Window, Window, MinEn, MaxEn, 0, 1);
			
			//Normalize Internal Energy 2
			Find_Max_and_Min(InEn2, Window, Window, &MinEn, &MaxEn);
			InEn2 = Normalize_Image(InEn2, Window, Window, MinEn, MaxEn, 0, 1);
			
			//Normalize External Energy
			Find_Max_and_Min(ExEn, Window, Window, &MinEn, &MaxEn);
			ExEn = Normalize_Image(ExEn, Window, Window, MinEn, MaxEn, 0, 1);
						
			/* Find the Minimum of the Sum of Energies */
			MinEnergySum = InEn1[0] + InEn2[0] + ExEn[0];			//Initial minimum energy sum
			MinEnergySumLocation = 0;
			for(j=1; j<(Window*Window); j++)
			{
				EnergySum = InEn1[j] + InEn2[j] + ExEn[j];
				if(MinEnergySum > EnergySum)
				{
					//Change Energy sum if current sum is less than the current Minimum
					MinEnergySum = EnergySum;
					MinEnergySumLocation = j;
				}
			}
			
			/* Find New Contour Points */
			ContRow[i] = (MinEnergySumLocation/Window) + ContRow[i] - W;
			ContCol[i] = (MinEnergySumLocation%Window) + ContCol[i] - W;
		}
	}
	
}


int main()
{
	unsigned char		*InputImage, *Init_Cont_Image, *SobelImageChar;
	float 					*SobelImage;
	FILE	*fpt;
//	unsigned char 	*NormalizedImage;
	char						header[320];
	int							ROWS,COLS,BYTES;						// The information for the input image
//	int							MSFROWS,MSFCOLS,MSFBYTES;		// The information for the MSF image
	int 						i=0, filesize=0;				// The boundaries of the range
	int 						*ContRow, *ContCol;
	char 						c;
	
	/* Read in the images */
	InputImage = Image_Read("hawk.ppm", header, &ROWS, &COLS, &BYTES);
	Init_Cont_Image = Image_Read("hawk.ppm", header, &ROWS, &COLS, &BYTES);
	
	/* read text file */
	if ((fpt=fopen("hawk_init.txt","r")) == NULL)
	{
		printf("Unable to open %s for reading\n", "hawk_init.txt");
		exit(0);
	}
	
	/* Find the number of lines in the text file of coordinates */
	for (c = getc(fpt); c != EOF; c = getc(fpt))
		if (c == '\n')
			filesize = filesize + 1; // Increment count if this character is newline

	rewind(fpt);		//Rewinds file to read the coordinates
	
	/* Allocate space for the row and column coordinates from the text file */
	ContCol = (int *)calloc(filesize, sizeof(int));
	ContRow = (int *)calloc(filesize, sizeof(int));
	
	/* Create arrays for the contour inputs */
	while (fscanf(fpt, "%d %d", &ContCol[i], &ContRow[i]) == 2) i++;		//Read in the coordinates
	
	fclose(fpt);							//Close text file after using
	
	/* Make image with initial contour coordinates */
	Draw_Contour_Points(Init_Cont_Image, COLS, ContRow, ContCol, filesize);
	Image_Write(Init_Cont_Image, "hawk_init.ppm", ROWS, COLS);					//Image with initial coords
	
	/* Create the sobel image */
	//SobelImage = (float *)calloc(ROWS*COLS, sizeof(float));
	//SobelImageInverted = (float *)calloc(ROWS*COLS, sizeof(float));
	SobelImageChar = (unsigned char *)calloc(ROWS*COLS, sizeof(unsigned char));
	
	SobelImage = Sobel_Edge_Gradient_Make(InputImage, ROWS, COLS);
	
	//Write the Sobel Image to a file
	Float_to_UnsignedChar(SobelImage, SobelImageChar, ROWS, COLS);	//Unsigned char image for writing
	Image_Write(SobelImageChar, "hawk_sobel.ppm", ROWS, COLS);			//Write image to file
	
	/* Active Contour */
	Active_Contour(SobelImage, ROWS, COLS, ContRow, ContCol, filesize);
	
	/* Make image with initial contour coordinates */
	Draw_Contour_Points(InputImage, COLS, ContRow, ContCol, filesize);
	Image_Write(InputImage, "hawk_final.ppm", ROWS, COLS);			//Image with initial coords
	
	
	/* Write new contour coordinates to text file */
	fpt=fopen("hawk_final.txt","w");
	for(i=0; i<filesize; i++)
		fprintf(fpt,"%d %d\n",ContCol[i],ContRow[i]);
	fclose(fpt);
}

	
	












