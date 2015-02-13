/*
 *  CrisprNorm.h
 *  Transform Crispr measure to MA plot and perform normalization
 *
 *  Created by Han Xu on 3/6/14.
 *  Copyright 2014 Dana Farber Cancer Institute. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#define NDEBUG
#include <assert.h>
#include <math.h>
#include "math_api.h"
#include "words.h"
#include "rvgs.h"
#include "rngs.h"

#define MAX_NAME_LEN 255           //maximum length of item name, group name or list name
#define MAX_WORD_IN_LINE 255	   //maximum number of words in a line

typedef struct
{
	char sgName[MAX_NAME_LEN];       //name of the sgRNA
	char geneName[MAX_NAME_LEN];	 //name of the gene
	double x1;                       //value of first measure
	double x2;						 //value of second measure
	double m;                        //log-mean
	double r;                        //log-ratio
	double adjustedR;                //adjusted log-ratio
} ITEM_STRUCT;


//Read input file. File Format: <sgRNA id> <gene id> <measure in library 1> <measure in library 2>. Return the number of items in the file
int ReadFile(char *fileName, ITEM_STRUCT **pItems);

//transform to log mean-ratio. m = x1'+x2', r = x2'-x1', x' = log2(x/median+pseudo-count), pseudo-count = <median of all values in a library>*0.01
int ComputeMR(ITEM_STRUCT *items, int itemNum);

//QuickSort items by log-mean
void QuickSortItemByM(ITEM_STRUCT *items, int lo, int hi);

//Adjust r using z-transform within a window sliding on items sorted by m
int AdjustMR(ITEM_STRUCT *items, int itemNum, int winSize);

//print the usage of Command
void PrintCommandUsage(const char *command);

//Read input file. File Format: <sgRNA id> <gene id> <measure in library 1> <measure in library 2>. Return the number of items in the file
int ReadFile(char *fileName, ITEM_STRUCT **pItems)
{
	FILE *fh;
	char **words, tmpS[MAX_WORD_IN_LINE*MAX_NAME_LEN];
	int wordNum;
	int totalItemNum;
	
	words = AllocWords(MAX_WORD_IN_LINE, MAX_NAME_LEN+1);
	
	assert(words!=NULL);
	
	if (words == NULL)
	{
		return -1;
	}
	
	fh = (FILE *)fopen(fileName, "r");
	
	if (!fh)
	{
		printf("Cannot open file %s\n", fileName);
		return -1;
	}
	
	//Read the header row to get the sample number
	fgets(tmpS, MAX_WORD_IN_LINE*MAX_NAME_LEN, fh);
	
	wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, MAX_WORD_IN_LINE, " \t\r\n\v\f");
	
	assert(wordNum == 4);
	
	if (wordNum != 4)
	{
		printf("Input file format: <sgRNA id> <gene id> <measure in library 1> <measure in library 2>.\n");
		return -1;
	}
	
	//read records of items
	
	totalItemNum = 0;
	
	fgets(tmpS, 255*(MAX_NAME_LEN+1)*sizeof(char), fh);
	wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, MAX_WORD_IN_LINE, " \t\r\n\v\f");
	
	while ((wordNum==4)&&(!feof(fh)))
	{		
		totalItemNum++;
		
		fgets(tmpS, 255*(MAX_NAME_LEN+1)*sizeof(char), fh);
		wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, MAX_WORD_IN_LINE, " \t\r\n\v\f");
	}
	
	fclose(fh);	
	
	*pItems = (ITEM_STRUCT *)malloc(totalItemNum*sizeof(ITEM_STRUCT));
	
	assert(pItems!=NULL);
	
	fh = (FILE *)fopen(fileName, "r");
	
	if (!fh)
	{
		printf("Cannot open file %s\n", fileName);
		return -1;
	}
	
	//Read the header row to get the sample number
	fgets(tmpS, MAX_WORD_IN_LINE*(MAX_NAME_LEN+1)*sizeof(char), fh);
	
	wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, MAX_WORD_IN_LINE, " \t\r\n\v\f");
	
	assert(wordNum == 4);
	
	//read records of items
	
	totalItemNum = 0;
	
	fgets(tmpS, MAX_WORD_IN_LINE*(MAX_NAME_LEN+1)*sizeof(char), fh);
	wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, MAX_WORD_IN_LINE, " \t\r\n\v\f");
	
	while ((wordNum==4)&&(!feof(fh)))
	{
		strcpy((*pItems)[totalItemNum].sgName, words[0]);
		strcpy((*pItems)[totalItemNum].geneName, words[1]);
		(*pItems)[totalItemNum].x1 = atof(words[2]);
		(*pItems)[totalItemNum].x2 = atof(words[3]);
		totalItemNum++;
		
		fgets(tmpS, 255*(MAX_NAME_LEN+1)*sizeof(char), fh);
		wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, MAX_WORD_IN_LINE, " \t\r\n\v\f");
	}
	
	fclose(fh);	
	
	FreeWords(words, MAX_WORD_IN_LINE);
	
	printf("%d sgRNAs read.\n", totalItemNum);
	
	return totalItemNum;
}

//transform to log mean-ratio. m = x1'+x2', r = x2'-x1', x' = log2(x/median+0.01), 0.01 is the pseudo-count
int ComputeMR(ITEM_STRUCT *items, int itemNum)
{
	double *tmpF;
	int i;
	double median1, median2,x1ba, x2ba;
	
	assert(itemNum>0);
	
	tmpF = (double *)malloc(itemNum*sizeof(double));
	
	assert(tmpF!=NULL);
	
	for (i=0;i<itemNum;i++)
	{
		tmpF[i] = items[i].x1;
	}
	
	QuicksortF(tmpF, 0, itemNum-1);
	
	median1 = tmpF[(itemNum+1)/2];
	
	for (i=0;i<itemNum; i++)
	{
		tmpF[i] =items[i].x2;
	}
	
	QuicksortF(tmpF, 0, itemNum-1);
	
	median2 = tmpF[(itemNum+1)/2];
	
	for (i=0;i<itemNum;i++)
	{
		x1ba = log2(items[i].x1/median1+0.01);
		x2ba = log2(items[i].x2/median2+0.01);
		
		items[i].m = x1ba + x2ba;
		items[i].r = x2ba - x1ba;
	}
	
	free(tmpF);
	
	return 1;
}

//QuickSort items by log-mean
void QuickSortItemByM(ITEM_STRUCT *items, int lo, int hi)
{
	int i=lo, j=hi;
	ITEM_STRUCT h;
	double x=items[(lo+hi)/2].m;
	
	if (hi<lo)
	{
		return;
	}
	
    //  partition
    while (i<=j)
    {    
		while ((items[i].m<x)&&(i<=j))
		{
			i++;
		}
		while ((items[j].m>x)&&(i<=j))
		{
			j--;
		}
        if (i<=j)
        {
			memcpy(&h,items+i,sizeof(ITEM_STRUCT));
			memcpy(items+i,items+j,sizeof(ITEM_STRUCT));
			memcpy(items+j,&h,sizeof(ITEM_STRUCT));
            i++; j--;
        }
    } 
	
    //  recursion
    if (lo<j) QuickSortItemByM(items, lo, j);
    if (i<hi) QuickSortItemByM(items, i, hi);	
}

//Adjust r using z-transform within a window sliding on items sorted by m
int AdjustMR(ITEM_STRUCT *items, int itemNum, int winSize)
{
	int i,j;
	double tmpMean, tmpStdev;
	int index1, index2, tmpRange;
	ITEM_STRUCT *tmpItems;
	double *tmpM;
	
	assert(itemNum>0);
	assert(winSize>0);
	
	tmpItems = (ITEM_STRUCT *)malloc(itemNum*sizeof(ITEM_STRUCT));
	tmpM = (double *)malloc(itemNum*sizeof(double));
	
	assert(tmpItems!=NULL);
	assert(tmpM!=NULL);
	
	memcpy(tmpItems, items, itemNum*sizeof(ITEM_STRUCT));
							
	QuickSortItemByM(items, 0, itemNum);
	
	for (i=0;i<itemNum;i++)
	{
		tmpM[i] = items[i].m;
	}
	
	for (i=0;i<itemNum;i++)
	{
		index1 = bTreeSearchingF(tmpItems[i].m-0.000000001, tmpM, 0, itemNum-1);
		index2 = bTreeSearchingF(tmpItems[i].m+0.000000001, tmpM, 0, itemNum-1);
		
		tmpRange = index2-index1+1;
		
		if (tmpRange<winSize)
		{
			index1 = index1-(winSize-tmpRange+1)/2;
			index2 = index2+(winSize-tmpRange+1)/2;
			
			index1 = index1>=0?index1:0;
			index2 = index2<itemNum?index2:itemNum-1;
		}
		
		index1 = bTreeSearchingF(tmpM[index1]-0.000000001, tmpM, 0, itemNum-1);
		index2 = bTreeSearchingF(tmpM[index2]+0.000000001, tmpM, 0, itemNum-1);
		
		assert(index1>=0);
		assert(index2<itemNum);
		
		tmpMean = 0.0;
		
		for (j=index1;j<=index2;j++)
		{
			tmpMean += items[j].r;
		}
		
		tmpMean = tmpMean/(index2-index1+1);
		
		tmpStdev = 0.0;
		
		for (j=index1;j<index2;j++)
		{
			tmpStdev += (items[j].r-tmpMean)*(items[j].r-tmpMean);
		}
		
		tmpStdev = sqrt(tmpStdev/(index2-index1+1));
		
		tmpItems[i].adjustedR = (tmpItems[i].r-tmpMean)/(tmpStdev+0.000000001);
	}
	
	memcpy(items, tmpItems, itemNum*sizeof(ITEM_STRUCT));
	
	free(tmpItems);
	free(tmpM);
	
	return 1;
}

//Save results to file. Format: <sgRNA id> <gene id> <measure in library 1> <measure in library 2> <normalized measure in library 1> <normalized measure in library 2> <mean> <ratio> <adjusted ratio>
int SaveToOuput(char *fileName, ITEM_STRUCT *items, int itemNum)
{
	FILE *fh;
	int i;
	
	fh = (FILE *)fopen(fileName, "w");
	
	if (!fh)
	{
		printf("Cannot open %s.\n", fileName);
		return -1;
	}
	
	fprintf(fh, "sgRNA_id\tgene_id\tmeasure_lib1\tmeasure_lib2\tnorm_measure_lib1\tnorm_measure_lib2\tmean\tratio\tadjusted_ratio\n");
	
	for (i=0;i<itemNum;i++)
	{
		fprintf(fh, "%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
				items[i].sgName,
				items[i].geneName,
				items[i].x1,
				items[i].x2,
				items[i].m-items[i].adjustedR/2,
				items[i].m+items[i].adjustedR/2,
				items[i].m,
				items[i].r,
				items[i].adjustedR);
	}
	
	fclose(fh);
	
	return 1;
}

int main (int argc, const char * argv[]) 
{
	int i, winSize;
	ITEM_STRUCT *items;
	int itemNum;
	char inputFileName[1000], outputFileName[1000];
	
	//Parse the command line
	if (argc == 1)
	{
		PrintCommandUsage(argv[0]);
		return -1;
	}
	
	inputFileName[0] = 0;
	outputFileName[0] = 0;
	winSize = 200;
	
	for (i=2;i<argc;i++)
	{
		if (strcmp(argv[i-1], "-i")==0)
		{
			strcpy(inputFileName, argv[i]);
		}
		if (strcmp(argv[i-1], "-o")==0)
		{
			strcpy(outputFileName, argv[i]);
		}
		if (strcmp(argv[i-1], "-w")==0)
		{
			winSize = atoi(argv[i]);
		}
	}
	
	if ((inputFileName[0]==0)||(outputFileName[0]==0))
	{
		printf("Command error!\n");
		PrintCommandUsage(argv[0]);
		return -1;
	}
	
	printf("read input file...");
	itemNum = ReadFile(inputFileName, &items);
	
	if (itemNum<=0)
	{
		printf("\nfailed.\n");
		printf("program exit!\n");
		
		return -1;
	}
	else
	{
		printf("done.\n");
	}
	
	printf("normalizing...");
	
	if (ComputeMR(items, itemNum)<=0)
	{
		printf("\nfailed.\n");
		printf("program exit!\n");
		
		return -1;
	}
	
	if (AdjustMR(items, itemNum, winSize)<=0)
	{
		printf("\nfailed.\n");
		printf("program exit!\n");
		
		return -1;
	}
	else
	{
		printf("done.\n");
	}
	
	printf("save to output file...");
	
	if (SaveToOuput(outputFileName, items, itemNum)<=0)
	{
		printf("\nfailed.\n");
		printf("program exit!\n");
		
		return -1;
	}
	else
	{
		printf("done.\n");
	}
	
	printf("finished.\n");
	
	free(items);
	
	return 0;
	
}

//print the usage of Command
void PrintCommandUsage(const char *command)
{
	//print the options of the command
	printf("%s - Crispr data normalization.\n", command);
	printf("usage:\n");
	printf("-i <input data file>. Format: <sgRNA id> <gene id> <measure in library 1> <measure in library 2>\n");
	printf("-o <output file>. Format: <sgRNA id> <gene id> <measure in library 1> <measure in library 2> <normalized measure in library 1> <normalized measure in library 2> <mean> <ratio> <adjusted ratio>\n");
	printf("-w <window size>. Default:200\n");
	printf("example:\n");
	printf("%s -i input.txt -o output.txt -w 200\n", command);
	
}
