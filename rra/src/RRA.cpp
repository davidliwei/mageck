/*
 *  RRA.c
 *  Implementation of Robust Rank Aggregation (RRA)
 *
 *  Created by Han Xu on 20/12/13.
 *  Modified by Wei Li.
 *  07/2014: updated the constants to accomodate pathway names
 *  Copyright 2013 Dana Farber Cancer Institute. All rights reserved.
 *
 */
#include "math_api.h"
#include "words.h"
#include "rvgs.h"
#include "rngs.h"
#include "classdef.h"
#include "fileio.h"

//C++ functions
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
using namespace std;

int PRINT_DEBUG=0;

//Global variables

//store control sequences
bool UseControlSeq=false;
map<string,int> ControlSeqMap; //save the control sequence name and their index in ControlSeqPercentile;
double* ControlSeqPercentile;

// used for calculation of lo-values
double* tmpLovarray=NULL;
int nLovarray=-1;


//Function declarations

//Process groups by computing percentiles for each item and lo-values for each group
int ProcessGroups(GROUP_STRUCT *groups, int groupNum, LIST_STRUCT *lists, int listNum, double maxPercentile);

//QuickSort groups by loValue
void QuickSortGroupByLoValue(GROUP_STRUCT *groups, int start, int end);

//Compute False Discovery Rate based on uniform distribution
int ComputeFDR(GROUP_STRUCT *groups, int groupNum, double maxPercentile, int numOfRandPass);

//print the usage of Command
void PrintCommandUsage(const char *command);


//Compute lo-value based on an array of percentiles
int ComputeLoValue(double *percentiles,     //array of percentiles
				   int num,                 //length of array
				   double &loValue,         //pointer to the output lo-value
				   double maxPercentile,   //maximum percentile, computation stops when maximum percentile is reached
           int &goodsgrna);// # of good sgRNAs

//WL: modification of lo_value computation
int ComputeLoValue_Prob(double *percentiles,     //array of percentiles
				   int num,                 //length of array
				   double &loValue,         //pointer to the output lo-value
				   double maxPercentile,   //maximum percentile, computation stops when maximum percentile is reached
				  double *probValue,// probability of each prob, must be equal to the size of percentiles
           int &goodsgrna);

//read control sequences
int loadControlSeq(const char* fname){
	ifstream fh;
  fh.open(fname);
  if(!fh.is_open()){
    cerr<<"Error opening "<<fname<<endl;
  }
  string oneline;
  int ncount=0;
  ControlSeqMap.clear();
  while(!fh.eof()){
    getline(fh,oneline);
    ControlSeqMap[oneline]=ncount;
    ncount++;
  }
  fh.close();
  ControlSeqPercentile = new double[ncount];
  for(int i=0;i<ncount;i++) ControlSeqPercentile[i]=-1.0;
  cout<<ControlSeqMap.size()<<" control sequences loaded.\n";
  return 0;
}

int main (int argc, const char * argv[]) {
	int i,flag;
	GROUP_STRUCT *groups;
	int groupNum;
	LIST_STRUCT *lists;
	int listNum;
	char inputFileName[1000], outputFileName[1000];
	double maxPercentile;
	
	//Parse the command line
	if (argc == 1)
	{
		PrintCommandUsage(argv[0]);
		return 0;
	}
	
	inputFileName[0] = 0;
	outputFileName[0] = 0;
	maxPercentile = 0.1;
	
	for (i=2;i<argc;i++)
	{
		if (strcmp(argv[i-1], "-i")==0){
			strcpy(inputFileName, argv[i]);
		}
		if (strcmp(argv[i-1], "-o")==0){
			strcpy(outputFileName, argv[i]);
		}
		if (strcmp(argv[i-1], "-p")==0){
			maxPercentile = atof(argv[i]);
		}
		if (strcmp(argv[i-1], "--control")==0){
       UseControlSeq=true;
       // load control sequences
       loadControlSeq(argv[i]);
    }
	}
	
	if ((inputFileName[0]==0)||(outputFileName[0]==0))
	{
		cerr<<"Error: input file or output file name not set.\n";
		PrintCommandUsage(argv[0]);
		return -1;
	}
	
	if ((maxPercentile>1.0)||(maxPercentile<0.0))
	{
		cerr<<("Error: maxPercentile should be within 0.0 and 1.0\n");
		return -1;
	}
	
	groups = (GROUP_STRUCT *)malloc(MAX_GROUP_NUM*sizeof(GROUP_STRUCT));
	lists = (LIST_STRUCT *)malloc(MAX_LIST_NUM*sizeof(LIST_STRUCT));
	assert(groups!=NULL);
	assert(lists!=NULL);
	
	printf("Reading input file...\n");
	
	flag = ReadFile(inputFileName, groups, MAX_GROUP_NUM, &groupNum, lists, MAX_LIST_NUM, &listNum);
	
	if (flag<=0){
	  cerr<<"\nError: reading ranking file ...\n";
		return -1;
	}
	
	cerr<<("Computing lo-values for each group...\n");
	
	if (ProcessGroups(groups, groupNum, lists, listNum, maxPercentile)<=0)
  {
		cerr<<("\nError: processing groups failed.\n");
		return -1;
	}
	
	cerr<<("Computing false discovery rate...\n");
	
	if (ComputeFDR(groups, groupNum, maxPercentile, RAND_PASS_NUM*groupNum)<=0)
	{
		cerr<<("\nError: computing FDR failed.\n");
		return -1;
	}
	
	cerr<<("Saving to output file...");
	
	if (SaveGroupInfo(outputFileName, groups, groupNum)<=0)
	{
		cerr<<("\nError: saving output file failed.\n");
		return -1;
	}
	
	cerr<<("RRA completed.\n");
	
  for(i=0;i<groupNum;i++)
    delete[] groups[i].items;
	free(groups);
	
	for (i=0;i<listNum;i++)
	{
		delete[] lists[i].values;
	}
	free(lists);
  if(UseControlSeq){
    delete[] ControlSeqPercentile;
  }
  if(nLovarray>0){
    delete []tmpLovarray;
  }

	return 0;

}

//print the usage of Command
void PrintCommandUsage(const char *command)
{
	//print the options of the command
	printf("%s - Robust Rank Aggreation.\n", command);
	printf("usage:\n");
	printf("-i <input data file>. Format: <item id> <group id> <list id> <value> [<probability>] [<chosen>]\n");
	printf("-o <output file>. Format: <group id> <number of items in the group> <lo-value> <false discovery rate>\n");
	printf("-p <maximum percentile>. RRA only consider the items with percentile smaller than this parameter. Default=0.1\n");
	printf("--control <control_sgrna list>. A list of control sgRNA names.\n");
	printf("example:\n");
	printf("%s -i input.txt -o output.txt -p 0.1 \n", command);
	
}




//Process groups by computing percentiles for each item and lo-values for each group
//groups: genes
//lists: a set of different groups. Comparison will be performed on individual list
int ProcessGroups(GROUP_STRUCT *groups, int groupNum, LIST_STRUCT *lists, int listNum, double maxPercentile)
{
	int i,j;
	int listIndex, index1, index2;
	int maxItemPerGroup;
	double *tmpF;
	double *tmpProb;
	bool isallone; // check if all the probs are 1; if yes, do not use accumulation of prob. scores
	
	maxItemPerGroup = 0;

  PRINT_DEBUG=1;
	
	for (i=0;i<groupNum;i++){
		if (groups[i].itemNum>maxItemPerGroup){
			maxItemPerGroup = groups[i].itemNum;
		}
	}
	
	assert(maxItemPerGroup>0);
	
	tmpF = new double[maxItemPerGroup];
	tmpProb = new double [maxItemPerGroup];
	
	for (i=0;i<listNum;i++){
		QuicksortF(lists[i].values, 0, lists[i].itemNum-1);
	}
	
	for (i=0;i<groupNum;i++){
		//Compute percentile for each item
		
		isallone=true;
    int validsgs=0;
		for (j=0;j<groups[i].itemNum;j++){
      if(groups[i].items[j].isChosen==0) continue;
			listIndex = groups[i].items[j].listIndex;
			
			index1 = bTreeSearchingF(groups[i].items[j].value-0.000000001, lists[listIndex].values, 0, lists[listIndex].itemNum-1);
			index2 = bTreeSearchingF(groups[i].items[j].value+0.000000001, lists[listIndex].values, 0, lists[listIndex].itemNum-1);
			
			groups[i].items[j].percentile = ((double)index1+index2+1)/(lists[listIndex].itemNum*2);
			tmpF[validsgs] = groups[i].items[j].percentile;
			tmpProb[validsgs]=groups[i].items[j].prob;
			if(tmpProb[validsgs]!=1.0){
				isallone=false;
			}
      // save for control sequences
      if(UseControlSeq){
        string sgname(groups[i].items[j].name);
        if(ControlSeqMap.count(sgname)>0){
          int sgindex=ControlSeqMap[sgname];
          ControlSeqPercentile[sgindex]=tmpF[validsgs];
        }
      }//end if
      validsgs++;
		}// end j
    if(validsgs<=1){
      isallone=true;
    }
    //printf("Gene: %s\n",groups[i].name);
		if(isallone){
			ComputeLoValue(tmpF, validsgs, groups[i].loValue, maxPercentile, groups[i].goodsgrnas);
		}
		else{
			ComputeLoValue_Prob(tmpF, validsgs, groups[i].loValue, maxPercentile, tmpProb,groups[i].goodsgrnas);
		}
    groups[i].isbad=0;
	}//end i loop

	delete[] tmpF;
	delete[] tmpProb;
  //check if all control sequences are properly assigned a value
  if(UseControlSeq){
    for(map<string,int>::iterator mit = ControlSeqMap.begin(); mit != ControlSeqMap.end(); mit++){
      if(ControlSeqPercentile[mit->second]<0){
        cerr<<"Warning: sgRNA "<<mit->first<<" not found in the ranked list. \n";
        //ControlSeqPercentile[mit->second]=0.5;
      }
    }
  }//end UseControlSeq
	
	return 1;
}

//Compute lo-value based on an array of percentiles. Return 1 if success, 0 if failure
int ComputeLoValue(double *percentiles,     //array of percentiles
				   int num,                 //length of array
				   double &loValue,         //pointer to the output lo-value
				   double maxPercentile,   //maximum percentile, computation stops when maximum percentile is reached
           int &goodsgrna){   // the number of " good" sgRNAs
	int i;
	double *tmpArray;
	double tmpLoValue, tmpF;
	
	if(num==0){
    loValue=1.00;
    goodsgrna=0;
    return 0;
  }
  if(num>nLovarray){
    delete[] tmpLovarray;
    tmpLovarray=new double[num];
    nLovarray=num;
  }
	//tmpArray = (double *)malloc(num*sizeof(double));
  tmpArray=tmpLovarray;
	
	if (!tmpArray){
		return -1;
	}
	
	memcpy(tmpArray, percentiles, num*sizeof(double));
	
	QuicksortF(tmpArray, 0, num-1);
	
	tmpLoValue = 1.0;
  goodsgrna=0;
	
	for (i=0;i<num;i++){
		// if ((tmpArray[i]>maxPercentile)&&(i>0)) //only calculate the value if at least 1 percentile is smaller than the cutoff
		if ((tmpArray[i]>maxPercentile)){
      if(i==0){}
      else{
			  break;
      }
		}else{
      goodsgrna++;
    }
		tmpF = BetaNoncentralCdf((double)(i+1),(double)(num-i),0.0,tmpArray[i],CDF_MAX_ERROR);
		if (tmpF<tmpLoValue){
			tmpLoValue = tmpF;
		}
	}
	
	loValue = tmpLoValue;

	//free(tmpArray);
	
	return 0;
	
}

//Compute lo-value based on an array of percentiles, by considering the probability of each sgRNAs. 
//Return 1 if success, 0 if failure
//Modified by Wei Li
int ComputeLoValue_Prob(double *percentiles,     //array of percentiles
				   int num,                 //length of array
				   double &loValue,         //pointer to the output lo-value
				   double maxPercentile,   //maximum percentile, computation stops when maximum percentile is reached
				  double *probValue,
          int &goodsgrna) {// probability of each prob, must be equal to the size of percentiles

	int i, pid;
	double *tmpArray;
	double tmpLoValue, tmpF;

	// used to calculate cumulative probability
	int real_i;
	int c_num;
	double c_prob;
	double accuLoValue;

	if(num==0){
    loValue=1.00;
    goodsgrna=0;
    return 0;
  }
		
  if(num>nLovarray){
    delete[] tmpLovarray;
    tmpLovarray=new double[num];
    nLovarray=num;
  }
	//tmpArray = (double *)malloc(num*sizeof(double));
  tmpArray=tmpLovarray;
	
	if (!tmpArray){
		return -1;
	}
	
	memcpy(tmpArray, percentiles, num*sizeof(double));
	
	QuicksortF(tmpArray, 0, num-1);
  
  goodsgrna=0;
  for (i=0;i<num;i++){
		if ((tmpArray[i]>maxPercentile)){
      if(i==0){}else{
			  break;
      }
		}else{
      goodsgrna++;
    }
  }

  
  if(PRINT_DEBUG) printf("probs:");
  for(i=0;i<num;i++)
  {
    if(PRINT_DEBUG) printf("%f,",probValue[i]);
  }
  if(PRINT_DEBUG) printf("\n");
	
	
	accuLoValue=0.0;
	for (pid=1;pid<(1<<num);pid++)
	{
		// decoding the selection
		tmpLoValue = 1.0;
		c_num=0;
		real_i=0;
		c_prob=1.0;
		for(i=0;i<num;i++)
		{
			if ( (1<<i)&pid )
			{
				c_num = c_num +1;
				c_prob= c_prob*probValue[i];
			}
			else
			{
				c_prob=c_prob*(1.0-probValue[i]);
			}
		}
		for (i=0;i<num;i++)
		{
			if ( ((1<<i)&pid) ==0)
			{// if this sgRNA is selected?
				continue;
			}
			if ((tmpArray[i]>maxPercentile)&&(i>0))
			{
				break;
			}
			// Beta (a,b, 0, frac)
			tmpF = BetaNoncentralCdf((double)(real_i+1),(double)(c_num-i),0.0,tmpArray[i],CDF_MAX_ERROR);
			if (tmpF<tmpLoValue)
			{
				tmpLoValue = tmpF;
			}
			real_i = real_i + 1;
		}
    if(PRINT_DEBUG) printf("pid: %d, prob:%e, score: %f\n",pid,c_prob,tmpLoValue);
		accuLoValue=accuLoValue+tmpLoValue*c_prob;
	}
  if(PRINT_DEBUG) printf("total: %f\n",accuLoValue);
	
	loValue = accuLoValue;

	//free(tmpArray);
	return 0;
	
}


//Compute False Discovery Rate based on uniform distribution
int ComputeFDR(GROUP_STRUCT *groups, int groupNum, double maxPercentile, int numOfRandPass)
{
	int i,j,k;
	double *tmpPercentile;
	int maxItemNum = 0;
	int scanPass = numOfRandPass/groupNum+1;
	double *randLoValue;
	int randLoValueNum;

	//WL
	double *tmpProb;
	double isallone;
	
	for (i=0;i<groupNum;i++){
		if (groups[i].itemNum>maxItemNum){
			maxItemNum = groups[i].itemNum;
		}
	}
	
	assert(maxItemNum>0);
	
	//tmpPercentile = (double *)malloc(maxItemNum*sizeof(double));
	//tmpProb= (double *)malloc(maxItemNum*sizeof(double));
  tmpPercentile=new double[maxItemNum];
  tmpProb=new double[maxItemNum];	

	randLoValueNum = groupNum*scanPass;
	
	assert(randLoValueNum>0);
	
	//randLoValue = (double *)malloc(randLoValueNum*sizeof(double));
  randLoValue=new double[randLoValueNum];
	
	randLoValueNum = 0;
	
	PlantSeeds(123456);
  
  PRINT_DEBUG=0;
  
  // set up control sequences
  int n_control=0;
  double* control_prob_array=NULL;
  double ufvalue=0.5;
  int rand_ctl_index=0;
  int tmp_int;
	if(UseControlSeq){
    for(map<string,int>::iterator mit = ControlSeqMap.begin(); mit != ControlSeqMap.end(); mit++){
      if(ControlSeqPercentile[mit->second]>=0){
        n_control++;
      }
    }
    control_prob_array=new double[n_control];
    n_control=0;
    for(map<string,int>::iterator mit = ControlSeqMap.begin(); mit != ControlSeqMap.end(); mit++){
      if(ControlSeqPercentile[mit->second]>=0){
        control_prob_array[n_control]=ControlSeqPercentile[mit->second];
        n_control++;
      }
    }
    cout<<"Total # control sgRNAs: "<<n_control<<endl;
  }
  
	for (i=0;i<scanPass;i++){
    for (j=0;j<groupNum;j++){
			isallone=true;
      int validsgs=0;
			for (k=0;k<groups[j].itemNum;k++)
			{
        if(groups[j].items[k].isChosen==0) continue;
        ufvalue=Uniform(0.0, 1.0);
        if(UseControlSeq){
          rand_ctl_index=(int)(n_control*ufvalue);
          if(rand_ctl_index>=n_control) rand_ctl_index=n_control-1;
          tmpPercentile[validsgs]=control_prob_array[rand_ctl_index];
        }else{
				  tmpPercentile[validsgs] = ufvalue;
        }
				tmpProb[validsgs]=groups[j].items[k].prob;
				if(tmpProb[validsgs]!=1.0)
				{
					isallone=false;
				}
        validsgs++;
			} //end for k
      if(validsgs<=1)
        isallone=true;
			
			if(isallone){
				ComputeLoValue(tmpPercentile, validsgs,randLoValue[randLoValueNum], maxPercentile, tmp_int);
			}
			else
			{
				ComputeLoValue_Prob(tmpPercentile, validsgs,randLoValue[randLoValueNum], maxPercentile,tmpProb,tmp_int);
			}
			
			randLoValueNum++;
		}// end for j
	}//end for i
	
	QuicksortF(randLoValue, 0, randLoValueNum-1);
						  
	QuickSortGroupByLoValue(groups, 0, groupNum-1);
	
  //FDR calcuation
  int goodGroupNum=0;
  for(i=0;i<groupNum;i++){
     //if(groups[i].isbad==0) 
    goodGroupNum+=1;
  }
  //save the index
  if(goodGroupNum==0) goodGroupNum=1;
  cout<<"Number of groups under permutation and FDR adjustment: "<<goodGroupNum<<endl;
  int* indexval=new int[goodGroupNum];
  int goodindex=0;
	for (i=0;i<groupNum;i++){
		//groups[i].fdr = (double)(bTreeSearchingF(groups[i].loValue-0.000000001, randLoValue, 0, randLoValueNum-1)
		//						 +bTreeSearchingF(groups[i].loValue+0.000000001, randLoValue, 0, randLoValueNum-1)+1)
		//						/2/randLoValueNum/((double)i+0.5)*groupNum;
    if(groups[i].isbad==0){
      groups[i].pvalue=(double)(bTreeSearchingF(groups[i].loValue-0.000000001, randLoValue, 0, randLoValueNum-1)
								 +bTreeSearchingF(groups[i].loValue+0.000000001, randLoValue, 0, randLoValueNum-1)+1)
								/2/randLoValueNum;
		  groups[i].fdr = groups[i].pvalue/((double)i+1.0)*goodGroupNum;
      indexval[goodindex]=i;
      goodindex++;
    }else{
      groups[i].pvalue=1.0;
      groups[i].fdr=1.0;
    }
	}
	
	if (groups[groupNum-1].fdr>1.0)
	{
		groups[groupNum-1].fdr = 1.0;
	}
	
	for (i=goodGroupNum-2;i>=0;i--){
    int g_i=indexval[i];
    int g_ip1=indexval[i+1];
		if (groups[g_i].fdr>groups[g_ip1].fdr){
			groups[g_i].fdr = groups[g_ip1].fdr;
		}
	}
	
	//free(tmpPercentile);
	//free(tmpProb);
	//free(randLoValue);
  delete []tmpPercentile;
  delete []tmpProb;
  delete []randLoValue;
  
  if(UseControlSeq){
    delete[] control_prob_array;
  }
  delete []indexval;
	
	return 1;
}

//QuickSort groups by loValue
void QuickSortGroupByLoValue(GROUP_STRUCT *groups, int lo, int hi)
{
	int i=lo, j=hi;
	GROUP_STRUCT tmpGroup;
	double x=groups[(lo+hi)/2].loValue;
	
	if (hi<lo)
	{
		return;
	}
	
    //  partition
    while (i<=j)
    {    
		while ((groups[i].loValue<x)&&(i<=j))
		{
			i++;
		}
		while ((groups[j].loValue>x)&&(i<=j))
		{
			j--;
		}
        if (i<=j)
        {
			memcpy(&tmpGroup,groups+i,sizeof(GROUP_STRUCT));
			memcpy(groups+i,groups+j,sizeof(GROUP_STRUCT));
			memcpy(groups+j,&tmpGroup,sizeof(GROUP_STRUCT));
            i++; j--;
        }
    } 
	
    //  recursion
    if (lo<j) QuickSortGroupByLoValue(groups, lo, j);
    if (i<hi) QuickSortGroupByLoValue(groups, i, hi);
	
}
