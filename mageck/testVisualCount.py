# 
#
from __future__ import print_function;

import sys;
import re;
import os;
import logging;
from fileOps import *
from mageckCount import *


class VisualRCount:
  '''
  Class for generating reports of count command 
  '''
  outprefix='sample1';
  
  # internal variable, for R file
  outrfh=None;  # file handle for R file
  
  # for Rnw file
  rnwtemplatestr='';
  outrnwfh=None;
  outrnwstring='';
  # for statistics of coutns
  fastqfile=[]; # fastq files
  fastqlabels=[]; # fastq labels
  reads=[]; # read counts
  mappedreads=[]; # 
  zerocounts=[];
  '''
  Member functions
  '''
  def startRTemplate(self):
    '''
    Open a template, create an R file
    
    '''
    # R files
    # Rnw files
    filename_rnw=os.path.join(os.path.dirname(__file__),'fastq_template.Rnw');
    if os.path.isfile(filename_rnw) and os.path.exists(filename_rnw):
      logging.info('Loading Rnw template file: '+filename_rnw+'.');
    else:
      logging.error('Cannot find template file: '+filename_rnw);
      return -1;
    logging.debug('Setting up the visualization module...');
    # 
    
    # R file
    outrfile=self.outprefix+'_countsummary.R';
    outrfh=open(outrfile,'w');
    self.outrfh=outrfh;
    # load Rnw file
    with open(filename_rnw,"r") as rtfile:
      rnw=rtfile.read();
      self.rnwtemplatestr=rnw;
      outrfile=self.outprefix+'_countsummary.Rnw';
      self.outrnwstring=self.rnwtemplatestr;
      outrfh=open(outrfile,'w');
      self.outrnwfh=outrfh;
      
    return 0;
  
  def closeRTemplate(self):
    '''
    Close the R file
    '''
    # write to R file
    #
    rnwfile=self.outprefix+'_countsummary.Rnw';
    rfile=self.outprefix+'_countsummary.R';
    summaryfile=self.outprefix+'_countsummary';
    latexfile=self.outprefix+'_countsummary.tex';
    (rnwfile_dir,rnwfile_base)=os.path.split(rnwfile);
    # write code in R file to generate PDF files
    print("Sweave(\""+rnwfile_base+"\");\nlibrary(tools);\n",file=self.outrfh);
    print("texi2dvi(\""+os.path.basename(latexfile)+"\",pdf=TRUE);\n",file=self.outrfh);
    # write to Rnw file
    print(self.outrnwstring,file=self.outrnwfh);
    
    self.outrnwfh.close();
    self.outrfh.close();
    



  def writeCountSummary(self):
    '''
    Write statistics from gene summary file to buffer
    '''
    # insert string
    insertstr='';
    insertstr+='filelist=c(' + ','.join(['"'+x+'"' for x in self.fastqfile])  +');\n';
    insertstr+='labellist=c('+ ','.join(['"'+x+'"' for x in self.fastqlabels]) +');\n';
    insertstr+='reads=c('+','.join([str(x) for x in self.reads])+');\n';
    insertstr+='mappedreads=c('+','.join([str(x) for x in self.mappedreads])+');\n';
    insertstr+='zerocounts=c('+','.join([str(x) for x in self.zerocounts])+');\n';
    # 
    nwktowrite=re.sub('#__COUNT_SUMMARY_STAT__',insertstr,self.outrnwstring);
    self.outrnwstring=nwktowrite;
  
  def insertReadCountBoxPlot(self,filename):
    '''
    Insert box plot of normalized read counts
    '''
    # inserted R code
    rtp="\ngenboxplot(\""+filename+"\");\n";
    #
    insertstr='';
    
    insertstr+=r"\n\\newpage\\section{Normalized read count distribution of all samples}\n";
    insertstr+="The following figure shows the distribution of median-normalized read counts in all samples.\n\n\n";
    
    insertstr+="<<fig=TRUE,echo=FALSE,width=4.5,height=4.5>>="+rtp+"@"+"\n";
    #
    insertstr+='\n%__INDIVIDUAL_PAGE__\n'
    nwktowrite=re.sub('%__INDIVIDUAL_PAGE__',insertstr,self.outrnwstring);
    self.outrnwstring=nwktowrite;
   
  def insertPCAPlot(self,filename):
    '''
    Insert box plot of PCA analysis 
    '''
    # inserted R code
    rtp="\ngenpcaplot(\""+filename+"\");\n";
    rtp2="\ngenpcavar("+");\n";
    #
    insertstr='';
    
    insertstr+=r"\n\\newpage\\section{Principle Component Analysis}\n";
    insertstr+="The following figure shows the first 2 principle components (PCs) from the Principle Component Analysis (PCA), and the percentage of variances explained by the top PCs.\n\n\n";
    
    insertstr+="\n<<fig=TRUE,echo=FALSE,width=4.5,height=4.5>>="+rtp+"@"+"\n";
    insertstr+="\n<<fig=TRUE,echo=FALSE,width=4.5,height=4.5>>="+rtp2+"@"+"\n";
    #
    insertstr+='\n%__INDIVIDUAL_PAGE__\n'
    nwktowrite=re.sub('%__INDIVIDUAL_PAGE__',insertstr,self.outrnwstring);
    self.outrnwstring=nwktowrite;
    
  
  
  def generatePDF(self,keeptmp=False):
    '''
    Call R and pdflatex
    '''
    rnwfile=self.outprefix+'_countsummary.Rnw';
    rfile=self.outprefix+'_countsummary.R';
    summaryfile=self.outprefix+'_countsummary';
    (rnwfile_dir,rnwfile_base)=os.path.split(rnwfile);
    if rnwfile_dir=='':
      rnwfile_dir='./';
    systemcall('cd '+rnwfile_dir+'; '+'Rscript '+os.path.basename(rfile));
    #systemcall('cd '+rnwfile_dir+'; '+ 'R CMD Sweave '+rnwfile_base);
    #systemcall('export SWEAVE_STYLEPATH_DEFAULT="TRUE";'+ 'cd '+rnwfile_dir+'; '+'pdflatex '+os.path.basename(summaryfile));
    # cleaning the fraction pdf
    if keeptmp==False:
      systemcall('cd '+rnwfile_dir+'; '+'rm -rf '+os.path.basename(summaryfile)+'-*.pdf');
    
    
    


