# 
#
from __future__ import print_function;

import sys;
import re;
import os;
import logging;


class VisualRValue:
  '''
  Class for visualization
  '''
  outprefix='sample1';
  genesummaryfile='';
  cpindex=[];
  targetgene=[];
  cplabel=''
  
  # internal variable
  rtemplatestr='';  # template string
  rtemplate_gene_str='';  # template string
  outrfh=None;  # file handle for R file
  '''
  Member functions
  '''
  def startRTemplate(self):
    '''
    Open a template, create an R file
    '''
    filename=os.path.join(os.path.dirname(__file__),'plot_template.RTemplate');
    if os.path.isfile(filename) and os.path.exists(filename):
      logging.info('Loading R template file: '+filename+'.');
    else:
      logging.error('Cannot find template file: '+filename);
      return -1;
    filename_indgene=os.path.join(os.path.dirname(__file__),'plot_template_indvgene.RTemplate');
    if os.path.isfile(filename_indgene) and os.path.exists(filename_indgene):
      logging.info('Loading R template file: '+filename_indgene+'.');
    else:
      logging.error('Cannot find template file: '+filename_indgene);
      return -1;
    logging.debug('Setting up the visualization module...');
    with open(filename, "r") as rtfile:
      rtp=rtfile.read();
      outpdffile=self.outprefix+'.pdf';
      # rtp=re.sub('__OUTPUT_FILE__',outpdffile,rtp); # pdf file
      #
      outrfile=self.outprefix+'.R';
      
      # write to R file
      outrfh=open(outrfile,'w');
      self.outrfh=outrfh;
      self.rtemplatestr=rtp;
      
      # write pdf loading
      pdfloadstr="pdf(file='"+outpdffile+"',width=6,height=6);";
      # write file reading
      # rtp=re.sub('__GENE_SUMMARY_FILE__',self.genesummaryfile,rtp); # gene_summary
      tableadstr="gstable=read.table('"+self.genesummaryfile+"',header=T)";
      print(pdfloadstr,file=outrfh);
      print(tableadstr,file=outrfh);
    # load individual gene code
    with open(filename_indgene, "r") as rtfile:
      rtp=rtfile.read();
      self.rtemplate_gene_str=rtp;
      
    return 0;
  
  def closeRTemplate(self):
    '''
    Close the R file
    '''
    print("dev.off()",file=self.outrfh);
    self.outrfh.close();
    
  def WriteRTemplate(self):
    '''
    Given a VisualRValue object, write an R file
    '''
    
    # load file
    rtp=self.rtemplatestr;
    # replace the variables
    
    
    indexchr=','.join([str(x) for x in self.cpindex]);
    rtp=re.sub('__INDEX_COLUMN__',indexchr,rtp); # index 
    logging.debug('Column index:'+indexchr);
    
    targetchr="c("+','.join(['"'+x+'"' for x in self.targetgene])+")";
    rtp=re.sub('__TARGET_GENE_LIST__',targetchr,rtp); # index 
    
    rtp=re.sub('__SAMPLE_LABEL__',self.cplabel,rtp);
    
    # write to R file
    print(rtp,file=self.outrfh);
  
  def loadTopK(self, filename, k=10):
    '''
    Load the top k gene names from the file
    '''
    n=0;
    self.targetgene=[];
    for line in open(filename):
      n+=1;
      if n==1:
        continue;
      if n<=k:
        field=line.strip().split();
        tgenename=field[0];
        self.targetgene+=[tgenename];
      else:
        break;
    # write to file?
    logging.info('Loading top '+str(k) +' genes from '+filename+': '+','.join(self.targetgene));
    self.WriteRTemplate();
    return 0;
  
  def loadTopKWithExp(self,filename,nttab,sgrna2genelist,collabels,k=10):
    self.loadTopK(filename,k);
    self.loadGeneExp(self.targetgene,nttab,sgrna2genelist,collabels);
  
  def loadGeneExp(self,genelist,nttab,sgrna2genelist,collabels):
    '''
    Load the sgRNA read counts of selected genes into file
    '''
    for gene in genelist:
      sglist=[ k for (k,v) in sgrna2genelist.iteritems() if v==gene];
      ntgene={k:v for (k,v) in nttab.iteritems() if k in sglist};
      # load to file
      valstring='list(';
      vstrlist=[];
      for (k,v) in ntgene.iteritems():
        vstr='c('+','.join([str(vv) for vv in v])+')';
        vstrlist+=[vstr];
      valstring+=','.join(vstrlist);
      valstring+=')';
      rtp=self.rtemplate_gene_str;
      rtp=re.sub('__TARGET_GENE__','"'+gene+'"',rtp);
      rtp=re.sub('__TARGET_MATRIX__',valstring,rtp);
      # labels
      clabstr='c(';
      clabstr+=','.join(['"'+ x+'"' for x in collabels]);
      clabstr+=')';
      rtp=re.sub('__COL_LABEL__',clabstr,rtp);
      # save to file
      print(rtp,file=self.outrfh);
    

