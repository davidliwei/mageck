#!/usr/bin/env python
"""MAGeCK test module
Copyright (c) 2014 Wei Li, Han Xu, Xiaole Liu lab 
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).
@status:  experimental
@version: $Revision$
@author:  Wei Li 
@contact: li.david.wei AT gmail.com
"""


from __future__ import print_function
import sys
import math
import types
import logging

from mageckCount import *

from fileOps import *
from testVisual import *

def getgeomean(v):
  meanval=sum([math.log(vx+0.1,2) for vx in v])/float(len(v));
  return 2**meanval-0.1;

def getMeans(matt):
  # arithmatic mean
  #meanvalue=[sum(v)/float(len(v)) for v in matt];
  # geometric mean
  meanvalue=[getgeomean(v) for v in matt];
  return meanvalue;

def getVars(matt):
  meanvalue=getMeans(matt);
  varvalue=[ sum([ (kj-meanvalue[i])*(kj-meanvalue[i])  for kj in matt[i]  ]  )/(float(len(matt[i]))-1)    for i in range(len(meanvalue))];
  #varvalue={k:sum([ (x-meanvalue[k])*(x-meanvalue[k]) for x in v])/(float(len(v))-1)  for (k,v) in ctable.iteritems()};
  return varvalue;

def leastsquare(x,y,weight=None):
  """
  least squares fitting
  coefficients from y= a+bx 
  return (b,a)
  reference: http://mathworld.wolfram.com/LeastSquaresFitting.html
  For weighted least square: http://goo.gl/pGpTZ6
  """
  n=len(x);
  if n != len(y):
    logging.error('Unequal length of vectors of x and y in least square');
    sys.exit(-1);
  if weight is None: 
    sy=sum(y);
    sx=sum(x);
    sx2=sum([t*t for t in x]);
    sxy=sum([x[i]*y[i] for i in range(n)]);
    a=(sy*sx2-sx*sxy)/(n*sx2-sx*sx);
    b=(n*sxy-sx*sy)/(n*sx2-sx*sx);
    return (b,a);
  else:
    nw=sum(weight);
    sy=sum([y[i]*weight[i] for i in range(n)]);
    sx=sum([x[i]*weight[i] for i in range(n)]);
    sx2=sum([x[i]*x[i]*weight[i] for i in range(n)]);
    sxy=sum([x[i]*y[i]*weight[i] for i in range(n)]);
    a=(sy*sx2-sx*sxy)/(nw*sx2-sx*sx);
    b=(nw*sxy-sx*sy)/(nw*sx2-sx*sx);
    return (b,a);

def modelmeanvar(ctable,method='edger'):
  """
  model the relation between mean and variance
  """
  # calculate the mean and variance
  tablemat=ctable.values();
  meanvalue=getMeans(tablemat);
  varvalue=getVars(tablemat);
  # choose values with variance greater than mean
  meangood=[meanvalue[i] for i in range(len(meanvalue)) if meanvalue[i]<varvalue[i]];
  vargood=[varvalue[i]-meanvalue[i] for i in range(len(varvalue)) if meanvalue[i]<varvalue[i]];
  # log
  meanglog=[math.log(x+1,2) for x in meangood];
  varglog=[math.log(x+1,2) for x in vargood];
  if method=='linear':
    # least square
    (k,b)=leastsquare(meanglog,varglog,meangood);
    if k<1:
      k=1;
    if b<0:
      b=0;
    return (k,b);
  elif method=='edger':
    dy=varglog;
    dx=[2*x for x in meanglog];
    ret=(sum(dy)-sum(dx))*1.0/len(dx);
    return ret;
  else:
    return 0;
  

def getadjustvar(coef,meanval,method='mixed'):
  """
  From the model, get the adjusted variance
  """
  if method=='linear':
    k=coef[0];b=coef[1];
    if type(meanval) is types.FloatType:
      return (meanval**k)*(2**b)+meanval;
    if type(meanval) is types.ListType:
      return [(z**k)*(2**b)+z for z in meanval];
  elif method=='edger':
    if type(meanval) is types.FloatType:
      return (meanval**2)*(2**coef)+meanval;
    if type(meanval) is types.ListType:
      return [(z**2)*(2**coef)+z for z in meanval];
  elif method=='mixed':
    var1=getadjustvar(coef,meanval,method='linear');
    var2=getadjustvar(coef[2],meanval,method='edger');
    return [ (lambda x,y: x if x>y else y)(var1[i],var2[i]) for i in range(len(var1))];
  else:
    return meanval;

def getnormcdf(x,lowertail=True):
  """
  Get the normal CDF function. used to calculate p-value
  """
  # ax=math.fabs(x);
  #axv=math.erfc(x/(2**0.5))/2; # higher tail
  if lowertail==False:
    #return axv;
    return math.erfc(x/(2**0.5))/2;
  else:
    #return 1-axv;
    return math.erfc(-x/(2**0.5))/2;
  #if (x>0 and lowertail==False) or (x<0 and lowertail==True):
  #  return axv;
  #else:
  #  return 1-axv;

def getNormalPValue(mean0,var0,mean1, lower=False):
  """
  Use truncated normal distribution to calculate the pvalue
  """
  # use ttmean to calculate the pvalue
  n=len(mean0);
  minmean1=min([x for x in mean1 if x>0]);
  mean1_adj=[(lambda x: x if x >minmean1 else minmean1)(t) for t in mean1];
  # first, convert to standard normal distribution values
  t_theta=[(mean1_adj[i]-mean0[i])/math.sqrt(var0[i]) for i in range(n)];
  t_theta_0=[(0.0-mean0[i])/math.sqrt(var0[i]) for i in range(n)];
  #
  t_p=[getnormcdf(x,lowertail=lower) for x in t_theta];
  t_p_0=[getnormcdf(x,lowertail=True) for x in t_theta_0];
  if lower==True:
    return [(t_p[i]-t_p_0[i])/(1-t_p_0[i]) for i in range(n)];
  else:
    return [t_p[i]/(1-t_p_0[i]) for i in range(n)];


def getNBPValue(mean0,var0,mean1, lower=False,log=False):
  """
  Use negative binomial to calculate p-value
  Reference:
  http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html#scipy.stats.nbinom
  """
  from scipy.stats import nbinom
  n=len(mean0);
  nb_p=[mean0[i]/var0[i] for i in range(n)]; # consisitent with R
  nb_n0=[mean0[i]*mean0[i]/(var0[i]-mean0[i]) for i in range(n)];
  nb_n=[ (lambda t: t if t>=1 else 1)(x) for x in nb_n0];
  #
  if lower==True:
    if log==False:
      nb_p_low=nbinom.cdf(mean1,nb_n,nb_p);
    else:
      nb_p_low=nbinom.logcdf(mean1,nb_n,nb_p);
    return list(nb_p_low);
  else:
    if log==False:
      nb_p_low=nbinom.sf(mean1,nb_n,nb_p);
    else:
      nb_p_low=nbinom.logsf(mean1,nb_n,nb_p);
    return list(nb_p_low);
 



# calculate the FDR of p values using false discovery rate
#def pFDR(x):
#    index = [i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True)]
#    rindex = [i[0] for i in sorted(enumerate(index), key=lambda v:v[1])]
#    qval = sum([1.0/i for i in xrange(1,len(x)+1)])
#    l = [qval*len(x)/i*x[j] for i,j in zip(reversed(xrange(1,len(x)+1)),index)]
#    l = [l[k] if l[k] < 1.0 else 1.0 for k in rindex]
#    return l

def pFDR(pvalues, method='fdr'):
  """
  calculate the FDR of p values using false discovery rate
  check http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python
  """
  n=len(pvalues);
  new_pvalues=[0]*n;
  
  if method=='holm':
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
    values.sort()
    for rank, vals in enumerate(values):                                                              
      pvalue, i = vals
      new_pvalues[i] = (n-rank) * pvalue                  
      if new_pvalues[i] > 1.0:
        new_pvalues[i]=1.0;
  elif method=='fdr':
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
    values.sort()
    values.reverse()                                                                                  
    new_values = []
    for i, vals in enumerate(values):                                                                 
      rank = n - i
      pvalue, index = vals                                                                          
      new_values.append((n/rank) * pvalue)                                                          
    for i in xrange(0, int(n)-1):  
      if new_values[i] < new_values[i+1]:                                                           
        new_values[i+1] = new_values[i]                                                           
    for i, vals in enumerate(values):
      pvalue, index = vals
      new_pvalues[index] = new_values[i]    
  else:
    logging.error('Unknown FDR mthod '+method);
    return pvalues;
  return new_pvalues

def crispr_test(tab,ctrlg,testg, destfile,sgrna2genelist,args):
  """
  main function of crispr test
  """
  n=len(tab);
  # control and test matrix
  tabctrl={k:[v[i] for i in range(len(v)) if i in ctrlg] for (k,v) in tab.iteritems()};
  tabtest={k:[v[i] for i in range(len(v)) if i in testg] for (k,v) in tab.iteritems()};
  # control matrix for mean-var estimation
  if len(ctrlg)>1 and args.variance_from_all_samples==False: # more than 1 controls
    tabctrlmod={k:[v[i] for i in range(len(v)) if i in ctrlg] for (k,v) in tab.iteritems()};
  else: # only 1 control: use all the samples for estimation
    tabctrlmod={k:[v[i] for i in range(len(v)) if i in (ctrlg+testg)] for (k,v) in tab.iteritems()};
  # training using control samples
  model1=modelmeanvar(tabctrlmod,method='linear');
  #model2=modelmeanvar(tabctrl,method='edger');
  model=[x for x in model1];#+[model2];
  if type(model) is types.ListType:
    logging.debug('Adjusted model: '+'\t'.join([str(x) for x in model]));
  else:
    logging.debug('Adjusted model: k='+str(model));
  
  tabctrl_mat=tabctrl.values();
  tabctrlmodel_mat=tabctrlmod.values();
  tabc_mean=getMeans(tabctrl_mat);
  tabcmodel_mean=getMeans(tabctrlmodel_mat);
  #
  # setup the valid sgRNA flag
  validsgrna=[1]*n;
  if hasattr(args,"remove_zero") and ( args.remove_zero=="control" or args.remove_zero=="both"):
    validsgrna=[ (lambda x: 1 if x>0 else 0)(t) for t in tabc_mean]; 
  # if mean of the control samples is 0: set it to greater than 0
  tabc_min=min([x for x in tabc_mean if x>0]);
  tabc_mean=[ (lambda x: x if x>tabc_min else tabc_min)(t) for t in tabc_mean];
  tabc_var=getVars(tabctrlmodel_mat);
  tabc_adjvar=getadjustvar(model,tabc_mean,method='linear');
  
  # testing using tebtest
  nt=tabtest[tabtest.keys()[0]];
  ttmat=tabtest.values();
  ttmean=getMeans(ttmat);
  # set up the valid sgRNA flag
  if hasattr(args,"remove_zero") and ( args.remove_zero=="treatment" or args.remove_zero=="both"):
    validsgrna2=[ (lambda x: 1 if x>0 else 0)(t) for t in ttmean]; 
    validsgrna=[validsgrna[t]*validsgrna2[t] for t in range(n)];
  # use ttmean to calculate the pvalue
  # first, convert to standard normal distribution values
  tt_theta=[(ttmean[i]-tabc_mean[i])/math.sqrt(tabc_adjvar[i]) for i in range(n)];
  tt_abstheta=[math.fabs(tt_theta[i]) for i in range(n)];
  #
  try:
    tt_p_lower=getNBPValue(tabc_mean,tabc_adjvar,ttmean,lower=True);
    tt_p_higher=getNBPValue(tabc_mean,tabc_adjvar,ttmean,lower=False);
    # tt_p_lower_score=getNBPValue(tabc_mean,tabc_adjvar,ttmean,lower=True,log=True);
    # tt_p_higher_score=getNBPValue(tabc_mean,tabc_adjvar,ttmean,lower=False,log=True);
  except ImportError:
    logging.warning('An error occurs while trying to compute p values using scipy. Will use normal model instead of Negative Binomial model, but please check with your scipy installation.');
    tt_p_lower=getNormalPValue(tabc_mean,tabc_adjvar,ttmean,lower=True);
    tt_p_higher=getNormalPValue(tabc_mean,tabc_adjvar,ttmean,lower=False);
  except:
    logging.error('An error occurs while trying to compute p values. Quit..');
    sys.exit(-1);
  # 
  # lower_score and higher_score are used to sort sgRNAs
  tt_p_lower_score=tt_theta;
  tt_p_higher_score=[-1*x for x in tt_theta];
  #
  tt_p_twosided=[ (lambda x,y: 2*x if x<y else 2*y)(tt_p_lower[i],tt_p_higher[i]) for i in range(n)];
  tt_p_fdr=pFDR(tt_p_twosided,method=args.adjust_method);
  
  # 
  # sort according to the score
  sort_id=[i[0] for i in sorted(enumerate(tt_abstheta), key=lambda x:x[1],reverse=True)];
  # 
  # write to file
  destfname=destfile+'.sgrna_summary.txt';
  destf=open(destfname,'w');
  destkeys=tabctrl.keys();
  dfmt="{:.5g}"
  
  # output to file
  print('\t'.join(['sgrna','Gene','control_count','treatment_count','control_mean','treat_mean', 'control_var','adj_var','score','p.low','p.high','p.twosided','FDR','high_in_treatment']),file=destf);
  for i in sort_id:
    # sgRNA mapping to genes?
    if sgrna2genelist is not None:
      destkeygene=sgrna2genelist[destkeys[i]];
    else:
      destkeygene='None';
    report=[destkeys[i], destkeygene, '/'.join([dfmt.format(x) for x in tabctrl_mat[i]]), '/'.join([dfmt.format(x) for x in ttmat[i]])];
    t_r=[tabc_mean[i],ttmean[i]];
    t_r+=[tabc_var[i],tabc_adjvar[i],tt_abstheta[i],tt_p_lower[i],tt_p_higher[i],tt_p_twosided[i],tt_p_fdr[i]];
    report+=[dfmt.format(x) for x in t_r];
    report+=[ttmean[i]>tabc_mean[i]];
    print('\t'.join([str(x) for x in report]),file=destf);
  destf.close();
  #
  # prepare files for gene test
  if sgrna2genelist is not None:
    destfname=destfile+'.plow.txt';
    destf=open(destfname,'w');
    destkeys=tabctrl.keys();
    sort_id=[i[0] for i in sorted(enumerate(tt_p_lower_score), key=lambda x:x[1],reverse=False)];
    # output to file
    print('\t'.join(['sgrna','symbol','pool','p.low','prob','chosen']),file=destf);
    for i in sort_id:
      report=[destkeys[i], sgrna2genelist[destkeys[i]],'list', tt_p_lower_score[i], '1', validsgrna[i]];
      print('\t'.join([str(x) for x in report]),file=destf);
    destf.close();
    tt_p_lower_fdr=pFDR(tt_p_lower,method=args.adjust_method);
    n_lower=sum([1 for x in tt_p_lower if x <= args.gene_test_fdr_threshold]);
    n_lower_p=n_lower*1.0/len(tt_p_lower);
    logging.debug('lower test FDR cutoff: '+str(n_lower_p));
    #
    destfname=destfile+'.phigh.txt';
    destf=open(destfname,'w');
    destkeys=tabctrl.keys();
    sort_id=[i[0] for i in sorted(enumerate(tt_p_higher_score), key=lambda x:x[1],reverse=False)];
    # output to file
    print('\t'.join(['sgrna','symbol','pool','p.high','prob','chosen']),file=destf);
    for i in sort_id:
      report=[destkeys[i], sgrna2genelist[destkeys[i]],'list', tt_p_higher_score[i], '1', validsgrna[i]];
      print('\t'.join([str(x) for x in report]),file=destf);
    destf.close();
    tt_p_higher_fdr=pFDR(tt_p_higher,method=args.adjust_method);
    n_higher=sum([1 for x in tt_p_higher if x <= args.gene_test_fdr_threshold]);
    n_higher_p=n_higher*1.0/len(tt_p_higher);
    logging.debug('higher test FDR cutoff: '+str(n_higher_p));
    # 
    return (n_lower_p,n_higher_p);
  else:
    return (None,None);

def rank_association_test(file,outfile,cutoff,args):
  if cutoff<0.05:
    cutoff=0.05;
  if cutoff>0.5:
    cutoff=0.5;
  #rrapath='/'.join(sys.argv[0].split('/')[:-1]+["../bin/RRA"])
  rrapath='RRA';
  command=rrapath+" -i "+file+" -o "+outfile+" -p "+str(cutoff);
  if hasattr(args,'control_sgrna') and args.control_sgrna != None :
    command+=" --control "+args.control_sgrna;
  systemcall(command);
  

def magecktest_removetmp(prefix):
  tmpfile=[prefix+'.plow.txt',prefix+'.phigh.txt',prefix+'.gene.low.txt',prefix+'.gene.high.txt'];
  for f in tmpfile:
    systemcall('rm '+f,cmsg=False);



  
def magecktest_main(args):
  """
  Main entry for MAGeCK test function
  """
  # stat test
  if args.subcmd == 'run' or args.subcmd == 'test':
    # read counts from file
    if args.subcmd == 'test':
      mapres=getcounttablefromfile(args.count_table);
    else:
      mapres=getcounttablefromfile(args.output_prefix+'.count.txt');
    cttab=mapres[0];
    sgrna2genelist=mapres[1];
    samplelabelindex=mapres[2];
    
    if len(cttab)==0:
      sys.exit(-1);
    nsample=len(cttab[cttab.keys()[0]]);
    
    # iterate control group and treatment group
    supergroup_control=args.control_id;
    supergroup_treat=args.treatment_id;
    # control group and treatment group labels
    labellist_control=[];
    labellist_treat=[];
    # R visualization init
    vrv=VisualRValue();
    vrv.outprefix=args.output_prefix;
    vrv.genesummaryfile=args.output_prefix+'.gene_summary.txt';
    vrv.startRTemplate();
    vrvrnwcplabel=[]; # labels to write in rnw
    
    # loop by comparisons
    for cpindex in range(len(supergroup_treat)):
      # convert the sample label to sample index
      if cpindex==0:
        cp_prefix=args.output_prefix;
      else:
        cp_prefix=args.output_prefix+'.'+str(cpindex);
      # labels
      (treatgroup,treatgrouplabellist)=parse_sampleids(supergroup_treat[cpindex],samplelabelindex);
      treatgroup_label=str(supergroup_treat[cpindex]);
      logging.info('Treatment samples:'+treatgroup_label);
      labellist_treat+=[treatgroup_label];
      if supergroup_control != None:
        (controlgroup,controlgrouplabellist)=parse_sampleids(supergroup_control[cpindex],samplelabelindex); 
        controlgroup_label=str(supergroup_control[cpindex]); # only for display
        logging.info('Control samples:'+controlgroup_label);
      else:
        #controlgroup=[x for x in range(nsample) if x not in treatgroup];
        #controlgrouplabellist=[samplelabelindex[x] for x in range(nsample) if x not in treatgroup];
        xls=[x for x in range(nsample) if x not in treatgroup];
        (controlgroup,controlgrouplabellist)=parse_sampleids(','.join([str(t) for t in xls]),samplelabelindex); 
        controlgroup_label='rest';
        logging.info('Control samples: the rest of the samples');
      labellist_control+=[controlgroup_label];
      # read the sgRNA-gene table for rank association
      # normalization
      cttab_sel={k:([v[i] for i in controlgroup + treatgroup]) for (k,v) in cttab.iteritems()}; # controlgroup do not overlap with treatgroup  
      if hasattr(args,'norm_method'):
        nttab=normalizeCounts(cttab_sel,method=args.norm_method);
      else:
        nttab=normalizeCounts(cttab_sel);
      # write normalized counts to file
      if hasattr(args,'normcounts_to_file'):
        if args.normcounts_to_file:
          # counts
          mageck_printdict(nttab,args,sgrna2genelist,samplelabelindex,controlgroup+treatgroup);
      
      controlgroup_ids=list(range(len(controlgroup)));
      treatgroup_ids=list(range(len(controlgroup),len(controlgroup+treatgroup)));
      # perform sgRNA test, and prepare files for gene test
      gene_as_cutoff=crispr_test(nttab, controlgroup_ids, treatgroup_ids, cp_prefix,sgrna2genelist,args);  
      #
      if gene_as_cutoff[0] is not None:
        rank_association_test(cp_prefix+'.plow.txt',cp_prefix+'.gene.low.txt',gene_as_cutoff[0],args);
      if gene_as_cutoff[1] is not None:
        rank_association_test(cp_prefix+'.phigh.txt',cp_prefix+'.gene.high.txt',gene_as_cutoff[1],args);
      # merge different files
      merge_rank_files(cp_prefix+'.gene.low.txt',cp_prefix+'.gene.high.txt',cp_prefix+'.gene_summary.txt',args);
      if cpindex>0:
        if cpindex>1:
          label1='';
        else:
          if len(labellist_treat)>0:
            label1=labellist_treat[0]+'_vs_'+labellist_control[0]+'.';
          else:
            label1='';
        label2=treatgroup_label+'_vs_'+controlgroup_label+'.';
        merge_rank_summary_files(args.output_prefix+'.gene_summary.txt',cp_prefix+'.gene_summary.txt',args.output_prefix+'.gene_summary.txt',args,lowfile_prefix=label1,highfile_prefix=label2);
      # visualization: load top k genes
      # print(str(samplelabelindex));
      vrv.cplabel=treatgroup_label+'_vs_'+controlgroup_label+' neg.';
      vrvrnwcplabel+=[vrv.cplabel];
      vrv.cpindex=[2+10*cpindex+1];
      vrv.loadTopKWithExp(cp_prefix+'.gene.low.txt',nttab,sgrna2genelist,controlgrouplabellist+treatgrouplabellist);
      vrv.cplabel=treatgroup_label+'_vs_'+controlgroup_label+' pos.';
      vrvrnwcplabel+=[vrv.cplabel];
      vrv.cpindex=[2+10*cpindex+5+1];
      vrv.loadTopKWithExp(cp_prefix+'.gene.high.txt',nttab,sgrna2genelist,controlgrouplabellist+treatgrouplabellist);
      
      # clean the file
      if args.keep_tmp==False:
        magecktest_removetmp(cp_prefix);
        if cpindex>0:
          systemcall('rm '+cp_prefix+'.gene_summary.txt',cmsg=False);
          systemcall('rm '+cp_prefix+'.sgrna_summary.txt',cmsg=False);
      # end cleaning
    # end cpindex loop
    
    # generate pdf file
    # write to rnw file buffer
    vrv.genesummaryfile=args.output_prefix+'.gene_summary.txt';
    vrv.getGeneSummaryStat(isplot=False);
    vrv.comparisonlabel=vrvrnwcplabel; # replace the label field
    vrv.writeGeneSummaryStatToBuffer();
    # write to rnw and R file
    vrv.closeRTemplate();
    if hasattr(args, "pdf_report") and args.pdf_report:
      vrv.generatePDF(args.keep_tmp);
  # end if        
 


