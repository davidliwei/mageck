#!/usr/bin/env python
""" MAGeCK count module
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

import sys;
import argparse;
import logging

def mageckcount_parseargs():
  """
  Parse arguments. Only used when mageckCount.py is executed directly.
  """
  parser=argparse.ArgumentParser(description='Collecting read counts for multiple samples.');
  
  parser.add_argument('-l','--list-seq',required=True,help='A file containing the list of sgRNA names, their sequences and associated genes. Support file format: csv and txt.');
  parser.add_argument('--sample-label',default='',help='Sample labels, separated by comma (,). Must be equal to the number of samples provided. Default "sample1,sample2,...".');
  parser.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.');
  parser.add_argument('--trim-5',type=int,default=0,help='Length of trimming the 5\' of the reads. Default 0');
  parser.add_argument('--sgrna-len',type=int,default=20,help='Length of the sgRNA. Default 20');
  parser.add_argument('--count-n',action='store_true',help='Count sgRNAs with Ns. By default, sgRNAs containing N will be discarded.');
  parser.add_argument('--fastq',nargs='+',help='Sample fastq files, separated by space; use comma (,) to indicate technical replicates of the same sample. For example, "--fastq sample1_replicate1.fastq,sample1_replicate2.fastq sample2_replicate1.fastq,sample2_replicate2.fastq" indicates two samples with 2 technical replicates for each sample.');

  
  args=parser.parse_args();
  
  
  return args;

def mageckcount_checkargs(args):
  """
  Check args
  """
  if args.sample_label!='':
    nlabel=args.sample_label.split(',');
    #nfq=args.fastq.split(',');
    nfq=(args.fastq);
    if len(nlabel)!=len(nfq):
      logging.error('The number of labels ('+str(nlabel)+') must be equal to the number of fastq files provided.');
      sys.exit(-1);
  return 0;

def normalizeCounts(ctable,method='median'):
  """
  Normalize read counts
  Return value: {sgRNA:[read counts]} 
  """
  # sums
  if len(ctable)==0:
    return ctable.copy();
  n=len(ctable[ctable.keys()[0]]); # samples
  m=len(ctable); # sgRNAs
  # calculate the sum
  sumsample=[0]*n;
  for (k,v) in ctable.iteritems():
    sumsample=[sumsample[i]+v[i] for i in range(n)];
  logging.info('Total read counts of each sample: '+' '.join([str(x) for x in sumsample]));
  logging.debug('Normalization method: '+method);
  # normalizing factor
  avgsample=sum(sumsample)/float(n);
  samplefactor=[avgsample/k for k in sumsample];
  logging.debug('Initial (total) size factor: '+' '.join([str(x) for x in samplefactor]));
  if method=='median':
    # calculate the average
    meanval={k:(sum(v)*1.0/n) for (k,v) in ctable.iteritems() if sum(v)>0};
    meanval={k:(lambda x: x if x>0 else 1)(v) for (k,v) in meanval.iteritems()};
    #samplefactor=[0]*n;
    usetotalnorm=False;
    medianfactor=[x for x in samplefactor];
    for ni in range(n):
      meanfactor=[ v[ni]/meanval[k] for (k,v) in ctable.iteritems() if k in meanval];
      #print(str(sorted(meanfactor)))
      xfactor=sorted(meanfactor)[m//2];
      if xfactor>0.0:
        medianfactor[ni]=1.0/xfactor;
        #logging.debug('xfactor:'+str(xfactor));
      else:
        logging.warning('Sample '+str(ni)+' has zero median count, so median normalization is not possible. Switch to total read count normalization.');
        usetotalnorm=True;
    # normalizing factor
    if usetotalnorm:
      pass;
    else:
      samplefactor=medianfactor;
      logging.debug('Median factor: '+' '.join([str(x) for x in samplefactor]));
  elif method=='none': 
    samplefactor=[1]*n;
  logging.debug('Final factor: '+' '.join([str(x) for x in samplefactor]));
    
  # normalize the table
  ntable={ k: [ samplefactor[i]*v[i] for i in range(n)] for (k,v) in ctable.iteritems()};
  
  return ntable;
  


def mageckcount_processonefile(filename,args,ctab,genedict,datastat):
  '''
  Go through one fastq file
  Parameters
  ----------
  filename
    Fastq filename to be sequence
  args
    Arguments
  ctab
    A dictionary of sgRNA sequence and count
  genedict
    {sequence:(sgRNA_id,gene_id)} dictionary
  datastat
    Statistics of datasets ({key:value})

  Return value
  ----------
  '''
  # ctab={};
  nline=0;
  logging.info('Parsing file '+filename+'...');
  nreadcount=0;
  for line in open(filename):
    nline=nline+1;
    if nline%1000000==1:
      logging.info('Processing '+str(round(nline/1000000))+ 'M lines..');
    if nline%4 == 2:
      nreadcount+=1;
      fseq=line.strip();
      if args.trim_5 >0:
        fseq=fseq[args.trim_5:];
      
      if len(fseq)<args.sgrna_len:
        continue;
      fseq=fseq[:args.sgrna_len];
      if fseq.count('N')>0 and args.count_n==False:
        continue;
      if fseq not in ctab:
        ctab[fseq]=0;
      ctab[fseq]=ctab[fseq]+1;
  # statistics
  datastat['reads']=nreadcount;
  # check if a library is provided
  if len(genedict)==0:
    datastat['mappedreads']=0;
    datastat['zerosgrnas']=0;
  else:
    nmapped=0;
    for (k,v) in ctab.iteritems():
      if k in genedict:
        nmapped+=v;
    nzerosg=0;
    for (k,v) in genedict.iteritems():
      if k not in ctab:
        nzerosg+=1;
    logging.info('mapped:'+str(nmapped));
    datastat['mappedreads']=nmapped;
    datastat['zerosgrnas']=nzerosg;
  #return ctab;
  return 0;

def mageckcount_mergedict(dict0,dict1):
  '''
  Merge all items in dict1 to dict0.
  '''
  nsample=0;
  if len(dict0)>0:
    nsample=len(dict0[dict0.keys()[0]]);
  for (k,v) in dict0.iteritems():
    if k in dict1:
      v+=[dict1[k]];
    else:
      v+=[0];
  for (k,v) in dict1.iteritems():
    if k not in dict0:
      if nsample>0:
        dict0[k]=[0]*nsample;
      else:
        dict0[k]=[];
      dict0[k]+=[v];
  # return dict0;

def mageckcount_printdict(dict0,args,ofile,sgdict,datastat,sep='\t'):
  '''
  Write the table count to file
  '''
  allfastq=args.fastq;
  nsample=len(allfastq);
  slabel=[datastat[f.split(',')[0]]['label'] for f in allfastq];
  # print header
  print('sgRNA'+sep+'Gene'+sep+sep.join(slabel),file=ofile);
  # print items
  if len(sgdict)==0:
    for (k,v) in dict0.iteritems():
      print(k+sep+'None'+sep+sep.join([str(x) for x in v]),file=ofile);
  else:
    for (k,v) in dict0.iteritems():
      if k not in sgdict: # only print those in the genedict
        continue;
      sx=sgdict[k];
      print(sep.join([sx[0],sx[1]])+sep+sep.join([str(x) for x in v]),file=ofile);
    # print the remaining counts, fill with 0
    for (k,v) in sgdict.iteritems():
      if k not in dict0:
        print(sep.join([v[0],v[1]])+sep+sep.join(["0"]*nsample),file=ofile);

def mageck_printdict(dict0,args,sgdict,sampledict,sampleids):
  """Write the normalized read counts to file
  
  Parameters
  ----------
  dict0 : dict
    a {sgRNA: [read counts]} structure
  args : class
    a argparse class
  sgdict: dict
    a {sgrna:gene} dictionary
  sampledict: dict
    a {sample name: index} dict
  sampleids: list
    a list of sample index. Should include control+treatment
  
  """
  # print header
  # print items
  dfmt="{:.5g}"
  ofile=open(args.output_prefix+'.normalized.txt','w');
  # headers
  mapres_list=['']*len(sampledict);
  for (k,v) in sampledict.iteritems():
    mapres_list[v]=k;
  if len(sampledict)>0:
    cntheader=[mapres_list[x] for x in sampleids]
  else:
    cntheader=None;
  logging.info('Writing normalized read counts to '+args.output_prefix+'.normalized.txt');
  if cntheader !=None:
    print('sgRNA\tGene\t'+'\t'.join(cntheader),file=ofile);
  if len(sgdict)==0:
    for (k,v) in dict0.iteritems():
      print(k+'\t'+'None'+'\t'+'\t'.join([str(x) for x in v]),file=ofile);
  else:
    for (k,v) in dict0.iteritems():
      if k not in sgdict: # only print those in the genedict
        logging.warning(k+' not in the sgRNA list');
        continue;
      print('\t'.join([k,sgdict[k]])+'\t'+'\t'.join([str(x) for x in v]),file=ofile);
    # print the remaining counts, fill with 0
  ofile.close();




def mageckcount_checklists(args):
  """
  Read sgRNAs and associated sequences and lists
  format: sgRNAid  seq  geneid
  """
  genedict={};
  hascsv=False;
  if args.list_seq.upper().endswith('CSV'):
    hascsv=True;
  n=0;
  seqdict={};
  for line in open(args.list_seq):
    if hascsv:
      field=line.strip().split(',');
    else:
      field=line.strip().split();
    n+=1;
    if field[0] in genedict:
      logging.warning('Duplicated sgRNA label '+field[0]+' in line '+str(n)+'. Skip this record.');
      continue;
    if len(field)<3:
      logging.warning('Not enough field in line '+str(n)+'. Skip this record.');
      continue;
    if field[1].upper() in seqdict:
      logging.warning('Duplicated sgRNA sequence '+field[1]+' in line '+str(n)+'. Skip this record.');
      continue;
    genedict[field[0]]=(field[1].upper(),field[2]);
  logging.info('Loading '+str(len(genedict))+' predefined sgRNAs.');
  return genedict;

def mageckcount_printstat(args,datastat):
  for (k,v) in datastat.iteritems():
    logging.info('Summary of file '+k+':');
    for (v1,v2) in v.iteritems():
      logging.info(str(v1)+'\t'+str(v2));

def mageckcount_main(args):
  """
  Main entry for mageck count module
  """
  # check arguments
  mageckcount_checkargs(args);
  # check the listed files
  # listfq=args.fastq.split(',');
  listfq=[[z for z in x.split(',')] for x in args.fastq];
  nsample=len(listfq);
  datastat={};
  # check labels
  alllabel=args.sample_label;
  if alllabel=='':
    slabel=['sample'+str(x) for x in range(1,nsample+1)];
  else:
    slabel=alllabel.split(',');
  for i in range(nsample):
    for fi in listfq[i]:
      datastat[fi]={};
      datastat[fi]['label']=slabel[i];
  # process gene dicts
  genedict={};
  if args.list_seq is not None:
    genedict=mageckcount_checklists(args);
  # save sgRNA ID and gene name
  sgdict={}; #
  for (k,v) in genedict.iteritems():
    sgdict[v[0]]=(k,v[1]);
  alldict={};
  # go through the fastq files
  for filenamelist in listfq:
    dict0={};
    for filename in filenamelist:
      mageckcount_processonefile(filename,args,dict0,sgdict,datastat[filename]);
    mageckcount_mergedict(alldict,dict0);
  # write to file
  ofilel=open(args.output_prefix+'.count.txt','w');
  mageckcount_printdict(alldict,args,ofilel,sgdict,datastat);
  ofilel.close();
  # write the median normalized read counts to csv file
  ofilel=open(args.output_prefix+'.count.median_normalized.csv','w');
  if len(sgdict)>0:
    allmappeddict={k:v for (k,v) in alldict.iteritems() if k in sgdict}; # only keep those with known sgRNAs
  else:
    allmappeddict=alldict;
  medalldict=normalizeCounts(allmappeddict);
  mageckcount_printdict(medalldict,args,ofilel,sgdict,datastat,sep=',');
  ofilel.close();
  # print statistics
  mageckcount_printstat(args,datastat);
  return 0;


def getcounttablefromfile(filename):
  """
  read count table from file
  Returns:
  ---------------
  x: dict
    {sgrna:[read counts]} 
  y: dict
    {sgrna:gene}
  z: dict
    z={sample_id:index}
  """
  gtab={};
  mapptab={};
  sampleids={};
  nline=0;
  nfield=-1;
  # if it is CSV file
  hascsv=False;
  if filename.upper().endswith('.CSV'):
    hascsv=True;
  logging.info('Loading count table from '+filename+' ');
  for line in open(filename):
    nline+=1;
    if nline % 100000 == 1:
      logging.info('Processing '+str(nline)+' lines..');
    try:
      if hascsv==False:
        field=line.strip().split();
      else:
        field=line.strip().split(',');
      sgid=field[0];
      geneid=field[1];
      # check if duplicate sgRNA IDs are detected
      if sgid in gtab:
        logging.warning('Duplicated sgRNA IDs: '+sgid+' in line '+str(nline)+'. Skip this record.');
        continue;
      mapptab[sgid]=geneid;
      sgrecs=[float(x) for x in field[2:]];
      # check the number of fields
      if nfield!=-1 and len(sgrecs)!=nfield:
        logging.error('Error: incorrect number of dimensions in line '+str(nline)+'. Please double-check your read count table file.');
        sys.exit(-1);
      nfield=len(sgrecs);
      gtab[sgid]=sgrecs;
    except ValueError:
      if nline!=1:
        logging.warning('Parsing error in line '+str(nline)+'. Skip this line.');
      else:
        logging.debug('Parsing error in line '+str(nline)+' (usually the header line). Skip this line.');
        ids=field[2:];
        for i in range(len(ids)):
          sampleids[ids[i]]=i;
      continue;
  logging.info('Loaded '+str(len(gtab))+' records.');
  return (gtab,mapptab,sampleids);



if __name__ == '__main__':
  try:
    args=mageckcount_parseargs();
    mageckcount_main(args);
  except KeyboardInterrupt:
    sys.stderr.write("Interrupted.\n")
    sys.exit(0)



#### These are old functions
# def getsgrnalist(filename):
#   sgrnalist=[];
#   if filename.upper().endswith('.FA') or filename.upper().endswith('.FASTA'):
#     for line in open(filename):
#       if line[0]=='>':
#         tf=line.strip();
#         tf=tf[1:];
#         sgrnalist+=[tf];
#   else:
#     for line in open(filename):
#       if line[0]=='#':
#         continue;
#       field=line.strip().split();
#       sgrnalist+=[field[0]];
#   logging.info('Loading '+str(len(sgrnalist))+' sgRNAs');
#   return sgrnalist;
# 
# def getsgrna2genelist(filename):
#   sgrna2gene={};
#   n=0;
#   for line in open(filename):
#     n=n+1;
#     if line[0]=='#':
#       continue;
#     field=line.strip().split();
#     if len(field)<2:
#       logging.warning('There are less than 2 fields in the sgRNA to gene mapping file (line '+str(n)+'). Skip the rank association test ...');
#       return None;
#     if field[0] in sgrna2gene:
#       logging.warning('Duplicated sgRNA name of '+field[0]+' in line '+str(n)+', skip this record..');
#       continue;
#     sgrna2gene[field[0]]=field[1];
#   return sgrna2gene;
# 
# 
# # retrive count tables
# def getcounttable(sgrnalist,samfiles,denovo=False):
#   sgrnatab={};
#   ns=len(samfiles);
#   for s in sgrnalist:
#     sgrnatab[s]=[0]*len(samfiles);
#   for si in range(ns):
#     samfilename=samfiles[si];
#     logging.info('Processing '+samfilename+' ... ');
#     nline=0;
#     for line in open(samfilename):
#       if line[0]=='@':
#         continue;
#       nline+=1;
#       if nline % 100000 ==1:
#         logging.info('Processing '+str(nline)+' lines..');
#       field=line.strip().split();
#       sid=field[2];
#       if denovo:
#         if sid not in sgrnatab:
#           sgrnatab[sid]=[0]*ns;
#         sgrnatab[sid][si]+=1;
#       else:
#         # pass;
#         if sid in sgrnatab:
#           sgrnatab[sid][si]+=1;
#   return sgrnatab;
# 
 
# # write count table
# def writecounttable(table,filename):
#   logging.info('writing to '+filename+' ...');
#   ofile=open(filename,'w');
#   for (ck, cv) in table.iteritems():
#     print(ck+'\t'+'\t'.join([str(x) for x in cv]),file=ofile);
#   ofile.close();
# 
# 
#### END of old functions

