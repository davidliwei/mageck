#!/usr/bin/env python
'''
MAGeCK set up script
'''


from __future__ import print_function;

import os
import sys
from distutils.core import setup, Extension
from subprocess import call as subpcall
from distutils.command.install import install as DistutilsInstall

def compile_rra():
  # 
  os.chdir('rra');
  subpcall('make',shell=True);
  rev=subpcall('../bin/RRA',shell=True);
  os.chdir('../');
  return rev;


class RRAInstall(DistutilsInstall):
  def run(self):
    # compile RRA
    if(compile_rra()!=0):
      print("CRITICAL: error compiling the RRA source code. Please check your c compilation environment.",file=sys.stderr);
      sys.exit(1);
    DistutilsInstall.run(self)



def main():
  # check python version
  if float(sys.version[:3])<2.7 or float(sys.version[:3])>=2.8:
    sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
    sys.exit(1);
  
  setup(name='mageck',
    version='0.5.0',
    description='Model-based Analysis of Genome-wide CRISPR-Cas9 Knockout',
    author='Wei Li, Han Xu',
    author_email='li.david.wei@gmail.com',
    url='http://mageck.sourceforge.net',
    packages=['mageck'],
    scripts=['bin/RRA','bin/mageck'],
    package_dir={'mageck':'mageck'},
    cmdclass={'install':RRAInstall},
    package_data={'mageck':['*.Rnw','*.RTemplate']}
    #package_data={'mageck':['mageck/Makefile','mageck/src/*.c','include/*','utils/*']}
    #data_files=[('',['Makefile','src/*.c','include/*','utils/*'])]
  );
    
    
if __name__ == '__main__':
  main();
    

