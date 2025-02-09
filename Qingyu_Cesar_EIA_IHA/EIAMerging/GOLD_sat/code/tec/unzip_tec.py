#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:24:58 2024

@author: qingyuzhu


Unzip tec file
"""

import gzip
import shutil
import os 
from glob import glob

fpath='/Volumes/Io2/proposal/2025/PostSunsetEIA/data/tec/temp/'
flist=sorted(glob(fpath+'*.gz'))

for ifile, fname in enumerate(flist[:]):
    
    print (fname)
    
    
    
    with gzip.open(fname,'rb') as f_in:
        
        with open(fname[:-11]+'.txt','wb') as f_out:
            
            shutil.copyfileobj(f_in, f_out)
            
            
            
    os.remove(fname)
        
        
    
