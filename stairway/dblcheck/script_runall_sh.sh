#!/bin/bash

cat all_sh_files.txt | while read line 

do 
bash ${line} 

done 

