#	MSConvert and MODa

##	20180706


1. Connect to ftp://massive.ucsd.edu/MSV000079053. There are over 100 bacteria in this folder. 
2. Download fasta files and raw files for each bacteria. 
    a. fasta files: Under sequence folder, there are many folders (one folder for one bacteria). In bacteria folder, download files ending with .fasta.
    b. RAW files: Under raw folder, there are many folders (one folder for one bacteria). In bacteria folder, download files ending with RAW.  
3. Convert RAW to mgf.
    a. Download ProteoWizard. http://proteowizard.sourceforge.net/download.html
        May have to download window version to convert raw files. 
    b. Convert RAW to mgf using msconvert command
    http://proteowizard.sourceforge.net/tools/msconvert.html
    ex) msconvert data.RAW --mgf --filter "msLevel 2" --filter "zeroSample removeExtra"
4. Run MODa. I attached an example parameter file (foo.txt) and a tutorial (README.pdf). For each mgf file, please change "Spectra=" and "Fasta=" and run MODa. 
5. Transfer results files to Azure storage. 

Please let me know if you have questions!
Thank you, Jake.





Seems msconvert.exe needs MSFileReader to read RAW files.

Still not working? Adding XcalDLL

convertion seems to be working.





C:\Users\jake>\Ruby25-x64\bin\ruby.exe msconvert.rb
















------




```

cp foo.txt base_config.txt


DELETE all of the Spectra


for i in *mgf; do
o=${i%.mgf}
echo $i $o
echo -e "Spectra=$i\r\r" > $o.config
cat base_config.txt >> $o.config
done




I created a .config file for each .mgf file, then I started a for loop processing each .config file producing a .config.out file.

I added a README file in both space/ and control/ describing what I did / am doing.


I copied moda_v1.51.jar file to the SpacePlant folder as just moda.jar. For some reason, using the full name wasn't working for me. I also wrote a little script called moda.bash


#!/usr/bin/env bash

echo "Running"

for i in ${1}*config; do
o=${i%.mgf}
echo $i $o
java -Xmx5000M -jar ../moda.jar -i $i -o $o.out
done


#!/usr/bin/env bash

echo "Running"

for i in ${1}*config; do
o=${i%.mgf}
echo $i $o
java -Xmx5000M -jar ../moda.jar -i $i -o $o.out
done
```
