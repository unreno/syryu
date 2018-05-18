#	DropSeq Processing

##	20180518

New sample prep.

Downloaded Project MK1779 fastq.gz files from Illumina for "Batch6" and "Batch7" (B6 and B7)

Copied all Batch6\*fastq.gz into B6 directory and ran convert\_fastq\_files\_to\_merged\_sorted\_bams.bash.

Copied all Batch7\*fastq.gz into B7 directory and ran convert\_fastq\_files\_to\_merged\_sorted\_bams.bash.





-r--r--r--  1 jakewendt 7934561487 May 18 16:11 B6.bam
-r--r--r--  1 jakewendt 4214557648 May 18 16:01 B7.bam

B6.bam is almost 8GB which is the size of 1, 2, 3 and 4 combined.
This may not be enough, but we shall see.

Restart Azure Linux Standard_E64-16s_v3 (16 vcpus, 432 GB memory)





###	UPDATE!


```BASH
ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@IP_ADDRESS
sudo apt update
sudo apt full-upgrade
sudo apt autoremove
sudo apt install gfortran

cd ~/syryu
git pull
make install


mkdir ~/working/dropseq
cd ~/working/dropseq
mkdir B6
mkdir B7
```


locally ...

```BASH
scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/20180518/*.bam jake@IP_ADDRESS:working/

```


Remotely ...
```BASH
ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@IP_ADDRESS
cd ~/working/dropseq
nohup drop_seq.bash ~/working/*.bam > drop_seq.log 2>&1 &
```





###	DOWNLOAD LOCALLY

```BASH
rsync --archive --verbose --compress --rsh "ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete jake@IP_ADDRESS:working/dropseq/ ~/github/unreno/syryu/drop_seq/20180518.drop_seq_alignment/
```


###	UPLOAD TO AZURE STORAGE

Browse to portal.azure.com, Storage Account -> ryulab -> Access Keys to find a key.

I copied mine into ~/dest-key

I tried the <( cat ~/dest-key ), but it errored

[ERROR] The syntax of the command is incorrect. The supplied storage key (dest-key) is not a valid Base64 string.


Cleanup and upload data to Azure Storage and prep to save VM image ...

Remotely ...

```BASH
mv ~/working/dropseq ~/working/20180518.drop_seq_alignment

azcopy --source ~/working/B6.bam --destination https://ryulab.blob.core.windows.net/ryulab/bams/ --dest-key YOUR_DEST_KEY
azcopy --source ~/working/B7.bam --destination https://ryulab.blob.core.windows.net/ryulab/bams/ --dest-key YOUR_DEST_KEY

azcopy --source ~/working/20180518.drop_seq_alignment --destination https://ryulab.blob.core.windows.net/ryulab/DropSeq/20180518.drop_seq_alignment/ --recursive --dest-key YOUR_DEST_KEY

sudo waagent -deprovision
```

Using the web portal GUI, save the image




##################################################


##	Post Run Notes

Probably should've stuck with a single name for the modified reference. Say mm10ryulab.


