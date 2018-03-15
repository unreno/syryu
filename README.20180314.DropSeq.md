#	DropSeq Processing

##	20180314

Restart Azure Linux Standard_E64-16s_v3 (16 vcpus, 432 GB memory)

Recreate empty ~/working/dropseq

Copy in error_detected.bam files to 1, 2, B3 and B4

find . -name error_detected.bam -execdir dge.bash \; > dge.log 2>&1



###	UPDATE!


```BASH
ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@52.226.33.193
sudo apt update
sudo apt full-upgrade
sudo apt autoremove

cd ~/syryu
git pull
make install


cd ~/working/dropseq
mkdir 1
mkdir 2
mkdir B3
mkdir B4
```



Locally ...

```BASH
scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/20180305.drop_seq_alignment/1/error_detected.bam jake@52.226.33.193:working/dropseq/1/

scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/20180305.drop_seq_alignment/2/error_detected.bam jake@52.226.33.193:working/dropseq/2/

scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/20180305.drop_seq_alignment/B3/error_detected.bam jake@52.226.33.193:working/dropseq/B3/

scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/20180305.drop_seq_alignment/B4/error_detected.bam jake@52.226.33.193:working/dropseq/B4/
```

Should've prepped a single directory locally and uploaded that. TomAto, TomAHto, though.


Remotely ...

```BASH
cd ~/working/dropseq

find . -name error_detected.bam -execdir pwd \;
#/home/jake/working/dropseq/B4
#/home/jake/working/dropseq/B3
#/home/jake/working/dropseq/2
#/home/jake/working/dropseq/1

nohup find . -name error_detected.bam -execdir dge.bash \; > dge.log 2>&1 &
```



###	DOWNLOAD

```BASH
rsync --archive --verbose --compress --rsh "ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete jake@52.226.33.193:working/dropseq/ ~/github/unreno/syryu/drop_seq/20180314.drop_seq_alignment/
```


##################################################


##	Post Run Notes

Probably should've stuck with a single name for the modified reference. Say mm10ryulab.



It seems that R doesn't clean up well.
Even after deleting object, it still takes up quite a bit of memory.
While the memory used doesn't show in mem\_used(), the system shows it as mostly taken.
Later in the script, another function can run out of memory and the script crashes out.
Rerunning seurat.R using my option --redo, loads the data from the stored data and never loads the dge file theu never fills the memory.
Success.
Perhaps future runs should have one script simply convert the dge to seurat and then quit leaving the analysis to load the saved seurat object and analyze.


