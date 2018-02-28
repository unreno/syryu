#	DropSeq Processing

##	20180228


###	Preprocess locally ...

* Download mm10c.fasta.gz, mm10c.refFlat.gz, mm10c.dict.gz, mm10c.gtf.gz from Azure
* Download picard.jar `wget https://github.com/broadinstitute/picard/releases/download/2.17.0/picard.jar`

* Download fastq data from Illumina
* Modify convert\_fastq\_files\_to\_bams.bash to be more flexible or change to convert\_fastq\_files\_to\_merged\_sorted\_bams.bash <----
* Convert fastq data to bams
* Merge bams and sort by name with picard so can correctly incorrectly sort

```
java -jar ~/picard.jar FastqToSam \
  F1=$fastq1 \
  F2=$fastq2 \
  O=$basename.bam \
  SM=$basename

java -jar ~/picard.jar MergeSamFiles \
	INPUT=1A_S4_L001.bam \
	INPUT=1A_S4_L002.bam \
	INPUT=1A_S4_L003.bam \
	INPUT=1A_S4_L004.bam \
	INPUT=1B_S3_L001.bam \
	INPUT=1B_S3_L002.bam \
	INPUT=1B_S3_L003.bam \
	INPUT=1B_S3_L004.bam \
	ASSUME_SORTED=true \
	SORT_ORDER=queryname \
	OUTPUT=1.bam
```


###	Start new AWS instance

```
create_ec2_instance.bash --profile syryu --key ~/.aws/JakeSYRyu.pem --instance-type x1e.xlarge --volume-size 100

ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
```

###	UPDATE!

```
sudo yum update
cd ~/syryu
git pull
make install
```

###	Destroy mm10/mm10a stuff and prep for mm10c

```
chmod -R +w ~/mm10a
/bin/rm -rf ~/mm10a
mkdir ~/mm10c
```

###	UPLOAD BAM FILES AND NEW FILES FOR MAKING REFERENCE

```
ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip




scp -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/20180227/?.bam ec2-user@$ip:working/





scp -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/mm10c.*.gz ec2-user@$ip:mm10c/
```


```
cd ~/mm10c/
gunzip mm10c.fasta.gz
gunzip mm10c.gtf.gz
gunzip mm10c.refFlat.gz
gunzip mm10c.dict.gz
```

###	create STAR reference for mm10c with just fasta and gtf (takes about 40-90 minutes depending on core count)


```
chmod -w ~/mm10c/mm10c*
cd ~/working/
mkdir -p ~/working/mm10c_star

gunzip ~/mm10c/mm10c.fasta.gz

nohup STAR --genomeFastaFiles ~/mm10c/mm10c.fasta --runMode genomeGenerate --genomeDir ~/working/mm10c_star --sjdbGTFfile ~/mm10c/mm10c.gtf --sjdbOverhang 100 --runThreadN 4 &
```

```
mkdir -p ~/working/dropseq
mv ~/working/Log.out ~/working/dropseq/star_mm10c_creation.log
chmod 444 ~/working/mm10c_star/*
chmod 444 ~/working/dropseq/star_mm10c_creation.log
rmdir ~/working/_STARtmp
```

###	RUN DATA (takes about 8 hours)

Drop-seq\_alignment.sh takes about an hour, then dge.bash/seurat.R takes about 2 hours or so.

```
ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@$ip
cd ~/working/dropseq
nohup drop_seq.bash ~/working/?.bam > drop_seq.log 2>&1 &
```


###	DOWNLOAD

```
ip=$( aws --profile syryu ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' | grep "\." | tr -d '"' | tr -d ' ' )
echo $ip

rsync --archive --verbose --compress --rsh "ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete ec2-user@$ip:working/dropseq/ ~/github/unreno/syryu/singlecell/20180228a.drop_seq_alignment/
```





##	Post Run Notes

Probably should've stuck with a single name for the modified reference. Say mm10ryulab.


