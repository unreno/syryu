#	DropSeq Processing

##	20180302


###	Preprocess locally ...

* Download mm10c.fasta.gz, mm10c.refFlat.gz, mm10c.dict.gz, mm10c.gtf.gz from Azure
* Download picard.jar `wget https://github.com/broadinstitute/picard/releases/download/2.17.0/picard.jar`
* Download fastq data from Illumina
* Convert B3 and B4 fastq data to bams
* Merge some of the previous sample data

```BASH
java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=1_L001.bam INPUT=1A_S4_L001.bam INPUT=1B_S3_L001.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=1_L002.bam INPUT=1A_S4_L002.bam INPUT=1B_S3_L002.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=1_L003.bam INPUT=1A_S4_L003.bam INPUT=1B_S3_L003.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=1_L004.bam INPUT=1A_S4_L004.bam INPUT=1B_S3_L004.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=2_L001.bam INPUT=2A_S2_L001.bam INPUT=2B_S1_L001.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=2_L002.bam INPUT=2A_S2_L002.bam INPUT=2B_S1_L002.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=2_L003.bam INPUT=2A_S2_L003.bam INPUT=2B_S1_L003.bam

java -jar ~/picard.jar MergeSamFiles ASSUME_SORTED=false SORT_ORDER=queryname OUTPUT=2_L004.bam INPUT=2A_S2_L004.bam INPUT=2B_S1_L004.bam
```



###	Gonna have a go at a Ubuntu Server on Azure


Started Standard_E16s_v3, 16CPU 128GB memory, 128GB disk via the web site.


ssh jake@52.168.35.1


###	UPDATE!

```BASH
sudo apt update
sudo apt full-upgrade
sudo apt install make

mkdir -p ~/working
mkdir -p ~/mm10c


git clone https://github.com/unreno/syryu
cd ~/syryu
ln -s Makefile.example Makefile
make install
```

```BASH
scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/bams/*.bam jake@52.168.35.1:working/

scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ~/github/unreno/syryu/drop_seq/mm10c.*.gz jake@52.168.35.1:mm10c/
```


```BASH


sudo apt install default-jre
sudo apt install unzip




wget -O Drop-seq_tools.zip http://mccarrolllab.com/download/1276/
unzip Drop-seq_tools.zip
/bin/rm -rf Drop-seq_tools.zip

wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz
tar xvfz 2.5.3a.tar.gz
/bin/rm -rf 2.5.3a.tar.gz 

wget https://github.com/broadinstitute/picard/releases/download/2.17.0/picard.jar

mkdir ~/bin
cp STAR-2.5.3a/bin/Linux_x86_64_static/STAR bin/

sudo apt install htop
sudo apt install r-base-core
sudo apt install openssl libssl-dev libcurl4-openssl-dev libssh2-1-dev


mkdir .R
cat > .Renviron <<EOF
R_LIBS="/home/jake/.R"
R_LIBS_USER="/home/jake/.R"
EOF

cat > .Rprofile <<EOF
local({
	r <- getOption("repos")
  r["CRAN"] <- "https://cloud.r-project.org/"
	options(repos = r)
})
EOF

R
install.packages("devtools")
library(devtools)
install.packages("httpuv")
#install.packages("igraph")
install_github("igraph/rigraph")
install.packages("Seurat")



cd ~/mm10c/
gunzip mm10c.fasta.gz
gunzip mm10c.gtf.gz
gunzip mm10c.refFlat.gz
gunzip mm10c.dict.gz


chmod -w ~/mm10c/mm10c*
cd ~/working/
mkdir -p ~/working/mm10c_star

nohup STAR --genomeFastaFiles ~/mm10c/mm10c.fasta --runMode genomeGenerate --genomeDir ~/working/mm10c_star --sjdbGTFfile ~/mm10c/mm10c.gtf --sjdbOverhang 100 --runThreadN 16 &

```








```BASH
mkdir -p ~/working/dropseq
mv ~/working/Log.out ~/working/dropseq/star_mm10c_creation.log
chmod 444 ~/working/mm10c_star/*
chmod 444 ~/working/dropseq/star_mm10c_creation.log
rmdir ~/working/_STARtmp
```

###	RUN DATA (individually, should take about 2 hours)

Drop-seq\_alignment.sh takes about 1 hours, then dge.bash/seurat.R takes about 1 hours or so.


```BASH
ssh -i /Users/jakewendt/.aws/JakeSYRyu.pem -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no jake@52.168.35.1
cd ~/working/dropseq
nohup drop_seq.bash ~/working/*.bam > drop_seq.log 2>&1 &
```


###	DOWNLOAD

```BASH
rsync --archive --verbose --compress --rsh "ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" --progress --delete jake@52.168.35.1:working/dropseq/ ~/github/unreno/syryu/drop_seq/20180302a.drop_seq_alignment/
```





##	Post Run Notes

Probably should've stuck with a single name for the modified reference. Say mm10ryulab.

