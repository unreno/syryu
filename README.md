# syryu

*	[AWS-CLI](#aws-cli)
*	[SSH PEM Keys](#ssh-pem-keys)
*	[Linux AMI](#linux-ami)
*	[Start Linux EC2 Instance](#start-linux-ec2-instance)





##	AWS-CLI

To access AWS resources via the terminal, `awscli` must be installed and configured.


###	Installation

pip is a python package manager and awscli is most easily distributed via it.

A thorough guide on installing the AWS CLI is available ...

http://docs.aws.amazon.com/cli/latest/userguide/installing.html

`pip install --user --upgrade awscli`


###	Configure

Visit ...

https://syryu.signin.aws.amazon.com/console

and login. Then go to ...

https://console.aws.amazon.com/iam/home#users/

Click on your user and then click on Security Credentials and then Create Access Key. Keep this open as you won't have access to it if you close it. Of course, you could create another set. Nevertheless.

At a terminal command line, type `aws configure`. It will prompt you for the Access Key and Secret Access Key from the web page, as well as region (us-west-2) and output format (json). This info is also kept in ~/.aws/config and ~/.aws/credentials . Once entered, your script should work.

The region is important as instances and AMIs are specific to each.


###	Usage

`aws s3 sync s3://somebucket/ ~/s3/somebucket`

Trailing slashes are important for folders.

`aws s3 ls s3://somebucket/somefolder`

... is different than ...

`aws s3 ls s3://somebucket/somefolder/`



## SSH PEM Keys

You will need to create a key pair in order to start and login to an EC2 linux instance.
There is no username or password, just this PEM file which is the key.
Without it you cannot gain access to a running instance.

### Creation

I recommend that you create a key of some simple name and save it in a file of the same name with a .pem extension.
I also recommend that you keep it in `~/.aws`.

`aws ec2 create-key-pair --key-name KEYNAME --query 'KeyMaterial' --output text > ~/.aws/KEYNAME.pem`

And lastly I REALLY recommend that you `chmod 400 ~/.aws/*pem` so that others cannot view it.




## Linux AMI

I will be creating and saving a Linux AMI based on Amazon Linux.
I will `yum update` to upgrade many of the built in packages.
I will then modify the base environment ...
*	Update ~/.bashrc with some aliases and PATH modifications
*	add ~/.inputrc
* edit ~/.vimrc and add ~/.vim/
*	install the latest gawk
*	edit ~/.Renviron
*	edit ~/.Rprofile
And then install 
*	Matlab
*	[Spectral Network](http://proteomics.ucsd.edu/software-tools/spectral-networks/)
*	[MODa](http://prix.hanyang.ac.kr/download/software_archive/release/moda_v1.51.zip)
*	[GutenTag](http://fields.scripps.edu/yates/wp/?page_id=17)
*	[InsPecT: Depreciated Tool](http://proteomics.ucsd.edu/Software/Inspect/)
*	[Specialize](http://proteomics.ucsd.edu/software-tools/specialize/)
*	[Skyline](https://skyline.gs.washington.edu/labkey/project/home/software/Skyline/begin.view)
*	[Proteowizard](http://proteowizard.sourceforge.net/downloads.shtml)
*	[R](https://www.r-project.org) (possibly through yum)







## Start Linux EC2 Instance





