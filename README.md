# syryu

*	[AWS-CLI](#aws-cli)
*	[SSH PEM Keys](#ssh-pem-keys)
*	[Linux AMI](#linux-ami)
*	[Start Linux EC2 Instance](#start-linux-ec2-instance)





##	AWS-CLI

To access AWS resources via the terminal, `awscli` must be installed and configured.




###	Installation

pip is a python package manager and awscli is most easily distributed via it,
although it can be installed on its own.

A thorough guide on installing the AWS CLI is available ...

http://docs.aws.amazon.com/cli/latest/userguide/installing.html

`pip install --user --upgrade awscli`

Once installed, confirm that the aws script is in your PATH by running `aws`.
A help page should be displayed.
If something like `-bash: aws: command not found` is displayed,
the `aws` scripts are not installed in your PATH.




###	Configure

Once the `aws` scripts are installed, they need to be configured for access to your account resources.
This will require the creation of an Access Key on the AWS website.
This key's values need to be copied to your own computer.

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
	*	http://proteomics.ucsd.edu/Software/sps/sps-linux-x64-static-latest.zip
*	[MODa](http://prix.hanyang.ac.kr/download/software_archive/release/moda_v1.51.zip)
*	[GutenTag](http://fields.scripps.edu/yates/wp/?page_id=17)
	*	http://fields.scripps.edu/download/GutenTag.zip
*	[InsPecT](http://proteomics.ucsd.edu/Software/Inspect/)
	*	http://proteomics.ucsd.edu/Software/Inspect/Inspect.20120109.zip
*	[Specialize](http://proteomics.ucsd.edu/software-tools/specialize/)
	*	http://proteomics.ucsd.edu/Software/Specialize/Specialize_v1.0.zip
*	[Skyline](https://skyline.gs.washington.edu/labkey/project/home/software/Skyline/begin.view)
	*	Windows Only??? ---------------------------------------
*	[Proteowizard](http://proteowizard.sourceforge.net/downloads.shtml)
	*	Web download followed by scp upload
*	[R](https://www.r-project.org) (sudo yum install R)
*	[MatLab](http://www.mathworks.com/downloads/web_downloads/download_release?release=R2016b)




I should probably NOT set up the vi-style command line editor.






## Start Linux EC2 Instance

I recommend that you downloand and install
[my AWS script repo](https://github.com/unreno/aws)

```BASH
git clone https://github.com/unreno/aws
cd aws
cp Makefile.example Makefile
make install
```

For the moment, it is just 2 scripts, only 1 of which is useful to you.


If you've made it this far, you have installed and configured the aws cli.
You have also created a ssh key pair.

Running the following with the appropriate KEYNAME ...

`create_ec2_instance.bash --key ~/.aws/KEYNAME.pem`

... should suggest success, among other things, with ...

`An error occurred (DryRunOperation) when calling the RunInstances operation: Request would have succeeded, but DryRun flag is set.`



If so, running ...

`create_ec2_instance.bash -h`

... should point out that adding the option `--NOT-DRY-RUN` will actually start an instance.


If you would like a different instance type, specify it with `--instance-type INSTANCE_TYPE`.

If you would like more disk space on the instance, specify the number of GB with `--volume-size ###`.


Running the following should create an t2.micro instance using our most recent AMI.

`create_ec2_instance.bash --key ~/.aws/KEYNAME.pem --NOT-DRY-RUN`

After the script runs, it should show how to get the instance's IP address
and an appropriate command to ssh to this instance.


When finished using this instance, running `sudo halt` will terminate it.
You could also stop it from the web console.
You will stop being billed and the end of the running hour.
All data on the instance will be gone forever.









