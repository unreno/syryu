# syryu

*	[s](#AWS-CLI)
*	[a](#SSH PEM Keys)
*	[x](#Start Linux EC2 Instance)





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







## Start Linux EC2 Instance





