# How to setup SSH key to avoid typing password

Author: Ji Huang

Date: 2023-11-17

Last update: 2023-11-17

1. Create a ssh key pair using `ssh-keygen -t rsa -f HP14_rsa`. `-f` is to specify the file name, change it to your preferred name.
This will create two files: public key `HP14_rsa.pub` and private key `HP14_ras`. **DO NOT share the private key. Keep it local!** On windows they will be in `\Users\username\.ssh\`. On Linux, they are in `~/.ssh`.

2. Upload the **public key** to server.
Just copy the content from the **public key** to the server `~/.ssh/authorized_keys`. On Linux, you can try  `ssh-copy-id -i ~/.ssh/HP14_rsa.pub username@greene.hpc.nyu.edu`.

3. Install VSCode-Remote SSH extension. Add location of the private key `IdentityFile C:\Users\YourUsername\.ssh\HP14_rsa` (change this to where your private key is). The config file should look like below:

```
Host greene.hpc.nyu.edu
  HostName greene.hpc.nyu.edu
  User username
  IdentityFile C:\Users\YourUsername\.ssh\HP14_rsa
```

Now you are able to connect to server without the need of typing password. Faster and safer.

You may need to use the `ssh -i` to specify the private key, for example `ssh -i ~/.ssh/desktop_greene username@username.hpc.nyu.edu`