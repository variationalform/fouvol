# `fouvol` - A Fourier proxy for a weakly singular Volterra kernel

>Home:
https://github.com/variationalform/fouvol

>Copyright (c) 2020, Simon Shaw
(https://github.com/variationalform, https://www.brunel.ac.uk/people/simon-shaw).

>The moral right of the author has been asserted.

>These codes are free software; you can redistribute them and/or
modify them under the terms of the GNU General Public License Version 3 - the terms of which should accompany this document.

This set of codes can be used to reproduce the results in the paper

>Approximate Fourier series recursion for problems
involving temporal fractional calculus
Simon Shaw and John R Whiteman
Brunel University London, 2022

To Appear in CMAME - **Computer Methods in Applied Mechanics and Engineering**

**_*TO DO:*_** Give DOI and link once the paper is in print.

The codes are in python  (with `python 2.7.17` used, as explained in the paper), with a bash script used to manage the batch solution. These scripts live in a folder called `fouvol`. To obtain this folder, `cd` to the directory you want to parent it and type...

```bash
git clone https://github.com/variationalform/fouvol.git
#git clone git@github.com:variationalform/fouvol.git
cd fouvol
./fouvol.py -h
./fouvol.py -v 0 -s 3 -a -0.5 -m 5 -T 10 --T1 0.5 -L 16 --Nt 500 -P
```

This should produce the  help page first, and then in the second command, `png` and `eps` versions of the graphics entitled `varphirep` and  `varphirepextended`. These correspond to the two Figures in Section 1 of the paper. All code corresponding to the creation and storing of `jpg` graphics has been commented out.

**_*TO DO:*_** Update with Figure numbers once the paper is in print.

Note that LaTeX is used as a `matplotlib` interpreter - if you  don't have LaTeX on the host machine you'll need to go through and comment/alter those lines. 

Next, try this:

```bash
./fouvol.py  -v 0 -s 1 -a -0.5 -T 10 -L 1 --Nt 512 -P
```
This should produce the error in the $N_t = 512$ in the results for the product-Rectangle method with $\alpha = −0.5$ and $T = 10$. Specifically,

```bash
FSerror =  0.0
|u(T)-U| =  0.00120137905322
```
Further,
```bash
./fouvol.py  -v 0 -s 3 -a -0.5 -m 5 -T 10 --T1 0.05 -L 32 --Nt 128
```

should give

```bash
FSerror =  0.123300371378
|u(T)-U| =  0.0253044042786
```
which corresponds to the error in the $N_t = 128$ and $L=32$ Fourier
method with $T = 10$, $\alpha = −0.5$, $T_1 = 0.05$ and a Hermite smoothness of $m = 5$.

The full set of results can be generated using the bash scripts `shortrun.sh` and `bigrun.sh`. These are designed to be run in a linux environment. In this environment execute this

```bash
THEN=`date`
./shortrun.sh | tee shortrun.out && ./compare.sh  | tee -a shortrun.out 
echo $THEN && date
```

to get a stripped-down set of the results in the paper. This should take less than an hour - perhaps around 20 minutes on a fast machine. The full set of results can be obtained by switching these commented lines over in `shortrun.sh`...

```bash
# JARGS="-J 5 16"
# LARGS="-L 3 11"
JARGS="-J 5 7"
LARGS="-L 3 7"
```
Beware though - with `16` and `11` determining the max number of time steps and Fourier components, it will take a while to run...

You can always tidy up with this...

```bash
rm -rf compare_?.sh compare.sh errortable.* *pyc results/ runout.txt timestable.* *.eps *.png *.txt *.out
```
This wont delete important files, but will get rid of everything that can be re-generated with other runs.

## Git management - some notes

Install, with an update first. For example, with Linux Mint:

```bash
sudo apt-get update
sudo apt install git
git config --global user.name "Your Name"
git config --global user.email your@email
git config --global core.editor vi
git config --global color.ui auto
git config --list
git --version
```

As usual, if you make changes, document them and then push as usual. For example,

```bash
git add *
git add .gitignore 
git status
git commit -m 'Initial commit of working code and explanatory README.md'
git push
git pull
git diff --staged
git fetch origin && git diff --name-only main origin/main # or master
# etc etc
```
The push may fail if the clone was made with `https`. If so,

```bash
# Paste ~/.ssh/id_rsa.pub into Git hub web page
# can't use http anymore - as in the clone above, so
ssh -T git@github.com
git remote -v
git remote set-url origin git@github.com:variationalform/fouvol.git
```
REF: <https://stackoverflow.com/questions/17659206/git-push-results-in-authentication-failed>
