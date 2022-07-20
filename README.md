# `fouvol` - A Fourier proxy for a weakly singular Volterra kernel

This code can be used to reproduce the results in the paper

>Approximate Fourier series recursion for problemsinvolving temporal fractional calculusSimon Shaw and John R Whiteman
Brunel University London, 2022

To Appear in CMAME - **Computer Methods in Applied Mechanics and Engineering**

The codes are in python, with a bash script used to manage the batch solution. These scripts live in a folder called `fouvol`. To obtain this folder, `cd` to the directory you want to parent it and type...

```bash
git clone https://github.com/variationalform/fouvol.git
#git clone git@github.com:variationalform/fouvol.git
cd fouvol
./fouvol.py -h
./fouvol.py -v 0 -s 3 -a -0.5 -m 5 -T 10 --T1 0.5 -L 16 --Nt 500 -P
```

This should produce the  help page first, and then in the second command, `png` and `eps` versions of the graphics entitled `varphirep` and  `varphirepextended`. These correspond to the two Figures in Section 1 of the paper. All code corresponding to the creation and storing of `jpg` graphics has been commented out.

**_*TO DO:*_** Update with Figure numbers once the paper is in print.

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
```bas
FSerror =  0.123300371378
|u(T)-U| =  0.0253044042786
```
which corresponds to the error in the $N_t = 128$ and $L=32$ Fourier
method with $T = 10$, $\alpha = −0.5$, $T_1 = 0.05$ and a Hermite smoothness of $m = 5$.

The full set of results can be generated using the bash scripts `shortrun.sh` and `bigrun.sh`. These are designed to be run in a linux environment. In this environment execute this

```bash
./shortrun.sh && ./compare.sh
```
to get a stripped-down set of the results in the paper. The full set of results can be obtained by switching these commented lines over in `shortrun.sh`...

```bash
# JARGS="-J 5 16"
# LARGS="-L 3 11"
JARGS="-J 5 6"
LARGS="-L 3 4"
```
Beware though - with `16` and `11` determing the max number of time steps and Fourier components, it will take a while to run...

You can always tidy up with this...

```bash
rm -rf compare_?.sh compare.sh errortable.* *pyc results/ runout.txt timestable.* *.eps *.png
```
This wont delete important files, but will get rid of everything that can be re-generated with other runs.

## Git management

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
```
The push failed because the clone was made with `https`.

```bash
Pasted ~/.ssh/id_rsa.pub into Git hub web page
can't use http anymore - as in the clone above, so
ssh -T git@github.com
git remote -v
git remote set-url origin git@github.com:variationalform/fouvol.git
```
