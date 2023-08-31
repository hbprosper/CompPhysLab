# Computational Physics Laboratory
Welcome to the GitHub repository for the Computational Physics Laboratory course (PHZ1140C) at Florida State University. In the spirit of open source this is a open repository which anyone can access. This course is designed for incoming physics majors. The premise of the course is that instructive physics investigations can be undertaken by dedicated students even if they have barely begun their formal university physics education. For example, a student who has made the effort to master an algorithm for numerically solving a set of coupled 1st-order ordinary differential equations (ODE) is fully equipped to investigate many interesting physical systems. The course assumes that you have a good command of high school algebra and trignometry, but not much else! But if you work through the examples and exercises diligently, you'll be able to solve problems from the relatively banal, projectile motion with friction, to the relatively exotic, computing the orbits of photons around a black hole! This is very much a hands-on course...so don't be afraid to get your hands dirty.

## Introduction
The notebooks in this repository depend on several well-known Python
modules, all open source, which are compatible with all well-known [operating systems](https://www.lifewire.com/operating-systems-2625912).

| __modules__   | __description__     |
| :---          | :---        |
| numpy         | array manipulation and numerical analysis      |
| matplotlib    | a widely used plotting module for producing high quality plots |
| vpython       | an excellent 3D animation module |
| sympy         | an excellent symbolic mathematics module |
| scipy         | scientific computing    |
| pandas        | data table manipulation, often with data loaded from csv files |
| iminuit       | an elegant wrapper around the venerable CERN minimizer Minuit |


The following modules are also of interest.

| __modules__   | __description__     |
| scikit-learn  | easy to use machine learning toolkit |
| pytorch       | a powerful, flexible, machine learning toolkit |
| imageio       | photo-quality image display module |
| emcee         | an Markov Chain Monte Carlo module |
| tqdm          | progress bar |
| joblib        | module to save and load Python object |
| importlib     | importing and re-importing modules |

##  Installation
The simplest way to install these Python modules is first to install __Anaconda__ on your laptop by following the instructions at

https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

for your operating system.

I recommend installing __miniconda3__, a slimmed down version of Anaconda, which comes pre-packaged with Python 3.

Software release systems such as Anaconda (__conda__ for short) make
it possible to have several separate self-consistent named software
*environments* on a single machine. For example, you
may need to use Python 3.7.5 and an associated set of compatible
Python software and at other times you may need to use Python 3.9.13 with
software that requires that version of Python.  If you install software without using *environments* there is
the danger that the software on your machine will eventually become
inconsistent. Anaconda and its lightweight companion miniconda, and other software environment systems,
provide a way, for example, to have a software *environment* on your machine that is
consistent with Python 3.7.5 and another that is consistent with
Python 3.9.13.  

Of course, like anything humans make, Anaconda or miniconda3 is not
perfect. There are times when the only solution is to remove an
environment using
```bash
conda env remove -n <name>
```
where \<name\> is the name of the environment and rebuild it by reinstalling the desired Python modules. Therefore, I recommend that you make a note of all the modules/packages you have installed and the order in which you installed.

In this course we shall be using the __Linux__ operating system. Therefore, the instructions below assume you are working within a __Terminal__ window on a Linux or Linux-like system such as Mac OSX.

### Miniconda3

After installing miniconda3, it is a good idea to update conda using the command
```bash
conda update conda
```
#### Step 1 
Assuming conda is properly installed and initialized on your laptop (or desktop), you can create an environment, here called *comphys* 
```bash
conda create --name comphys
```
and activate it using the command
```bash
conda activate comphys
```
The environment need be created only once, but you must activate it whenever you launch a new terminal window.

#### Step 2 

Make sure you have activated the __comphys__ environment. Then execute the commands below in the order given. 
```bash
  conda install -c conda-forge numpy
	conda install –c conda-forge matplotlib
  conda install -c conda-forge vpython
  conda install -c conda-forge sympy
  conda install -c conda-forge scipy
	conda install –c conda-forge pandas
  conda install -c conda-forge iminuit
```
If all goes well, these commands will have installed the Python modules numpy, matplotlib, vpython, etc. The command that installs vpython also installs __jupyterlab__, which includes jupyter notebook. However, verify this by looking at the printout to the screen. The main tool we'll use in the course is the jupyter notebook.

#### Step 3
If time permits, you'll be introduced to the basics of machine learning-based artificial intelligence (ML/AI). This will require installing the __pytorch__ and __scikit-learn__ modules: 

```bash
	conda install –c conda-forge pytorch
	conda install –c conda-forge scikit-learn
```

#### Step 4
The command __git__ is needed to download the __CompPhysLab__ package from GitHub. If __git__ is not on your machine, it can be installed using the command
```bash
	conda install –c conda-forge git
```
To install __CompPhysLab__, whether on your laptop, or on your classroom RCC virtual desktop, do
```bash
  	cd 
	git clone https://github.com/hbprosper/CompPhysLab
```
This should clone, that is, replicate, a directory (aka folder) in your home directory called __CompPhysLab__. That folder will contain all of the examples and exercises we shall work through during the semester.
