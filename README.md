# Rye
Rye is a large scale global ancestry inference algorithm that works from principal component analysis (PCA) data.  The PCA data (eigenvector and eigenvalue) reduces the massive genomic scale comparison to a much smaller matrix solving problem which is solved by using non-negative least square (NNLS) optimization.  
  
Codebase stage: development  
Developers and maintainers: Andrew Conley, Lavanya Rishishwar, Shivam Sharma  
Testers: Lavanya Rishishwar, Shivam Sharma, Emily Norris  

## Installation
We are working on improving the user experience with Rye.  This section will be updated in near future with easier installation procedure for pastrami and it's dependencies.  

### Dependencies
Pastrami requires the following OS/programs/modules to run:
* Linux/Mac OS - tested on RedHat OS, expected to run on Mac environments, will possibly work on WSL and Windows; though Mac, WSL, and Windows haven't been tested so far
* R - tested on several versions between 3.5.1 and 4.X.X
* R libraries = nnls, Hmisc, parallel, optparse, and crayon
  
#### git/R based installation
This is a straightforward implmenetation procedure with git and R:
```
# Install R on linux, skip this step if already installed
sudo apt-get install r-base

# Alternatively, install R through conda
conda install -c r r

# Install R libraries
# NOTE: If you maintain several environments, we recommend creating a dedicated environment for Rye, 
# though its dependencies shouldn't conflict with other dependencies
Rscript -e 'install.packages(c('nnls','Hmisc','parallel', 'optparse', 'crayon'))'

# Download the repository
git clone https://github.com/healthdisparities/rye

# Ensure permissions are correct
cd rye
chmod +x rye.R

# Check if the install happened correctly
./rye.R -h
```

## Quickstart guide
This section will be populated in near future with small example datasets and commands for analyzing them.

## Basic usage
### General program structure overview
Rye is a three-step process:
1. Computing eigenvalue/eigenvectors from genomic data (VCF/BED/PED)
2. Creating a population to group mapping file
3. Estimating ancestry through Rye

These steps are described in more details in the following sections

### General usage of Rye
```
Usage: ./rye.R [options]


Options:
        --eigenval=<EVAL_FILE>
                Eigenvalue file [REQUIRED]

        --eigenvec=<EVEC_FILE>
                Eigenvector file [REQUIRED]

        --pop2group=<P2G_FILE>
                Population-to-group mapping file [REQUIRED]

        --output=<OUTPUT_PREFIX>
                Output prefix (Default = output)

        --threads=<THREADS>
                Number of threads to use (Default = 4)

        --pcs=<#PCS>
                Number of PCs to use (Default = 20)

        --rounds=<OPTIM-ROUNDS>
                Number of rounds to use for optimization (higher number = more accurate but slower; Default=200)

        --iter=<OPTIM-ITERS>
                Number of iterations to use for optimization (higher number = more accurate but slower; Default=100)

        -h, --help
                Show this help message and exit
```

### Step 1: Computing eigenvalue/eigenvectors
The following example is developed with 1000 Genomes Project data as example dataset.  Run PLINK command on vcf.gz file to generate the eigenvector/eigenvalue files
```
PLINK command - to be filled
```

### Step 2: Creating a population to group mapping file
A population to group mapping file simply provides instructions on how to aggregate populations into groups (e.g., GBR and CEU into European or Western European groups).  
  
This is a tab-separated file with headers and two mandatory columns - Pop and Group - rest of the columns are ignored:
```
Pop Group
CEU	European
GBR	European
YRI African
GWD African
CHB Asian
JPT Asian
...
```
  
The name of the file is inconsequential and the user can decide a name that is intuitive to them.  For our example, we will call this file ***pop2group.txt***.  

### Step 3: Estimating ancestry through Rye
Finally, run the following command:

```
./rye.R --eigenval=input.eigenval --eigenvec=input.eigenvec --pop2group=pop2group.txt --output=rye_output
```
  
This will create two output files:
1. rye_outputM.N.Q: An admixture style file with `M = number of PCs` and `N = number of groups in the pop2group.txt file`
2. rye_outputM.fam: A FAM file containing the mapping of the individuals

## Fine tuning
Rye currently supports the following options for fine tuning of the results:
1. `pcs`
2. `rounds`
3. `iter`

*Describe how changing `pcs`, `rounds`, and, `iter` changes.*
