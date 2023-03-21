<h1 align="center">
    <img src="logo/logo.png" height="100">Rye
</h1>
Rye (<b>R</b>apid ancestr<b>Y</b> <b>E</b>stimation) is a large scale global ancestry inference algorithm that works from principal component analysis (PCA) data.  The PCA data (eigenvector and eigenvalue) reduces the massive genomic scale comparison to a much smaller matrix solving problem which is solved by using non-negative least square (NNLS) optimization.  
  
Codebase stage: development  
Developers and maintainers: Andrew Conley, Lavanya Rishishwar, Shivam Sharma  
Testers: Lavanya Rishishwar, Shivam Sharma, Emily Norris, Maria Ahmad

## Installation
We are working on improving the user experience with Rye.  This section will be updated in near future with easier installation procedure for pastrami and it's dependencies.  

### Dependencies
Rye requires the following OS/programs/modules to run:
* Linux/Mac OS - tested on RedHat OS, expected to run on Mac environments, will possibly work on WSL and Windows; though Mac, WSL, and Windows haven't been tested so far
* R - tested on several versions between 3.5.1 and 4.1.2
* R libraries = nnls, Hmisc, parallel, optparse, and crayon
  
#### git/R based installation
This is a straightforward installation procedure with git and R:
```
# Install R on linux, skip this step if already installed
sudo apt-get install r-base

# Alternatively, install R through conda
# we recommend switching from conda to mamba
conda install -c r r

# Install R libraries
# NOTE: If you maintain several environments, we recommend creating a dedicated environment for Rye, 
# though its dependencies shouldn't conflict with other dependencies
Rscript -e 'install.packages(c('nnls','Hmisc','parallel', 'optparse', 'crayon'))'

# Download the repository
git clone https://github.com/healthdisparities/rye

# Ensure permissions are correct
cd rye
# the uploader (me) wasn't smart enough to change the permissions the first time and refuses to reupload the file with right permissions
chmod +x rye.R

# Check if the install happened correctly
./rye.R -h
```

## Quickstart guide
Example files are provided to help user understand the input and test Rye easily.  If the environment is setup correctly, the user should be able to run the following command to run Rye:

```
./rye.R --eigenvec=examples/example.eigenvec --eigenval=examples/example.eigenval --pop2group=examples/pop2group.txt --rounds=5 --threads=2 --iter=5 --out=out
```

If the above command run successfully, the user should expect to see the following output on the screen:
```
[ Mar 03 2022 - 04:14:46 PM ] Parsing user supplied arguments...
[ Mar 03 2022 - 04:14:46 PM ] Arguments passed validation
[ Mar 03 2022 - 04:14:46 PM ] Running core rye with 2 threads
[ Mar 03 2022 - 04:14:46 PM ] Reading in Eigenvector file
[ Mar 03 2022 - 04:14:46 PM ] Reading in Eigenvalue file
[ Mar 03 2022 - 04:14:46 PM ] Reading in pop2group file
[ Mar 03 2022 - 04:14:46 PM ] Creating individual mapping
[ Mar 03 2022 - 04:14:46 PM ] Scaling PCs
[ Mar 03 2022 - 04:14:46 PM ] Weighting PCs
[ Mar 03 2022 - 04:14:46 PM ] Aggregating individuals to population groups
[ Mar 03 2022 - 04:14:46 PM ] Optimizing estimates using NNLS
Round 1/5 Mean error: 0.021117, Best error: 0.021098
Round 2/5 Mean error: 0.021097, Best error: 0.021092
Round 3/5 Mean error: 0.021101, Best error: 0.021100
Round 4/5 Mean error: 0.021100, Best error: 0.021098
Round 5/5 Mean error: 0.021100, Best error: 0.021096
[ Mar 03 2022 - 04:14:58 PM ] Calculate per-individual ancestry estimates
[ Mar 03 2022 - 04:14:58 PM ] Create output files
[ Mar 03 2022 - 04:14:59 PM ] Process completed
[ Mar 03 2022 - 04:14:59 PM ] 0.683 [ Mar 03 2022 - 04:14:59 PM ] 0.083 [ Mar 03 2022 - 04:14:59 PM ] 13.018 [ Mar 03 2022 - 04:14:59 PM ] 21.972 [ Mar 03 2022 - 04:14:59 PM ] 1.086
[ Mar 03 2022 - 04:14:59 PM ] The process took 0.68 seconds
```
The time calculation is off due to some silliness with R's time function and will be fixed one day.  
  
This command will produce the following files:
* out-20.7.Q
* out-20.fam
  
Here are the first 10 lines of the out-20.7.Q file:
```
European        Asian   Amerindian      SouthAsian      African NorthAfrican    Persian
ACB1    0.145498855162626       0       0.0129968919507834      0       0.84150425288659        0       0
ACB2    0.0369714135850651      0.0610180538310907      0.0183089317328864      0.326429353526824       0.557272247324134        0       0
ACB4    0.0361467110577818      0.00339369533154827     0       0       0.96045959361067        0       0
ACB5    0.126426064996172       0       0.00902862209257718     0       0.86454531291125        0       0
ACB7    0.0413888454515936      0       0.00476606015740593     0       0.953845094391001       0       0
ACB8    0.103024574359851       0       0       0       0.896975425640149       0       0
ACB11   0.0175766318211168      0       0       0       0.982423368178883       0       0
ACB12   0.124852156493598       0       0.0132053883019439      0       0.861942455204458       0       0
ACB14   0.0708302420422627      0       0       0       0.929169757957737       0       0
```

The first column is the user supplied individual ID (IID in the input file), the next 7 columns (in this example) are the groups that we provided in the pop2group.txt file.  If these grouping are changed, the resulting output groups will also change.  The \*.Q and \*.fam format are analogous to the output files produced by ADMIXTURE.

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
        
        --attempts=<ATTEMPTS>
                Number of attempts to find the optimum values (Default = 4)

        -h, --help
                Show this help message and exit
```

### Step 1: Computing eigenvalue/eigenvectors
The following example is developed with 1000 Genomes Project data as example dataset.  Run PLINK command on vcf.gz file to generate the eigenvector/eigenvalue files.

For plink v1.9:
```
plink --vcf ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --pca --out chr22_pca
```

For plink v2.0
```
plink2 --vcf ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --pca 20 --out chr22_pca
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
4. `attempts`

In our experience, the top 20 PCs are usually sufficient for assessing accurate results.  The more the number of PCs, the better the result.  There is a point of diminishing returns as later PCs explain very little of the total variance - we haven't tested this out thoroughly.

Rounds, iterations, and attempts control how many times Rye will try to identify the most optimal parameters sets (i.e., the weights and shrinkage).  Generally speaking, the longer you let this process run, the more likely we are to stumble upon the most optimal value.  We recommend users to start with a short run (i.e., small rounds, iterations and attempts) to get a sense of the ancestry prediction data and runtime of the software and then increase these values to higher numbers (e.g., 100-200 range) for a more accurate run.

## Publications

Conley AB, Rishishwar L, Ahmad M, Sharma S, Norris ET, Jordan I, Mariño-Ramírez L. Rye: genetic ancestry inference at biobank scale. Nucleic Acids Res. 2023. gkad149. doi: 10.1093/nar/gkad149 [[PubMed]](https://pubmed.ncbi.nlm.nih.gov/36928108/) [[Article]] (https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkad149/7079633)
