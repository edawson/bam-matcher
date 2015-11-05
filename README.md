# BAM-match #

A tool for determining whether two BAM files were sequenced from the same sample or patient. 

## Installation ##

```
cd /directory/path/where/bam-matcher/is/to/be/installed/
git clone https://bitbucket.org/sacgf/bam-matcher.git
```

Either include the path to BAM-matcher to the environment variable PATH, or move ```bam-matcher.py```, ```bam-matcher.conf``` and ```bam_matcher_html_template``` to a directory that is in the PATH already.


## Dependencies ##

**Python** 

(version 2.7)

**Python libraries**

* PyVCF
* HTSeq
* ConfigParser
* Cheetah

**Variant Callers**

(Require at least one)

* GATK (requires Java)
* VarScan2 (requires Java and Samtools)
* Freebayes


## Configuration ##




## Running tests ##


* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact