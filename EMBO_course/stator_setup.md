---
output:
  pdf_document: default
  html_document: default
---

# Stator --- Nextflow + conda setup

## Prerequisites

-   A machine where you can module add or otherwise load `Java` (Example for Slurm).
-   `conda` / `mamba` available for creating environments. curl, git, bash available.
-   Paths in the examples assume `/home/yy605` (adjust to your user paths).

## 1. Add Java & install Nextflow

If you're on Eddie you may need to load Java differently — example for slurm HPC:

``` bash
# Example (Cambridge MBU cluster)
module add java/17
```

Download Nextflow locally (pipe to bash):

``` bash
# download nextflow
curl -s https://get.nextflow.io | bash
```

Place the nextflow binary somewhere accessible (e.g. /home/yy605/nextflow) and make it executable.

## 2. Clone the Stator repo

``` bash
git clone -b develop https://github.com/AJnsm/Stator.git
```

Assume the repo is at `/home/yy605/Stator`.

## 3. Create Python conda environment (`Stator_py_env`)

Create and activate the env:

``` bash
conda create -n Stator_py_env
conda activate Stator_py_env
```

Install packages (as provided):

``` bash
conda install python==3.7.12
conda install numpy=1.19.4
conda install -c anaconda scipy=1.5.3
conda install -c conda-forge matplotlib=3.3.3
conda install -c anaconda pandas=1.1.5
conda install -c bioconda scanpy=1.7.2
conda install -c bioconda scrublet=0.2.1
conda install -c anaconda mkl=2020.4
conda install -c conda-forge seaborn=0.11.0
conda install -c conda-forge upsetplot=0.8.0
```

Pip installs

``` bash
pip install hypernetx==1.2.5
pip install pillow==9.4.0
pip install python-igraph==0.10.8
```

## 4. Create R conda environment (`rEnv`)

``` bash
conda create -n rEnv r-base=4.0.3 \
  bioconductor-graph=1.68.0 \
  bioconductor-rbgl=1.66.0 \
  r-ggm=2.3 \
  r-abind=1.4_5 \
  r-pcalg=2.6_12 \
  bioconductor-rgraphviz=2.34.0 \
  r-stringr=1.4.0 \
  r-devtools=2.3.2 \
  r-dagitty=0.3_0
```

Within R (inside `rEnv` `conda` environment), install:

``` r
install.packages("BiDAG")
```

## 5. Edit `./Stator/scripts/iterMCMC.R`

Because of package updates, convert objects to matrices / data.frames before writing to CSV. Replace or ensure the following code exists:

``` r
result_iterMCMC_endspace <- as.matrix(result_iterMCMC$endspace)
result_iterMCMC_endspace <- as.data.frame(result_iterMCMC_endspace)

# due to the updates of packages, we can’t directly save result_iterMCMC$CPDAG to csv now
result_iterMCMC_CPDAG <- as.matrix(result_iterMCMC$CPDAG)
result_iterMCMC_CPDAG <- as.data.frame(result_iterMCMC_CPDAG)

write.csv(result_iterMCMC_endspace, paste('MCMCgraph_', DSname, '.csv', sep=''), row.names=TRUE)
write.csv(result_iterMCMC_CPDAG, paste('CPDAGgraph_', DSname, '.csv', sep=''), row.names=TRUE)
```

## 6. Enable conda & add `conda_slurm` profile

In `./Stator/nextflow.config`,profiles block `profiles {`, add `conda.enabled = true` within conda block:

``` groovy
conda.enabled = true
```

In profiles block, add `conda_slurm` block, so that jobs can be submitted to Slurm (example block):

``` groovy
conda_slurm {
    conda.enabled = true
    includeConfig 'configs/conda_slurm.config'
}
```

The final `./Stator/nextflow.config` looks like:

``` groovy
manifest {
    name = 'AJnsm/Stator'
    author = 'Abel Jansma'
    homePage = 'https://github.com/AJnsm/Stator'
    description = 'Stator: Inferring cell states from gene expression data.'
    mainScript = 'Stator.nf'
    nextflowVersion = '>=23.04'
    version = '1.1'
}

params {
    report_dir = "${launchDir}/reports"
    dataType = "agnostic"
    boundBool = 0
    asympBool = 0
    dataDups = 0
    estimationMethod = "expectations"
    fracMito = 1
    minGenes = 0
    minCells = 1
    PCalpha = 0.05
    estimationMode = "MFI"
    bsResamps = 1000
    nRandomHOIs = 1000
    plotPairwiseUpsets = 0
    sigHOIthreshold = 0.05
    minStateDeviation = 0
    stateDevAlpha = 1.0
    doubletFile = false
    genesToOne = false
}


profiles {
    eddie_singularity {
        includeConfig 'configs/eddie_singularity.config'
    }

    singularity {
        includeConfig 'configs/singularity.config'
    }
    
    docker {
        includeConfig 'configs/docker.config'
    }

    conda {
        conda.enabled = true
    includeConfig 'configs/conda.config'
    }

    conda_slurm {
        conda.enabled = true
    includeConfig 'configs/conda_slurm.config'
    }
}


report {
    enabled = true
    overwrite = true
    file = "${params.report_dir}/report.html"
}

timeline {
    enabled = true
    overwrite = true
    file = "${params.report_dir}/timeline.html"
}

trace {
    enabled = true
    overwrite = true
    file = "${params.report_dir}/trace.txt"
}

dag {
    enabled = true
    overwrite = true
    file = "${params.report_dir}/pipeline.html"
}
```

## 7. Edit `./Stator/configs/conda.config` — per-process conda envs

Give the right path to conda environment.

The final `./Stator/configs/conda.config` looks like:

``` groovy
process {


    withName: makeData {
        conda = "/home/yy605/.conda/envs/Stator_py_env"
    }

    withName: estimatePCgraph {
        conda = "/home/yy605/.conda/envs/rEnv"
    }

    withName: iterMCMCscheme {
        conda = "/home/yy605/.conda/envs/rEnv"
    }

    withName: estimateCoups_2345pts_WithinMB {
        conda = "/home/yy605/.conda/envs/Stator_py_env"
    }

    withName: estimateCoups_6n7pts {
        conda = "/home/yy605/.conda/envs/Stator_py_env"
    }
    
    withName: identifyDTuples {
        conda = "/home/yy605/.conda/envs/Stator_py_env"
    }

}
```

## 8. Create `./Stator/configs/conda_slurm.config`

Create `./Stator/configs/conda_slurm.config` so that we can also submit jobs through slurm.

It is referenced by the conda_slurm profile, see that in `./Stator/nextflow.config`). This file should contain Slurm-specific executor settings (partition, cpus, memory, etc.).

``` groovy
process {
    executor="slurm"   

    withName: makeData {
        conda = "/home/yy605/.conda/envs/Stator_py_env"
        cpus = "${->params.cores_makeData}"
        memory = "${->params.mem_makeData}"
        time = "${->params.time_makeData}"
    }

    withName: estimatePCgraph {
        conda = "/home/yy605/.conda/envs/rEnv"
        cpus = "${->params.cores_PC}"
        memory = "${->params.mem_PC}"
        time = "${->params.time_PC}"
    }

    withName: iterMCMCscheme {
        conda = "/home/yy605/.conda/envs/rEnv"
        cpus = "${->params.cores_MCMC}"
        memory = "${->params.mem_MCMC}"
        time = "${->params.time_MCMC}"
    }

    withName: estimateCoups_2345pts_WithinMB {
        conda = "/home/yy605/.conda/envs/Stator_py_env"
        time = "${->params.time_HOIs_MB}"
        memory = "${->params.mem_HOIs_MB}"
        cpus = "${->params.cores_HOIs_MB}"
    }

    withName: estimateCoups_6n7pts {
        conda = "/home/yy605/.conda/envs/Stator_py_env"
        time = "${->params.time_HOIs_MB}"
        memory = "${->params.mem_HOIs_MB}"
        cpus = "${->params.cores_HOIs_MB}"
    }
    
    withName: identifyDTuples {
        conda = "/home/yy605/.conda/envs/Stator_py_env"
    }


}

executor {
    queueSize = "${->params.maxQueueSize}"
```

## 9. Prepare vignette test run

-   Edit `./Stator/vignette/params.json` to give the correct paths to raw data and any other parameters.

-   Ensure the data files are present at:

```         
./Stator/vignette/astrocytesVignette.csv
./Stator/vignette/userGenes.csv
```

Also set "executor" : "slurm", inside params.json if you want to use Slurm.

The final json file:

``` groovy
{
      "dataType"    : "expression",
      "rawDataPath" : "/home/yy605/Stator/vignette/astrocytesVignette.csv",
      "userGenes"   : "/home/yy605/Stator/vignette/userGenes.csv",      
      "nGenes"      : 50,
      "fracMito"    : 0.1,
      "minGenes"    : 1,
      "minCells"    : 1,
      "nCells"      : 4000,
      "PCalpha"     : 0.1,
      "bsResamps"   : 1000,

      "nRandomHOIs" : 10,
      "sigHOIthreshold":0.1,

      "asympBool"   : 0,
      "estimationMode" : "MFI",
    
      "executor"    : "slurm",
      "maxQueueSize"  : 25,
      "cores_makeData": 1,
      "cores_PC"      : 6,
      "cores_MCMC"    : 2,
      "cores_1pt"     : 4,
      "cores_HOIs_MB" : 6,
      "cores_HOIs_6n7" : 6,
      "cores_HOIs_plots": 6,

      "mem_makeData"  : "16G",
      "mem_PC"        : "32G",
      "mem_MCMC"      : "16G",
      "mem_1pt"       : "16G",
      "mem_HOIs_MB"   : "16G",
      "mem_HOIs_6n7"   : "16G",
      "mem_HOIs_plots" : "16G",

      "time_makeData" : "1h",
      "time_PC"       : "1h",
      "time_MCMC"     : "1h",
      "time_1pt"      : "1h",
      "time_HOIs_MB"  : "1h",
      "time_HOIs_6n7"  : "1h",
      "time_HOIs_plots"  : "1h"
  }
```

## 10. Paths (example)

-   Nextflow binary path: `/home/yy605/nextflow`

-   Stator repo path: `/home/yy605/Stator`

-   Conda environments root: `/home/yy605/.conda/envs/` (update if different)

## 11. Run Stator

### Run locally (use conda profile which uses `./Stator/configs/conda.config`), that does not submit jobs

``` bash
NXF_VER=23.04.4 /home/yy605/nextflow run /home/yy605/Stator  -profile conda   -params-file /home/yy605/Stator/vignette/params.json
```

`-profile conda` should make Nextflow use the `./Stator/configs/conda.config` profile which maps process -\> conda env.

### Submit to Slurm using `conda_slurm` profile

`-profile conda_slurm`

``` bash
NXF_VER=23.04.4 /home/yy605/nextflow run /home/yy605/Stator  -profile conda_slurm   -params-file /home/yy605/Stator/vignette/params.json
```
