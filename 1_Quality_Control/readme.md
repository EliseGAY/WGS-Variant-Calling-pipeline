# Run FastQC  

## Author  
**Elise GAY**  
📅 *02/2022*  

> ⚠️ **Please inform the authors before sharing.**  

---

## 📌 Aim  
Run FastQC on a list of FASTQ files.  

## 📁 Input  
- FASTQ file (zipped or unzipped) from a directory.  

## 🛠 Methods  
- Uses **FASTQC** tool with default parameters.  
- Adapted for **Genotoul** cluster with **SLURM** command.  
- To run the script on the cluster:  
`sh fastqc.sh`
- note that is included inside `fastqc.sh` itself (line 38) to run one SLURM script per sample.

## 📤 Output
Standard FASTQC output files.
