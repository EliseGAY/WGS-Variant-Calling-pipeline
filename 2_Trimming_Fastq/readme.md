# Run TRIMMOMATIC  

## Author  
**Elise GAY**  
📅 *02/2022*  

---  

## 📌 Aim  
Trim Illumina sequences.  

## 📂 Input  
- **R1 and R2 FASTQ.gz files**  
- **Adapter file**: A FASTA file containing adapter sequences.  
  - `Adapter.fasta` should contain the adapter sequences.  
  - These sequences depend on the sequencing technology.  
  - For **Illumina-PE**, you can use the FASTA file provided by Trimmomatic or create your own.  
  - **Tip:** Look at the FASTQC results to detect the type of adapters in your FASTQ files.  

## 🛠 Methods  
- Uses **Trimmomatic** with standard parameters for paired-end Illumina sequences:  
  ```plaintext
  ILLUMINACLIP:${Adapter}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:100 LEADING:3 TRAILING:3```
Adapted for Genotoul cluster with SLURM command.
To run the script on the cluster:

```sh trim.sh```

Note: The command
```sbatch script.sh```

is included inside trim.sh itself (line 53) to run one SLURM script per sample.

## 📤 Output
Four FASTQ.gz files:

- paired_R1.fastq.gz
- paired_R2.fastq.gz
- unpaired_R1.fastq.gz
- unpaired_R2.fastq.gz

One trimming report per sample

