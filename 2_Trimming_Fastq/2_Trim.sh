#!/usr/bin/bash

#=================#
# Load Directories
#=================#

# Absolute Path of your local directory (where the script is launch)
Local_PATH="/work/egay/white_shark_project/samples_mapping/TRIMMING_DATA/"

# INPUT FASTQ DIRECTORY
DIR_samples="/work/egay/white_shark_project/data/raw_fastq/samples/"

# Absolute Path of the trimmomatic ".jar" file. 
# Look at the documentation of your cluster. On genotoul, you need to call the absolute path of the '.jar' file. It depends on the soft management.
Trimmomatic="/usr/local/bioinfo/src/Trimmomatic/Trimmomatic-0.39/trimmomatic.jar"

# Absolute Path of adapter.fasta file. 
Adapter="/work/egay/white_shark_project/samples_mapping/TRIMMING_DATA_test2/adapt.fasta"

#=============#
# Trimming
#=============#

# write your samples prefix in the list separated by space character. Put the variable part of your fastq.gz file (usually it is the sample name)
Basename_samples="GN17608_S40_L002 GN17608_S52_L003 GN17768_S41_L002 GN1434_S33_L002 GN17768_S53_L003 GN1434_S45_L003 GN17774_S42_L002 GN1436_S34_L002 GN17774_S54_L003 GN1436_S46_L003 GN17775_S43_L002 GN17525_S35_L002 GN17775_S55_L003 GN17525_S47_L003 GN18404_S36_L002 GN17582_S38_L002 GN17582_S50_L003 GN18404_S48_L003 GN18411_S37_L002 GN17597_S39_L002 GN18411_S49_L003"

# LOOP ON SAMPLES NAME 
for name in $Basename_samples
do
    name_fastq_R1=${name}"_R1_001.fastq.gz" 	# get fastqR1 and fastqR2 file name of the sample x
    name_fastq_R2=${name}"_R2_001.fastq.gz"
    mkdir ${name}_trim 							# create folder for the trimming output for the sample x
    cd ${name}_trim 							# get into the folder for a specific samples x and create a script to run the trimming 
    cat > ${name}_TRIM.sh << EOF
#!/usr/bin/bash
#SBATCH -V
#SBATCH -o trim2.out
#SBATCH -e trim2.err
#SBATCH -J trimm
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --node=1

# IMPORT MODULE
module load bioinfo/Trimmomatic-0.38

# If your seqs have already good qualities, don't change the following parameters: ILLUMINACLIP / SLIDINGWINDOW / MINLEN / LEADING / TRAILING.
# If note look at the trimmomatic documentation
java -jar ${Trimmomatic}  PE -threads 20 -phred33 $DIR_samples${name_fastq_R1} $DIR_samples${name_fastq_R2} ${name}_R1.trim.paired.fastq.gz ${name}_R1.trim.unpaired.fastq.gz ${name}_R2.trim.paired.fastq.gz ${name}_R2.trim.unpaired.fastq.gz ILLUMINACLIP:${Adapter}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:100 LEADING:3 TRAILING:3

EOF
    sbatch ${name}_TRIM.sh # launch the trimming for the sample x
    cd ${Local_PATH} # return to the directory where the script is launch
done
