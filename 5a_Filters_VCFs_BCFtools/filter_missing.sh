#============================#
# Filter VCF for Mising data
#============================#

#===========================================#
# Fill the Pop_List and Max_missing arg
#===========================================#

# arguments infos:
#----------------#

# Each population in your gVCF file is comprised in quotes ""
# In Each population fill the column index of genotype of each ind

Pop_List=$4

# Int, number max of individue with missing genotype "./." per pop
# $2 is the second arg of the script, as to be set during the call of the script
Max_missing_in_pop=$2


# Int, Number max allowed of population which not pass the Max_missing_in_pop filter
# set to zero, because usually none population are allowed to have more missing data than antoher
Max_pop=0


# function help
#----------------#
if [ "$1" == "-h" ] ; then
    echo "Usage: `basename $0` [-h]
        This function select position in vcf file which correspond the chosen rate of Missing genotype per population

        Fill the following list in the script
        Pop_List=("col_ind1-pop1 col_ind2-pop1" "col_ind1-pop2 col_ind2-pop2" )
                Each population in your gVCF file is comprised in quotes ""
                In Each population fill the column index of genotype of each ind
        Max_missing_in_pop (int) : Number max of individue with missing genotype "./." per pop
        Max_pop (int) : Number max of population which not pass the Max_missing_in_pop filter

        Run the filter_vcf.sh with the command :
        sh filter_vcf.sh my_vcf_not_compressed.sh Max_missing_in_pop output_name
        "
    exit 0
fi

# function
#----------#

Max_missing_pop() {
        len=${#Pop_List[@]}
        #echo "nb pop= "$len
        nb_pop_miss=0
        for pop_number in `seq $len`
        do
                #echo "begin with pop number = "$pop_number
                pop_i=${Pop_List[$pop_number - 1 ]}
                #echo "ind in pop are "$pop_i
                count_missing=0
                Nb_Missing_Pop=0

                for ind_col in $pop_i
                do
                        Position=`echo $POS | awk '{print $1,$2}'`
                        #echo "current ind treated is :"$ind_col" at position : "$Position

                        genotype=`echo $POS | awk -v col=$ind_col -F" " '{print $col}' | awk -F":" '{print $1}'`
                        #echo "genotype of current ind is :"$genotype
                        if [ "${genotype}" == "./." ]
                        then
                                count_missing=$(($count_missing+1))

                        else
                         continue
                        fi
                done
                #echo "total number of missing data in pop i :"$count_missing
                if [ ${count_missing} -gt ${Max_missing_in_pop} ]
                then

                        nb_pop_miss=$(($nb_pop_miss + 1))
                        #echo ${nb_pop_miss}

                else
                    continue
                fi
        done
        #echo "current nb of pop with >2 missing data in pop i :"${nb_pop_miss}
}


# Reads each line (position) of the vcf file and run the Max_missing_pop function for each line
#-----------------------------------------------------------------------------------------------#
while read POS;
do
        Max_missing_pop
        if [ ${nb_pop_miss} -eq ${Max_pop} ]
                then
                echo "${POS}" >> $3
                else
                continue
        fi
done < $1