#!/bin/bash

# This script executes all mapping steps consecutive, according to nick lomans pipeline (with modifications)
# with some modifications and additions

#the path to all related scripts and executables
execDir=$(dirname `readlink -f $0`)
Rcmd="R --vanilla -q --slave"

usage () {
    echo -e "\n  Usage: $0 -i <input.fas> -r <ref.fas> [ -o <outpath> -p <prefix>]"
    echo -e "-i  input.fas: path to sequencing output file (multiple fasta)"
    echo -e "-r  ref.fas: path to reference genome fasta file"
    echo -e "-o  outpath: path where the output is stored (default .) "
    echo -e "-p  prefix: if given this prefix will be used for all output data [default=input]"
    exit 1
}

while getopts 'i:o:r:p:h' opt
do
   case $opt in
        i) input=$OPTARG;;
        o) outpath=$OPTARG;;
        r) ref=$OPTARG;;
        p) prefix=$OPTARG;;
        h) usage;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

if [[ -z $input ]] || [[ -z $ref ]]; then
    usage
fi

refStrain=`basename ${ref}`
refStrain=`echo -e ${refStrain} | sed 's/\.fasta$//'` #just in case
input=`readlink -f ${input}`

if [ ! -r ${input} ]; then
    echo -e "Cannot read sequencing file: ${input}."
    exit 1
fi

if [ ! -r ${ref} ]; then
    echo -e "Cannot read reference file ${ref}."
    exit 1
fi

if [ -z "${outpath}" ]; then
    outpath=.
    echo -e "No output folder defined, using: ./"
fi

if [ -z "${prefix}" ]; then
    suffix=`echo "${input}" | awk -F . '{print $NF}'`
    prefix=`basename ${input} .${suffix}`
    echo -e "No prefix specified, using: ${prefix}"
fi

if [ ! -r ${outpath}/${prefix}_mapping ]; then
    echo -e "Creating output base directory"
    mkdir -p ${outpath}/${prefix}_mapping
fi

outpath=${outpath}/${prefix}_mapping


####start the mapping
#first link the reference and reads into the outdir for easier processing
cwd=`echo $PWD`
refAbs=`readlink -f ${ref}`
inAbs=`readlink -f ${input}`
cd ${outpath}
ln -s ${refAbs} REFERENCE
ln -s ${inAbs} INPUT
cd ${cwd}
mapOut=${outpath}/${prefix}
echo -e "\n[1]Start mapping\n"

# echo -e "   ${execDir}/gotmap.sh ${mapOut} ${input} ${ref} \n"
# ${execDir}/gotmap.sh ${mapOut} ${input} ${ref}

echo -e "   ${execDir}/gobwa.sh ${mapOut} ${input} ${ref} \n"
${execDir}/gobwa.sh ${mapOut} ${input} ${ref}

if [ $? -ne 0 ]; then exit $err ; fi
echo -e "\n   Mapping completed\n"
prefixUni=${prefix}.uni
mapOutUni=${mapOut}.uni

#####read bam and summarizes the indels/subs 
aliTABUni=${outpath}/${prefixUni}.alignment.tab
echo -e "[2]Creating alignment based error summary file (tsv table) ${aliTABUni} \n"
echo -e "   ${execDir}/read_bam ${mapOutUni}.sorted.bam mapping ${refStrain} > ${aliTABUni}\n" 
${execDir}/read_bam ${mapOutUni}.sorted.bam raw ${refStrain} > ${aliTABUni} 
if [ $? -ne 0 ]; then exit $err ; fi
echo -e "   done.\n"

#### create indel summary table
if [ $? -ne 0 ]; then exit $err ; fi
indelTABUni=${aliTABUni}.res
echo -e "[3]Creating R error summary Table ${aliTABUni}.res \n"
echo -e "   cat ${execDir}/indel_summary_table.R | ${Rcmd} --args ${aliTABUni}\n"  
cat ${execDir}/indel_summary_table.R | ${Rcmd} --args ${aliTABUni}  
if [ $? -ne 0 ]; then exit $err ; fi
echo -e "   done.\n"

##### calculate mean coverage
echo -e "[4]Calculating mean coverage  ${mapOutUni}.sorted.bam \n"
${execDir}/count_mean_coverage ${mapOutUni}.sorted.bam > ${outpath}/${prefixUni}.meanCovs
if [ $? -ne 0 ]; then exit $err ; fi
echo -e "   done.\n"

echo -e "[5]Calculating error rate  ${mapOutUni}.sorted.bam \n"
${execDir}/error_rate ${mapOutUni}.sorted.bam > ${outpath}/${prefixUni}.errorRate
if [ $? -ne 0 ]; then exit $err ; fi
echo -e "   done.\n"

echo -e "\n MAPPING pipeline finished.\n"
