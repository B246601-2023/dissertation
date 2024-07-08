#!/bin/bash

helpFunction()
{
   echo "Align a set of sequences in parallel to a specified reference sequence"
   echo "Usage: $0 -i unaligned_fasta_path -o output_file -t threads -r reference_seq_id"
   echo "\t-i Full path to unaligned fasta file of sequences"
   echo "\t-o Output file path for aligned fasta file. New sequences will be *APPENDED* to this file"
   echo "\t-t Number of threads to use"
   echo "\t-r Reference sequence ID to use for alignment"
   exit 1 # Exit script after printing help
}

while getopts "i:o:t:r:" opt
do
   case "$opt" in
      i ) inputfasta="$OPTARG" ;;
      o ) outputfasta="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      r ) refseqid="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inputfasta" ] || [ -z "$outputfasta" ] || [ -z "$threads" ] || [ -z "$refseqid" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Extract reference sequence
refaln=reference.fa
awk -v RS=">" -v FS="\n" -v seqid="$refseqid" '$1 == seqid { print ">" $0 }' $inputfasta > $refaln
#grep -A 1 "^>$refseqid" $inputfasta | sed '/^--$/d' > $refaln

# so we can get at it as a global variable
export REFERENCE_ALN=$refaln

echo "Aligning all sequences w.r.t. the reference ID: $refseqid"
echo ""
echo "Splitting unaligned input sequences into individual files"
N=$(grep ">" $inputfasta | wc -l)
N=$(( N*2 )) # doubling it is to ensure each record goes to one file

inputdir=$(dirname $inputfasta)
splitdir="$inputdir/split_files"  
mkdir -p $splitdir  
export splitdir
echo "temporary work dir is $splitdir.."
faSplit sequence $inputfasta $N $splitdir/individual_seq

echo ""
echo "Profile aligning each new sequence to the target alignment"
echo "This can take a while, be patient"
echo ""

profile_align()
{
   seqfile="$splitdir/$1"
   alfile=$seqfile"_profile_aligned.fa"
   final=$seqfile"_ind_aligned.fa"

   mafft --thread 1 --keeplength --quiet --add $seqfile "$REFERENCE_ALN" > $alfile

   name=$(grep ">" $seqfile | tr -d ">")
   echo "$name" | faSomeRecords $alfile /dev/stdin $final

   rm $seqfile
   rm $alfile
}

export -f profile_align

ls $splitdir | grep individual_seq | parallel -j $threads --bar "profile_align {}" > /dev/null

# join it all together and clean up
# note that here we *APPEND* to the global alignment, which allows us to add small batches of new stuff whenever we like
find $splitdir -name \*.fa_ind_aligned.fa -exec cat {} \; >> $outputfasta
find $splitdir -maxdepth 1 -name "individual_seq*" -delete
rm $refaln
