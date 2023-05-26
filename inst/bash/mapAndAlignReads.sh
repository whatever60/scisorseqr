#!/bin/bash
# by Hagen Tilgner for mapping of pacBio reads : December 31st, 2011; updated July 2012
# Edited AJ 04/2020

{
	echo "#################################"
	echo "# RUNNING [$0]"
	echo "# Current date: $(date)"
	echo "# Current dir: $(pwd)"
} 1>&2

###############
# 0. the arguments
echo "+++++++++++++++++ 1. arguments"

fastqGuide=${1}
echo "fastqGuide="$fastqGuide

outdir=${2}
echo "outdir="$outdir

tmpDir=${3}
echo "tmpDir="$tmpDir

numThreads=${4}
echo "numThreads="$numThreads

seqDirectory=${5}
echo "seqDirectory="$seqDirectory

annoGZ=${6}
echo "annoGZ="$annoGZ

mappingBAMGuide=${7}
echo "mappingBAMGuide="$mappingBAMGuide

scriptDir=${8}
echo "scriptDir="$scriptDir

pyScriptDir=${9}

genV=${10}
echo "genome + version"=$genV

echo "+++++++++++++++++ 1b. deduced from arguments"
tmpdir1=$tmpDir
mkdir -p $tmpdir1

echo "+++++++++++++++++ 1c. input check"
# if [ ! -f $outdir"/mapping.bestperRead.bam" ]
# then
mappingBAM=$tmpdir1"/mapping.bam"
samtools merge -f $mappingBAM $(cat $mappingBAMGuide)

echo "+++++++++++++++++ 2b. sorting annotation and projection"

echo "+++++++++++++++++ 2b1 sorting annotation"
zcat $annoGZ | awk '$3=="exon"' | sort -gk4,4 >$tmpdir1"/sortedAnno"

echo "+++++++++++++++++ 2c. how many were well mapped"
mappingWellMappedSAM=$tmpdir1"/mapping.bestperRead.sam"
mappingWellMappedBAM=$tmpdir1"/mapping.bestperRead.bam"

echo "+++++++++++++++++ 2c1 finding best hits"

samtools view -H $mappingBAM >$mappingWellMappedSAM

for i in $(cat $mappingBAMGuide); do
	samtools view $i | awk '$3!="*"' | awk -v file=$i '
	BEGIN{
		u1=0;
		u2=0;
		d1=0;
		d2=0;
	}
	{
		n[$1]++;
		v[$1"\t"n[$1]]=$5;
		line[$1"\t"n[$1]]=$0;
	}
	END{
		for(k in n){
			if(n[k]==1){
				if(v[k"\t1"]>=20){
					print line[k"\t1"];
					u1++;
					continue;
				}else{
					d1++;
					continue;
				}
			} 
			max1=0;
			max2=0;
			argmax1=0;
			argmax2=0;
			for(i=1;i<=n[k];i++){
				if(v[k"\t"i]>=max1){
					max2=max1;
					argmax2=argmax1;
					max1=v[k"\t"i];
					argmax1=i;
					continue;
				}
				if(v[k"\t"i]>=max2){
					max2=v[k"\t"i];
					argmax2=i;
					continue;
				}
			}
			if(max1>=20 && max1>=max2+20 && max2<=30){
				print line[k"\t"argmax1];
				u2++;
			}else{
				d2++;
			}
		} 
		x=(u1+u2)/(u1+u2+d1+d2);  
		print "wellMapped:\t"file"\t"u1"\t"u2"\t"d1"\t"d2"\t"x > "/dev/stderr";
	}'
done >>$mappingWellMappedSAM

echo "+++++++++++++++++ 2c2 sam to bam conversion"
samtools view -b -S $mappingWellMappedSAM >$mappingWellMappedBAM

# else
# mappingWellMappedBAM=$outdir"/mapping.bestperRead.bam";
# zcat $annoGZ | awk '$3=="exon"' | sort -gk4,4 > $tmpdir1"/sortedAnno";
# echo "---> Well mapped bam already exists, moving on to next step"
# fi

echo "+++++++++++++++++ 2d. removing ribosomal RNA hits"
# if [ ! -f $outdir"/mapping.bestperRead.noRiboRNA.bam" ]
# then
echo "+++++++++++++++++ 2d1. finding ribosomal RNA hits"
zcat $annoGZ | awk '$3=="exon" && ($14~/rRNA/ || $14~/rRNA_pseudogene/ || $14~/Mt_rRNA/)' >$tmpdir1"/tmp.AnnoRibo.gff"
bedtools intersect -abam $mappingWellMappedBAM -b $tmpdir1"/tmp.AnnoRibo.gff" -split -bed >$tmpdir1"/reads2Remove.bed"

echo "+++++++++++++++++ 2d2. removing ribosomal RNA hits"
samtools view -h $mappingWellMappedBAM | awk -v hits=$tmpdir1"/reads2Remove.bed" \
	'BEGIN{comm="cat "hits; while(comm|getline){rm[$4]=1;}}{if($1 in rm){next;} print;}' |
	samtools view -b -S - >$outdir"/mapping.bestperRead.noRiboRNA.bam"
# else
# echo "---> ribo-depleted bam already exists, moving on to next step"
# fi

echo "+++++++++++++++++ 3. further analysis"

echo "+++++++++++++++++ 3a. checking consensus"

# if [ ! -f $outdir"/mapping.bestperRead.noRiboRNA.introns.gff.gz" ]
# then
echo "+++++++++++++++++ 3a1. getting introns and exons in gff format"
# Convert BAM file to BED format
bedtools bamtobed -i $outdir"/mapping.bestperRead.noRiboRNA.bam" -bed12 | awk -v genV=$genV 'BEGIN {
    # Set the output field separator to a tab character
    OFS = "\t";
}

{
    # Increment the count for the read ID
    seen[$4]++;

    # Split the block sizes and block starts into arrays
    n = split($11, blockLen, ",");
    m = split($12, offset, ",");

    # If the number of block sizes and block starts do not match, print an error and exit
    if (m != n) {
        print "ERROR after bamtobed: m=" m "!=n=" n > "/dev/stderr";
        exit(0);
    }

    # If the read is not spliced (has less than 2 exons), skip to the next line
    if (n < 2) {
        next;
    }

    # Calculate the end of the first exon
    exEnd = $2 + blockLen[1];

    # Loop over the remaining exons
    for (i = 2; i <= n; i++) {
        # Calculate the start of the next exon
        exStart = $2 + offset[i] + 1;

        # Print the GFF line for the intron (between the previous exon and the next exon)
        print $1, genV, "intron", exEnd + 1, exStart - 1, ".", $6, ".", "transcript_id_with_chr=" $4 ".path" seen[$4] "@" $1;

        # Update the end of the "previous" exon for the next iteration
        exEnd = exStart + blockLen[i] - 1;
    }
}' >$outdir"/mapping.bestperRead.noRiboRNA.introns.gff"

# Use bedtools to convert the BAM file to BED format
bedtools bamtobed -i $outdir"/mapping.bestperRead.noRiboRNA.bam" -bed12 | awk -v genV=$genV \
	'BEGIN {
    # Set the output field separator to a tab character
    OFS = "\t";
}

{
    # Increment the count for the read ID
    seen[$4]++;

    # Split the block sizes and block starts into arrays
    n = split($11, blockLen, ",");
    m = split($12, offset, ",");

    # If the number of block sizes and block starts do not match, print an error
    if (m != n) {
        print "ERROR after bamtobed: m=" m "!=n=" n > "/dev/stderr";
        exit(0);
    }

    # Loop over the blocks
    for (i = 1; i <= n; i++) {
        # Calculate the start and end of the exon
        exStart = $2 + offset[i] + 1;
        exEnd = exStart + blockLen[i] - 1;

        # Print the GFF line
        print $1, genV, "cDNA_match", exStart, exEnd, ".", $6, ".", "read_id \"" $4 ".path" seen[$4] "\";";
    }
}' >$outdir"/mapping.bestperRead.noRiboRNA.gff"

# else
# echo "---> intron.gff file already exists, moving on"
# fi

echo "+++++++++++++++++ 3a2. getting dinucleotides at intron borders"
echo "+++++++++++++++++ 3a2.a. preparing commands for parallel execution"
# if [ ! -f $outdir"/parallel.comp.anno.guide.3.siteSeq" ]
# then

collectCommand="cat "
# Loop over chromosomes contained in the GFF file
for i in $(cat $outdir"/mapping.bestperRead.noRiboRNA.introns.gff" | 
    awk '{if(!($1 in s))print $1; s[$1]=1;}' | sort -u |
    awk '{if($1~/chr|\./) {print } else {print "chr"$1}}'); do
    
    seqfile=$seqDirectory"/"$i".fa.gz"  # Define the path to the sequence file for this chromosome
    f2tFunc=$scriptDir/"f2t.sh"  # Define the path to a script
    
    # Build a command string
    str="source "$f2tFunc"; zcat "$seqfile" | FastaToTbl | awk -v chr="$i" -v file="$outdir"/mapping.bestperRead.noRiboRNA.introns.gff -f "$scriptDir"/v0.1.getSpliceSiteSequence.awk > "$tmpdir1"/site_sequence."$i
    
    echo $str >> $outdir"/parallel.comp.anno.guide.3.siteSeq"  # Append the command to a file
    
    collectCommand=$collectCommand" "$tmpdir1"/site_sequence."$i  # Build a collection command
    
    # rmCommand=$rmCommand" "$tmpdir1"/site_sequence."$i  # Commented out: would build a remove command
done

cat $outdir"/parallel.comp.anno.guide.3.siteSeq" | shuf >$outdir"/tmp"
mv $outdir"/tmp" $outdir"/parallel.comp.anno.guide.3.siteSeq"

echo "+++++++++++++++++ 3a2b execution "
python3 $pyScriptDir/v0.2.executeInParallel.py --commandFile $outdir"/parallel.comp.anno.guide.3.siteSeq" --n $numThreads

echo "+++++++++++++++++ 3a2c collecting results and removing temporary files"
# echo $collectCommand;
# echo $rmCommand;
`$collectCommand > $outdir"/mapping.bestperRead.RNAdirection.introns.type.gff"`;
# `$rmCommand`;

gzip -c -f $outdir"/mapping.bestperRead.RNAdirection.introns.type.gff" >$outdir"/mapping.bestperRead.RNAdirection.introns.type.gff.gz"

# else
# echo "---> parallel process already done, moving on"
# fi

echo "+++++++++++++++++ 3b. strand correction and follow up analysis"

echo "+++++++++++++++++ 3b1. getting mapping.bestperRead.noRiboRNA.bam in gff format "

# if [ ! -f $outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz" ]
# then

awk -v bestMatchFile=$outdir"/mapping.bestperRead.noRiboRNA.gff" \
	-v intronTypeFile=$outdir"/mapping.bestperRead.RNAdirection.introns.type.gff.gz" \
	-v outputCorrectFile=$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.gff" \
	-v outputUnClassifiableFile=$outdir"/mapping.bestperRead.RNAdirection.withWeirdIntron.gff" \
	-f $scriptDir/v1.0b.correct.strand.bySpliceSiteConsensus.awk

gzip -c -f $outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.gff" > \
	$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.gff.gz"

samtools view -h $outdir"/mapping.bestperRead.noRiboRNA.bam" |
	awk -v CSMM=$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.gff.gz" 'BEGIN {
    comm="zcat " CSMM;
    while (comm | getline) {
        split($10, a, /.path|"/);
        use[a[2]] = 1;
    }
    use["@PG"] = 1;
    use["@HD"] = 1;
    use["@SQ"] = 1;
}
{
    if (!($1 in use)) {
        next;
    }
    print;
}' |
	samtools view -b -S - >$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.bam"

echo "+++++++++++++++++ 3b3. getting introns"

zcat $outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.gff.gz" |
	awk -v introns=$outdir"/mapping.bestperRead.noRiboRNA.introns.gff" '{
    split($10, a, /"/);
    use[a[2]] = $7;
}
END {
    comm="cat " introns;
    while (comm | getline) {
        split($9, a, /(transcript_id_with_chr=|@)/);
        if (a[2] in use) {
            s = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" use[a[2]];
            for (i = 8; i <= NF; i++) {
                s = s "\t" $i;
            }
            print s;
        }
    }
}' >$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.introns.gff"

transIdColumn=14

echo "+++++++++++++++++ 3b4. getting genes"
awk -v sortedAnno=$tmpdir1"/sortedAnno" \
    -v intronGFF=$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.introns.gff" \
	-v feature=exon \
    -v transIdColumn=$transIdColumn \
    -v geneIdColumn=10 \
	-f $scriptDir"/"v1.0b.getGenes_andHighlightProblematicReadsAndGenes.awk |
	gzip -c >$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.transcriptWise.genes.gz"

echo "+++++++++++++++++ 4. zipping what has not be zipped before"
gzip -c -f $outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.introns.gff" >$outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.introns.gff.gz"

# else
# echo "---> transcriptWise.genes.gz already exists, skipping this step and moving on"
# fi

echo "+++++++++++++++++ 5. new isoforms"
mkdir -p $outdir"/newIsoforms_vs_Anno_ignoreAnno"

# if [ ! -f $outdir"/newIsoforms_vs_Anno_ignoreAnno/stretches.gz" ]
# then
echo "+++++++++++++++++ 5a. annotation"
zcat $annoGZ | awk -v feature=exon -v transcol=$transIdColumn -f $scriptDir"/"exon.gtf.to.intron.gff3.awk |
	sort -k1,1 -gk4,4 >$tmpdir1"/anno.introns.gff"
gzip -cf $tmpdir1"/anno.introns.gff" > $tmpdir1"/anno.introns.gff.gz"

echo "+++++++++++++++++ 5b. RNAseq"
zcat $outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.introns.gff.gz" |
	sort -k1,1 -gk4,4 | gzip -c >$tmpdir1"/rnaSeq.introns.gff.gz"

zcat $outdir"/mapping.bestperRead.RNAdirection.withConsensIntrons.gff.gz" |
	awk '{exon=$1"_"$4"_"$5"_"$7; if(!a[$10]) {a[$10]=$10"\t;%;"exon} else {a[$10]=a[$10]";%;"exon} }END \
{for(i in a) {print a[i]}}' |
	gzip -c >$outdir"/newIsoforms_vs_Anno_ignoreAnno/exonStretches.gz"

echo "+++++++++++++++++ 5c. comparison"
awk -v keyword=intron -v anno=$tmpdir1"/anno.introns.gff.gz" -v pred=$tmpdir1"/rnaSeq.introns.gff.gz" \
	-v outfile=$outdir"/newIsoforms_vs_Anno_ignoreAnno/stretches" -v trIDcolumn=9 \
	-f $scriptDir/v0.1.findNotAnnotatedStretches.awk

gzip -cf $outdir"/newIsoforms_vs_Anno_ignoreAnno/stretches" >$outdir"/newIsoforms_vs_Anno_ignoreAnno/stretches.gz"

# else
# echo "---> stretches file already exists, finishing up."
# fi

echo "+++++++++++++++++ 6. done"
if [ ! -f $outdir"/mapping.bestperRead.bam" ]; then
	mv $tmpdir1"/mapping.bestperRead.bam" $outdir
fi
# rm $tmpdir1"/"*
# rmdir $tmpdir1
