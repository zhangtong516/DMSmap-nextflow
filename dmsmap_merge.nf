/*
 * File: dmsmap.nf
 * Created: Wednesday, 2nd December 2020 1:56:21 pm
 * Author: Zhang Tong (zhangtong516@gmail.com)
 * -----
 * Last Modified: Wednesday, 2nd December 2020 1:56:21 pm
 * Modified By: Zhang Tong (zhangtong516@gmail.com)
 * -----
 * Copyright (c) 2020 GIS
 * 
 */

params.help = false
def helpMessage() {
    log.info"""
    The typical command for running the pipeline is as follows:
    nextflow run ${baseDir}/icshape_RT.nf --sampleinfo <> --designinfo <> --batchName <>
    example:
        cd WORKINF_DIR
        nextflow run ${baseDir}/icshape_RT.nf --batchName batch1 --sampleInfo RHN1259.sample.tsv 

    By default, the pipeline will write the Output files into "results/" folder. It can be changed if user specify another directory.

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

// Run ID
if (params.containsKey('batchName') && params.batchName){
    params.outDir = "results/"
    // the reftype could be genome, transcriptome or longesttranscriptome and so on.
    params.collapsedSeqDir = params.outDir + "/" + params.genomeType + "/" + params.batchName + "/" + "collapsedFastq/" 
    params.trimSeqDir = params.outDir + "/" + params.genomeType + "/" + params.batchName + "/" + "trimmedFastq/" 
    params.mapResultDir = params.outDir + "/" + params.genomeType + "/" + params.batchName + "/" +  "mappedResult/" 
    params.qcDir = params.outDir + "/" + params.genomeType + "/" + params.batchName + "/" + "qc/"
    // params.coverageDir = params.outDir + "/" + params.genomeType + "/" + params.batchName + "/" + "coverageInfo/"
    params.mutationDir = params.outDir + "/" + params.genomeType + "/" + params.batchName + "/" + "mutation/"
    params.rpkmDir = params.outDir + "/" + params.genomeType + "/" + params.batchName + "/" + "rpkm/"
    params.reactivityDir = params.outDir + "/" + params.genomeType + "/" + params.batchName + "/" + "reactivity/"
} else {
    exit 1, "[Pipeline error] Please specify your design info file using `--batchName`!\n"
}

// sample info file
assert (params.sampleInfo != true) : "[Pipeline error] Please specify your sample info file using `--sampleInfo`!\n"
// assert (params.designInfo != true) : "[Pipeline error] Please specify your sample info file using `--designInfo`!\n"

// single end reads as in 2020 nov
ch_reads = process_sampleinfo(file(params.sampleInfo), params.batchName, "reads")

process trim_reads {
    tag "${name}"
    publishDir "${params.collapsedSeqDir}", mode: 'copy'

    input:
    set val(name), file(fastq1), file(fastq2) from ch_reads

    output:
    set val(name), file("*.R1.trimmed.fq.gz"), file("*.R2.trimmed.fq.gz") into trim_5end_reads

    script:
    """
    /mnt/software/bin/seqtk trimfq  -b ${params.SSII_trim_length} ${fastq1} | pigz -c -p ${task.cpus} > ${name}.R1.trimmed.fq.gz
    /mnt/software/bin/seqtk trimfq  -b ${params.SSII_trim_length} ${fastq2} | pigz -c -p ${task.cpus} > ${name}.R2.trimmed.fq.gz

    """
}

process collapse_reads {
    tag "${name}"
    publishDir "${params.collapsedSeqDir}", mode: 'copy'

    input:
    set val(name), file(trim_5end_read1), file(trim_5end_read2) from trim_5end_reads

    output:
    set val(name), file("*.R1.trimmed.collapsed.fq.gz"), file("*.R2.trimmed.collapsed.fq.gz") into collapsed_reads
    set val(name), file("*.collapsed.read_counts.txt") into collapsed_read_count

    script:
    """
    /mnt/projects/wenm/rnaStructure/ming/miniconda3/bin/clumpify.sh \
        in=${trim_5end_read1} in2=${trim_5end_read2} \
        out=${name}.R1.trimmed.collapsed.fq.gz out2=${name}.R2.trimmed.collapsed.fq.gz  \
        dedupe optical spany adjacent


    zcat ${trim_5end_read1} | wc -l | awk '{print "raw_reads\t" \$1/4}' > ${name}.collapsed.read_counts.txt
    zcat ${name}.R1.trimmed.collapsed.fq.gz | wc -l | awk '{print "collapsed_reads\t"\$1/4}' >> ${name}.collapsed.read_counts.txt

    """
}

process trim_adaptor {
    tag "${name}"
    publishDir "${params.trimSeqDir}", mode: 'copy'
    input:
    set val(name), file(collapsed_read1), file(collapsed_read2) from collapsed_reads  

    
    output:
    set val(name), file("*.R1.collapsed.trimmed.fastq.gz"), file("*.R2.collapsed.trimmed.fastq.gz") into trimmed_reads
    set val(name), file("*.trimmed.read_counts.txt") into trimmed_read_count
    file("*.read.trim.log.sum")
    file("*.trimAdaptor.log")

    script:
    """
    /mnt/projects/wenm/rnaStructure/ming/miniconda3/bin/cutadapt -j ${task.cpus} -n 2 -a ${params.nextera_adp1} -g ${params.nextera_adp2} \
        -m ${params.minLength} -q 20 -O 1 \
        -o ${name}.R1.collapsed.trimmed.fastq.gz -p ${name}.R2.collapsed.trimmed.fastq.gz \
        ${collapsed_read1} ${collapsed_read2} > ${name}.trimAdaptor.log
    
    perl ${baseDir}/scripts/parse_trim_log.pl ${name}.trimAdaptor.log  > ${name}.read.trim.log.sum

    zcat ${name}.R1.collapsed.trimmed.fastq.gz | wc -l | awk '{print "trimmed_reads\t"\$1/4}' > ${name}.trimmed.read_counts.txt
    """
}

process bowtie2_mappings {
    tag "${name}"
    publishDir "${params.mapResultDir}", mode: 'copy'

    input:
    set val(name), file(trimmed_reads1), file(trimmed_reads2) from trimmed_reads

    output:
    set val(name), file("*.mapped.best.sorted.bam") into mapped_sam_for_rpkm
    set val(name), file("*.mapped.best.sorted.bam") into mapped_sam_for_mutation
    set val(name), file("*.mapped.best.sorted.bam.bai") into mapped_sam_idx
    // set val(name), file("*.mapped.best.bam") into mapped_sam_for_coverage
    set val(name), file("*.mapped.read_counts.txt") into mapped_read_count

    file("*.mapping.log")


    script:
    """
    bowtie2 --local --no-unal -k 30 --un ${name}.trim.fastq.unAlign \
        -D ${params.tryTimes} -R 10 -L ${params.seedLength} -N ${params.mismatch} \
        -p ${task.cpus} --mp 3  -x ${params.genome} \
        -1 ${trimmed_reads1} -2 ${trimmed_reads2} \
        -S ${name}.trim.fastq.raw_mapping.sam 2> ${name}.mapping.log

    perl ${baseDir}/scripts/parse_bam_best_random.pl ${name}.trim.fastq.raw_mapping.sam > ${name}.mapped.best.sam
    
    

    cat ${name}.mapped.best.sam | /mnt/projects/wenm/rnaStructure/ming/miniconda3/bin/samtools view -Sb -o ${name}.mapped.best.bam
    /mnt/software//bin/sambamba sort  ${name}.mapped.best.bam
    
    awk '{if(NR==1){tot=\$1};if(NR==3){unmapped=\$1};if(NR==4){unique=\$1}}END{print tot"\t"tot-unmapped"\t"unique}' ${name}.mapping.log > ${name}.mapped.read_counts.txt

    rm -f ${name}.trim.fastq.raw_mapping.sam ${name}.mapped.best.sam ${name}.R1.collapsed.trimmed.fastq ${name}.R2.collapsed.trimmed.fastq ${name}.mapped.best.bam

    """
}


ch_for_qc = collapsed_read_count.join(trimmed_read_count).join(mapped_read_count)

process reads_qc{
    tag "${name}"
    publishDir "${params.qcDir}", mode: 'copy'
    input:
    set val(name), file(collapsed_read_count_file), file(trimmed_read_count_file), file(mapped_read_count_file) from ch_for_qc
    
    output:
    set val(name), file("*.reads_qc.txt") into read_qc

    script:
    """
    python ${baseDir}/scripts/read_stats_qc.py ${collapsed_read_count_file} ${trimmed_read_count_file} ${mapped_read_count_file} ${name} > ${name}.reads_qc.txt
    """
}

process estimate_rpkm{
    tag "${name}"
    publishDir "${params.rpkmDir}", mode: 'copy'
    input:
    set val(name), file(mapped_sam_file) from mapped_sam_for_rpkm
    
    output:
    set val(name), file("*.rpkm") into raw_rpkm_for_rt
    // set val(name), file("*.rpkm") into raw_rpkm_for_shape

    script:
    """
    perl ${baseDir}/scripts/estimate_RPKM.pl -i ${mapped_sam_file} -o ${name}.rpkm 
    """
}


/* 
merge all bam files from all replicates
*/ 

ch_bam_to_merge = mapped_sam_for_mutation.map{ [remove_last_item(it[0]), it[1] ] }.groupTuple().map{
    [ it[0], wrap_items(it[1], sep=" ")]
}

process merge_bam {
    tag "${name}"
    publishDir "${params.mapResultDir}", mode: 'copy'
    
    input:
    set val(name), val(bams_to_merge) from ch_bam_to_merge
    
    output:
    set val(name), file("*.merged.sorted.bam") into merged_bam_for_mutation
    file("*.merged.sorted.bam.bai")

    script:
    """
    samtools merge -@ ${task.cpus} ${name}.merged.bam ${bams_to_merge}
    /mnt/software/bin/sambamba sort ${name}.merged.bam 
    """

}


process mutation_calling {
    tag "${name}"
    publishDir "${params.mutationDir}", mode: 'copy'
    input:
    set val(name), file(mapped_sam_file) from merged_bam_for_mutation
    
    output:
    set val(name), file("*.mutation.tsv.gz") into mutation_count_file
    set val(name), file("*.rc") into mutation_count_file_binary 
    set val(name), file("*.rci") into mutation_count_file_binary_index  
    // set val(name), file("*.rpkm") into raw_rpkm_for_shape

    script:
    """
    ~/tools/RNAFramework/rf-count -t ${task.cpus} -wt ${task.cpus} \
        -r -m -f ${params.genome} ${mapped_sam_file} 

    mv ./rf_count/${name}.merged.rc ${name}.rc
    mv ./rf_count/index.rci ${name}.rci 
    
    ~/tools/RNAFramework/rf-rctools view  ${name}.rc -t  |\
        awk '{if(NF==1){gene=\$1};if(NF==3){ print gene"\t"\$0 } }' |\
        pigz -c -p ${task.cpus} > ${name}.mutation.tsv.gz
    """
}

process normalize {
    tag "${name}"
    publishDir "${params.reactivityDir}", mode: 'copy'
    input:
    set val(name), file(mutation_count) from mutation_count_file_binary
    set val(name), file(mutation_count_index) from mutation_count_file_binary_index  
    
    output:
    set val(name), file("*.complete") into normalize_done

    script:
    """
    ~/tools/RNAFramework/rf-norm -p ${task.cpus} -t ${mutation_count} -i ${mutation_count_index} \
        -sm 4 -nm 2 -rb AC -o ${params.reactivityDir}/${name}/
    touch ${name}.complete
    """
}

// extra functions to parse my sample.tsv input file
def wrap_items(input_files, sep=";") {
    // example input: [file1, file2]
    // output: file1+sep+file2
    def result =  input_files instanceof Path ? input_files.toString() : (input_files as List).join(sep)
    return result.toString()
}

def take_first_item(condition, sep="__") {
    // remove the batch name at the end
    def result = condition.toString().split(sep)[0]
    return result
}


// to remove conditions from name string
def remove_last_item(condition, sep="__") {
    // remove the batch name at the end
    def splitted = condition.toString().split(sep)
    def result = Arrays.copyOf(splitted, splitted.length - 1).join(sep)
    return result
}


def process_sampleinfo(tsvFile, batchName, mode) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            def sample      = row.stage // D0
            def libID       = row.libid // example: RHN1259Lib1
            def compBatch   = row.comparisonBatch //example: batch2
            def libType     = row.treatment // NAIN3 or DMSO
            def replicate   = row.rep // rep1
            def fastq1      = row.r1 // file with full path
            def fastq2      = row.r2
            def condition   = row.stage + "__" + row.treatment
            def unique_name = condition + "__" + replicate  

            if (mode == "reads" && compBatch == batchName) return [ unique_name, file(fastq1), file(fastq2) ]
        }
        .unique()
    }

