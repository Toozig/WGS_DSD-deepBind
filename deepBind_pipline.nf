params.fastaFile = ''
params.outputFolder = ''
// params.upload = 'false'

nextflow.enable.dsl=2

process cleanData {
    // there is a problem to work it with another venv when using the HUJI slurm
    label "medium_slurm"
   tag "process_${chrom}"

   publishDir params.outputDir , mode: 'copy'
    
    input:
    path chrom_tsv
    path regionFile
    path sampleFile
    
    output:
    path "${regionFile.simpleName}.${sampleFile.simpleName}.parquet"
    
    script:
    outputName = "${regionFile.simpleName}.${sampleFile.simpleName}"  
    """
    source ${params.VENV}
    clean_tsv.py ${outputName} ${chrom_tsv.join(' ')}
    """
}

process uploadData {
  label "small_slurm"
   tag 'upload'
    
    input:
    path parquetFile
    path sampleFile
    
    script:
    if (params.upload){
        println "upload files"
        """
        current_date=`date +'%d_%m_%y_%H'`

        ${params.DBXCLI} put $parquetFile $params.uploadDir/\$current_date/$parquetFile
        ${params.DBXCLI} put $sampleFile $params.uploadDir/\$current_date/$sampleFile
        """
    }
    else{
        println "not uploading"
        """
            echo "not uploading"
        """
    }
}

process createRegionFile{

    input:
        val chrom
        path regionFile
    output:
        path "${chrom}.${regionFile.simpleName}.bed*"
        
    script:
    output="${chrom}.${regionFile.simpleName}.bed"
    """
    grep -E '^${chrom}\\b' ${regionFile} > "${output}" 
    line_count=\$(wc -l < "$output")
    # Check if the line count is greater than 1000
    # if [ "\$line_count" -gt ${params.MAX_REGIONS} ]; then
    if [ "\$line_count" -gt ${params.REGION_SPLIT_SIZE} ]; then # fOR DEBUGGING!!!!
    split -l 1000 "${output}" "${output}_"
    fi
    """
}


process getSamples{
   publishDir params.sampleRawDir , mode: 'copy'
   label "medium_slurm"
   tag "process_${pathFile.split("/")[-1].replace('.vcf.gz','')}"

    input:
        tuple  val(pathFile), path(regionFile)
        path data
    output:
         path output
        // stdout
    script:
    output="${pathFile.split("/")[-1].replace('.vcf.gz','')}.tsv" 
    """
    # Loading necessary modules
    module load hurcs bcftools
    getSamples.sh ${pathFile} ${regionFile} ${output}
    """
}

def checkIfExistsResult(regionFile, sampleFile){

        if (file(params.outputDir).exists()){
            println "The file $regionFile.simpleName already processed with the samples in $sampleFile.simpleName \
            \n find the result in $params.outputDir "
         //   exit 1
        }

}

process createSegmentFile {
    label "small_slurm"
    tag "segmentFile_${id}"
    
    input:
    tuple val(id), val(seq)
    val window_size

    output:
    path "${chrom}.tsv"
    
    script:
    outputName = "${chrom}"  
    """
    source ${params.VENV}
    merge_chrom.py ${chrom} ${gnomAD} ${outputName} ${file_list.join(' ')}
    """
}



workflow {

    log.info """
        V C F - T S V   P I P E L I N E 
         fastaFile: ${params.fastaFile}
         outputFolder: ${params.outputFolder} 
         """
         .stripIndent()
    

    sequence = Channel.fromPath(params.fastaFile).splitFasta( record: [id: true, seqString: true ])
    sequence.view()

}
