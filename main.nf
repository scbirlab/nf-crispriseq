#!/usr/bin/env nextflow

/*
========================================================================================
   INSPECT Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-bccount
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2
pipeline_title = """\
                 S C B I R   C R I S P R i   P O O L E D   F I T N E S S   P I P E L I N E
                 =========================================================================
                 """
                 .stripIndent()

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println pipeline_title + """\
         Nextflow pipeline to count guides from SRA files and calculate
         fitness changes.

         Usage:
            nextflow run sbcirlab/nf-inspect --help
            nextflow run sbcirlab/nf-inspect --sample_sheet <csv> --fastq-dir <dir>
            nextflow run sbcirlab/nf-inspect -c <config-file>

         Required parameters:
            sample_sheet      Path to SRA run table identifying the SRA IDs to download, and with columns 
                              corresponding to conditions and replicates.
            ---
            conditions
            sequencing_group
            expansion_group
            reference


         Optional parameters (with defaults):
            sample_names = "sample_id"      Column heading from `sample_sheet` containing the names of the samples.
            name_column = ${params.name_column}         Which column from the guide table to use as the guide name.
            sequence_column = ${params.sequence_column}         Which column from the guide table to use as the guide sequence.
            trim_qual = ${params.trim_qual}             For `cutadapt`, the minimum Phred score for trimming 3' calls
            min_length = ${params.min_length}           For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded.

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   System.exit(0)
}

/*
========================================================================================
   Check parameters
========================================================================================
*/

required_params = [
   'sample_sheet', 
   'sample_names', 
   'guides'
]

for ( p in params ) {
   if ( p.value == null && p.key in required_params ) {
      throw new Exception("ERROR: required parameter ${p.key} is not set. Check your config file.")
   }
}

guide_name = "source_name"

working_dir = params.outputs

features_o = "${params.outputs}/features"
processed_o = "${params.outputs}/processed"
mapped_o = "${params.outputs}/mapped"
counts_o = "${params.outputs}/counts"
model_o = "${params.outputs}/model"
multiqc_o = "${params.outputs}/multi_qc"

sample_names = params.from_sra ? params.sample_names : "sample_id"

log.info pipeline_title + """\
   inputs
      input dir.     : ${params.inputs}
      FASTQ dir.     : ${params.fastq_dir}
      sample sheet   : ${params.sample_sheet}
      sample names   : ${sample_names}
   UMI mode          : ${params.use_umis}
   SRA options
      SRA mode       : ${params.from_sra}
   table to FASTA
      name columns   : ${params.name_column}
      seq columns    : ${params.sequence_column}
   trimming 
      quality        : ${params.trim_qual}
      minimum length : ${params.min_length}
   output
      Features       : ${features_o}
      Processed      : ${processed_o}
      Mapped         : ${mapped_o}
      Counts         : ${counts_o}
      MultiQC        : ${multiqc_o}
   """
   .stripIndent()

dirs_to_make = [features_o, processed_o, 
                counts_o, mapped_o, 
                model_o, multiqc_o]

log.info  """
          Making directories: 
          """.stripIndent()

dirs_to_make.each {
   log.info "${it}: " 
   log.info file(it).mkdirs() ? "OK" : "Cannot create directory: ${it}"
}

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

   Channel.fromPath( "${params.inputs}/${params.sample_sheet}", 
                     checkIfExists: true )
      .splitCsv( header: true, quote: '"', strip: true )
      .set { csv_ch }

   csv_ch
      .map { tuple( it[sample_names],
                     it.adapter_5prime,
                     it.adapter_3prime ) }
      .set { adapter_ch }  // sample_name, adapt5, adapt3 
   csv_ch
      .map { tuple( it[sample_names],
                     it.genome,
                     it.pam,
                     it.scaffold ) }
      .set { genome_pam_ch }

   if ( params.from_sra ) {
     Channel.of( params.ncbi_api_key ).set { ncbi_api_key }
      csv_ch
         .map { tuple( it[sample_names], 
                       it.Run ) }
         .combine( ncbi_api_key )  // sample_id, run_id, api_key
         | PULL_FASTQ_FROM_SRA
      PULL_FASTQ_FROM_SRA.out
         .set { reads_ch }  // sample_id, reads
   } else {
      csv_ch
         .map { tuple( it[sample_names], 
                       file( "${params.fastq_dir}/${it.reads}*.fastq.gz", 
                              checkIfExists: true ) ) }
         .set { reads_ch }  // sample_id, reads
   }

   reads_ch | FASTQC
   
   Channel.value( 
      tuple( params.trim_qual, params.min_length )
   ).set { trim_params } 
   TRIM_CUTADAPT(
      reads_ch
         .combine( adapter_ch, by: 0 ),  // sample_id, reads, adapt5, adapt3
       trim_params
   )

   csv_ch
      .map { it.genome }  // genome_acc
      .unique()
      | DOWNLOAD_GENOME   // genome_acc, genome, gff
   
   if ( params.guides ) {
      csv_ch
         .map { tuple( 
            it.sample_id, 
            it.guides_filename,
            file( "${params.inputs}/${it.guides_filename}", 
                  checkIfExists: true ) 
         ) }
         .set { guide_csv }  // sample_id, guide_filename, guide_file

      Channel.value( 
         tuple( params.sequence_column, params.name_column )
      ).set { table2fasta_params }
      guide_csv
         .map { it[1..2] }  // guide_filename, guide_file
         .unique()
         .combine( table2fasta_params )  // guide_filename, guide_file, seq_col, name_col
         | TABLE2FASTA  // guide_filename, guide_fasta

      guide_csv
         .map { it[1..0] }  // guide_filename, sample_id
         .combine( TABLE2FASTA.out, by: 0 )  // guide_filename, sample_id, guide_fasta
         .map { it[1..-1] }  // sample_id, guide_fasta      
         .set { guide_fasta0 }

      guide_fasta0
         .combine( genome_pam_ch, by: 0 )  // sample_id, guide_fasta, genome_acc, pam, scaffold
         .map { it[2..4] + [ it[1] ] }  // genome_acc, pam, scaffold, guide_fasta
         .unique()
         .combine( DOWNLOAD_GENOME.out, by: 0)  // genome_acc, pam, scaffold, guide_fasta, genome_fasta, genome_gff
         | MAP_GUIDES_TO_FEATURES 
      MAP_GUIDES_TO_FEATURES.out.main
         .set { guide_gff }  // genome_acc, pam, guide_gff
   } else {
      genome_pam_ch
         .map { it[1..-1] }  // genome_acc, pam, scaffold
         .unique()
         .combine( DOWNLOAD_GENOME.out, by: 0 )  // genome_acc, pam, genome_fasta, genome_gff
         .combine( Channel.of( params.guide_length ) )  // genome_acc, pam, genome_fasta, genome_gff, guide_length
         | DESIGN_GUIDES  
      DESIGN_GUIDES.out.main   
         .set { guide_gff }  // genome_acc, pam, guide_gff
   }

   guide_gff | GFF2TABLE  // genome_acc, pam, guide_tsv

   if ( ! params.guides ) {
      GFF2TABLE.out
         .map { tuple( it[-1].getSimpleName(), it[-1], "guide_sequence", "source_name" ) } 
         | TABLE2FASTA  // guide_id, guide_fasta
      GFF2TABLE.out
         .map { [ it[-1].getSimpleName() ] + it[0..1] }  // guide_id, genome_acc, pam
         .combine( TABLE2FASTA.out, by: 0 )  // guide_id, genome_acc, pam, guide_fasta
         .map { it[1..-1] }  // genome_acc, pam, guide_fasta
         .combine( 
            genome_pam_ch
               .map { it[1..2] + [ it[0] ] },  // genome_acc, pam, sample_id
            by: [0,1] 
         )  // genome_acc, pam, guide_fasta, sample_id
         .map { it[-1..-2] }  // sample_id, guide_fasta
         .set { guide_fasta0 }
   }

   if ( params.rc ) {  // reverse complement
      guide_fasta0 | REVERSE_COMPLEMENT
      REVERSE_COMPLEMENT.out.set { guide_fasta }
   }
   else {  // pass through
      guide_fasta0.set { guide_fasta }
   }

   if ( params.use_umis ) {
      csv_ch
         .map { tuple( it[sample_names], it.umi_pattern ) }
         .combine( TRIM_CUTADAPT.out.main, by: 0 )  // sample_id, umi_pattern, reads 
         .set { pre_umi }
      pre_umi | UMITOOLS_EXTRACT  // sample_id, reads
      UMITOOLS_EXTRACT.out.main.set { pre_demux }
   } else {
      TRIM_CUTADAPT.out.set { pre_demux }
   }

   pre_demux  // sample_id, reads
      .combine( guide_fasta, by: 0 )   // sample_id, reads, guide_fasta
      | CUTADAPT_DEMUX  // sample_id, reads
   
   if ( params.use_umis ) {
      CUTADAPT_DEMUX.out.main 
      | FASTQ2TAB  // sample_id, tab
      | UMITOOLS_COUNT_TAB  // sample_id, counts
      UMITOOLS_COUNT_TAB.out.main
      | PLOT_READS_VS_UMIS
      FASTQ2TAB.out 
      | READS_PER_UMI_AND_PER_GUIDE  // sample_id, per_umi, per_guide
      UMITOOLS_COUNT_TAB.out.main
         .set { guide_counts }
   } else {
      CUTADAPT_DEMUX.out.main
      .combine( guide_fasta, by: 0 )  // sample_id, reads, guide_fasta
      | COUNTS_PER_GUIDE  // sample_id, per_guide
      COUNTS_PER_GUIDE.out
         .set { guide_counts }
   }

   GFF2TABLE.out  // genome_acc, pam, guide_tsv
      .map { tuple( it[0], it[2] ) }  // genome_acc, guide_tsv
      .combine( genome_pam_ch.map { it[1..0] }, by: 0 )  // genome_acc, guide_tsv, sample_id
      .map { it[2..1] }  // sample_id, guide_tsv
      .combine( guide_counts, by: 0 )  // sample_id, guide_tsv, counts
      | ANNOTATE_COUNTS_WITH_GENOME_FEATURES

   if ( params.do_fitness ) {
      csv_ch
         .map { 
            [ it[sample_names], it.ref_guide, it.ref_timepoint ] +
            it.findAll { k, v -> k.startsWith( "condition_" ) }
         }
         .set { conditions_ch }
      STACK_JOIN_CONDITIONS( 
         guide_counts
            .map { [ it[1] ] }
            .collect(),
         conditions_ch 
      )  
      | FITNESS
      JOIN_GFF(FITNESS.out, GFF2TABLE.out)
      PLOT_FITNESS(JOIN_GFF.out, essential_ch)
   } 

   TRIM_CUTADAPT.out.multiqc_logs
   .concat(
      CUTADAPT_DEMUX.out.multiqc_logs,
      FASTQC.out.multiqc_logs 
   )
   .flatten()
   .unique()
   .collect() 
   | MULTIQC

}

process DOWNLOAD_GENOME {

   tag "${accession}"
   label 'some_mem'

   input:
   val accession

   output:
   tuple val( accession ), path( "ncbi_dataset/data/*/${accession}_*_genomic.fna" ), path( "ncbi_dataset/data/*/*.gff" )

   script:
   """
   wget "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${accession}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED" -O genome-out
   unzip -o genome-out ncbi_dataset/data/${accession}/{${accession}_*_genomic.fna,*.gff}
   """
}


// Do quality control checks
process FASTQC {

   label 'med_mem'

   tag "${sample_id}"

   input:
   tuple val( sample_id ), path ( reads )

   output:
   tuple val( sample_id ), path ( "*.zip" ), emit: logs
   path "*.zip", emit: multiqc_logs

   script:
   """
   zcat ${reads} > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads ${task.cpus} ${sample_id}.fastq
   rm ${sample_id}.fastq
   """
   stub:
   """
   zcat ${reads} | head -n1000 > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads ${task.cpus} ${sample_id}.fastq
   rm ${sample_id}.fastq
   """

}


// Convert a table of barcodes to a FASTA file for mapping.
process TABLE2FASTA {

   tag "${table}"

   publishDir( features_o, 
               mode: 'copy' )

   input:
   tuple val( table_filename ), path( table ), val( sequence_column ), val( name_column )

   output:
   tuple val( table_filename ), path( "*.fasta" )

   script:
   if ( table.getExtension() == "csv" )
      """
      #cat ${table} | tr \$'\\t' , > ${table.getSimpleName()}.csv
      bioino table2fasta ${table.getSimpleName()}.csv \
         --sequence ${sequence_column} \
         --format CSV \
         --name ${name_column} \
         --output ${table.getSimpleName()}.fasta
      """
   else
      """
      cp ${table} __copied_fasta__.fasta
      """
}

// Reverse complement the FASTA file.
process REVERSE_COMPLEMENT {

   tag "${fasta}"

   publishDir( features_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( fasta )

   output:
   tuple val( sample_id ), path( "*.fasta" )

   script:
   """
   seqtk seq -r ${fasta} > ${fasta.getBaseName()}.rc.fasta
   """
}

// Design guides from scratch.
process DESIGN_GUIDES {
   
   tag "${genome_acc}-${pam}"
   label 'big_time'

   publishDir( features_o, 
               mode: 'copy' )

   input:
   tuple val( genome_acc ), val( pam ), path( genome ), path( gff ), val( guide_length )

   output:
   tuple val( genome_acc ), val( pam ), path( "*.gff" ), emit: main
   tuple val( genome_acc ), val( pam ), path( "*.log" ), emit: logs

   script:
   """
   crispio generate ${fasta} \
      --genome ${genome} \
      --annotations ${gff} \
      --pam ${pam} \
      -o ${genome.getBaseName()}-${pam}-l=${guide_length}.gff
      2> ${genome.getBaseName()}-${pam}-l=${guide_length}.map.log
   """
}


// Map a FASTA of guides to a genome and annotate.
process MAP_GUIDES_TO_FEATURES {
   
   tag "${genome_acc}-${pam}-${scaffold}"
   label 'big_time'

   publishDir( features_o, 
               mode: 'copy' )

   input:
   tuple val( genome_acc ), val( pam ), val( scaffold ), path( guide_fasta ), path( genome ), path( gff )

   output:
   tuple val( genome_acc ), val( pam ), path( "*_mapped.gff" ), emit: main
   tuple val( genome_acc ), val( pam ), path( "*.log" ), emit: logs

   script:
   """
   crispio map ${guide_fasta} \
      --genome ${genome} \
      --annotations ${gff} \
      --pam ${pam} \
      2> ${guide_fasta.getBaseName()}.map.log \
   | crispio featurize \
      --scaffold ${scaffold} \
      > ${guide_fasta.getBaseName()}_mapped.gff
   """
}


// Convert a GFF to a TSV table
process GFF2TABLE {
   
   tag "${genome_acc}-${pam}"

   publishDir( features_o, 
               mode: 'copy' )

   input:
   tuple val( genome_acc ), val( pam ), path( gff )

   output:
   tuple val( genome_acc ), val( pam ), path( "*.tsv" )

   script:
   """
   bioino gff2table ${gff} > ${gff.getSimpleName()}.tsv
   """
}


// Get FASTQ
process PULL_FASTQ_FROM_SRA {

   tag "${sample_id}-${sra_run_id}" 

   label 'big_mem'
   time '24 h'

   input:
   tuple val( sample_id ), val( sra_run_id ), val( ncbi_api_key )

   output:
   tuple val( sample_id ), path( "*.fastq.gz" )

   script:
   """
   NCBI_API_KEY=${ncbi_api_key} \
   fastq-dump \
      --stdout \
      --read-filter pass \
      --split-files ${sra_run_id} \
      | gzip -v --best \
      > ${sample_id}.fastq.gz
   """

   stub:
   """
   NCBI_API_KEY=${ncbi_api_key} \
   fastq-dump \
      --stdout \
      --read-filter pass \
      --split-files ${sra_run_id} \
      | head -n100000 \
      | gzip -v --best \
      > ${sample_name}.fastq.gz
   """
}

// Trim adapters from reads
process TRIM_CUTADAPT {
   
   tag "${sample_id}"

   label 'med_mem'
   time '12h'

   publishDir( processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( reads, stageAs: "???/*" ), val( adapt5 ), val( adapt3 )
   tuple val( trim_qual ), val( min_length )

   output:
   tuple val( sample_id ), path( "*.trimmed.fastq.gz" ), emit: main
   tuple val( sample_id ), path( "*.log" ), emit: logs
   path "*.log", emit: multiqc_logs

   script:
   """
   zcat */*.fastq.gz | gzip --best > ${sample_id}_3prime_R1.fastq.gz

   cutadapt \
		-a '${adapt3}' \
      --no-indels \
      -q ${trim_qual} \
      --nextseq-trim ${trim_qual} \
		--minimum-length ${min_length} \
		--report full \
      --action trim \
		--discard-untrimmed \
		-o ${sample_id}_5prime_R1.fastq.gz \
		${sample_id}_3prime_R?.fastq.gz \
      > ${sample_id}.3prime.cutadapt.log

   cutadapt \
		-g '${adapt5}' \
      --no-indels \
		--report full \
      --action retain \
		--discard-untrimmed \
		-o ${sample_id}_R1.trimmed.fastq.gz \
		${sample_id}_5prime_R?.fastq.gz \
      > ${sample_id}.5prime.cutadapt.log

   rm *_?prime_R?.fastq.gz
   
   """
}

// Extract UMIs
process UMITOOLS_EXTRACT {

   time '6h'

   tag "${sample_id}"

   publishDir( processed_o, 
               mode: 'copy', 
               pattern: "*.extract.log" )
   publishDir( processed_o, 
               mode: 'copy', 
               pattern: "*.extracted.fastq.gz" )

   input:
   tuple val( sample_id ), val( umis ), path( reads )
   
   output:
   tuple val( sample_id ), path( "*.extracted.fastq.gz" ), emit: main
   tuple val( sample_id ), path( "*.log" ), emit: logs

   script:
   """
   umi_tools extract \
		--bc-pattern "${umis}" \
      --extract-method regex \
      --quality-filter-mask ${params.trim_qual} \
      --quality-encoding phred33 \
      --log ${sample_id}.extract.log \
		--stdin ${reads} \
		--stdout ${sample_id}_R1.extracted.fastq.gz 

   """
}

// Trim adapters from reads
process CUTADAPT_DEMUX {

   tag "${sample_id}" 

   label 'big_mem'
   time '6h'

   publishDir( path: mapped_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( reads ), path( fastas )

   output:
   tuple val( sample_id ), path( "*.matched.fastq.gz" ), emit: main
   tuple val( sample_id ), path( "*.log" ), emit: logs
   path "*.log", emit: multiqc_logs
   tuple val( sample_id ), path( "*.unmatched.fastq.gz" ), emit: unmatched

   script:
   def seq_to_append = "TGGG"  // to anchor guides of different lengths
   def qual_to_append = "FFFF" * 2
   """
   APPEND=${seq_to_append}
   RC_APPEND=\$(echo \$APPEND | tr ACGTacgt TGCAtgca | rev)
   PAL_APP=\$APPEND\$RC_APPEND

   awk '/^>/; /^[ATCG]/ { print "'\$PAL_APP'"\$0"'\$PAL_APP'" }' ${fastas} \
      > ${fastas.getSimpleName()}.appended.fasta

   # de-duplicate
   awk '/^[ATCG]/' ${fastas.getSimpleName()}.appended.fasta \
      | sort \
      | uniq -c \
      | awk -F' ' '\$1 > 1 { print \$2 }' \
      > duplicate-seqs.txt

   if [ \$(cat duplicate-seqs.txt | wc -l) -gt 0 ] 
   then
      grep -Fx -B1 \
         --no-group-separator \
         -f duplicate-seqs.txt \
         ${fastas.getSimpleName()}.appended.fasta \
         > duplicate-seqs2.txt
   
      grep -Fvx \
         --no-group-separator \
         -f duplicate-seqs2.txt \
         ${fastas.getSimpleName()}.appended.fasta \
         > ${fastas.getSimpleName()}.appended-dedup.fasta
   else
      ln -s ${fastas.getSimpleName()}.appended.fasta \
         ${fastas.getSimpleName()}.appended-dedup.fasta
   fi

   zcat ${reads} \
      | awk '((NR + 3) % 4 == 0 || (NR + 3) % 4 == 2); (NR + 3) % 4 == 1 { print "'\$PAL_APP'"\$0"'\$PAL_APP'" }; (NR + 3) % 4 == 3 { print "${qual_to_append}"\$0"${qual_to_append}" }' \
      | gzip --best \
      > ${reads.getSimpleName()}.appended.fastq.gz

   cutadapt \
		-g '^file:${fastas.getSimpleName()}.appended-dedup.fasta' \
      -e 1 \
      -j 1 \
      --no-indels \
		--report full \
      --action lowercase \
      --rename '{id} {adapter_name} {comment}' \
		--untrimmed-output ${sample_id}.unmatched.fastq.gz \
		-o ${sample_id}.matched.fastq.gz \
		${reads.getSimpleName()}.appended.fastq.gz \
      > ${sample_id}.matched.cutadapt.log

   """
}

process FASTQ2TAB {

   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( fastqs )

   output:
   tuple val( sample_id ), path( "*.tab.tsv" )

   script:
   """
   zcat ${fastqs[0]} \
      | awk '(NR + 3) % 4 == 0' \
      | tr ' ' \$'\t' \
      | cut -f1-2 \
      | sort -k2 \
      > ${sample_id}.tab.tsv
   
   ## Hack because of bug in `umitools count_tab`. It expects read_id_UMI_CB 
   ## when `umitools extract` makes read_id_CB_UMI (!!!)
   #f=${sample_id}.tab0.tsv
   #paste <(paste -d_ <(cut -d_ -f1 \$f) <(cut -f1 \$f | cut -d_ -f3) <(cut -d_ -f2 \$f)) \
   #   <(cut -f2 \$f) \
   #   > ${sample_id}.tab.tsv
   #rm \$f
   """
}


// Count unique UMIs per cell per gene
process UMITOOLS_COUNT_TAB {

   tag "${sample_id}"

   label 'big_mem'
   time '48h'

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( tabfile )

   output:
   tuple val( sample_id ), path( "*.umitools_count.tsv" ), emit: main
   tuple val( sample_id ), path( "*.umitools_count.log" ), emit: logs

   script:
   """
   umi_tools count_tab \
      --method unique \
		--stdin ${tabfile} \
      --stdout ${sample_id}.umitools_count0.tsv \
      --log ${sample_id}.umitools_count.log

   ## Hack for older versions of UMI-tools
   #tail -n+2 ${sample_id}.umitools_count0.tsv \
   #   | sed 's/^b'\\''//;s/'\\''\\t/\\t/' \
   #   > ${sample_id}.umitools_count0.tsv.tail
   NLINES=\$(tail -n+2 ${sample_id}.umitools_count0.tsv | wc -l)

   printf 'guide_name\\tumi_count\\tsample_id\\n' \
      > ${sample_id}.umitools_count-a.tsv
   paste \
      <(tail -n+2 ${sample_id}.umitools_count0.tsv) \
      <(yes ${sample_id} | head -n \$NLINES) \
      | sort -k1 \
      >> ${sample_id}.umitools_count-a.tsv

   cut -f2 ${tabfile} \
      > ${sample_id}.read_count0.tsv

   printf 'guide_name\\tread_count\\n' \
      > ${sample_id}.read_count.tsv
   sort ${sample_id}.read_count0.tsv \
      | uniq -c \
      | awk -F' ' -v OFS=\$'\\t' '{ print \$2,\$1 }' \
      | sort -k1 \
      >> ${sample_id}.read_count.tsv

   join --header ${sample_id}.umitools_count-a.tsv ${sample_id}.read_count.tsv \
      | awk -F' ' -v OFS=\$'\\t' '{ print \$3,\$1,\$2,\$4 }' \
      > ${sample_id}.umitools_count.tsv

   rm ${sample_id}.read_count0.tsv
   """
}

process PLOT_READS_VS_UMIS {

   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( guide_umi_counts )

   output:
   tuple val( sample_id ), path( "*.png" )

   script:
   """
   #!/usr/bin/env python

   from carabiner.mpl import figsaver, scattergrid
   import pandas as pd
   
   df = pd.read_csv(
      "${guide_umi_counts}", 
      sep='\\t',
   ) 

   fig, axes = scattergrid(
      df,
      grid_columns=["read_count", "umi_count"],
      log=["read_count", "umi_count"]
   )
   figsaver()(
      fig=fig,
      name='${sample_id}.umi-vs-reads',
   )
   
   """
}


process READS_PER_UMI_AND_PER_GUIDE {

   tag "${sample_id}"

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( tabfile )

   output:
   tuple val( sample_id ), path( "*.umi.tsv" ), path( "*.guide_name.tsv" )

   script:
   """
   cut -f1 ${tabfile} | cut -d _ -f2 > umi.tsv
   cut -f2 ${tabfile} > guide_name.tsv

   for f in umi.tsv guide_name.tsv
   do
      BASENAME=\$(basename \$f .tsv)
      NLINES=\$(cat \$f | wc -l)

      printf 'sample_id\\t'\$BASENAME'\\t'\$BASENAME'_read_count\\n' \
         > ${sample_id}.\$BASENAME.tsv

      sort \$f | uniq -c \
         | awk -F' ' -v OFS=\$'\\t' '{ print "${sample_id}", \$2, \$1 }' \
         | sort -k3 -n \
         >> ${sample_id}.\$BASENAME.tsv
   done

   """
}


process COUNTS_PER_GUIDE {

   tag "${sample_id}"

   // publishDir( counts_o, 
   //             mode: 'copy' )

   input:
   tuple val( sample_id ), path( fastqs ), path( fastas )

   output:
   tuple val( sample_id ), path( "*.counts.tsv" )

   script:
   """
   printf '${guide_name}\\tsample_id\\tguide_count\\n' \
      > ${sample_id}.counts.tsv

   zcat ${fastqs} \
      | awk '(NR + 3) % 4 == 0' \
      | tr ' ' \$'\\t' \
      | cut -f2 \
      | sort -k1 \
      | uniq -c \
      | awk -F' ' -v OFS=\$'\\t' '{ print \$2,"${sample_name}",\$1 }' \
      | sort -k1 \
      > ${sample_id}.counts0.tsv

   grep '^>' ${fastas} \
      | cut -d'>' -f2 \
      | tr -d ' ' \
      | sort -u \
      > guide-names.txt

   join -j 1 -a 1 -t\$'\\t' \
      guide-names.txt ${sample_id}.counts0.tsv  \
      | tr ' ' \$'\\t' \
      | awk -F\$'\\t' -v OFS=\$'\\t' 'NF==1 { print \$0,"${sample_id}",0 }; NF>1' \
      | sort -k4 -n \
      >> ${sample_id}.counts.tsv

   """
}

// stack count TSV files and merge the condition table
process STACK_JOIN_CONDITIONS {

   tag "${conditions}"

   label 'med_mem'
   time '24h'

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   path 'counts*.in.tsv'
   path conditions

   output:
   path "counts.tsv" 

   script:
   """
   for f in counts*.in.tsv
   do
      cat \$f | python ${projectDir}/bin/join.py ${conditions} "," > \$f.joined.tsv
   done

   cat <(head -n 1 counts1.in.tsv.joined.tsv) <(tail -q -n +2 counts*.joined.tsv) > counts.tsv
   """
}

// merge the guide table
process ANNOTATE_COUNTS_WITH_GENOME_FEATURES {
   
   tag "${sample_id}"
   label 'med_mem'

   publishDir( counts_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( guide_tsv ), path( counts )

   output:
   tuple val( sample_id ), path( "*.tsv" )

   script:
   """
   GUIDE_NAME="source_name"
   GUIDE_NAME_COL=\$(head -n1 ${guide_tsv} | tr \$'\\t' \$'\\n' | grep -n \$GUIDE_NAME  | cut -d: -f 1)

   cat ${guide_tsv} \
   | cut -f10-21,\$GUIDE_NAME_COL \
   | sed 's/'\$GUIDE_NAME'/guide_name/' \
   > ${guide_tsv}.mini
   cat ${counts} \
   | python ${projectDir}/bin/join.py ${guide_tsv}.mini \
   > ${counts.getSimpleName()}-annotated.tsv
   """
}

// merge the guide table
process JOIN_GFF {
   tag "${gfftable}"

   label 'med_mem'

   publishDir(counts_o, mode: 'copy')

   input:
   path fit_params 
   path fitted 
   path gfftable 

   output:
   path( "*params*.tsv", includeInputs: true )
   path "*-fit-annotated.tsv"

   script:
   """
   cat fitness_params-guide_name.tsv \
      | python ${projectDir}/bin/join.py ${gfftable} \
      > fitness_params-guide_name-annotated.tsv
   cat ${fitted} \
      | python ${projectDir}/bin/join.py ${gfftable} \
      > ${fitted.getBaseName()}-annotated.tsv
   """
}


// Use `crispio` to calculate fitness 
process FITNESS {

   tag "${counts}" 

   label 'big_gpu'
   time '24h'

   input:
   path counts 

   output:
   path "fitness_params-*.tsv"
   path "fitness_*-fit.tsv"

   script:
   """
   guidefitness ${counts} \
      --sequencing_group "${params.sequencing_group}" \
      --expansion_group "${params.expansion_group}" \
      --culture "${params.culture_group}" \
      --reference "${params.reference}" \
      --initial "${params.initial}" \
      --name ${guide_name} \
      --count guide_count \
      --format TSV \
      -o fitness
   """
}


// Use `crispin` to plot fitness 
process PLOT_FITNESS {
   tag{"${fitted}"}

   label 'big_mem'

   publishDir( model_o, 
               mode: 'copy' )

   input:
   path fit_params 
   path fitted 
   path essentials 

   output:
   path "*.png"

   script:
   """
   guideplot \
      --fitness fitness_params-guide_name-annotated.tsv \
      --expansion fitness_params-exp_group.tsv \
      --essentials ${essentials} \
      --essential_calls ${params.essential_call} \
      --essential_scores ${params.essential_score} \
      --fitted ${fitted} \
      --reference "${params.reference}" \
      --initial "${params.initial}" \
      --control_column "${params.control_column}" \
      --negative ${params.negative} \
      --count guide_count \
      --format TSV \
      -o fitness
   """
}

// Make log report
process MULTIQC {

   publishDir( multiqc_o, 
               mode: 'copy' )

   input:
   path '*'

   output:
   tuple path( "*.html" ), path( "multiqc_data" )

   script:
   """
   multiqc .
   """
}

/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/