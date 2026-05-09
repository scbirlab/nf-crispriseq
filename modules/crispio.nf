// Design guides from scratch.
process design_guides_with_crispio {
   
   tag "${id}:${pam}"
   label 'big_time'

   publishDir( 
      "${params.outputs}/guides", 
      mode: 'copy',
      saveAs: { "${id}.${pam}-l=${guide_length}.${it}" },
   )

   input:
   tuple val( id ), val( pam ), path( genome ), path( gff ), val( guide_length )

   output:
   tuple val( id ), val( pam ), path( "guide-design.gff" ), emit: main
   path "guide-design.log", emit: logs

   script:
   """
   crispio generate "${genome}" \
      --annotations "${gff}" \
      --pam ${pam} \
      -o guide-design.gff \
      2> guide-design.log

   """
}


// Map a FASTA of guides to a genome and annotate.
process map_guides_to_genome_features {
   
   tag "${id}:${pam}:${scaffold}"
   label 'big_time'

   publishDir( 
      "${params.outputs}/guides", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), val( pam ), val( scaffold ), path( guide_fasta ), path( genome_fasta ), path( gff )

   output:
   tuple val( id ), val( pam ), path( "mapped.gff" ), emit: main
   path "map.log", emit: logs

   script:
   """
   crispio map "${guide_fasta}" \
      --genome "${genome_fasta}" \
      --annotations "${gff}" \
      --pam "${pam}" \
   2> map.log \
   | crispio featurize \
      --scaffold "${scaffold}" \
   > mapped.gff

   """
}
