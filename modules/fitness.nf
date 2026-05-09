process Bartab_fit {

   tag "${id}"
   label 'big_time'

   publishDir(
      "${params.outputs}/fitness", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( counts )
   path sample_sheet
   path growth
   val guide_name
   val reference_guide
   val use_umis
   val use_spike
   val timepoint_column
   val concentration_column

   output:
   tuple val( id ), path( "bartab.h5ad" ), emit: h5ad
   tuple val( id ), path( "bartab-fit.csv" ), emit: table
   path "*.log", emit: logs

   script:
   """
   bartab fit "${counts}" \
      --sample-sheet "${sample_sheet}" \
      --barcode-sheet strain_meta.csv \
      --reference "${reference_guide}" \
      --spike-name spike \
      ${use_spike ? "--use-spike" : "--growth ${growth}"} \
      --barcode-column "${guide_name}" \
      --sample-column sample_id \
      --culture-column culture_id \
      --count-column ${use_umis ? "umi_count" : "read_count"} \
      --timepoint-column "${timepoint_column}" \
      ${concentration_column ? "--concentration-column ${concentration_column} --model-type HillFitnessModel" : "--model-type WLS"} \
      --output results.h5ad \
   2> fitness.log

   """
}

process Bartab_plot {

   tag "${id}"

   publishDir(
      "${params.outputs}/fitness/plots", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( results )
   val concentration_column
   val control_guides

   output:
   tuple val( id ), path( "*.png" ), emit: plots
   path "*.log", emit: logs

   script:
   """
   bartab plot "${results}"
      --output bartab \
      ${control_guides ? "--highlight control_guides" : ""} \
      --model-type ${concentration_column ? "HillFitnessModel" : "WLS"}  \
      --plot-format png \
   2> bartab-plot.log

   """
}
