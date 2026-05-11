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
   val growth
   val guide_name
   val reference_guide
   val use_umis
   val use_spike
   val timepoint_column
   val concentration_column
   val culture_column
   val growth_type

   output:
   tuple val( id ), path( "bartab.h5ad" ), emit: h5ad
   tuple val( id ), path( "bartab.csv" ), emit: table
   path "*.log", emit: logs

   script:
   """
   python -c "
   import pandas as pd
   df = (
      pd.read_csv('${counts}', sep='\\t')
      .assign(**{
         '${guide_name}': lambda x: x['${guide_name}'].fillna('unmapped-' + x['guide_name'])
      })
   )
   df[['${guide_name}']].drop_duplicates().to_csv('strain-meta.csv', index=False)
   (
      df
      #.groupby(['${guide_name}', 'sample_id'], as_index=False)
      #['${use_umis ? "umi_count" : "read_count"}']
      #.mean()
      .to_csv('counts_deduped.tsv.gz', sep='\\t', index=False)
   )
   "

   bartab fit "counts_deduped.tsv.gz" \
      --sample-sheet "${sample_sheet}" \
      --barcode-sheet strain-meta.csv \
      --reference "${reference_guide}" \
      --spike-name spike \
      ${use_spike ? "--use-spike" : "--growth-column ${growth} --growth-type ${growth_type}"} \
      --barcode-column "${guide_name}" \
      --sample-column sample_id \
      --culture-column ${culture_column} \
      --count-column ${use_umis ? "umi_count" : "read_count"} \
      --timepoint-column "${timepoint_column}" \
      ${concentration_column ? "--concentration-column ${concentration_column} --model-type HillFitnessModel" : "--model-type WLS"} \
      --output bartab.h5ad \
   2> >(tee fitness.log >&2)

   """
}

process Bartab_plot {

   tag "${id}"
   label 'big_time'

   publishDir(
      "${params.outputs}/fitness/plots", 
      mode: 'copy',
      saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( results )
   val concentration_column
   val control_guides
   val highlight_guides

   output:
   tuple val( id ), path( "*.png" ), emit: plots
   path "*.log", emit: logs

   script:
   """
   bartab plot "${results}" \
      --output bartab \
      ${control_guides ? "--control ${control_guides}" : ""} \
      ${highlight_guides ? "--highlight ${highlight_guides}" : ""} \
      --model-type ${concentration_column ? "HillFitnessModel" : "WLS"}  \
      --plot-format png \
   2> >(tee bartab-plot.log >&2)
   """
}
