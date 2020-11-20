// nextflow run -c nf-AnalysesOnCoordinates.config  --name "Testing" nf-AnalysesOnCoordinates.nf 
Channel
   .fromPath(params.input_design)
   .splitCsv(header:true, sep:';')
   .map { row -> [ row.LibName,  
                    file("$row.LibBam", checkIfExists: true),
                    file("${row.LibBam}.bai", checkIfExists: true),
                    file("$row.LibBW", checkIfExists: true),
                    row.LibSequenced,
                    row.LibMapped,
                    row.LibUnique,
                    row.LibInsertSize,
                    row.LibQpcrNorm,
                    row.LibType,
                    row.LibProj,
                    row.LibExp,
                    row.LibCondition,
                    row.LibOrder,
                    row.LibIsControl,
                    row.LibControl ]
                    }
   .into { design_bigwig_csv; testbw_ch; ch_control_bam; ch_sample_bam }

Channel
    .fromPath(params.bed_design)
    .splitCsv(header:true, sep:";")
    .map { row -> [row.BedName,
		file("$row.BedFile", checkIfExists: true),
		row.BedPref,
		row.BedFls,
		row.BedExts,
		row.BedExtls,
		row.BedExtvs ]
        }
    .into { design_bed_csv; testbed_ch }

/* Combine works to have a channel with both Bed & Bw
testbed_ch
    .combine(testbw_ch)
    .view()

*/

/* Filtering files that are not control datasets
    applying array modification (.map)
    grouping per "LibIsControl"
    setting into new channel */
ch_control_bam
    .filter { it[14] != "" && it[15] == "" }
    .map { it -> [ it[14], it[0],it[1], it[2]]}
    .groupTuple(by: 0)
    .set { ch_macs2_control_bam } 

/* Filtering out files that are control datasets
    applying same organization as ch_macs2_control_bam
    setting into new channel */
ch_sample_bam
    .filter { it[14] == "" && it[15] != "" } // && it[15].ifEmpty()
    .map { it -> [ it[15], it[0], it[1], it[2] ] } //.view()
    .set { ch_macs2_sample_bam }


/* Crossing  control_bam channel with sample_bam channel
    outputs one group of value per association
    - CTRL1 with ChIP1a
    - CRTL1 with ChIP1b
    - CTRL2 with ChIP2*/
/*
Cross     sample.control => 1 output
Cross     control.sample => 2 output [array control, array sample] GOOD !
Combine   control.sample => 4 output not matching
Join      control.sample => 1 output */

 ch_macs2_control_bam
   .cross(ch_macs2_sample_bam)
   .set { ch_macs2_run }
   //.map {it.flatten()}
   //.view()
   /*.transpose().map{ it -> [it[1]]}.view()
  .map { it.flatten() }
   .map { it -> [ it[1], it[2], it[-1] ] }.view()
   */


process macs2_run {
    label 'local'
    echo true
    input:
    tuple CtrlArray, SampleArray from ch_macs2_run
    output:
    stdout into toto_ch
    script:
    def CtrlBam=CtrlArray[2]
    def SampleBam=SampleArray[2]
    def SampleName=SampleArray[1]
    """
    echo "macs2 callpeak \\
        -t ${SampleBam} \\
        -c ${CtrlBam.join(' ')} \\
        -g '${params.macs2_genome}' \\
        -f BAMPE \\
        --outdir $
        -n ${SampleName}"
    """
}
/*
testbw_ch
    .map { it -> [ name:it[0], name:it[1] ]}
    .collect()
    .view()

*/