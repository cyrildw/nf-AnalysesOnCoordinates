// nextflow run -c nf-AnalysesOnCoordinates.config  --name "Testing" nf-AnalysesOnCoordinates.nf 
Channel
   .fromPath(params.input_design)
   .splitCsv(header:true, sep:';')
   .map { row -> [ row.LibName,  
                    file("$params.input_dir/$row.LibBam", checkIfExists: true),
                    file("$params.input_dir/${row.LibBam}.bai", checkIfExists: true),
                    file("$params.input_dir/$row.LibBW", checkIfExists: true),
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
   .into { design_bigwig_csv; ch_before_dt_lib}

Channel
    .fromPath(params.bed_design)
    .splitCsv(header:true, sep:";")
    .map { row -> [row.BedName,
		file("$params.input_dir/$row.BedFile", checkIfExists: true),
		row.BedPref,
		row.BedFls,
		row.BedExts,
		row.BedExtls,
		row.BedExtvs ]
        }
    .into { design_bed_csv; ch_before_dt_bed }
/* Macs analyses contains : 
    -Channel split as control or sample
    -Channel cross to have sample & control on same channel emission
    -Macs2 peak calling
    -Add new peaks to the bed files to analyse
    -Merge peaks within a w-nt window size
    -Group peaks with
        - 10% most enriched
        - 50% most enriched

if(params.macs2_analyses){

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
}
*/
/*Deeptools analyses contains :
- MultiBigwigSymmary & PlotCorrelation
    -Should plot correlation & Scatterplot
    -see for --removeOutliers
- ComputeMatrix & PlotHeatmap
    -For ComputeMatrix
        -generate bedfiles from groups
        -take scale-region/reference-point into account.

    -For PlotHeatmap: 
        -account for plotType (lines, fill, se, std)
        --linesAtTickMarks ?

*/


/*process toto {
    tag "$BedName"
    echo true
    input:
    val Labels from ch_dt_labels
    val Files from ch_dt_files
    tuple BedName, file(BedFile), BedPref, BedFls, BedExts, BedExtls, BedExtvs from design_bed_csv.take(1)
    """
    echo "${BedName} \n ${BedFile} \n ${Labels.join(' ')} \n ${Files.join(' ')}
    """
}
*/
if(params.deeptools_analyses){

/* Deeptools process requires
    -bw & Bed files for computation
    -labels & BedName for plotting
*/
    ch_before_dt_bed
        .into{ ch_dt_bed_multiBWsummary; ch_dt_bed_computeMatrix}

    ch_before_dt_lib.map {it -> [ it[0], it[3]]}
        .multiMap { it ->
                    labels: it[0]
                    files: it[1]
                }
        .set{ch_dt_input}

    ch_dt_input.labels.collect()
        .into {ch_dt_labels_plotCor; ch_dt_labels_plotHeatmap}

    ch_dt_input.files.collect()
        .into{ch_dt_files_multiBWsummary; ch_dt_files_computeMatrix}
ch_dt_files_computeMatrix.view()
    process dt_MultiBWsummary {
        tag "$BedName"
        label "multiCpu"
        publishDir "${params.outdir}/DeeptoolsData", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.npz')) "./$filename"
            else null
        }
        input:
        tuple BedName, file(BedFile), BedPref, BedFls, BedExts, BedExtls, BedExtvs from ch_dt_bed_multiBWsummary.take(1)
        file(Files) from ch_dt_files_multiBWsummary
        output:
        "dt_MultiBWsummary.Matix.$BedName.npz" into ch_multibw_matrix //the computed matrix
        val(BedName) into ch_multibw_bedname
        script:
        """
        multiBigwigSummary BED-file \
        -b ${Files.join(' ')} \
        -o dt_MultiBWsummary.Matrix.${BedName}.npz \
        --BED ${BedFile} \
        --numberOfProcessors ${task.cpus}
        """
    }
    process dt_PlotCorrelation {
        tag "$BedName"
        publishDir "${params.outdir}/", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.tab')) "./DeeptoolsData/$filename"
            else if (filename.endsWith('.pdf')) "./DeeptoolsFigures/$filename"
            else null
        }
        input:
        file(Matrix) from ch_multibw_matrix
        val Labels from ch_dt_labels_plotCor
        val(BedName) from ch_multibw_bedname
        output:
        "dt_MultiBWsummary.CorTable.${BedName}.tab"
        "Heatmap.dt_MultiBWsummary.${BedName}.pdf"
        script:
        """
        plotCorrelation \
        --corData ${Matrix} \
        --corMethod spearman \
        --whatToPlot heatmap \
        --plotNumbers \
        --outFileCorMatrix dt_MultiBWsummary.CorTable.${BedName}.tab \
        -o Heatmap.dt_MultiBWsummary.${BedName}.pdf \
        --plotTitle ${BedName} \
        --labels ${Labels.join(' ')}
        """
    }}/*
    // Adujst for reference point
    
    process dt_ComputeMatrix {
        tag "$BedName"
        label "multiCpu"
        publishDir "${params.outdir}/DeeptoolsData", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.gz')) "./$filename"
            else null
        }
        input:
        output:
        script:
        """
        computeMatrix scale-regions \
        -S Array with bigwig files \
        -R ${BedFile} \
        -b amount of nt before \
        -a amount of nt after \
        -m size of the middle region \
        --skipZeros \
        -p ${task.cpus} \
        -o dt_ComputeMatrix.${BedName}.gz
        """
    }
    // Adjust for :
    //    - groups
    //    - kmean
    
    
    process dt_PlotHeatmap {
        tag "$BedName"
        publishDir "${params.outdir}/", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.pdf')) "./DeeptoolsFigures/$filename"
            else null
        }
        input:
        matrix from computMatrix_ch
        output:
        "Heatmap.dt_PlotHeatmap.${BedName}.pdf"
        script:
        """
        plotHeatmap \
        --matrixFile ${matrix} \
        -o Heatmap.dt_PlotHeatmap.${BedName}.pdf \
        --startLabel '1' \
        --endLabel '0'
        --yMin 0 \
        --xAxisLabel ${BedName} \
        --samplesLabel Array with Labels
        """
    }
}
*/

/*R analyses includes : 
- Combine all libraries with all bedfiles
- Getting Tag Density over the bed files
- Converting all density tables to R object with scaling
- Producing combined graphics for
    -all elements
    -grouped elements
    -quantiles
- Outputing R objects (and R scripts ?)

if(params.r_analyses){
    process tag_density {
        tag "$LibName"
        input:
        combinedLib_and_Bed from some_channel
        file(r_scaling_function.R) from ${params.r_scaling}
        output:
        file(temp_file)
        file("r_file_2_run.R")
        file("${LibName}.${BedName}.R")
        
        script:
        """
        get_tag_density -f ${BwFile} ${BedFile} | awk '{print \$4"\\t"\$6"\\t"\$7}' - > temp_file
        echo "#!/usr/bin/env Rscript
        source ${r_scaling_function.R}
        finalL=${BedFinalLength}
        ext=${BedExtension}
        extLL=${BedExtLengthLeft};extLR=${BedExtLengthRight};
        extVL=${BedExtValLeft};extVR=${BedExtValRight};
        t=as.list(read.table(temp_file, stringsAsFactor=FALSE, header=FALSE))
        names(t)=c('Q_id', 'Strand', 'PerBP')
        t[['Splt_PerBP']]=lapply(as.character(t[['PerBP']]), function(x) as.numeric(strsplit(x,';')[[1]]))
        for(k in 1:length(t[['Strand']])){ 
            if(as.character(t[['Strand']][k])=='-'){t[['Splt_PerBP']][[k]]=rev(t[["Splt_PerBP"]][[k]])            }
        } # Reversing order for minus strand
        t_scaled=c()
        t_scaled=rbind(t_scaled, sapply(t[['Splt_PerBP']], function(y) Scale_Vector(Data=y,FinalLength=finalL, Extention=ext, Ext_length=c(extLL, extLR), Ext_value=c(extVL, extVR))))
        colnames(t_scaled)=t[['Q_id']]
        save(x=t_scaled, file='${LibName}.${BedName}.R')" > r_file_2_run.R
        R --vanilla r_file_2_run.R
        """
    }

    /* Must get all bed-associated R file (containing data)
    process combine_R {

    }

}
*/
/*
testbw_ch
    .map { it -> [ name:it[0], name:it[1] ]}
    .collect()
    .view()

*/

/*
if(testing_combine_cross_etc){
 Combine works to have a channel with both Bed & Bw
testbed_ch
    .combine(testbw_ch)
    .view()

*/

/* Filtering files that are not control datasets
    applying array modification (.map)
    grouping per "LibIsControl"
    setting into new channel */
/*ch_control_bam
    .filter { it[14] != "" && it[15] == "" }
    .map { it -> [ it[14], it[0],it[1], it[2]]}
    .groupTuple(by: 0)
    .set { ch_macs2_control_bam } 
*/
/* Filtering out files that are control datasets
    applying same organization as ch_macs2_control_bam
    setting into new channel */
/*ch_sample_bam
    .filter { it[14] == "" && it[15] != "" } // && it[15].ifEmpty()
    .map { it -> [ it[15], it[0], it[1], it[2] ] } //.view()
    .set { ch_macs2_sample_bam }
*/

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

 /*ch_macs2_control_bam
   .cross(ch_macs2_sample_bam)
   .set { ch_macs2_run }
 */
   //.map {it.flatten()}
   //.view()
   /*.transpose().map{ it -> [it[1]]}.view()
  .map { it.flatten() }
   .map { it -> [ it[1], it[2], it[-1] ] }.view()
 }*/
