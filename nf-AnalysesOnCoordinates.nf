// nextflow run -c nf-AnalysesOnCoordinates.config  --name "Testing" nf-AnalysesOnCoordinates.nf 
Channel
   .fromPath(params.input_design)
   .splitCsv(header:true, sep:';')
   .map { row -> [ row.LibName,
                    file("$params.input_dir/$row.LibBam", checkIfExists: false),
                    file("$params.input_dir/${row.LibBam}.bai", checkIfExists: false),
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
   .into { design_bigwig_csv; ch_before_dt_lib; ch_before_R_lib}

Channel
    .fromPath(params.bed_design)
    .splitCsv(header:true, sep:";")
    .map { row -> [row.BedName,
		file("$params.input_dir/$row.BedFile", checkIfExists: true),
        file("$params.input_dir/$row.BedGroupFile"),
		row.BedReferencePoint, //Used by Deeptools in opposition to "scale-region"
		row.BedExtLengthLeft,  //Used by Deeptools and R to extend the bed coordinates upstream
		row.BedExtLengthRight, //Used by Deeptools and R to extend the bed coordinates downstream
		row.BedFinalLength,    //Used by Deeptools and R to set the bed coordinates final length
		row.BedExtension,      // Should the bed coordinates be extended
		row.BedExtValLeft,     // Used by R, how much of the FinalLength should the upstream extension represent
        row.BedExtValRight]     // Used by R, how much of the FinalLength should the upstream extension represent
        }
    .into { design_bed_csv; ch_before_dt_bed; ch_before_R_bed }
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

/*Making groups
 - from bed channel + groupFile
 - groupFile = 1 line per group; 
    GroupName_1;ID;ID;ID;ID
    GroupName_2;ID;ID;ID;ID
 - Creating 1 file per group
    - using IDs to grep in the bed file
    - 1 subBED per group
 
rm groupenames.txt
cat testgroup.txt | while read line; do 
    echo -e ${line//;/"\t"} | \
    awk '{ print "#"$0 > $1".txt"; print $1".txt" >> "groupenames.txt" } ';
done

cat groupenames.txt | while read line;do
    sed "s/\t/\n/g" $line | grep -v "#" | 
    while read id; do 
        grep $id TSS_TES_steinmetz_jacquier.mRNA.bed >> "$line.bed";
    done;
done
*/

if(params.deeptools_analyses){
/* Deeptools process requires :
    -bw & Bed files for computation
    -labels & BedName for plotting
*/
    process bedGroups {
        tag "$BedName"
        container=''
        publishDir "${params.outdir}/${params.name}/", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.bed')) "./BedFiles/$filename"
            else null
        }

        input:
        tuple BedName, file(BedFile), file(BedGroupFile),BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedFinalLength, BedExtension,BedExtValLeft, BedExtValRight from ch_before_dt_bed
        output:
        files "*.txt"
        files "*.bed"
        tuple BedName, file(BedFile), file("${BedName}.GrpFiles.txt"), file("${BedName}.*.bed"),BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedFinalLength into (ch_test,ch_dt_bedGroup_computeMatrix)
        tuple BedName, file(BedFile), BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedFinalLength into (ch_dt_bed_multiBWsummary, ch_dt_bed_computeMatrix)

        script:
        if(BedGroupFile.isFile() && BedGroupFile.size()!=0 ){
            //Creating 1 file per group + 1 file with goupefile-names.
            //Then greping ids from each group file into the BED file to produce 1 GroupBedFile/group.
            """
            cat ${BedGroupFile} | while read line; do 
                echo -e \${line//;/"\\t"} | \
                awk '{ print "#"\$0 > \$1".txt"; print \$1".txt" >> "${BedName}.GrpFiles.txt" } ';
            done

            cat ${BedName}.GrpFiles.txt | while read line;do
                sed "s/\\t/\\n/g" \$line | grep -v "#" | 
                while read id; do 
                    grep \$id ${BedFile} >> "${BedName}.\$line.bed";
                done;
            done
            """
        }
        else{
            """
            touch ${BedName}.GrpFiles.txt ${BedName}.nogroup.bed 
            """
            
        }
    }
    //ch_test.view()

    //ch_before_dt_bed
    //    .into{ ch_dt_bed_multiBWsummary; ch_dt_bed_computeMatrix}

    ch_before_dt_lib.map {it -> [ it[0], it[3]]} // From the bigwig channel : sending names into a labels-subchannel and sending files into a files-subchannel
        .multiMap { it ->
                    labels: it[0]
                    files: it[1]
                }
        .set{ch_dt_input}

    ch_dt_input.labels.collect() //Using the collect() method to have a single list of values  which will be used in the next processes
        .into {ch_dt_labels_plotCor; ch_dt_labels_plotHeatmap; ch_dt_labels_groupHeatmap; test3_ch}
    ch_dt_input.files.collect() //Using the collect() method to have a single list of values which will be used in the next processes
        .into{ch_dt_files_multiBWsummary; ch_dt_files_computeMatrix; ch_dt_files_groupcomputeMatrix; test4_ch}


    process dt_MultiBWsummary {
        // Computes the npz Matrix file required for the dt_PlotCorrelation process.
        // The npz file is saved in the DeeptoolsData directory
        tag "$BedName"
        label "multiCpu"
        publishDir "${params.outdir}/${params.name}/DeeptoolsData", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.npz')) "./$filename"
            else null
        }
        input:
        tuple BedName, file(BedFile), BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedFinalLength from ch_dt_bed_multiBWsummary
        file(Files) from ch_dt_files_multiBWsummary
        output:
        file("dt_MultiBWsummary.Matrix.${BedName}.npz") into ch_multibw_matrix //the computed matrix
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
        // Produce the Correlation matrix graphical output.
        tag "$BedName"
        publishDir "${params.outdir}/${params.name}/", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.tab')) "./DeeptoolsData/$filename"
            else if (filename.endsWith('.pdf')) "./DeeptoolsFigures/$filename"
            else null
        }
        input:
        file(Matrix) from ch_multibw_matrix
        val(Labels) from ch_dt_labels_plotCor
        val(BedName) from ch_multibw_bedname
        output:
        file("dt_MultiBWsummary.CorTable.${BedName}.tab")
        file("Heatmap.dt_MultiBWsummary.${BedName}.pdf")
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
    }

    
    process dt_ComputeMatrix {
        // Compute the matrix file (.gz)
        tag "$BedName"
        label "multiCpu"
        publishDir "${params.outdir}/${params.name}/DeeptoolsData", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.gz')) "./$filename"
            else null
        }

        input:
        tuple BedName, file(BedFile), BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedFinalLength  from ch_dt_bed_computeMatrix
        file(Files) from ch_dt_files_computeMatrix
        output:
        file("dt_ComputeMatrix.${BedName}.gz") into ch_computeMatrix_matrix //the computed matrix
        val(BedName) into ch_computeMatrix_bedname
        script:
        
        if(BedReferencePoint=='false')
            """
            computeMatrix scale-regions \
            -S ${Files.join(' ')} \
            -R ${BedFile} \
            -b ${BedExtLengthLeft} \
            -a ${BedExtLengthRight} \
            -m ${BedFinalLength} \
            --skipZeros \
            -p ${task.cpus} \
            -o dt_ComputeMatrix.${BedName}.gz
            """
                
        else
            """
            computeMatrix reference-point \
            -S ${Files.join(' ')} \
            -R ${BedFile} \
            -b ${BedExtLengthLeft} \
            -a ${BedExtLengthRight} \
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
        publishDir "${params.outdir}/${params.name}/", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.pdf')) "./DeeptoolsFigures/$filename"
            else null
        }
        input:
        file(Matrix) from ch_computeMatrix_matrix
        val(Labels) from ch_dt_labels_plotHeatmap
        val(BedName) from ch_computeMatrix_bedname
        output:
        file("Heatmap.dt_PlotHeatmap.${BedName}.pdf")
        script:
        """
        plotHeatmap \
        --matrixFile ${Matrix} \
        -o Heatmap.dt_PlotHeatmap.${BedName}.pdf \
        --startLabel 'st' \
        --endLabel 'end' \
        --yMin 0 \
        --xAxisLabel ${BedName} \
        --samplesLabel ${Labels.join(' ')}
        """
    }
    
    //If groups have been mentionned then produce heatmaps per groups

    process dt_Group_ComputeMatrix {
        tag "$BedName"
        label "multiCpu"
        echo true
        publishDir "${params.outdir}/${params.name}/", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.gz')) "./DeeptoolsData/$filename"
            else if (filename.endsWith('.pdf')) "./DeeptoolsFigures/$filename"
            else null
        }
        input:
        tuple BedName, file(BedFile),file(BedGrpFile), file(BedGrpBedFiles), BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedFinalLength from ch_dt_bedGroup_computeMatrix
        file(Files) from ch_dt_files_groupcomputeMatrix
        val(Labels) from ch_dt_labels_groupHeatmap
        output:
        file("dt_ComputeMatrix.Group.${BedName}.gz")//the computed matrix
        file("Heatmap.dt_PlotHeatmap.Group.${BedName}.pdf")//the heatmap
        //val(BedName) into ch_computeMatrix_bedname
        
        when:
        BedGrpFile.size() != 0
        script:
        
        if(BedReferencePoint=='false')
            """
            computeMatrix scale-regions \
            -S ${Files.join(' ')} \
            -R ${BedGrpBedFiles.join(' ')} \
            -b ${BedExtLengthLeft} \
            -a ${BedExtLengthRight} \
            -m ${BedFinalLength} \
            --skipZeros \
            -p ${task.cpus} \
            -o dt_ComputeMatrix.Group.${BedName}.gz
            """
        
        
        else
            """
            computeMatrix reference-point \
            -S ${Files.join(' ')} \
            -R ${BedGrpBedFiles.join(' ')} \
            -b ${BedExtLengthLeft} \
            -a ${BedExtLengthRight} \
            --skipZeros \
            -p ${task.cpus} \
            -o dt_ComputeMatrix.Group.${BedName}.gz
            """
        

        """
        plotHeatmap \
        --matrixFile dt_ComputeMatrix.Group.${BedName}.gz \
        -o Heatmap.dt_PlotHeatmap.Group.${BedName}.pdf \
        --startLabel '-${BedExtLengthLeft}' \
        --endLabel '${BedExtLengthRight}' \
        --yMin 0 \
        --xAxisLabel ${BedName} \
        --samplesLabel ${Labels.join(' ')}
        """
       
    }
}


/*R analyses includes : 
TEST    - Create the bed files if it requires some bp extension.
OK      - Combine all libraries with all bedfiles
OK      - Getting Tag Density over the bed files
TEST    - Converting all density tables to R object with scaling
- Producing combined graphics for
    -all elements
    -grouped elements
    -quantiles
- Outputing R objects (and R scripts ?)*/

//r_func=Channel.fromPath(params.r_scaling)





if(params.r_analyses){

    /*process create_bed_with_ext {
        tag "$BedName:$BedExtension-$BedExtLengthLeft:$BedExtLengthRight"
        input:
        tuple BedName, file(BedFile),file(BedGrpFile), BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedFinalLength, BedExtension, BedExtValLeft,BedExtValRight from ch_before_R_bed
        output:
        tuple BedName, file("${BedFile.baseName}.ext.bed"),file(BedGrpFile), BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedFinalLength, BedExtension, BedExtValLeft,BedExtValRight into ch_for_R_ext_Bed
        file("${BedFile.baseName}.ext.bed")

        script:
        if(BedExtension=='true')
            """
            awk '{print \$1"\\t"\$2-${BedExtLengthLeft}"\\t"\$3+${BedExtLengthRight}"\\t"\$4"\\t"\$5"\\t"\$6}' ${BedFile} > ${BedFile.baseName}.ext.bed
            """
        else
            """
            cp ${BedFile}  ${BedFile.baseName}.ext.bed
            """
        
    }


    ch_for_R_ext_Bed.combine(ch_before_R_lib) // This combines the Bed channel with the lib channel.
        .set{ch_R_TD}

    //ch_R_test.view()
    */
     ch_before_R_bed.combine(ch_before_R_lib) // This combines the Bed channel with the lib channel.
        .set{ch_R_TD}

    process tag_density {
        tag "$LibName - $BedName"
        input:
        tuple BedName, file(BedFile),file(BedGrpFile), BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedFinalLength, BedExtension, BedExtValLeft, BedExtValRight,
            LibName, file(LibBam), file(LibBai), file(LibBW), LibSequenced, LibMapped, LibUnique, LibInsertSize, LibQpcrNorm, LibType, LibProj, LibExp, LibCondition, LibOrder, LibIsControl, LibControl   from ch_R_TD
        file(r_function) from Channel.fromPath(params.r_scaling)
        output:
        file(temp_file)
        //file("r_file_2_run.R")
        //file("${LibName}.${BedName}.R")
        
        script:
        """
        get_tag_density -f ${LibBW} ${BedFile} | awk '{print \$4"\\t"\$6"\\t"\$7}' - > temp_file
        """
        /*echo "#!/usr/bin/env Rscript
        source('${r_function}')
        finalL=${BedFinalLength}
        ext='${BedExtension}'
        if(ext=='false'){ext=FALSE}
        if(ext=='true'){ext=TRUE}
        extLL=${BedExtLengthLeft};extLR=${BedExtLengthRight};
        extVL=${BedExtValLeft};extVR=${BedExtValRight};
        t=as.list(read.table('temp_file', stringsAsFactor=FALSE, header=TRUE))
        names(t)=c('Q_id', 'Strand', 'PerBP')
        t[['Splt_PerBP']]=lapply(as.character(t[['PerBP']]), function(x) as.numeric(strsplit(x,';')[[1]]))
        for(k in 1:length(t[['Strand']])){ 
            if(as.character(t[['Strand']][k])=='-'){t[['Splt_PerBP']][[k]]=rev(t[['Splt_PerBP']][[k]])            }
        } # Reversing order for minus strand
        t_scaled=c()
        t_scaled=rbind(t_scaled, sapply(t[['Splt_PerBP']], function(y) Scale_Vector(Data=y,FinalLength=finalL, Extention=ext, Ext_length=c(extLL, extLR), Ext_value=c(extVL, extVR))))
        colnames(t_scaled)=t[['Q_id']]
        save(x=t_scaled, file='${LibName}.${BedName}.R')" > r_file_2_run.R
        Rscript r_file_2_run.R
        """*/
    }

    /* Must get all bed-associated R file (containing data)
    process combine_R {

    }*/

}
/**/


/*process toto {
    tag "$BedName"
    echo true
    input:
    val Labels from ch_dt_labels
    val Files from ch_dt_files
    tuple BedName, file(BedFile), BedReferencePoint, BedFls, BedExts, BedExtls, BedExtvs from design_bed_csv.take(1)
    """
    echo "${BedName} \n ${BedFile} \n ${Labels.join(' ')} \n ${Files.join(' ')}
    """
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
