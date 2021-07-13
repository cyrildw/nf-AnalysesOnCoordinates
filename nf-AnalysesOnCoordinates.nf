/* nf-AnalysesOnCoordinates.nf
    Pipeline using, bam files, bigwig file and bedfiles to produce
        - peak calling (macs2)
        - heatmaps and correlations (deeptools)
            - options for groups of regions using their ID
        - average density per regions (get_tag_density, R)
        - scaled tag density for each region (get_tag_density, R)


*/

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
		row.BedDTlength, //Used by Deeptools for -m option in plot heatmap
        row.BedReferencePoint, //Used by Deeptools in opposition to "scale-region"
		row.BedExtLengthLeft,  //Used by Deeptools and R to extend the bed coordinates upstream
		row.BedExtLengthRight, //Used by Deeptools and R to extend the bed coordinates downstream
		row.BedRFinalLength,    //Used by R to set the final length of the vector
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
        /* Creates sub bed file from the bedgroupe file
        input   :   ch_before_dt_bed
        output  :   ch_dt_bedGroup_computeMatrix
                    ch_dt_bed_multiBWsummary, ch_dt_bed_computeMatrix
                    bedfiles, saved in /BedFiles/
        */
        tag "$BedName"
        container=''
        publishDir "${params.outdir}/${params.name}/", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.bed')) "./BedFiles/$filename"
            else null
        }

        input:
        tuple BedName, file(BedFile), file(BedGroupFile),BedDTlength,BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedRFinalLength, BedExtension,BedExtValLeft, BedExtValRight from ch_before_dt_bed
        output:
        files "*.txt"
        files "*.bed"
        tuple BedName, file(BedFile), NbGroup, file("${BedName}.GrpFiles.txt"), file("${BedName}.*.bed"),BedDTlength,BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedRFinalLength into (ch_test,ch_dt_bedGroup_computeMatrix)
        tuple BedName, file(BedFile),BedDTlength, BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedRFinalLength into (ch_dt_bed_multiBWsummary, ch_dt_bed_computeMatrix)

        script:
        if(BedGroupFile.isFile() && BedGroupFile.size()!=0 ){
            NbGroup=1
            //Creating 1 file per group + 1 file with goupefile-names.
            //Then greping ids from each group file into the BED file to produce 1 GroupBedFile/group.
            // In bash it takes too much time to iterate over all IDs, I'll try in R
 /*           """
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
*/
                   """
        echo "R --no-save --no-restore --slave <<RSCRIPT
        R.Version()
        grp=readLines('${BedGroupFile}')
        bed=read.table('${BedFile}', stringsAsFactor=FALSE, head=FALSE)
        a=unname(sapply(grp, function(x) strsplit(x, split=';')[[1]])) #spliting the line
        if(is.list(a) && length(a) != 1){ #In case their is more than 1 group
            names(a)=sapply(a, function(x) x[1]) #using first element as name for the group
            a=sapply(a, function(x) x[-1]) #removing the first element

            for(i in 1:length(a)){
                filename=paste0('${BedName}', '.', names(a)[i], '.bed')
                write.table(x=bed[bed[[4]] %in% a[[i]],], file=filename,quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\\t')
                cat(x=paste(names(a)[i], '\\t', filename), append=TRUE, file='${BedName}.GrpFiles.txt', fill=TRUE)
            }
        } else{
            filename=paste0('${BedName}', '.', a[1], '.bed')
            write.table(x=bed[bed[[4]] %in% a[-1],], file=filename,quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\\t')
            cat(x=paste(a[1], '\\t', filename), append=TRUE, file='${BedName}.GrpFiles.txt', fill=TRUE)
        }
        " > r_file_2_run.R
        bash r_file_2_run.R
        """
        }
        else{
            NbGroup=0
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
        /* Computes the npz Matrix file required for the dt_PlotCorrelation process.
            The npz file is saved in the DeeptoolsData directory
        input   :   ch_dt_bed_multiBWsummary
                    ch_dt_files_multiBWsummary
        ouput   :   ch_multibw_matrix
                    ch_multibw_bedname
                    npz file saved in /DeeptoolsData/
        */
        tag "$BedName"
        label "multiCpu"
        publishDir "${params.outdir}/${params.name}/DeeptoolsData", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.npz')) "./$filename"
            else null
        }
        input:
        tuple BedName, file(BedFile),BedDTlength, BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedRFinalLength from ch_dt_bed_multiBWsummary
        file(Files) from ch_dt_files_multiBWsummary
        output:
        file("${params.deeptools_PC_prefix}.Matrix.${BedName}.npz") into ch_multibw_matrix //the computed matrix
        val(BedName) into ch_multibw_bedname
        script:
        """
        multiBigwigSummary BED-file \
        -b ${Files.join(' ')} \
        -o ${params.deeptools_PC_prefix}.Matrix.${BedName}.npz \
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
        file("${params.deeptools_PC_prefix}.CorTable.${BedName}.tab")
        file("${params.deeptools_PC_prefix}.HeatMap.${BedName}.pdf")
        script:
        """
        plotCorrelation \
        --corData ${Matrix} ${params.deeptools_PC_options} \
        --outFileCorMatrix ${params.deeptools_PC_prefix}.CorTable.${BedName}.tab \
        -o ${params.deeptools_PC_prefix}.HeatMap.${BedName}.pdf \
        --plotTitle ${BedName} \
        --labels ${Labels.join(' ')}
        """
    }

    
    process dt_ComputeMatrix {
        // Compute the matrix file (.gz) for the plotHeatmap
        tag "$BedName"
        label "multiCpu"
        publishDir "${params.outdir}/${params.name}/DeeptoolsData", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.gz')) "./$filename"
            else null
        }

        input:
        tuple BedName, file(BedFile), BedDTlength, BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedRFinalLength  from ch_dt_bed_computeMatrix
        file(Files) from ch_dt_files_computeMatrix
        output:
        file("${params.deeptools_HM_prefix}.Matrix.${BedName}.gz") into ch_computeMatrix_matrix //the computed matrix
        val(BedName) into ch_computeMatrix_bedname
        script:
        
        if(BedReferencePoint=='false')
            """
            computeMatrix scale-regions \
            -S ${Files.join(' ')} \
            -R ${BedFile} \
            -b ${BedExtLengthLeft} \
            -a ${BedExtLengthRight} \
            -m ${BedDTlength} \
            --skipZeros \
            -p ${task.cpus} \
            -o ${params.deeptools_HM_prefix}.Matrix.${BedName}.gz
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
            -o ${params.deeptools_HM_prefix}.Matrix.${BedName}.gz
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
        file("${params.deeptools_HM_prefix}.${BedName}.pdf")
        
        script:
        """
        plotHeatmap \
        --matrixFile ${Matrix} \
        -o ${params.deeptools_HM_prefix}.${BedName}.pdf \
        --startLabel 'st' \
        --endLabel 'end' \
        --refPointLabel 0 \
        --regionsLabel '' \
        --heatmapHeight ${params.deeptools_HM_heatmapHeight} ${params.deeptools_HM_options} \
        --averageTypeSummaryPlot ${params.deeptools_HM_TypeSummaryPlot} \
        --labelRotation ${params.deeptools_HM_labelRotation} \
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
        tuple BedName, file(BedFile), NbGroup, file(BedGrpFile), file(BedGrpBedFiles),BedDTlength, BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedRFinalLength from ch_dt_bedGroup_computeMatrix
        file(Files) from ch_dt_files_groupcomputeMatrix
        val(Labels) from ch_dt_labels_groupHeatmap
        output:
        file("${params.deeptools_HM_prefix}.Matrix.Group.${BedName}.gz")//the computed matrix
        file("${params.deeptools_HM_prefix}.Group.${BedName}.pdf")//the heatmap
        //val(BedName) into ch_computeMatrix_bedname
        
        when:
        NbGroup == 1
        script:        
        if(BedReferencePoint=='false'){
            """
            computeMatrix scale-regions \
            -S ${Files.join(' ')} \
            -R ${BedGrpBedFiles.join(' ')} \
            -b ${BedExtLengthLeft} \
            -a ${BedExtLengthRight} \
            -m ${BedDTlength} \
            --skipZeros \
            -p ${task.cpus} \
            -o ${params.deeptools_HM_prefix}.Matrix.Group.${BedName}.gz

            plotHeatmap \
            --matrixFile ${params.deeptools_HM_prefix}.Matrix.Group.${BedName}.gz \
            -o ${params.deeptools_HM_prefix}.Group.${BedName}.pdf \
            --startLabel 'start' \
            --endLabel 'end' \
            --refPointLabel 0 \
            --heatmapHeight ${params.deeptools_HM_heatmapHeight} ${params.deeptools_HM_options} \
            --averageTypeSummaryPlot ${params.deeptools_HM_TypeSummaryPlot} \
            --labelRotation ${params.deeptools_HM_labelRotation} \
            --xAxisLabel ${BedName} \
            --samplesLabel ${Labels.join(' ')}
            """
        }
        else{
            """
            computeMatrix reference-point \
            -S ${Files.join(' ')} \
            -R ${BedGrpBedFiles.join(' ')} \
            -b ${BedExtLengthLeft} \
            -a ${BedExtLengthRight} \
            --skipZeros \
            -p ${task.cpus} \
            -o ${params.deeptools_HM_prefix}.Matrix.Group.${BedName}.gz
            
            plotHeatmap \
            --matrixFile ${params.deeptools_HM_prefix}.Matrix.Group.${BedName}.gz \
            -o ${params.deeptools_HM_prefix}.Group.${BedName}.pdf \
            --startLabel '-${BedExtLengthLeft}' \
            --endLabel '${BedExtLengthRight}' \
            --refPointLabel 0 \
            --heatmapHeight ${params.deeptools_HM_heatmapHeight} ${params.deeptools_HM_options} \
            --averageTypeSummaryPlot ${params.deeptools_HM_TypeSummaryPlot} \
            --labelRotation ${params.deeptools_HM_labelRotation} \
            --xAxisLabel ${BedName} \
            --samplesLabel ${Labels.join(' ')}
            """
        }
        /*"""
        plotHeatmap \
        --matrixFile ${params.deeptools_HM_prefix}.Matrix.Group.${BedName}.gz \
        -o ${params.deeptools_HM_prefix}.Group.${BedName}.pdf \
        --startLabel '-${BedExtLengthLeft}' \
        --endLabel '${BedExtLengthRight}' \
        --refPointLabel 0 \
        --labelRotation ${params.deeptools_HM_labelRotation} ${params.deeptools_HM_options} \
        --xAxisLabel ${BedName} \
        --samplesLabel ${Labels.join(' ')}
        """*/
    
    }
}


/*R analyses includes : 
VERIF    - Create the bed files if it requires some bp extension.
OK      - Combine all libraries with all bedfiles
OK      - Getting Tag Density over the bed files
OK    - Converting all density tables to R object with scaling
- Producing combined graphics for
    -all elements
    -grouped elements
    -quantiles
- Outputing R objects (and R scripts ?)*/

//r_func=Channel.fromPath(params.r_scaling)



if(params.r_analyses){
    //process check_bed_format {
        /* The bed file for the r_analyses should be a 6column format.
TODO    - Count the number of columns "awk '{print NF}' file | sort -nu | tail -n 1"
TODO    - If < 6 columns : add name (peak_$i); score (0); strand (.)
TODO    - Gather all avg density in one file per bed
TODO    - report file with file location per bigwig file
*/

    //}

    process create_bed_with_ext {
    /*In case the bed conf requires extension on the flanking regions, this process modifies the bed file accordingly.
ok    - save the bedfile
TODO    - remove the unnecessary fields from input.
*/
        tag "$BedName:$BedExtension-$BedExtLengthLeft:$BedExtLengthRight"
        publishDir "${params.outdir}/${params.name}/", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.bed')) "./BedFiles/$filename"
            else null
        }

        input:
        tuple BedName, file(BedFile),file(BedGrpFile),BedDTlength, BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedRFinalLength, BedExtension, BedExtValLeft,BedExtValRight from ch_before_R_bed
        output:
        tuple BedName, file("${BedFile.baseName}.ext.bed"),file(BedGrpFile), BedDTlength, BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedRFinalLength, BedExtension, BedExtValLeft,BedExtValRight into ch_for_R_ext_Bed
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
        .set{ch_R_bed_n_lib}
        
    /*Channel.fromPath(params.r_function_file) // Requires to combine the r_function file with the bed n lib channel
        .combine(ch_R_bed_n_lib).view()
        .set{ch_R_rfunc_bed_lib}
*/

    process tag_density {
    /* Get the read density (from bw file) on coordinates (bed file) using get_tag_density script.
OK      - get_tag_density only getting the first 8 columns
OK    - use 1, 2, 3, 4, 8, 6 columns to produce bed file with average tag density per coordinates
OK    - send the initial file to the ch_ToScale channel
OK    - Gather all avg density in one file per bed
?TODO?  - Output everything according to initial order.


    */
        tag "$LibName - $BedName"
        /*publishDir "${params.outdir}/${params.name}/", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.avgdensity.bed')) "./RData/$filename"
            else null
        }*/

        input:
        tuple BedName, file(BedFile), file(BedGrpFile), BedDTlength, BedReferencePoint, BedExtLengthLeft, BedExtLengthRight, BedRFinalLength, BedExtension, BedExtValLeft, BedExtValRight,
                LibName, file(LibBam), file(LibBai), file(LibBW), LibSequenced, LibMapped, LibUnique, LibInsertSize, LibQpcrNorm, LibType, LibProj, LibExp, LibCondition, LibOrder, LibIsControl, LibControl   from ch_R_bed_n_lib
        
        output:
        tuple BedName, LibName, file("${LibName}.${BedName}.avgdensity.bed") into ch_avgTD
        tuple file("${LibName}.${BedName}.tagdensity_output"), LibName, BedName, BedExtLengthLeft, BedExtLengthRight, BedRFinalLength, BedExtension, BedExtValLeft, BedExtValRight into ch_ToScale
        
        script:
        """
        get_tag_density -f ${LibBW} ${BedFile} | awk '{for(i=1;i<=8;i++) printf \$i"\\t"; print ""}' - > ${LibName}.${BedName}.tagdensity_output
        grep -v "#" ${LibName}.${BedName}.tagdensity_output | awk '{ print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$8"\\t"\$6 }' > ${LibName}.${BedName}.avgdensity.bed
        """
    }
    ch_avgTD.groupTuple(by: 0).set{ch_grouped_avgTD}

    process combine_avgTD_per_bed {
        tag "$BedName"
        publishDir "${params.outdir}/${params.name}/", mode: 'copy', //params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.tsv')) "./RData/$filename"
            else null
        }

        input: 
        tuple BedName, LibNames, path(R_files) from ch_grouped_avgTD
        output:
        file("r_file_2_run.R")
        file("${BedName}.avg_TagDensity.tsv")
        script:
        """
        echo "R --no-save --no-restore --slave <<RSCRIPT
        R.Version()
        bedName='${BedName}'
        libNames=c('${LibNames.join('\',\'')}')
        libFiles=c('${R_files.join('\',\'')}')
        sortby=order(libNames);
        libNames=libNames[sortby]; libFiles=libFiles[sortby]
        resTable=read.table(libFiles[1], head=FALSE, stringsAsFactor=FALSE)[,c(1, 2, 3, 6)]
        colnames(resTable)=c(\'chr\', \'start\', \'end\', \'strand\')
        for( i in 1:length(libNames)){
            resTable=cbind(resTable, read.table(libFiles[i], head=FALSE, stringsAsFactor=FALSE)[,c(5)])
            colnames(resTable)[dim(resTable)[2]]=libNames[i]}
        write.table(resTable, file='${BedName}.avg_TagDensity.tsv', quote=FALSE, row.names=FALSE, col.names=TRUE, sep=\'\t\')
        " > r_file_2_run.R
        bash r_file_2_run.R
        """

    }



    if(params.r_scaling){

        Channel.fromPath(params.r_function_file) // Requires to combine the r_function file with the bed n lib channel
            .combine(ch_ToScale)
            .set{ch_R_rfunc_toScale}
        
        process density_R_scaling {
        /* From the result of get_tag_density script, use R to rescale regions to a fixed number of values.
    OK      - get only 3 columns from get_tag_density for each feature : ID, strand, density for each base pair (#4, #6, #7).
    OK      - create a R-script allowing to scale each feature to a fixed length set by BedFinalLength
    OK      - execute the R-script to only save the final table of dimension nb_of_bed_coordinates x BedFinalLength
    OK      - save the R_table
    OK      - output a channel with the BedName, LibName, r_table

        */
            tag "$LibName - $BedName"
            input:
            tuple file(R_function), file(TagDensity) ,LibName, BedName, BedExtLengthLeft, BedExtLengthRight, BedRFinalLength, BedExtension, BedExtValLeft, BedExtValRight from ch_R_rfunc_toScale
            
            output:
            file(temp_file)
            file("r_file_2_run.R")
            tuple BedName, LibName, file("${LibName}.${BedName}.R") into ch_scaled_R
            
            script:
            """
            grep -v "#" ${TagDensity} | awk '{ print \$4"\\t"\$6"\\t"\$7 }' > temp_file
            
            echo "R --no-save --no-restore --slave <<RSCRIPT
            R.Version()
            source('${R_function}')
            finalL=${BedRFinalLength}
            ext='${BedExtension}'
            if(ext=='false'){ext=FALSE}
            if(ext=='true'){ext=TRUE}
            extLL=${BedExtLengthLeft};extLR=${BedExtLengthRight};
            extVL=${BedExtValLeft};extVR=${BedExtValRight}; 
            t=as.list(read.table('temp_file', stringsAsFactor=FALSE, header=TRUE))
            names(t)=c('Q_id', 'Strand', 'PerBP')
            t[['Splt_PerBP']]=lapply(as.character(t[['PerBP']]), function(x) as.numeric(strsplit(x,';')[[1]]))
            for(k in 1:length(t[['Strand']])){ 
                if(as.character(t[['Strand']][k])=='-'){t[['Splt_PerBP']][[k]]=rev(t[['Splt_PerBP']][[k]])}
            } # Reversing order for minus strand
            t_scaled=c()
            t_scaled=rbind(t_scaled, sapply(t[['Splt_PerBP']], function(y) Scale_Vector(Data=y,FinalLength=finalL, Extention=ext, Ext_length=c(extLL, extLR), Ext_value=c(extVL, extVR))))
            colnames(t_scaled)=t[['Q_id']]
            save(x=t_scaled, file='${LibName}.${BedName}.R')
            #write.table(x=t_scaled, file='${LibName}.${BedName}.R', quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
            RSCRIPT
            " > r_file_2_run.R
            bash r_file_2_run.R
            """
        }
        ch_scaled_R.groupTuple(by: 0).set{ch_grouped_scaled_R}

        /* For each bed, get the R_table files and LibName*/
        process combine_R_per_bed {
            tag "$BedName"
            publishDir "${params.outdir}/${params.name}/", mode: 'copy', //params.publish_dir_mode,
            saveAs: { filename ->
                if (filename.endsWith('.scaledData.R')) "./RData/$filename"
                else null
            }

            input: 
            tuple BedName, LibNames, path(R_files) from ch_grouped_scaled_R
            output:
            file("r_file_2_run.R")
            file("${BedName}.scaledData.R")
            script:
            """
            echo "R --no-save --no-restore --slave <<RSCRIPT
            R.Version()
            bedName='${BedName}'
            libNames=c('${LibNames.join('\',\'')}')
            libFiles=c('${R_files.join('\',\'')}')
            sortby=order(libNames);
            libNames=libNames[sortby]; libFiles=libFiles[sortby]
            bedData=list();for( i in 1:length(libNames)){load(libFiles[i]); bedData[[libNames[i]]]=t_scaled}
            save(x=bedData, file='${BedName}.scaledData.R')
            " > r_file_2_run.R
            bash r_file_2_run.R
            """

        }
    }

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
