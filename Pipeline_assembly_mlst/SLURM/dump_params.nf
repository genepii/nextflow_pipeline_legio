#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// -----------------------------------------------------------------------------
// WORFLOW FOR GETTING PATH FROM CONFIG FILE
// -----------------------------------------------------------------------------

workflow {
    println "input_dir=${params.input_dir}"
    println "output_dir=${params.output_dir}"
    println "save_dir=${params.save_dir}"
    println "tmp_dir=${params.tmp_dir}"
    println "work_dir=${workflow.workDir}"
    println "result=${params.result}"
    println "analyse_id=${params.analyse_id}"
    println "partition=${params.rep_partition}"
    println "metadata=${params.rep_metadata}"
}