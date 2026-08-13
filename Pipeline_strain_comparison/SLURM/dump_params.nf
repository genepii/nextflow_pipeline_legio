#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// -----------------------------------------------------------------------------
// WORFLOW FOR GETTING PATH FROM CONFIG FILE
// -----------------------------------------------------------------------------

workflow {
    println "output_dir=${params.output_dir}"
    println "tmp_dir=${params.tmp_dir}"
    println "work_dir=${workflow.workDir}"
    println "result=${params.result}"
    println "analyse_id=${params.analyse_id}"
    println "metadata=${params.metadata}"
    println "db_snippy=${params.snippy_dir}"
}