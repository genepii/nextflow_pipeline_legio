#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl=2

// -----------------------------------------------------------------------------
// WORFLOW FOR GETTING PATH FROM CONFIG FILE
// -----------------------------------------------------------------------------

workflow {
    println params.input_dir
    println params.output_dir
    println params.save_dir
    println params.tmp_dir
    println workflow.workDir
    println params.result
}