#!/usr/bin/env nextflow

/*
 * Enables modules
 */
nextflow.enable.dsl = 2

/*
 * Validate required parameters up front. Without this, a missing --input/--covarfile/--phenofile
 * only surfaces much later as a generic "process input channel evaluates to null" error deep in
 * the DAG (e.g. at MAKEANALYSISSETS) -- costly to hit on a cloud/Google Batch run.
 */
def missingParams = []
if (!params.input)     missingParams << '--input (path to genotype VCF/PLINK files)'
if (!params.covarfile) missingParams << '--covarfile (path to covariates file)'
if (!params.phenofile) missingParams << '--phenofile (path to phenotype file)'
if (missingParams) {
    error("Missing required parameter(s):\n" +
          missingParams.collect { "  - ${it}" }.join('\n') +
          "\nProvide them via -params-file or --param value. See conf/examples/ for templates.")
}

/*
 * Validate that requested phenotype/covariate columns actually exist in
 * --phenofile / --covarfile before running anything. Without this, a typo'd
 * column name only surfaces at EXPORT_PLINK/GWASGALLOP/GWASCPH -- after
 * GENETICQCPLINK/MERGER_CHRS/GWASQC/MAKEANALYSISSETS/COMPUTE_PCA have all
 * already run on a real dataset -- very costly to hit that late on a cloud run.
 */
def readHeaderColumns(path, label) {
    def line
    try {
        line = file(path).withReader { it.readLine() }
    } catch (Exception e) {
        error("Could not read ${label} '${path}': ${e.message}")
    }
    if (!line) {
        error("${label} '${path}' appears to be empty (no header line found).")
    }
    return line.split('\t').collect { it.trim() } as Set
}

def phenoHeader = readHeaderColumns(params.phenofile, '--phenofile')
def covarHeader = readHeaderColumns(params.covarfile, '--covarfile')

def missingColumns = []

(params.pheno_name ?: '').split(/\s+/).findAll { it }.each { col ->
    if (!(col in phenoHeader)) missingColumns << "--pheno_name '${col}' not found in --phenofile columns"
}

def numericList = (params.covar_numeric ?: '').split(/\s+/).findAll { it }
numericList.each { col ->
    // PC1/PC2/... are computed by COMPUTE_PCA and merged in later -- not
    // expected to already exist in the raw --covarfile, so skip those.
    if (!(col ==~ /(?i)^PC\d+$/) && !(col in covarHeader)) {
        missingColumns << "--covar_numeric '${col}' not found in --covarfile columns"
    }
}

(params.covar_categorical ?: '').split(/\s+/).findAll { it }.each { col ->
    if (!(col in covarHeader)) missingColumns << "--covar_categorical '${col}' not found in --covarfile columns"
}

if (params.study_arm_col && !(params.study_arm_col in covarHeader)) {
    missingColumns << "--study_arm_col '${params.study_arm_col}' not found in --covarfile columns"
}

// time_col is used unconditionally by TABLEONE (modules/results.nf), not just
// longitudinal/survival analyses, so it's always checked.
if (params.time_col && !(params.time_col in phenoHeader)) {
    missingColumns << "--time_col '${params.time_col}' not found in --phenofile columns"
}

if (params.covar_interact && !(params.covar_interact in numericList)) {
    missingColumns << "--covar_interact '${params.covar_interact}' must also be listed in --covar_numeric"
}

if (missingColumns) {
    error("Phenotype/covariate column validation failed:\n" +
          missingColumns.collect { "  - ${it}" }.join('\n') +
          "\n\nphenofile ('${params.phenofile}') columns found: ${phenoHeader.sort()}" +
          "\ncovarfile ('${params.covarfile}') columns found: ${covarHeader.sort()}")
}

/*
 * Main workflow log
 */
if (params.longitudinal_flag) {
    MODEL = "lmm_gallop"
} 
else if (params.survival_flag) {
    MODEL = "cph"
}
else {
    MODEL = "glm"
}

log.info """\
 LONG-GWAS - GWAS P I P E L I N E
 ======================================
 Chunk size for genetic processing        : ${params.chunk_size}
 Kinship matrix threshold                 : ${params.kinship}
 R2 threshold                             : ${params.r2thres}
 MAF threshold                            : ${params.minor_allele_freq}
 data ancestry                            : ${params.ancestry}
 genetic data assemble                    : ${params.assembly}
 phenotype name                           : ${params.pheno_name}
 numeric covariates                       : ${params.covar_numeric}
 categorical covariates                   : ${params.covar_categorical}
 interaction covariate                    : ${params.covar_interact}
 analysis                                 : ${MODEL}
 project directory                        : ${params.project_dir}
 analysis name                            : ${params.analysis_name}
 genetic cache key                        : ${params.genetic_cache_key}
 """

/*
 * Datetime
 */
datetime = new java.util.Date()
params.datetime = new java.text.SimpleDateFormat("YYYY-MM-dd'T'HHMMSS").format(datetime)

/* 
 * Import consolidated modules
 */
include { CHECK_REFERENCES; SPLIT_VCF; GENETICQC; GENETICQCPLINK; MERGER_CHUNKS; LD_PRUNE_CHR; MERGER_CHRS; SIMPLE_QC; GWASQC } from './modules/qc.nf'
include { MAKEANALYSISSETS; COMPUTE_PCA; MERGE_PCA; HARMONIZE_CATEGORICAL_COVARS; RAWFILE_EXPORT; EXPORT_PLINK } from './modules/dataprep.nf'
include { GWASGLM; GWASGALLOP; GWASCPH } from './modules/gwas.nf'
include { SAVEGWAS; MANHATTAN; TABLEONE } from './modules/results.nf'

/*
 * Get the cache and the input check channels
 */
Channel
  .fromPath("${params.project_dir}/genotypes/${params.genetic_cache_key}/chromosomes/*/*.{pgen,pvar,psam,log}", checkIfExists: false)
  .map{ f -> tuple(f.getSimpleName(), f) }
  .set{ cache_raw }

Channel
   .fromPath(params.input)
   // Fail fast if the glob matched nothing. Without this, an empty match here
   // silently propagates to zero tasks in every downstream genotype-processing
   // process, and the pipeline reports "succeeded" having done nothing -- e.g.
   // a VWB-mounted resource path (/home/jupyter/...) used instead of the
   // actual gs:// URI, or a typo in the glob.
   .ifEmpty { error("No genotype files matched --input: '${params.input}'. Check the path/glob is correct and reachable -- on Google Batch this must be a real gs:// URI, not a local VM-mounted resource path.") }
   .map{ f -> tuple(f.getSimpleName(), f) }
   .set{ input_check_ch }

/*
 * Guard against genetic_cache_key/genetic_data_id being reused across incompatible
 * input file sets (e.g. two analyses that each pre-filter the same "cohort" to a
 * different, non-overlapping sample subset, but share a genetic_data_id).
 * cache_raw above lists everything ever cached under this cache key on disk,
 * regardless of whether it belongs to THIS run's --input -- restrict it to only
 * fileTags this run's --input actually matches, and warn if stale entries were
 * found and excluded. Without this, MERGER_CHRS would silently merge chromosome
 * files from an unrelated, non-overlapping sample set into this run's dataset.
 * (List values are wrapped in an extra list -- Channel.combine() auto-flattens a
 * bare List value into separate tuple elements, so wrapping is what keeps it as
 * one list-typed value in the combined tuple.)
 */
input_check_ch.map{ fileTag, f -> fileTag }.unique().toList()
    .map{ tags -> [tags] }
    .set{ validTagsWrapped }

cache_raw.map{ fileTag, f -> fileTag }.unique().toList()
    .map{ tags -> [tags] }
    .combine(validTagsWrapped)
    .subscribe{ cacheTags, validTags ->
        def stale = cacheTags - validTags
        if (stale) {
            log.warn "genetic_cache_key '${params.genetic_cache_key}' has ${stale.size()} cached chromosome fileTag(s) that don't match this run's --input (e.g. ${stale.take(3).join(', ')}${stale.size() > 3 ? ', ...' : ''}) -- excluding them from this run. If this genetic_data_id was previously used for a different input file set, give each distinct input set its own genetic_data_id instead of relying on this filter."
        }
    }

cache_raw
    .combine(validTagsWrapped)
    .filter{ fileTag, fCache, validTags -> fileTag in validTags }
    .map{ fileTag, fCache, validTags -> tuple(fileTag, fCache) }
    .set{ cache }

/* 
 * Get the phenotypes arg on a channel
 */
Channel
    .of(params.pheno_name)
    .splitCsv(header: false)
    .collect()
    .set{ phenonames }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {
    // ==================================================================================
    // PROCESS 0: CHECK REFERENCE GENOMES (runs once)
    // ==================================================================================
    CHECK_REFERENCES()
    
    // Prepare reference files channel
    def refDir = params.reference_dir
    reference_files = Channel.fromPath([
        "${refDir}/Genome/hg38.fa.gz",
        "${refDir}/Genome/hg38.fa.gz.fai",
        "${refDir}/Genome/hg38.fa.gz.gzi"
    ] + (params.assembly != 'hg38' ? [
        "${refDir}/Genome/${params.assembly}.fa.gz",
        "${refDir}/Genome/${params.assembly}.fa.gz.fai",
        "${refDir}/Genome/${params.assembly}.fa.gz.gzi",
        "${refDir}/liftOver/${params.assembly}ToHg38.over.chain.gz"
    ] : []), checkIfExists: true)
    .collect()
    .view{ "Reference files: ${it}" }
    
    // ==================================================================================
    // QUALITY CONTROL (QC) PHASE
    // ==================================================================================
    input_check_ch
        .join(cache, remainder: true)
        .filter{ fileTag, fOrig, fCache -> fCache == null }
        .map{ fileTag, fOrig, fCache -> tuple(fileTag, fOrig) }
        .set{ chrvcf }

    // Determine input format from params.input pattern
    def isPlink = params.input =~ /\.(bed|pgen)$/
    
    if (isPlink) {
        // ============================================================
        // PLINK INPUT PATHWAY: Direct cache, no chunking
        // ============================================================
        
        // Gather all companion files (.pgen, .pvar, .psam or .bed, .bim, .fam)
        chrvcf
        .map{ fileTag, fOrig ->
            // Use toUri() to preserve full path (works for both GCS and local files)
            // For GCS: gs://bucket/path/file.pgen
            // For local: file:///path/to/file.pgen
            def fullPath = fOrig.toUri().toString()
            def basePath = fullPath.replaceFirst(/\.(bed|pgen)$/, '')
            def ext = fOrig.name =~ /\.bed$/ ? ['bed', 'bim', 'fam'] : ['pgen', 'pvar', 'psam']
            def files = ext.collect{ file(basePath + '.' + it) }
            tuple(fileTag, files)
        }
        .combine(CHECK_REFERENCES.out.references_flag)
        .map{ fileTag, chr_pfiles, references_flag -> tuple(fileTag, chr_pfiles) }
        .set{ plink_input_ch }

        // Process PLINK files directly to cache
        GENETICQCPLINK(plink_input_ch, reference_files)
        
        // Collect processing status for tracking
        // NOTE: uses params.analyses_dir (conf/params.config), not a profile-level env{} var --
        // see conf/params.config for why this was made an explicit param.
        GENETICQCPLINK.out.chunk_status
            .map{ fileTag, statusFile -> statusFile.text }
            .collectFile(name: "geneticqc_chunk_status_${params.datetime}.tsv",
                         storeDir: "${params.analyses_dir}/${params.genetic_cache_key}/genetic_qc/logs/",
                         seed: "fileTag\tchunkId\tinput\tstart_time\tend_time\texit_code\tstatus\tvariants\n",
                         newLine: false)
        
        // PLINK output goes directly to chrsqced (already in pgen format, no merge needed)
        GENETICQCPLINK.out.plink_qc_cached
            .collect()
            .flatten()
            .map{ fn -> tuple(fn.getSimpleName(), fn) }
            .concat(cache)
            .set{ chrsqced }
            
    } else {
        // ============================================================
        // VCF INPUT PATHWAY: Chunk, process, merge
        // ============================================================
        
        // Split VCF files into chunks using a process (faster on cloud)
        SPLIT_VCF(chrvcf)
        
        // Flatten chunks: [fileTag, fOrig, [chunk1, chunk2, ...]] → multiple [fileTag, fOrig, chunk]
        SPLIT_VCF.out.vcf_chunks
        .transpose()
        .map{ fileTag, fOrig, fChunk -> tuple(fileTag, fOrig, fChunk) }
        .combine(CHECK_REFERENCES.out.references_flag)
        .map{ fileTag, fOrig, fChunk, references_flag -> tuple(fileTag, fOrig, fChunk) }
        .set{ vcf_chunks_ch }

        // Process VCF chunks (adds headers internally)
        GENETICQC(vcf_chunks_ch, reference_files)
        
        // Collect processing status for tracking
        GENETICQC.out.chunk_status
            .map{ fileTag, chunkId, statusFile -> statusFile.text }
            .collectFile(name: "geneticqc_chunk_status_${params.datetime}.tsv",
                         storeDir: "${params.analyses_dir}/${params.genetic_cache_key}/genetic_qc/logs/",
                         seed: "fileTag\tchunkId\tinput\tstart_time\tend_time\texit_code\tstatus\tvariants\n",
                         newLine: false)

        // Merge VCF chunks per chromosome
        GENETICQC.out.snpchunks_names
            .collectFile(newLine: true) 
                            { fileTag, chunkId -> ["${fileTag}.mergelist.txt", chunkId] }
            .set{ chunknames }

        MERGER_CHUNKS(chunknames, GENETICQC.out.snpchunks_merge.collect())
        
        // VCF merged output goes to chrsqced
        MERGER_CHUNKS.out
            .collect()
            .flatten()
            .map{ fn -> tuple(fn.getSimpleName(), fn) }
            .concat(cache)
            .set{ chrsqced }
    }

    // Branch based on skip_pop_split mode
    if (params.skip_pop_split) {
        // Skip population splitting mode: LD prune per chromosome before merging
        LD_PRUNE_CHR(chrsqced.groupTuple(by: 0).map{ fileTag, files -> files })
        
        LD_PRUNE_CHR.out
            .flatten()
            .map{ fn -> tuple(fn.getSimpleName(), fn) }
            .set{ chrsqced_pruned }
        
        // For GWAS: use unpruned chromosome-level data
        chrsqced
            .groupTuple(by: 0)
            .set{ gallop_plink_input }

        // For QC/PCA: merge pruned chromosomes
        chrsqced_pruned
            .map{ fileTag, f -> fileTag }
            // f contains .log, .pgen, .pvar, .psam for each fileTag. Reduce to one per fileTag.
            .unique()
            .collectFile() { fileTag ->
                ["allchr.mergelist.txt", fileTag + '\n'] }
            .set{ list_files_merge }
        chrsqced_pruned
            .map{ fileTag, f -> file(f) }
            .set{ chrfiles }

        MERGER_CHRS(list_files_merge, chrfiles.collect())
        MERGER_CHRS.out
            .flatten()
            .filter{ fName -> ["pgen", "pvar", "psam"].contains(fName.getExtension()) }
            .collect()
            .set{ input_compute_pca }

        // Run simplified QC (no ancestry inference)
        SIMPLE_QC(MERGER_CHRS.out)
        qc_h5_file = SIMPLE_QC.out.simpleqc_h5_file

    } else {
        // Standard mode: merge first, then full QC with ancestry inference
        
        // Prepare channels for downstream analysis
        chrsqced
            .groupTuple(by: 0)
            .set{ gallop_plink_input }

        // Merge all chromosomes
        chrsqced
            .map{ fileTag, f -> fileTag }
            .unique()
            .collectFile() { fileTag ->
                ["allchr.mergelist.txt", fileTag + '\n'] }
            .set{ list_files_merge }
        chrsqced
            .map{ fileTag, f -> file(f) }
            .set{ chrfiles }

        MERGER_CHRS(list_files_merge, chrfiles.collect())
        MERGER_CHRS.out
            .flatten()
            .filter{ fName -> ["pgen", "pvar", "psam"].contains(fName.getExtension()) }
            .collect()
            .set{ input_compute_pca }

        // Run GWAS QC
        GWASQC(MERGER_CHRS.out)
        qc_h5_file = GWASQC.out.gwasqc_h5_file
    }

    // ==================================================================================
    // DATA PREPARATION PHASE
    // ==================================================================================
    MAKEANALYSISSETS(qc_h5_file, params.covarfile)
    COMPUTE_PCA(MAKEANALYSISSETS.out.study_arm_files.flatten(), input_compute_pca)
    MERGE_PCA(COMPUTE_PCA.out.eigenvec)
    HARMONIZE_CATEGORICAL_COVARS(MERGE_PCA.out.flatten())

    // Branch based on analysis type
    if (params.longitudinal_flag | params.survival_flag) {
        // For longitudinal/survival: chunk variants and export to raw format
        // RAWFILE_EXPORT now handles both chunking and export internally
        RAWFILE_EXPORT(gallop_plink_input, HARMONIZE_CATEGORICAL_COVARS.out)
        
        // Flatten to process each raw file individually
        RAWFILE_EXPORT.out.gwas_rawfile
            .transpose()
            .set{ CHUNKS }
        
        PLINK_SAMPLE_LIST = Channel.empty()

    } else {
        // For cross-sectional: use PLINK binary directly (no chunking, no raw export)
        EXPORT_PLINK(HARMONIZE_CATEGORICAL_COVARS.out.flatten(), params.phenofile)
        
        // Collect outputs from EXPORT_PLINK: pheno.tsv, covar_names.txt, n_covar.txt
        // Log files (output[3]) are published automatically via publishDir
        EXPORT_PLINK.out[0]
            .mix(EXPORT_PLINK.out[1], EXPORT_PLINK.out[2])
            .flatten()
            .filter{ it != null }
            .map{ file ->
                // Extract study arm from filename
                def matcher = file.name =~ /(.+)_filtered\.pca\.pheno\.tsv/
                if (matcher.find()) {
                    return [matcher[0][1], file, 'pheno']
                }
                matcher = file.name =~ /(.+)_covar_names\.txt/
                if (matcher.find()) {
                    return [matcher[0][1], file, 'covar_names']
                }
                matcher = file.name =~ /(.+)_n_covar\.txt/
                if (matcher.find()) {
                    return [matcher[0][1], file, 'n_covar']
                }
                return null
            }
            .filter{ it != null }
            .groupTuple(by: 0)
            .map{ study_arm, files, types ->
                // Return all three files grouped by study arm
                def pheno_file = files[types.indexOf('pheno')]
                def covar_names = files[types.indexOf('covar_names')]
                def n_covar = files[types.indexOf('n_covar')]
                return tuple(study_arm, pheno_file, covar_names, n_covar)
            }
            .set{ PLINK_SAMPLE_LIST }
        
        // For GLM: use gallop_plink_input (already grouped per chromosome)
        // Unpack PLINK files: convert from [fileTag, [files]] to [fileTag, log, pgen, psam, pvar]
        // Files are selected by extension for robustness (not positional indexing)
        // Then combine each chunk with PLINK_SAMPLE_LIST (1 sample list applies to all 22 chromosomes)
        gallop_plink_input
            .map{ fileTag, plinkFiles ->
                tuple(
                    fileTag,
                    plinkFiles.find { it.extension == 'log' },
                    plinkFiles.find { it.extension == 'pgen' },
                    plinkFiles.find { it.extension == 'psam' },
                    plinkFiles.find { it.extension == 'pvar' }
                )
            }
            .combine(PLINK_SAMPLE_LIST)
            .set{ CHUNKS }
    }

    // ==================================================================================
    // GWAS ANALYSIS PHASE
    // ==================================================================================
    if (params.longitudinal_flag) {
        GWASGALLOP(CHUNKS, params.phenofile, phenonames)
        GWASRES = GWASGALLOP.out
    }
    else if (params.survival_flag) {
        GWASCPH(CHUNKS, params.phenofile, phenonames)
        GWASRES = GWASCPH.out
    } else {
        GWASGLM(CHUNKS, phenonames)
        
        // Use manifest to create proper tuples
        // GWASGLM.out[0] = result files, GWASGLM.out[1] = manifest files
        
        // Parse manifest: tuple(filename, key)
        GWASGLM.out[1]
            .splitCsv(header: true, sep: '\t')
            .map{ row -> tuple(row.filename, row.key) }
            .set{ manifest_ch }
        
        // Flatten result files and map to tuple(filename, file)
        GWASGLM.out[0]
            .flatten()
            .map{ file -> tuple(file.name, file) }
            .set{ results_ch }
        
        // Join by filename, then remap to (key, file)
        manifest_ch
            .join(results_ch)
            .map{ filename, key, file -> tuple(key, file) }
            .set{ GWASRES }
    }

    GWASRES
        .groupTuple(sort: true)
        .set{ GROUP_RESULTS }

    // ==================================================================================
    // RESULTS MANAGEMENT PHASE
    // ==================================================================================
    SAVEGWAS(GROUP_RESULTS, MODEL)
    if (params.mh_plot) {
        MANHATTAN(SAVEGWAS.out.res_all.collect(), MODEL)
    }
    
    // ==================================================================================
    // TABLE 1 AND DESCRIPTIVE STATISTICS
    // ==================================================================================
    TABLEONE(MAKEANALYSISSETS.out.analytical_set, file(params.covarfile), file(params.phenofile))
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
