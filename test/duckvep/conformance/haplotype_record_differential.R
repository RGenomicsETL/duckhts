#!/usr/bin/env Rscript
# Source-record geometry audit against the pinned, unmodified Haplosaurus runner.
# Complete sequence/count/provenance differences remain failures, including conflicts.

canonical <- function(rows, provenance = TRUE) {
  if (!length(rows)) return(list())
  keys <- vapply(rows, function(x) jsonlite::toJSON(list(cds=x$cds, protein=x$protein),
    auto_unbox=TRUE, na='null'), '')
  lapply(split(rows, keys), function(group) {
    value <- list(cds=group[[1L]]$cds, protein=group[[1L]]$protein,
      count=sum(vapply(group, function(x) as.numeric(x$count), 0)))
    if (provenance) value$contributors <- sort(unique(unlist(lapply(group, `[[`, 'contributors'),
      use.names=FALSE)))
    value
  })
}

main <- function() {
  opt <- optparse::parse_args(optparse::OptionParser(option_list=list(
    optparse::make_option('--seed', type='integer', default=173L),
    optparse::make_option('--random-cases', dest='random_cases', type='integer', default=512L),
    optparse::make_option('--rare-per-stratum', dest='rare_per_stratum', type='integer', default=0L),
    optparse::make_option('--extension-receipt', dest='extension_receipt', default=NULL),
    optparse::make_option('--vep-prefix', dest='vep_prefix', default='/root/miniconda3/envs/vep')
  )))
  stopifnot(!is.na(opt$seed), opt$random_cases >= 0L, opt$random_cases <= 65536L-144L)
  stopifnot(!is.na(opt$rare_per_stratum), opt$rare_per_stratum >= 0L)
  source('scripts/duckvep_evidence.R', local=TRUE)
  root <- normalizePath('.')
  revision <- duckvep_evidence_revision(root)
  extension <- normalizePath('build/release/duckhts.duckdb_extension')
  binding <- 'diagnostic_unbound'
  if (!is.null(opt$extension_receipt)) {
    duckvep_evidence_assert_checkout(root, revision)
    binding <- duckvep_evidence_read_extension_receipt(opt$extension_receipt, root, extension, revision)$binding
  }
  out <- tempfile(paste0('haplotype_records_seed',opt$seed,'_'), tmpdir='test/duckvep/conformance/results')
  dir.create(out)
  message('Record artifacts: ', out)
  pins <- c(vep='57ea5c52340acc1f156267f810ad162e26597082',
    variation='2fb834b987ede3824e200197a838ce11e91aeb4b')
  mirrors <- normalizePath(c('.sync/ensembl-vep','.sync/ensembl-variation'))
  for (i in seq_along(pins)) stopifnot(
    identical(duckvep_evidence_command('git', c('-C',mirrors[i],'rev-parse','HEAD'), 'oracle revision'),
      unname(pins[i])),
    !length(duckvep_evidence_command('git', c('-C',mirrors[i],'status','--porcelain'), 'oracle checkout')))
  prefix <- normalizePath(opt$vep_prefix)
  environment <- duckvep_evidence_command('micromamba', c('list','-p',prefix,'--explicit'), 'oracle environment')
  writeLines(environment, file.path(out,'environment.txt'))
  stopifnot(identical(duckvep_evidence_explicit_packages(environment),
    duckvep_evidence_explicit_packages(readLines('test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt'))))
  Sys.setenv(DUCKHTSBENCH_REGISTRY=file.path(root,'r/duckhtsbench/inst/benchmark_registry.tsv'))
  paths <- duckhtsbench::duckhts_bench_stage_repository_fixtures(root,'duckvep-haplotypes')
  paths <- paths['haplotype_benchmark_reference']
  cds <- readLines(paths[['haplotype_benchmark_reference']])[2L]
  stopifnot(nchar(cds) == 180L)
  complement <- function(x) paste(rev(strsplit(chartr('ACGT','TGCA',x),'',fixed=TRUE)[[1L]]),collapse='')
  change <- function(x) chartr('ACGT','CGTA',x)
  dna <- function(n) paste(sample(c('A','C','G','T'),n,replace=TRUE),collapse='')
  shapes <- c('disjoint','adjacent','same_start_snv','same_start_duplicate',
    'mnv_retained_middle','containing_deletion','containing_insertion','partial_overlap',
    'shared_deletion_anchor','same_end','same_start_length_change','duplicate_deletion')
  genotypes <- data.frame(phase=c('cis','trans','homozygous','missing','prefixed','unphased'),
    a=c('1|0','1|0','1|1','.|1','|0|1','0/1'), b=c('1|0','0|1','1|1','1|0','1|1','1/0'))
  cases <- expand.grid(shape=shapes, phase=genotypes$phase, strand=c(1L,-1L), stringsAsFactors=FALSE)
  cases$generated <- FALSE
  set.seed(opt$seed)
  if (opt$random_cases) cases <- rbind(cases, data.frame(shape='random_overlap',
    phase=sample(genotypes$phase,opt$random_cases,replace=TRUE),
    strand=sample(c(1L,-1L),opt$random_cases,replace=TRUE),generated=TRUE))
  cases$source_ploidy <- NA_integer_
  cases$rare <- FALSE
  rare_gt <- expand.grid(phase=c('called_pipe','called_slash','mixed','leading_pipe',
    'leading_slash','missing_first','missing_last','all_missing','late_alt'),
    source_ploidy=c(1L,2L,4L,8L,16L,64L),stringsAsFactors=FALSE)
  rare_gt <- subset(rare_gt,
    !(source_ploidy == 1L & phase %in% c('called_slash','missing_first','missing_last')) &
    !(source_ploidy < 4L & phase %in% c('mixed','late_alt')))
  strata <- merge(expand.grid(shape=shapes,strand=c(1L,-1L),stringsAsFactors=FALSE),rare_gt)
  stopifnot(nrow(cases) + as.double(nrow(strata))*opt$rare_per_stratum <= 65536L)
  if (opt$rare_per_stratum) {
    rare <- strata[rep(seq_len(nrow(strata)),each=opt$rare_per_stratum),]
    rare$generated <- rare$rare <- TRUE
    cases <- rbind(cases,rare[names(cases)])
  }
  raw_gt <- function(kind,n) {
    alleles <- sample(c('0','1'),n,replace=TRUE)
    alleles[sample.int(n,1L)] <- '1'
    separators <- rep('|',n-1L)
    if (kind == 'called_slash') separators[] <- '/'
    if (kind == 'mixed') {
      separators <- sample(c('/','|'),n-1L,replace=TRUE)
      separators[1:2] <- c('/','|')
    }
    if (kind == 'missing_first') { alleles[1L] <- '.'; alleles[n] <- '1' }
    if (kind == 'missing_last') { alleles[n] <- '.'; alleles[1L] <- '1' }
    if (kind == 'all_missing') alleles[] <- '.'
    if (kind == 'late_alt') { alleles[] <- '0'; alleles[n] <- '1' }
    prefix <- switch(kind,leading_pipe='|',leading_slash='/', '')
    paste0(prefix,paste0(alleles,c(separators,''),collapse=''))
  }
  cases$transcript_index <- cases$seq_region <- seq_len(nrow(cases))-1L
  cases$chrom <- sprintf('chrR%06d',seq_len(nrow(cases)))
  cases$transcript <- sprintf('HR%06d',seq_len(nrow(cases)))
  cases$cds <- cds
  fasta <- gff <- records <- vector('list',nrow(cases))
  for (i in seq_len(nrow(cases))) {
    genome <- paste0(strrep('A',10L),if (cases$strand[i] == 1L) cds else complement(cds),strrep('A',10L))
    fasta[[i]] <- c(paste0('>',cases$chrom[i]),genome)
    id <- cases$transcript[i]
    attrs <- c(paste0('ID=gene:',id,';biotype=protein_coding'),
      paste0('ID=transcript:',id,';Parent=gene:',id,';biotype=protein_coding'),
      paste0('ID=exon:',id,';Parent=transcript:',id),paste0('Parent=transcript:',id))
    gff[[i]] <- paste(cases$chrom[i],'records',c('gene','mRNA','exon','CDS'),11,190,'.',
      if (cases$strand[i] == 1L) '+' else '-',c('.','.','.','0'),attrs,sep='\t')
    geometry <- switch(cases$shape[i],
      disjoint=c(40L,1L,47L,1L), adjacent=c(40L,1L,41L,1L),
      same_start_snv=c(40L,1L,40L,1L), same_start_duplicate=c(40L,1L,40L,1L),
      mnv_retained_middle=c(40L,7L,43L,1L), containing_deletion=c(40L,7L,43L,1L),
      containing_insertion=c(40L,4L,42L,1L), partial_overlap=c(40L,6L,43L,5L),
      shared_deletion_anchor=c(40L,5L,44L,3L), same_end=c(40L,7L,42L,5L),
      same_start_length_change=c(40L,4L,40L,1L), duplicate_deletion=c(40L,5L,40L,5L),
      random_overlap={a <- sample(35:120,1L); n <- sample(2:10,1L); c(a,n,a+sample(0:(n-1L),1L),sample(1:10,1L))})
    if (cases$rare[i]) {
      start <- sample(15:105,1L)
      n <- sample(3:25,1L)
      inside <- sample.int(n-2L,1L)
      geometry <- switch(cases$shape[i],
        disjoint=c(start,1L,start+sample(2:25,1L),1L), adjacent=c(start,1L,start+1L,1L),
        same_start_snv=c(start,1L,start,1L), same_start_duplicate=c(start,1L,start,1L),
        mnv_retained_middle=c(start,n,start+inside,1L),
        containing_deletion=c(start,n,start+inside,1L),
        containing_insertion=c(start,n,start+inside,1L),
        partial_overlap=c(start,n,start+inside,n-inside+sample(1:10,1L)),
        shared_deletion_anchor=c(start,n,start+n-1L,sample(2:10,1L)),
        same_end=c(start,n,start+inside,n-inside),
        same_start_length_change=c(start,n,start,1L), duplicate_deletion=c(start,n,start,n))
    }
    positions <- c(geometry[c(1L,3L)],165L)
    refs <- substring(genome,positions,positions+c(geometry[c(2L,4L)],1L)-1L)
    alts <- vapply(refs,change,'')
    shape <- cases$shape[i]
    if (shape == 'same_start_snv') alts[2L] <- change(alts[1L])
    if (shape == 'mnv_retained_middle') alts[1L] <- paste0(change(substr(refs[1L],1L,1L)),
      substr(refs[1L],2L,nchar(refs[1L])-1L),change(substring(refs[1L],nchar(refs[1L]))))
    if (shape %in% c('containing_deletion','duplicate_deletion','shared_deletion_anchor'))
      alts[1L] <- substr(refs[1L],1L,1L)
    if (shape %in% c('duplicate_deletion','shared_deletion_anchor')) alts[2L] <- substr(refs[2L],1L,1L)
    if (shape %in% c('containing_insertion','same_start_length_change'))
      alts[1L] <- paste0(substr(refs[1L],1L,1L),
        if (cases$rare[i]) dna(sample(c(1:8,15L,16L,31L),1L)) else 'AC',substring(refs[1L],2L))
    if (shape == 'partial_overlap') alts[1:2] <-
      if (cases$rare[i]) vapply(sample(1:20,2L,replace=TRUE),dna,'') else c('ACGT','GTC')
    if (shape == 'random_overlap') {
      alts[1:2] <- vapply(sample(1:10,2L,replace=TRUE),dna,'')
      for (j in 1:2) if (alts[j] == refs[j]) alts[j] <- change(alts[j])
    }
    if (cases$rare[i]) for (j in 1:2) if (alts[j] == refs[j]) alts[j] <- change(alts[j])
    gt <- if (cases$rare[i]) list(a=raw_gt(cases$phase[i],cases$source_ploidy[i]),
      b=raw_gt(cases$phase[i],cases$source_ploidy[i])) else genotypes[match(cases$phase[i],genotypes$phase),]
    records[[i]] <- data.frame(event_index=3L*(i-1L)+1:3,seq_region=cases$seq_region[i],
      chrom=cases$chrom[i],transcript_index=cases$transcript_index[i],position=positions,
      reference=refs,alternate=alts,source_id=c('a','b','anchor'),sample_index=0L,gt=c(gt$a,gt$b,'1|1'))
  }
  records <- do.call(rbind,records)
  records <- records[order(records$seq_region,records$position,records$event_index),]
  stopifnot(!anyDuplicated(records$event_index),nrow(records) == 3L*nrow(cases),
    all(records$reference != records$alternate))
  pairs <- which(records$source_id != 'anchor')
  pairs <- split(pairs,records$transcript_index[pairs])
  pairs <- pairs[match(as.character(cases$transcript_index),names(pairs))]
  for (i in which(cases$rare)) {
    gt <- records$gt[pairs[[i]]]
    atoms <- strsplit(sub('^[|/]','',gt),'[|/]')
    stopifnot(length(gt) == 2L,all(lengths(atoms) == cases$source_ploidy[i]))
    kind <- cases$phase[i]
    if (kind == 'mixed') stopifnot(all(grepl('|',gt,fixed=TRUE) & grepl('/',gt,fixed=TRUE)))
    if (kind == 'leading_pipe') stopifnot(all(startsWith(gt,'|')))
    if (kind == 'leading_slash') stopifnot(all(startsWith(gt,'/')))
    if (kind == 'all_missing') stopifnot(all(unlist(atoms) == '.'))
    if (kind %in% c('missing_first','missing_last')) stopifnot(all(vapply(atoms,
      function(x) sum(x == '.') == 1L && any(x == '1'),TRUE)))
    if (kind == 'late_alt') stopifnot(all(vapply(atoms,
      function(x) all(head(x,-1L) == '0') && tail(x,1L) == '1',TRUE)))
  }
  saveRDS(list(cases=cases,records=records),file.path(out,'inputs.rds'))
  writeLines(unlist(fasta),file.path(out,'reference.fa'))
  writeLines(c('##gff-version 3',unlist(gff)),file.path(out,'model.gff3'))
  writeLines(c('##fileformat=VCFv4.4',paste0('##contig=<ID=',cases$chrom,',length=200>'),
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample',
    with(records,paste(chrom,position,source_id,reference,alternate,'.','PASS','.','GT',gt,sep='\t'))),
    file.path(out,'calls.vcf'))
  run <- function(command,args,name) {
    status <- system2(command,shQuote(args),stdout=file.path(out,paste0(name,'.stdout')),
      stderr=file.path(out,paste0(name,'.stderr')))
    if (status != 0L) stop(name,' failed; artifacts retained: ',out)
  }
  run('samtools',c('faidx',file.path(out,'reference.fa')),'faidx')
  run('bcftools',c('norm','-c','e','-f',file.path(out,'reference.fa'),'-o',
    file.path(out,'ref_checked.vcf'),file.path(out,'calls.vcf')),'reference_check')
  run('bgzip',file.path(out,'model.gff3'),'bgzip')
  run('tabix',c('-p','gff',file.path(out,'model.gff3.gz')),'tabix')
  libs <- paste(c(file.path(mirrors,'modules'),file.path(prefix,'share/ensembl-vep-116.0-0')),collapse=':')
  run('micromamba',c('run','--clean-env','--env',paste0('PERL5LIB=',libs),'-p',prefix,
    'perl','test/duckvep/conformance/haplotype_oracle.pl',file.path(out,'calls.vcf'),
    file.path(out,'reference.fa'),file.path(out,'model.gff3.gz')),'oracle')
  oracle <- lapply(readLines(file.path(out,'oracle.stdout')),jsonlite::fromJSON,simplifyVector=FALSE)
  names(oracle) <- vapply(oracle,`[[`,'','transcript')
  stopifnot(!anyDuplicated(names(oracle)),setequal(names(oracle),cases$transcript))
  oracle <- oracle[match(cases$transcript,names(oracle))]
  con <- DBI::dbConnect(duckdb::duckdb(config=list(allow_unsigned_extensions='true')))
  on.exit(DBI::dbDisconnect(con,shutdown=TRUE),add=TRUE)
  q <- function(x) as.character(DBI::dbQuoteString(con,x))
  DBI::dbExecute(con,paste('LOAD',q(extension)))
  DBI::dbExecute(con,'SET threads=4')
  DBI::dbWriteTable(con,'models',cases)
  DBI::dbWriteTable(con,'records',records)
  queries <- c('SELECT seq_region::UINTEGER seq_region FROM models ORDER BY seq_region',
    'SELECT transcript_index::UINTEGER transcript_index,seq_region::UINTEGER seq_region,
     11::UBIGINT transcript_start,190::UBIGINT transcript_end,strand::TINYINT strand,
     transcript_index::UINTEGER gene_index,3::UBIGINT transcript_flags,11::UBIGINT cds_start,
     190::UBIGINT cds_end,cds::BLOB cds_sequence,1::UTINYINT codon_table FROM models ORDER BY transcript_index',
    'SELECT transcript_index::UINTEGER transcript_index,11::UBIGINT exon_start,190::UBIGINT exon_end,
     1::UBIGINT exon_cdna_start,180::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase
     FROM models ORDER BY transcript_index')
  stopifnot(DBI::dbGetQuery(con,paste0("SELECT loaded FROM duckvep_model_load('records',",
    paste(q(queries),collapse=','),')'))$loaded)
  actual <- DBI::dbGetQuery(con,"SELECT * FROM duckvep_haplotypes(
    'SELECT *,[alternate] alternates FROM records','records',input_mode:='source_records',phase_policy:='vep116_compat')")
  saveRDS(actual,file.path(out,'actual.rds'))
  stopifnot(all(actual$carrier_count == vapply(actual$carriers,nrow,1L)))
  rows_by_tx <- split(seq_len(nrow(actual)),actual$transcript_index)
  rows_by_tx <- rows_by_tx[match(as.character(cases$transcript_index),names(rows_by_tx))]
  record_index <- match(seq_len(nrow(records)),records$event_index)
  stopifnot(!anyNA(record_index))
  comparisons <- lapply(seq_len(nrow(cases)),function(i) {
    a <- actual[rows_by_tx[[i]],,drop=FALSE]
    observed <- lapply(seq_len(nrow(a)),function(j) {
      ids <- suppressWarnings(as.numeric(unlist(a$coding_blocks[[j]]$event_indices,use.names=FALSE)))
      stopifnot(all(is.finite(ids) & ids >= 1 & ids <= nrow(records) & ids == floor(ids)))
      matched <- record_index[ids]
      stopifnot(!anyNA(matched),all(records$seq_region[matched] == cases$seq_region[i]))
      list(cds=a$cds[j],protein=a$protein[j],count=a$carrier_count[j],contributors=records$source_id[matched])
    })
    expected <- oracle[[i]]$haplotypes
    list(expected=canonical(expected),observed=canonical(observed),
      equal=identical(canonical(expected),canonical(observed)),
      sequences_equal=identical(canonical(expected,FALSE),canonical(observed,FALSE)),
      counts_equal=sum(a$carrier_count) == oracle[[i]]$total_haplotype_count,
      unavailable_carriers=sum(a$carrier_count[is.na(a$cds)]))
  })
  saveRDS(comparisons,file.path(out,'comparisons.rds'))
  summary <- cases[setdiff(names(cases),'cds')]
  for (field in c('equal','sequences_equal','counts_equal','unavailable_carriers'))
    summary[[field]] <- vapply(comparisons,`[[`,if(field == 'unavailable_carriers') 0 else TRUE,field)
  write.csv(summary,file.path(out,'summary.csv'),row.names=FALSE)
  witness <- oracle[[1L]]$haplotypes
  corrupt <- list(duplicate=c(witness,witness[1L]),cds=witness,protein=witness,contributor=witness)
  corrupt$cds[[1L]]$cds <- paste0('C',substring(witness[[1L]]$cds,2L))
  corrupt$protein[[1L]]$protein <- paste0('X',substring(witness[[1L]]$protein,2L))
  corrupt$contributor[[1L]]$contributors <- c(witness[[1L]]$contributors,'deliberate_extra_record')
  rejected <- vapply(corrupt,function(x) !identical(canonical(x),canonical(witness)),TRUE)
  stopifnot(all(rejected),all(summary$equal[!summary$rare & summary$shape %in% c('disjoint','adjacent')]))
  write.csv(data.frame(control=names(rejected),rejected),file.path(out,'controls.csv'),row.names=FALSE)
  coverage <- strata
  stratum_key <- function(x) do.call(paste,c(x[c('shape','strand','phase','source_ploidy')],sep=':'))
  coverage$required <- opt$rare_per_stratum
  coverage$observed <- tabulate(match(stratum_key(summary[summary$rare,]),stratum_key(strata)),nrow(strata))
  stopifnot(all(coverage$observed == coverage$required))
  write.csv(coverage,file.path(out,'coverage.csv'),row.names=FALSE)
  identities <- unique(c(paths,extension,'test/duckvep/conformance/haplotype_record_differential.R',
    'test/duckvep/conformance/haplotype_oracle.pl','r/duckhtsbench/inst/benchmark_registry.tsv',
    list.files(out,full.names=TRUE)))
  jsonlite::write_json(list(source_revision=revision,extension_build_binding=binding,
    scope='source_record_geometry_not_complete_phased_conformance',
    source_artifact='haplotype_benchmark_reference',
    generator='test/duckvep/conformance/haplotype_record_differential.R',
    oracle_revisions=as.list(pins),seed=opt$seed,random_cases=opt$random_cases,
    rare_per_stratum=opt$rare_per_stratum,rare_strata=nrow(strata),rare_profiles=sum(summary$rare),
    input_records=nrow(records),profiles=nrow(cases),failures=sum(!summary$equal),
    threads=4L,output_leaves=nrow(actual),observed_carriers=sum(actual$carrier_count),
    oracle_leaves=sum(vapply(oracle,function(x) length(x$haplotypes),1L)),
    oracle_lanes=sum(vapply(oracle,`[[`,0,'total_haplotype_count')),
    sequence_failures=sum(!summary$sequences_equal),count_failures=sum(!summary$counts_equal),
    available_sequence_failures=sum(!summary$sequences_equal & summary$unavailable_carriers == 0L),
    profiles_with_unavailable=sum(summary$unavailable_carriers > 0L),
    controls_rejected=sum(rejected),sha256=as.list(vapply(identities,duckvep_evidence_sha256,''))),
    file.path(out,'receipt.json'),pretty=TRUE,auto_unbox=TRUE)
  if (!is.null(opt$extension_receipt)) duckvep_evidence_assert_checkout(root,revision)
  print(aggregate(cbind(profiles=rep(1L,nrow(summary)),failures=as.integer(!summary$equal),
    sequence_failures=as.integer(!summary$sequences_equal),unavailable_carriers) ~ shape,
    data=summary,FUN=sum),row.names=FALSE)
  if (any(!summary$equal)) stop('Source-record replay differences retained: ',out,call.=FALSE)
}
main()
