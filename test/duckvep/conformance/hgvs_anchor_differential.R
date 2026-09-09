#!/usr/bin/env Rscript
# Independent source-record HGVSp across stop codons, downstream codons and anchors.
suppressPackageStartupMessages({ library(DBI); library(duckdb); library(optparse) })
op <- OptionParser()
op <- add_option(op, '--extension', default = 'build/release/duckhts.duckdb_extension')
op <- add_option(op, '--vep-prefix', default = Sys.getenv('VEP_PREFIX'))
op <- add_option(op, '--tail-codons', type = 'integer', default = 64L)
op <- add_option(op, '--out', default = '')
opt <- parse_args(op, convert_hyphens_to_underscores = TRUE)
stopifnot(opt$tail_codons %in% 1:64, nzchar(opt$vep_prefix))

compare_hgvsp <- function(actual, expected) {
  stopifnot(identical(names(actual), c('id','hgvsp')),identical(names(expected),names(actual)),
    nrow(actual)==nrow(expected),!anyNA(actual$id),!anyNA(expected$id),
    !anyDuplicated(actual$id),!anyDuplicated(expected$id),setequal(actual$id,expected$id))
  expected <- expected[match(actual$id,expected$id), ]
  same <- (is.na(actual$hgvsp) & is.na(expected$hgvsp)) |
    (!is.na(actual$hgvsp) & !is.na(expected$hgvsp) & actual$hgvsp==expected$hgvsp)
  same[is.na(same)] <- FALSE
  same
}

run <- function() {
  extension <- normalizePath(opt$extension, mustWork = TRUE)
  prefix <- normalizePath(opt$vep_prefix, mustWork = TRUE)
  out <- if (nzchar(opt$out)) opt$out else tempfile('hgvs_anchor_',
    'test/duckvep/conformance/results')
  stopifnot(!dir.exists(out), dir.create(out, recursive = TRUE))
  out <- normalizePath(out)
  sha <- function(p) digest::digest(file = p, algo = 'sha256', serialize = FALSE)
  source_files <- c('test/duckvep/conformance/hgvs_anchor_differential.R',
    'src/duckvep/kernel/src/duckvep_hgvs.c','src/duckvep/kernel/src/duckvep_delta.c',
    'src/duckvep/kernel/src/duckvep_hgvs.h','src/duckvep/duckvep_annotate.c',
    'src/duckvep/duckvep_haplotype_sql.c', 'src/duckvep/duckvep_reference.c',
    'src/duckvep/duckvep_reference.h', 'src/duckvep/kernel/src/duckvep_kernel.c',
    'src/duckvep/kernel/src/duckvep_annotation_internal.h',
    'third_party/patches/htslib/0003-add-caller-buffer-faidx-fetch.patch',
    'third_party/patches/htslib/0004-check-faidx-coordinate-and-seek-arithmetic.patch',
    'third_party/htslib/faidx.c', 'third_party/htslib/htslib/faidx.h')
  source_hashes <- vapply(source_files,sha,'')
  extension_hash <- sha(extension)
  pins <- c(vep = '57ea5c52340acc1f156267f810ad162e26597082',
    variation = '2fb834b987ede3824e200197a838ce11e91aeb4b')
  for (name in names(pins)) {
    repo <- file.path('.sync', paste0('ensembl-', name))
    stopifnot(identical(system2('git', c('-C', repo, 'rev-parse', 'HEAD'), stdout = TRUE),
      unname(pins[name])), !length(system2('git', c('-C', repo, 'status', '--porcelain'), stdout = TRUE)))
  }
  command <- function(exe, args, label) {
    status <- system2(exe, shQuote(args), stdout = file.path(out, paste0(label, '.stdout')),
      stderr = file.path(out, paste0(label, '.stderr')))
    stopifnot(status == 0L)
  }
  rc <- function(x) paste0(rev(strsplit(chartr('ACGT', 'TGCA', x), '', fixed = TRUE)[[1L]]),
                          collapse = '')
  codons <- do.call(paste0, expand.grid(rep(list(c('A','C','G','T')), 3L)))
  models <- expand.grid(stop = c('TAA','TAG','TGA'), tail_codon = head(codons, opt$tail_codons),
    strand = c(1L,-1L), stringsAsFactors = FALSE)
  models$transcript_index <- models$seq_region <- seq_len(nrow(models)) - 1L
  models$chrom <- sprintf('chrA%04d', seq_len(nrow(models)))
  models$cds <- paste0('ATGGGTCCT', models$stop)
  models$post_cds <- paste0(models$tail_codon, 'GAACAATAATAACTAGCTGA')
  models$transcript <- paste0(models$cds, models$post_cds)
  models$tx_start <- 11L
  models$tx_end <- 10L + nchar(models$transcript)
  models$cds_start <- ifelse(models$strand == 1L, 11L, models$tx_end - 11L)
  models$cds_end <- models$cds_start + 11L
  genomes <- fasta <- gff <- records <- vector('list', nrow(models))
  for (i in seq_len(nrow(models))) {
    m <- models[i, ]
    genomes[[i]] <- paste0(strrep('A', 10L),
      if (m$strand == 1L) m$transcript else rc(m$transcript), strrep('A', 10L))
    fasta[[i]] <- c(paste0('>', m$chrom), genomes[[i]])
    id <- paste0('ANCHOR', i)
    attrs <- c(paste0('ID=gene:',id,';biotype=protein_coding'),
      paste0('ID=transcript:',id,';Parent=gene:',id,';biotype=protein_coding'),
      paste0('ID=exon:',id,';Parent=transcript:',id), paste0('Parent=transcript:',id))
    gff[[i]] <- paste(m$chrom,'anchor',c('gene','mRNA','exon','CDS'),
      c(rep(m$tx_start,3L),m$cds_start), c(rep(m$tx_end,3L),m$cds_end),'.',
      if (m$strand == 1L) '+' else '-',c('.','.','.','0'),attrs,sep='\t')
    r <- expand.grid(position1 = 9:12, base = c('A','C','G','T'), right = c(FALSE,TRUE),
      stringsAsFactors = FALSE)
    r$scene <- seq_len(nrow(r))
    p <- r$position1 - as.integer(!r$right)
    anchor <- substring(m$cds, p, p)
    alt <- ifelse(r$right, paste0(r$base,anchor), paste0(anchor,r$base))
    r$position <- if (m$strand == 1L) m$cds_start+p-1L else m$cds_end-p+1L
    r$reference <- if (m$strand == 1L) anchor else vapply(anchor,rc,'')
    r$alternate <- if (m$strand == 1L) alt else vapply(alt,rc,'')
    r$chrom <- m$chrom
    r$seq_region <- r$transcript_index <- m$transcript_index
    r$id <- paste0(id,'_',r$position1,'_',r$base,'_',as.integer(r$right))
    r$expected_cds <- vapply(seq_len(nrow(r)), function(j)
      paste0(substr(m$cds,1L,r$position1[j]-1L),r$base[j],substring(m$cds,r$position1[j])), '')
    for (j in seq_len(nrow(r))) {
      at <- r$position[j]
      stopifnot(substr(genomes[[i]],at,at) == r$reference[j])
      rebuilt <- paste0(substr(genomes[[i]],1L,at-1L),r$alternate[j],substring(genomes[[i]],at+1L))
      tx <- substr(rebuilt,m$tx_start,m$tx_end+1L)
      if (m$strand == -1L) tx <- rc(tx)
      stopifnot(tx == paste0(r$expected_cds[j],m$post_cds))
    }
    records[[i]] <- r
  }
  records <- do.call(rbind, records)
  stopifnot(nrow(models)==6L*opt$tail_codons,nrow(records)==192L*opt$tail_codons,
    !anyDuplicated(records$id))
  records$event_index <- seq_len(nrow(records))
  records <- records[order(records$chrom,records$position,records$id), ]
  rownames(records) <- NULL
  writeLines(unlist(fasta),file.path(out,'reference.fa'))
  writeLines(c('##gff-version 3',unlist(gff)),file.path(out,'model.gff3'))
  writeLines(c('##fileformat=VCFv4.4',paste0('##contig=<ID=',models$chrom,',length=',
    nchar(unlist(genomes)),'>'),
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    with(records,paste(chrom,position,id,reference,alternate,'.','PASS','.',sep='\t'))),
    file.path(out,'input.vcf'))
  saveRDS(list(models=models,records=records),file.path(out,'inputs.rds'))
  command('samtools',c('faidx',file.path(out,'reference.fa')),'faidx')
  command('bgzip',file.path(out,'model.gff3'),'bgzip')
  command('tabix',c('-p','gff',file.path(out,'model.gff3.gz')),'tabix')
  command('micromamba',c('list','--explicit','-p',prefix),'environment')
  libs <- paste(normalizePath(c('.sync/ensembl-vep/modules','.sync/ensembl-variation/modules',
    file.path(prefix,'share/ensembl-vep-116.0-0'))),collapse=':')
  oracle <- list()
  for (buffer in c(1L,5000L)) {
    label <- paste0('vep',buffer)
    command('micromamba',c('run','--clean-env','--env',paste0('PERL5LIB=',libs),'-p',prefix,
      'perl',normalizePath('.sync/ensembl-vep/vep'),'--gff',file.path(out,'model.gff3.gz'),
      '--fasta',file.path(out,'reference.fa'),'--format','vcf','--vcf','--hgvs',
      '--buffer_size',buffer,'--no_stats','--force_overwrite','--input_file',file.path(out,'input.vcf'),
      '--output_file',file.path(out,paste0(label,'.vcf'))),label)
    lines <- readLines(file.path(out,paste0(label,'.vcf')))
    stopifnot(any(startsWith(lines,'##VEP="v116.0" API="v116"')))
    columns <- strsplit(sub('.*Format: ([^"]+).*','\\1',
      lines[startsWith(lines,'##INFO=<ID=CSQ,')]),'|',fixed=TRUE)[[1L]]
    rows <- strsplit(lines[!startsWith(lines,'#')],'\t',fixed=TRUE)
    stopifnot(length(rows)==nrow(records))
    fields <- do.call(rbind,lapply(rows,function(x) x[1:7]))
    input_fields <- cbind(as.matrix(records[c('chrom','position','id','reference','alternate')]),'.','PASS')
    stopifnot(identical(unname(fields),unname(input_fields)))
    csq <- lapply(rows,function(x) {
      stopifnot(startsWith(x[8L],'CSQ='),!grepl(',',x[8L],fixed=TRUE))
      a <- strsplit(sub('^CSQ=','',x[8L]),'|',fixed=TRUE)[[1L]]
      length(a) <- length(columns)
      a[!is.na(a) & a==''] <- NA_character_
      # VEP percent-escapes VCF INFO strings, including equality's %3D.
      a[!is.na(a)] <- vapply(a[!is.na(a)],utils::URLdecode,'')
      setNames(a,columns)
    })
    oracle[[label]] <- do.call(rbind,csq)
  }
  stopifnot(identical(oracle[[1L]],oracle[[2L]]))
  expected <- sub('^.*:p\\.','p.',oracle[[1L]][,'HGVSp'])
  witness <- data.frame(id=records$id,hgvsp=expected)
  stopifnot(all(compare_hgvsp(witness,witness)),anyNA(expected),any(!is.na(expected)))
  duplicate <- changed <- missing <- extra <- witness
  duplicate$id[2L] <- duplicate$id[1L]
  changed$hgvsp[which(!is.na(expected))[1L]] <- 'p.DELIBERATE_CORRUPTION'
  missing$hgvsp[which(!is.na(expected))[1L]] <- NA_character_
  extra$hgvsp[which(is.na(expected))[1L]] <- 'p.(=)'
  controls <- vapply(list(drop=witness[-1L,],duplicate=duplicate,changed=changed,
    missing_hgvsp=missing,extra_hgvsp=extra),function(x)
      !isTRUE(tryCatch(all(compare_hgvsp(x,witness)),error=function(e) FALSE)),TRUE)
  stopifnot(all(controls))
  con <- dbConnect(duckdb(config=list(allow_unsigned_extensions='true')))
  on.exit(dbDisconnect(con,shutdown=TRUE))
  q <- function(x) as.character(dbQuoteString(con,x))
  dbExecute(con,paste('LOAD',q(extension)))
  dbWriteTable(con,'models',models)
  dbWriteTable(con,'events',records)
  tx <- paste('SELECT transcript_index::UINTEGER transcript_index,seq_region::UINTEGER seq_region,',
    'tx_start::UBIGINT transcript_start,tx_end::UBIGINT transcript_end,strand::TINYINT strand,',
    'transcript_index::UINTEGER gene_index,3::UBIGINT transcript_flags,',
    'cds_start::UBIGINT cds_start,cds_end::UBIGINT cds_end,cds::BLOB cds_sequence,',
    "1::UTINYINT codon_table,''::BLOB pre_cds_sequence,post_cds::BLOB post_cds_sequence",
    'FROM models ORDER BY transcript_index')
  ex <- paste('SELECT transcript_index::UINTEGER transcript_index,tx_start::UBIGINT exon_start,',
    'tx_end::UBIGINT exon_end,1::UBIGINT exon_cdna_start,',
    '(tx_end-tx_start+1)::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase',
    'FROM models ORDER BY transcript_index')
  stopifnot(dbGetQuery(con,paste0("SELECT loaded FROM duckvep_model_load('anchor',",
    q(paste('SELECT seq_region::UINTEGER seq_region,(tx_end+10)::UBIGINT sequence_length,',
      'chrom::VARCHAR seq_region_name FROM models ORDER BY seq_region')),',',q(tx),',',q(ex),
    ',reference_fasta:=',q(file.path(out,'reference.fa')),')'))$loaded)
  dbExecute(con,'SET threads=1')
  dbExecute(con,paste('CREATE TABLE independent_events AS SELECT event_index,seq_region,',
    'position,reference,alternate,NULL::UBIGINT end_position,NULL::VARCHAR structural_type,',
    'NULL::VARCHAR copy_change,NULL::UINTEGER mate_seq_region,NULL::UBIGINT mate_position',
    'FROM events ORDER BY seq_region,position,event_index'))
  dbExecute(con,paste("CREATE TABLE independent AS SELECT * FROM duckvep_annotate(",
    "'independent_events','anchor',hgvs:=true,upstream_distance:=0,downstream_distance:=0)"))
  independent <- dbGetQuery(con,'SELECT event_index,protein_hgvs FROM independent ORDER BY event_index')
  ids <- order(records$event_index)
  stopifnot(nrow(independent)==nrow(records),
    identical(as.numeric(independent$event_index),as.numeric(records$event_index[ids])))
  independent_equal <- compare_hgvsp(data.frame(id=records$id[ids],hgvsp=independent$protein_hgvs),witness)
  write.csv(data.frame(id=records$id[ids],independent,vep_hgvsp=expected[ids],
    hgvsp_equal=independent_equal),file.path(out,'independent_comparisons.csv'),row.names=FALSE,na='')
  dbExecute(con,paste('COPY independent TO',q(file.path(out,'independent.parquet')),'(FORMAT PARQUET)'))
  comparisons <- list()
  for (threads in c(1L,4L)) for (route in c('strict','source_records')) {
    dbExecute(con,paste('SET threads=',threads))
    table <- paste0('result_',route,'_',threads)
    for (scene in 1:32) {
      # One source per transcript/query prevents raw-file duplicate retention
      # from changing an independent-event oracle comparison.
      calls <- paste('SELECT event_index,seq_region,position,reference,transcript_index,',
        'event_index AS sample_index,',if (route=='strict')
          'alternate,1 alt_index,[1,1] alleles,[true,true] phase_before,NULL::BIGINT phase_set' else
          "[alternate] alternates,'1|1' gt",'FROM events WHERE scene=',scene)
      select <- paste0('SELECT * FROM duckvep_haplotypes(',q(calls),",'anchor',hgvs:=true,phase_policy:=",
        q(if (route=='strict') 'strict' else 'vep116_compat'),',input_mode:=',
        q(if (route=='strict') 'alt_events' else route),')')
      dbExecute(con,paste(if (scene==1L) paste('CREATE TABLE',table,'AS') else paste('INSERT INTO',table),select))
    }
    actual <- dbGetQuery(con,paste('SELECT contributors[1].event_index event_index,cds,protein,',
      'hgvsp,hgvsp_status,carrier_count,len(contributors) contributor_count FROM',table,'ORDER BY event_index'))
    ids <- order(records$event_index)
    stopifnot(nrow(actual)==nrow(records),identical(as.numeric(actual$event_index),as.numeric(records$event_index[ids])),
      all(actual$contributor_count==1L),all(actual$carrier_count==2L))
    same <- compare_hgvsp(data.frame(id=records$id[ids],hgvsp=gsub('[()]','',actual$hgvsp)),witness)
    comparisons[[table]] <- data.frame(threads=threads,route=route,id=records$id[ids],
      expected_cds=records$expected_cds[ids],actual,vep_hgvsp=expected[ids],hgvsp_equal=same,
      cds_equal=actual$cds==records$expected_cds[ids])
    dbExecute(con,paste('COPY',table,'TO',q(file.path(out,paste0(table,'.parquet'))),'(FORMAT PARQUET)'))
  }
  thread_differences <- setNames(vapply(c('strict','source_records'),function(route) {
    first <- paste0('SELECT * FROM result_',route,'_1')
    last <- paste0('SELECT * FROM result_',route,'_4')
    as.numeric(dbGetQuery(con,paste('SELECT count(*) n FROM ((',first,'EXCEPT ALL',last,') UNION ALL (',
      last,'EXCEPT ALL',first,'))'))$n)
  },0),c('strict','source_records'))
  comparisons <- do.call(rbind,comparisons)
  comparisons$hgvsp_equal[is.na(comparisons$hgvsp_equal)] <- FALSE
  comparisons$cds_equal[is.na(comparisons$cds_equal)] <- FALSE
  write.csv(comparisons,file.path(out,'comparisons.csv'),row.names=FALSE,na='')
  saveRDS(oracle,file.path(out,'oracle.rds'))
  stopifnot(identical(source_hashes,vapply(source_files,sha,'')),sha(extension)==extension_hash)
  for (path in source_files) {
    destination <- file.path(out,'source',path)
    dir.create(dirname(destination),recursive=TRUE,showWarnings=FALSE)
    stopifnot(file.copy(path,destination),sha(destination)==source_hashes[[path]])
  }
  files <- list.files(out,full.names=TRUE,recursive=TRUE)
  jsonlite::write_json(list(scope='independent source-record phased HGVSp versus pinned VEP 116',
    oracle_revisions=as.list(pins),source_revision=system2('git',c('rev-parse','HEAD'),stdout=TRUE),
    source_dirty=length(system2('git',c('status','--porcelain'),stdout=TRUE))>0L,
    build_binding='diagnostic_unbound',extension_sha256=extension_hash,models=nrow(models),
    source_records=nrow(records),tail_codons=opt$tail_codons,buffers=c(1L,5000L),
    comparisons=nrow(comparisons),hgvsp_mismatches=sum(!comparisons$hgvsp_equal),
    independent_hgvsp_comparisons=length(independent_equal),
    independent_hgvsp_mismatches=sum(!independent_equal),
    cds_mismatches=sum(!comparisons$cds_equal),failures_waived=0L,
    corruption_controls=as.list(controls),complete_output_thread_differences=as.list(thread_differences),
    source_sha256=as.list(source_hashes),
    input_scope='Literal source alleles, including right-retained replacements; no normalization or VCF padding repair.',
    compound_hgvs_certified=FALSE,
    sha256=as.list(setNames(vapply(files,sha,''),substring(files,nchar(out)+2L)))),
    file.path(out,'receipt.json'),pretty=TRUE,auto_unbox=TRUE)
  cat(nrow(comparisons),'comparisons;',sum(!comparisons$hgvsp_equal),'HGVSp mismatches;',
      sum(!comparisons$cds_equal),'CDS mismatches;',sum(!independent_equal),
      'independent-event HGVSp mismatches;',out,'\n')
  stopifnot(all(comparisons$hgvsp_equal),all(comparisons$cds_equal),all(thread_differences==0),
    all(independent_equal))
}
run()
