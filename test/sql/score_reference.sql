-- Test-only relational reference for DuckHTS bcftools_score (src/score_udf.c).
-- Upstream score algorithm: Giulio Genovese, bcftools +score 1.24 (MIT).
-- Required variables: score_path VARCHAR, score_summaries VARCHAR[].
-- Optional variables (native defaults unless stated): score_use (default 'GT'; 'GT' or 'DS'),
-- score_columns (default 'PLINK'), score_use_variant_id, score_q_score_thr DOUBLE[],
-- score_counts, score_samples, score_regions, score_targets, score_regions_file,
-- score_targets_file, score_regions_overlap, score_targets_overlap, score_apply_filters.
-- score_bounds adds <metric>_BOUND columns for the float32 summation error bound.
-- TSV summaries only; no include/exclude, GWAS-VCF, summary-list files or auto tag selection.
-- END/symbolic reference spans are rejected for overlap modes 1/2.
-- Custom column files, force_samples, log output, and caller-supplied PRS names are not supported.
-- Summary rows are scanned with read_tabix, which reports truncated compressed input.
-- Counts are exact. Scores sum float32 products in double precision. For n terms,
-- |native - SQL| <= (gamma(n, 2^-24) + gamma(n, 2^-53)) * A / (1-gamma(n, 2^-53)),
-- where A is SQL sum(abs(term)) and gamma(n,u) = n*u/(1-n*u).
-- The denominator accounts for rounding in A. This assumes finite normal arithmetic
-- and n*2^-24 < 1; n rather than n-1 also leaves slack for bound evaluation.
-- A weight converted from an odds ratio adds C = sum(|dose| * beta_tol + ulp(term)) twice:
-- native uses float32 logf, whose last bits depend on the platform libm.
-- query() binds a projection generated solely from the input PRS names/thresholds.
FROM query(CASE
    WHEN getvariable('score_regions') IS NOT NULL AND getvariable('score_regions_file') IS NOT NULL
        THEN error('score_reference: regions and regions_file cannot be used together')
    WHEN getvariable('score_targets') IS NOT NULL AND getvariable('score_targets_file') IS NOT NULL
        THEN error('score_reference: targets and targets_file cannot be used together')
    WHEN getvariable('score_include') IS NOT NULL OR getvariable('score_exclude') IS NOT NULL
        THEN error('score_reference: include/exclude expressions are unsupported')
    WHEN upper(coalesce(getvariable('score_columns')::VARCHAR, 'PLINK')) NOT IN
         ('PLINK', 'PLINK2', 'REGENIE', 'SAIGE', 'BOLT', 'BOLT_LMM', 'METAL', 'PGS', 'SSF', 'GWAS-SSF')
        THEN error('score_reference: unsupported columns preset')
    WHEN len(list_distinct(list_transform(getvariable('score_summaries'), lambda p:
         lower(regexp_replace(regexp_extract(p, '[^/]+$'), '(\.(gz|txt|tsv|vcf|bcf))+$', '')))))
         <> len(getvariable('score_summaries'))
        THEN error('score_reference: duplicate PRS names')
    ELSE $score$
WITH
options AS NOT MATERIALIZED (
    SELECT coalesce(getvariable('score_use'), 'GT') AS mode,
           replace(replace(upper(coalesce(getvariable('score_columns')::VARCHAR, 'PLINK')),
               'GWAS-SSF', 'SSF'), 'BOLT_LMM', 'BOLT') AS preset,
           coalesce(getvariable('score_use_variant_id'), false) AS by_id
),
summaries AS (
    SELECT path, regexp_replace(regexp_extract(path, '[^/]+$'),
               '(\.(gz|txt|tsv|vcf|bcf))+$', '') AS prs
    FROM unnest(getvariable('score_summaries')) AS files(path)
),
-- The header is the first line that is not a '#' comment; as in native, a line starting
-- with #CHR or #ID is the header with its '#' removed. read_tabix skips every '#' line,
-- so the header line is read as whole lines through read_csv.
summary_headers AS MATERIALIZED (
    SELECT path, starts_with(line, '#') AS hashed,
           string_split(CASE WHEN starts_with(line, '#') THEN substr(line, 2) ELSE line END, chr(9)) AS header_fields
    FROM query(array_to_string(list_transform(getvariable('score_summaries'),
        lambda p: 'SELECT ' || chr(39) || replace(p, chr(39), chr(39) || chr(39)) || chr(39) ||
             ' AS path, line FROM (SELECT line FROM read_csv(' ||
             chr(39) || replace(p, chr(39), chr(39) || chr(39)) || chr(39) ||
             ', columns := {''line'': ''VARCHAR''}, header := false, delim := chr(1), quote := '''', ' ||
             'escape := '''', comment := '''', auto_detect := false, compression := ' ||
             CASE WHEN regexp_matches(p, '\.(gz|bgz)$') THEN '''gzip''' ELSE '''none''' END ||
             ') WHERE NOT starts_with(line, ''#'') OR starts_with(line, ''#CHR'') OR starts_with(line, ''#ID'') LIMIT 1)'),
        ' UNION ALL '))
),
-- Each scan has its own ordinal: UNION scheduling cannot decide which duplicate wins.
-- Fields stay VARCHAR (auto_detect := false) so parsing starts from the original text, as native does.
-- read_tabix, rather than read_csv, reads the data rows so truncated compressed input is an error.
summary_rows AS MATERIALIZED (
    SELECT * FROM query(array_to_string(list_transform(getvariable('score_summaries'),
        lambda p: 'SELECT ' || chr(39) || replace(p, chr(39), chr(39) || chr(39)) || chr(39) ||
             ' AS path, row_number() OVER () AS ordinal, list_value(*COLUMNS(*)) AS field_values FROM read_tabix(' ||
             chr(39) || replace(p, chr(39), chr(39) || chr(39)) || chr(39) ||
             ', header := false, auto_detect := false)'), ' UNION ALL '))
),
-- Alias order follows the pinned native presets.
aliases(preset, role, names) AS (
    VALUES
        ('PLINK', 'CHR', ['CHR']), ('PLINK', 'BP', ['BP']),
        ('PLINK', 'SNP', ['SNP']), ('PLINK', 'A1', ['A1']),
        ('PLINK', 'BETA', ['BETA']), ('PLINK', 'OR', ['OR']), ('PLINK', 'P', ['P']),
        ('PLINK2', 'LP', ['LOG10_P']),
        ('PLINK2', 'CHR', ['CHROM']), ('PLINK2', 'BP', ['POS']),
        ('PLINK2', 'SNP', ['ID']), ('PLINK2', 'A1', ['A1']),
        ('PLINK2', 'BETA', ['BETA']), ('PLINK2', 'OR', ['OR']), ('PLINK2', 'P', ['P']),
        ('REGENIE', 'CHR', ['CHROM']), ('REGENIE', 'BP', ['GENPOS']),
        ('REGENIE', 'SNP', ['ID']), ('REGENIE', 'A1', ['ALLELE1']),
        ('REGENIE', 'BETA', ['BETA']), ('REGENIE', 'LP', ['LOG10P']),
        ('SAIGE', 'CHR', ['CHR']), ('SAIGE', 'BP', ['POS']),
        ('SAIGE', 'SNP', ['SNPID', 'markerID']), ('SAIGE', 'A1', ['Allele2']),
        ('SAIGE', 'BETA', ['BETA']), ('SAIGE', 'P', ['p.value']),
        ('BOLT', 'CHR', ['CHR']), ('BOLT', 'BP', ['BP']),
        ('BOLT', 'SNP', ['SNP']), ('BOLT', 'A1', ['ALLELE1']),
        ('BOLT', 'BETA', ['BETA']), ('BOLT', 'P', ['P_LINREG', 'P_BOLT_LMM_INF', 'P_BOLT_LMM']),
        ('METAL', 'CHR', ['Chromosome']), ('METAL', 'BP', ['Position']),
        ('METAL', 'SNP', ['MarkerName']), ('METAL', 'A1', ['Allele1']),
        ('METAL', 'BETA', ['Effect']), ('METAL', 'P', ['P-value']), ('METAL', 'LP', ['log(P)']),
        ('SSF', 'CHR', ['chromosome']), ('SSF', 'BP', ['base_pair_location']),
        ('SSF', 'SNP', ['variant_id', 'rsid', 'rs_id']), ('SSF', 'A1', ['effect_allele']),
        ('SSF', 'BETA', ['beta']), ('SSF', 'OR', ['odds_ratio', 'hazard_ratio']), ('SSF', 'P', ['p_value']),
        ('PGS', 'OR', ['OR', 'HR']),
        ('PGS', 'CHR', ['chr_name']), ('PGS', 'BP', ['chr_position']),
        ('PGS', 'SNP', ['rsID']), ('PGS', 'A1', ['effect_allele']),
        ('PGS', 'BETA', ['effect_weight'])
),
selected_aliases AS (
    SELECT role, names FROM aliases, options WHERE aliases.preset = options.preset
),
-- Each role maps to the 1-based position of its header field; the last preset alias present wins.
header_mapping AS (
    SELECT path, map_from_entries(list({'key': role, 'value':
        list_position(list_transform(header_fields, lambda k: lower(k)), lower(
            list_last(list_filter(names, lambda a: list_contains(list_transform(header_fields, lambda k: lower(k)), lower(a))))))})) AS mapping,
           bool_or(hashed) AS hashed, any_value(header_fields) AS header_fields
    FROM summary_headers CROSS JOIN selected_aliases
    GROUP BY path
),
fields AS (
    SELECT path, ordinal, map_from_entries(list_transform(map_entries(mapping), lambda m:
        {'key': m['key'], 'value': field_values[m['value']]})) AS f,
        mapping['CHR'] IS NULL OR mapping['BP'] IS NULL AS lack_position
    FROM summary_rows JOIN header_mapping USING (path)
    -- Without a '#' header, read_tabix returns the header line as its first row, unless a
    -- tabix index with line_skip already skipped it. This filter stays here: joining
    -- summary_headers into summary_rows changed DuckDB's estimates so the marker hash join
    -- built over the V x S calls instead of the markers (9 GiB at 4 PRS).
    WHERE NOT (ordinal = 1 AND NOT hashed AND field_values = header_fields)
      AND CASE WHEN mapping['A1'] IS NULL OR (mapping['BETA'] IS NULL AND mapping['OR'] IS NULL)
               THEN error('score_reference: missing effect allele/weight columns in summary header: ' || path)
               WHEN (coalesce(getvariable('score_use_variant_id'), false) OR lack_position) AND mapping['SNP'] IS NULL
               THEN error('score_reference: missing marker name column in summary header: ' || path)
               ELSE true END
),
parsed_markers AS (
    SELECT path, ordinal, f['CHR'] AS chrom, try_cast(f['BP'] AS BIGINT) AS pos,
           f['SNP'] AS id, upper(f['A1']) AS effect,
           coalesce(nullif(try_cast(f['BETA'] AS FLOAT), 'NaN'::FLOAT),
               CASE WHEN try_cast(f['OR'] AS FLOAT) > 0 THEN ln(try_cast(f['OR'] AS FLOAT))::FLOAT END) AS beta,
           -- Two float32 ulps of a weight converted from an odds ratio (native logf); 0 for BETA.
           CASE WHEN nullif(try_cast(f['BETA'] AS FLOAT), 'NaN'::FLOAT) IS NOT NULL THEN 0
                ELSE 2 * pow(2, floor(log2(greatest(abs(beta), pow(2, -126)))) - 23) END AS beta_tol,
           -- Native parses LP first and falls back to P when LP is missing, unparsable or NaN.
           CASE WHEN NOT isnan(coalesce(try_cast(f['LP'] AS FLOAT), 'NaN'::FLOAT)) THEN try_cast(f['LP'] AS FLOAT)
                WHEN try_cast(split_part(lower(f['P']), 'e', 1) AS FLOAT) > 0
                THEN (-log10(try_cast(split_part(lower(f['P']), 'e', 1) AS FLOAT))::FLOAT -
                      coalesce(try_cast(nullif(split_part(lower(f['P']), 'e', 2), '') AS FLOAT), 0))::FLOAT
                ELSE NULL END AS lp,
           -- Native computes -log10f(mantissa) in float32. Its last bit depends on the platform
           -- libm (glibc 2.39 log10f is not correctly rounded), which SQL cannot reproduce.
           -- lp_exact is the double-precision value, and lp_tol allows two float32 ulps of the
           -- logarithm and two of the result. Both are NULL when the LP column was read directly.
           CASE WHEN NOT isnan(coalesce(try_cast(f['LP'] AS FLOAT), 'NaN'::FLOAT)) THEN NULL
                WHEN try_cast(split_part(lower(f['P']), 'e', 1) AS FLOAT) > 0
                THEN -log10(try_cast(split_part(lower(f['P']), 'e', 1) AS FLOAT)::DOUBLE) -
                     coalesce(try_cast(nullif(split_part(lower(f['P']), 'e', 2), '') AS FLOAT), 0)
                ELSE NULL END AS lp_exact,
           2 * (pow(2, floor(log2(greatest(abs(lp_exact + coalesce(try_cast(nullif(split_part(lower(f['P']), 'e', 2), '') AS FLOAT), 0)), pow(2, -126)))) - 23) +
                pow(2, floor(log2(greatest(abs(lp_exact), pow(2, -126)))) - 23)) AS lp_tol,
           by_id OR lack_position AS match_id
    FROM fields CROSS JOIN options
),
contigs AS MATERIALIZED (
    SELECT id FROM read_hts_header(getvariable('score_path')) WHERE record_type = 'contig'
),
chromosome_aliases(name, candidates) AS (
    VALUES ('23', ['X', 'chrX']), ('25', ['X', 'chrX']), ('XY', ['X', 'chrX']),
           ('XX', ['X', 'chrX']), ('PAR1', ['X', 'chrX']), ('PAR2', ['X', 'chrX']),
           ('24', ['Y', 'chrY']), ('26', ['MT', 'chrM']), ('MT', ['MT', 'chrM']), ('chrM', ['MT', 'chrM'])
),
-- Native keeps a literal '.' summary ID and, matching by ID, pairs it with VCF records whose
-- ID is '.'. read_tabix returns NULL both for '.' and for an empty field, which native skips,
-- and read_geno returns NULL for a '.' ID. The reference rejects the one ambiguous case (in
-- markers): a by-ID summary row without an ID while the VCF has records without one.
vcf_missing_id AS (
    SELECT EXISTS (SELECT 1 FROM read_bcf(getvariable('score_path')) WHERE ID IS NULL) AS present
),
markers AS (
    SELECT m.*, CASE WHEN match_id THEN 'I:' || m.id ELSE
        'P:' || coalesce(exact.id, stripped.id, prefixed.id, special.id, fallback.id) || ':' || pos END AS marker_key
    FROM parsed_markers m
    LEFT JOIN contigs exact ON m.chrom = exact.id
    LEFT JOIN contigs stripped ON regexp_replace(m.chrom, '^chr', '') = stripped.id
    LEFT JOIN contigs prefixed ON len(m.chrom) <= 2 AND 'chr' || m.chrom = prefixed.id
    LEFT JOIN chromosome_aliases a ON m.chrom = a.name
    LEFT JOIN contigs special ON a.candidates[1] = special.id
    LEFT JOIN contigs fallback ON a.candidates[2] = fallback.id
    WHERE beta IS NOT NULL AND NOT isnan(beta) AND effect IS NOT NULL
      AND CASE WHEN NOT match_id THEN pos IS NOT NULL
               WHEN m.id IS NOT NULL THEN true
               WHEN (SELECT present FROM vcf_missing_id)
               THEN error('score_reference: summary rows and VCF records without IDs under ID matching; ' ||
                          'native matches a literal ''.'' ID')
               ELSE false END
    QUALIFY row_number() OVER (PARTITION BY path, marker_key ORDER BY ordinal) = 1
),
thresholds AS (
    SELECT threshold,
           CASE WHEN getvariable('score_q_score_thr') IS NULL THEN ''
                -- Native names the column from pow(10, -lp) of the float32-parsed threshold.
                ELSE '_p' || printf('%g', pow(10, log10(CASE WHEN threshold::FLOAT = 0 THEN threshold
                                                         ELSE threshold::FLOAT::DOUBLE END))) END AS suffix,
           CASE WHEN threshold <= 0 OR threshold > 1 THEN error('score_reference: invalid q_score_thr')
                ELSE -log10(CASE WHEN threshold::FLOAT = 0 THEN threshold ELSE threshold::FLOAT::DOUBLE END) END AS lp
    FROM unnest(coalesce(getvariable('score_q_score_thr')::DOUBLE[], [1.0])) AS q(threshold)
),
samples AS (
    SELECT s.* FROM read_bcf_samples(getvariable('score_path'), samples := getvariable('score_samples')) s
    WHERE CASE WHEN EXISTS (SELECT 1 FROM read_hts_header(getvariable('score_path'))
                            WHERE record_type = 'FORMAT' AND id = coalesce(getvariable('score_use'), 'GT'))
               THEN true ELSE error('score_reference: required FORMAT/' || coalesce(getvariable('score_use'), 'GT') ||
                                    ' is absent from header: ' || getvariable('score_path')) END
),
-- GT and DS use a common normalized call layout without decoding unrelated FORMAT tags.
-- The per-record widths are computed before unnest: beside unnest, a list lambda over
-- calls is evaluated once per call, which makes the work grow with samples squared.
calls AS (
    SELECT * FROM query(
        'SELECT record_index, CHROM, POS, ID, REF, ALT, FILTER, gt_width, ds_width, c.sample_index, c.alleles, ' ||
        CASE coalesce(getvariable('score_use'), 'GT')
            WHEN 'GT' THEN 'NULL::FLOAT[]'
            WHEN 'DS' THEN 'CASE WHEN typeof(c.format.DS) LIKE ''%[]'' THEN c.format.DS::FLOAT[] ELSE [c.format.DS::FLOAT] END'
            ELSE error('score_reference: supported use values are GT and DS') END ||
        ' AS dosage FROM (SELECT * EXCLUDE (calls), unnest(calls) AS c FROM (SELECT *, ' ||
        'list_max(list_transform(calls, lambda c: len(c.alleles))) AS gt_width, ' ||
        CASE coalesce(getvariable('score_use'), 'GT') WHEN 'DS' THEN
            'list_max(list_transform(calls, lambda c: CASE WHEN typeof(c.format.DS) LIKE ''%[]'' THEN len(c.format.DS::FLOAT[]) ELSE 1 END))'
            ELSE 'NULL::BIGINT' END || ' AS ds_width FROM read_geno(' ||
        chr(39) || replace(getvariable('score_path'), chr(39), chr(39) || chr(39)) || chr(39) ||
        ', scan_mode := ''sequential'', include_filter := true, samples := ' ||
        CASE WHEN getvariable('score_samples') IS NULL THEN 'NULL' ELSE
            chr(39) || replace(getvariable('score_samples')::VARCHAR, chr(39), chr(39) || chr(39)) || chr(39) END ||
        ', format_fields := ' || CASE WHEN coalesce(getvariable('score_use'), 'GT') = 'DS'
            THEN '[''DS'']' ELSE '[]' END || ')))')
),
filtered_calls AS (
    SELECT * FROM calls
    WHERE coalesce(getvariable('score_apply_filters')::VARCHAR, '') = ''
       -- htslib's synced reader (has_filter): '.' selects records without FILTER entries;
       -- PASS selects only records whose FILTER is PASS. Filter names are header IDs.
       OR len(list_intersect(coalesce(FILTER, ['.']), string_split(getvariable('score_apply_filters')::VARCHAR, ','))) > 0
),
selection(kind, spec, path, mode) AS (
    VALUES ('regions', getvariable('score_regions')::VARCHAR, getvariable('score_regions_file')::VARCHAR,
            coalesce(getvariable('score_regions_overlap'), 1)),
           ('targets', getvariable('score_targets')::VARCHAR, getvariable('score_targets_file')::VARCHAR,
            coalesce(getvariable('score_targets_overlap'), 0))
),
-- Region and target files are read as whole lines, as native reads them with hts_getline
-- (score_regions_from_file): leading whitespace is skipped, then empty and '#' lines.
-- score_split_fields splits on runs of tabs and ends a field at '\r', so empty fields vanish
-- and a literal '.' stays a field. BED is a case-insensitive .bed or .bed.gz name.
file_rows AS (
    SELECT selection.kind, selection.mode, f.path, f.data_row, f.fields, is_bed
    FROM query(coalesce(nullif(array_to_string(list_transform(
        list_filter([getvariable('score_regions_file')::VARCHAR, getvariable('score_targets_file')::VARCHAR], lambda p: p IS NOT NULL),
        lambda p: 'SELECT ' || chr(39) || replace(p, chr(39), chr(39) || chr(39)) || chr(39) ||
             ' AS path, row_number() OVER () AS data_row, ' ||
             'list_filter(regexp_split_to_array(trimmed, ''[\t\r]''), lambda x: x <> '''') AS fields ' ||
             'FROM (SELECT ltrim(line, '' '' || chr(9) || chr(11) || chr(12) || chr(13)) AS trimmed FROM read_csv(' ||
             chr(39) || replace(p, chr(39), chr(39) || chr(39)) || chr(39) ||
             ', columns := {''line'': ''VARCHAR''}, header := false, delim := chr(1), quote := '''', ' ||
             'escape := '''', comment := '''', auto_detect := false, compression := ' ||
             CASE WHEN regexp_matches(p, '\.(gz|bgz)$') THEN '''gzip''' ELSE '''none''' END ||
             ')) WHERE coalesce(trimmed, '''') <> '''' AND NOT starts_with(trimmed, ''#'')'),
        ' UNION ALL '), ''), 'SELECT NULL::VARCHAR AS path, NULL::BIGINT AS data_row, []::VARCHAR[] AS fields WHERE false')) AS f
    JOIN selection USING (path),
    LATERAL (SELECT regexp_matches(lower(path), '\.bed(\.gz)?$') AS is_bed)
),
-- A one-field row of a non-BED file is a region expression, parsed like the string option.
interval_items AS (
    SELECT kind, mode, item, coalesce(contigs.id, regexp_extract(item, '^(\{[^}]+\}|[^:]+)', 1)) AS contig
    FROM (SELECT kind, mode, spec AS expression FROM selection WHERE coalesce(spec, '') <> ''
          UNION ALL
          SELECT kind, mode, fields[1] FROM file_rows WHERE len(fields) = 1 AND NOT is_bed) AS expressions
    CROSS JOIN unnest(regexp_extract_all(expression, '(\{[^}]*\}[^,]*|[^,]+)')) AS items(item)
    LEFT JOIN contigs ON item = contigs.id
),
-- Region coordinates follow htslib hts_parse_region with HTS_PARSE_LIST | HTS_PARSE_ONE_COORD
-- (third_party/htslib/hts.c), as native score_parse_region_list uses it. Each number is
-- read as hts_parse_decimal reads it: leading whitespace, a sign, digits, an optional
-- fraction, then an exponent or a k/M/G suffix; without digits the text is left unread.
-- The same pattern parses the start and, after a '-', the end.
region_starts AS (
    SELECT kind, mode, contig, trim(contig, '{}') AS chrom,
           substr(item, len(contig) + 1) AS after_contig,
           CASE WHEN starts_with(after_contig, ':') THEN substr(after_contig, 2) END AS coordinates,
           regexp_extract(coordinates,
               '^[ \t\n\v\f\r]*([+-]?)([0-9]*)(?:\.([0-9]*))?(?:[eE]([+-]?[0-9]*)|([kKmMgG]))?(.*)$',
               ['sign', 'whole', 'fraction', 'exponent', 'suffix', 'rest']) AS start_number,
           CASE WHEN start_number['whole'] || start_number['fraction'] = '' THEN coordinates
                ELSE start_number['rest'] END AS after_start
    FROM interval_items
),
region_numbers AS (
    SELECT *,
           CASE WHEN starts_with(after_start, '-') THEN substr(after_start, 2) END AS end_text,
           regexp_extract(end_text,
               '^[ \t\n\v\f\r]*([+-]?)([0-9]*)(?:\.([0-9]*))?(?:[eE]([+-]?[0-9]*)|([kKmMgG]))?(.*)$',
               ['sign', 'whole', 'fraction', 'exponent', 'suffix', 'rest']) AS end_number,
           CASE WHEN end_number['whole'] || end_number['fraction'] = '' THEN end_text
                ELSE end_number['rest'] END AS after_end,
           -- hts_parse_decimal: the digits scaled by 10^(exponent - fraction digits),
           -- truncated toward zero. Values beyond int64 are outside the reference.
           list_transform([start_number, end_number], lambda g: CASE WHEN g IS NULL THEN NULL ELSE
               list_transform([{'digits': ltrim(g['whole'] || g['fraction'], '0'),
                                'scale': CASE lower(g['suffix']) WHEN 'k' THEN 3 WHEN 'm' THEN 6 WHEN 'g' THEN 9
                                         WHEN '' THEN CASE WHEN len(g['exponent']) > 6
                                             THEN error('score_reference: region coordinate exponent is outside the reference')
                                             ELSE coalesce(try_cast(g['exponent'] AS INTEGER), 0) END END
                                         - len(g['fraction'])}],
                   lambda n: CASE WHEN n['digits'] = '' THEN 0::HUGEINT
                                  WHEN len(n['digits']) + greatest(n['scale'], 0) > 18
                                      THEN error('score_reference: region coordinate is outside the int64 range of the reference')
                                  WHEN n['scale'] >= 0 THEN (n['digits'] || repeat('0', n['scale']))::HUGEINT
                                  WHEN len(n['digits']) + n['scale'] <= 0 THEN 0::HUGEINT
                                  ELSE left(n['digits'], len(n['digits']) + n['scale'])::HUGEINT END
                             * CASE WHEN g['sign'] = '-' THEN -1 ELSE 1 END)[1] END) AS number_values
    FROM region_starts
),
text_intervals AS (
    SELECT kind, mode, chrom,
           number_values[1] - 1 AS start0,
           CASE WHEN number_values[2] = 0 THEN 9223372036854775807 ELSE number_values[2] END AS end_value,
           CASE WHEN contains(contig, '{') AND NOT contains(contig, '}')
                     OR after_contig <> '' AND coordinates IS NULL
                    THEN error('score_reference: invalid region/target expression')
                WHEN coordinates IS NULL THEN [1, 9223372036854775807]
                -- chr:-N is chr:1-N, and chr:0 or chr: is the whole contig.
                WHEN start0 < 0 AND start0 <> -1 AND starts_with(after_start, '-') AND coordinates <> ''
                    THEN error('score_reference: invalid region/target expression')
                WHEN start0 < 0 AND (after_start = '' OR regexp_matches(after_start, '^[0-9]'))
                    THEN [1, CASE WHEN start0 = -1 THEN 9223372036854775807 ELSE -(start0 + 1) END]
                WHEN start0 < -1 THEN error('score_reference: invalid region/target expression')
                WHEN after_start = '' THEN [start0 + 1, start0 + 1]
                WHEN NOT starts_with(after_start, '-') OR after_end <> '' OR start0 >= end_value
                    THEN error('score_reference: invalid region/target expression')
                ELSE [start0 + 1, end_value] END::BIGINT[] AS bounds
    FROM region_numbers
),
file_intervals AS (
    SELECT kind, mode, fields[1] AS chrom,
           CASE WHEN try_cast(fields[2] AS BIGINT) IS NULL
                     OR try_cast(coalesce(fields[3], fields[2]) AS BIGINT) IS NULL
                     OR try_cast(fields[2] AS BIGINT) < CASE WHEN is_bed THEN 0 ELSE 1 END
                     OR try_cast(coalesce(fields[3], fields[2]) AS BIGINT) <
                        try_cast(fields[2] AS BIGINT) + CASE WHEN is_bed THEN 1 ELSE 0 END
                THEN error('score_reference: invalid interval in ' || path || ' at data row ' || data_row)
                ELSE fields[2]::BIGINT + CASE WHEN is_bed THEN 1 ELSE 0 END END AS lo,
           try_cast(coalesce(fields[3], fields[2]) AS BIGINT) AS hi
    FROM file_rows
    WHERE len(fields) > 1 OR is_bed
),
intervals AS (
    SELECT kind, mode, chrom, bounds[1] AS lo, bounds[2] AS hi
    FROM text_intervals
    UNION ALL
    SELECT * FROM file_intervals
),
selected_calls AS (
    SELECT c.* FROM filtered_calls c
    WHERE
$score$ || CASE WHEN coalesce(getvariable('score_regions')::VARCHAR, '') = ''
                    AND coalesce(getvariable('score_targets')::VARCHAR, '') = ''
                    AND getvariable('score_regions_file') IS NULL AND getvariable('score_targets_file') IS NULL
               THEN 'true' ELSE $selection$
       NOT EXISTS (
        SELECT 1 FROM selection s
        -- A file without data rows gives native an empty region string, and
        -- score_parse_region_list leaves an empty string inactive: every record stays.
        WHERE (coalesce(s.spec, '') <> '' OR EXISTS (SELECT 1 FROM file_rows f WHERE f.kind = s.kind))
          AND (EXISTS (SELECT 1 FROM intervals i WHERE i.kind = s.kind AND i.chrom = c.CHROM
                       -- Native modes 1 and 2 both use the record reference span, including its anchor.
                       AND c.POS <= i.hi
                       AND c.POS + CASE WHEN i.mode = 0 THEN 0
                           WHEN EXISTS (SELECT 1 FROM read_hts_header(getvariable('score_path'))
                                        WHERE record_type = 'INFO' AND id = 'END')
                                OR list_any_value(list_filter(c.ALT, lambda a: regexp_matches(a, '[<>\[\]]'))) IS NOT NULL
                           THEN error('score_reference: END/symbolic span overlap is unsupported: ' || getvariable('score_path'))
                           ELSE len(c.REF) - 1 END >= i.lo) = false)
    )
$selection$ END || $score$
),
matched AS (
    SELECT c.*, m.path, m.beta, m.beta_tol, m.lp, m.lp_exact, m.lp_tol,
           list_position(list_prepend(REF, ALT), effect) - 1 AS allele,
           CASE WHEN mode = 'GT' THEN gt_width >= 2
                WHEN ds_width IS NULL THEN false
                WHEN ds_width <> len(ALT) THEN error('score_reference: invalid DS FORMAT cardinality at ' ||
                    getvariable('score_path') || ':' || c.CHROM || ':' || c.POS)
                ELSE true END AS available,
           CASE WHEN mode = 'GT' THEN alleles[1] IS NULL OR (len(alleles) >= 2 AND alleles[2] IS NULL)
                -- Native score_is_missing() also treats an ordinary NaN dosage as missing.
                ELSE len(list_filter(dosage, lambda x: x IS NOT NULL AND NOT isnan(x))) <> len(ALT) END AS missing,
           CASE WHEN mode = 'GT' THEN len(list_filter(alleles[1:2], lambda x: x = allele))::FLOAT
                WHEN allele > 0 THEN dosage[allele]
                ELSE list_reduce(list_prepend(2::FLOAT, dosage), lambda a, b: (a - b)::FLOAT) END AS dose
    FROM (SELECT *, unnest(['I:' || ID, 'P:' || CHROM || ':' || POS]) AS marker_key FROM selected_calls) c
    JOIN markers m USING (marker_key)
    CROSS JOIN options
),
terms AS (
    SELECT sample_index, path, t.suffix, (beta * dose)::FLOAT::DOUBLE AS term,
           CASE WHEN beta_tol = 0 THEN 0
                ELSE abs(dose) * beta_tol + pow(2, floor(log2(greatest(abs(term), pow(2, -126)))) - 23) END AS conversion
    FROM matched CROSS JOIN thresholds t
    WHERE available AND NOT missing AND allele >= 0
      AND (getvariable('score_q_score_thr') IS NULL OR matched.lp IS NULL OR isnan(matched.lp) OR
           CASE WHEN abs(matched.lp_exact - t.lp) <= matched.lp_tol
                THEN error('score_reference: a p-value is within float32 log10 rounding of a q_score_thr ' ||
                           'threshold; the native selection depends on the platform log10f')
                ELSE matched.lp >= t.lp END)
),
reductions AS (
    SELECT sample_index, path, suffix, sum(term) AS score, count(*)::BIGINT AS n, sum(abs(term)) AS magnitude,
           sum(conversion) AS conversion
    FROM terms GROUP BY sample_index, path, suffix
),
results AS (
    SELECT sample_name AS SAMPLE, prs || t.suffix AS metric,
           coalesce(score, 0) AS score, coalesce(n, 0)::BIGINT AS n,
           prs || '_CNT' || t.suffix AS count_metric,
           CASE WHEN n * pow(2, -24) >= 1 THEN 'Infinity'::DOUBLE
                ELSE coalesce(magnitude * (n * pow(2, -24) / (1 - n * pow(2, -24)) +
                                           n * pow(2, -53) / (1 - n * pow(2, -53))) /
                    (1 - n * pow(2, -53) / (1 - n * pow(2, -53))) + 2 * conversion, 0) END AS bound
    FROM samples CROSS JOIN summaries CROSS JOIN thresholds t
    LEFT JOIN reductions r ON samples.sample_index = r.sample_index AND summaries.path = r.path AND t.suffix = r.suffix
)
SELECT SAMPLE,
$score$ ||
-- The one-element list binds the generated metric names once for the name check and the
-- projection. Native rejects any repeated generated column name, ignoring case
-- (score_validate_prs_names); here a repeated threshold would count each term twice.
list_transform([
    flatten(list_transform(
        list_transform(getvariable('score_summaries'), lambda p: regexp_replace(regexp_extract(p, '[^/]+$'),
            '(\.(gz|txt|tsv|vcf|bcf))+$', '')),
        lambda p: list_transform(coalesce(getvariable('score_q_score_thr')::DOUBLE[], [1.0]), lambda q:
            {'metric': p || CASE WHEN getvariable('score_q_score_thr') IS NULL THEN '' ELSE '_p' || printf('%g', pow(10, log10(CASE WHEN q::FLOAT = 0 THEN q ELSE q::FLOAT::DOUBLE END))) END,
             'count_metric': p || '_CNT' || CASE WHEN getvariable('score_q_score_thr') IS NULL THEN '' ELSE '_p' || printf('%g', pow(10, log10(CASE WHEN q::FLOAT = 0 THEN q ELSE q::FLOAT::DOUBLE END))) END})))],
  lambda metrics: CASE
    WHEN len(list_distinct(list_transform(
             ['SAMPLE'] || list_transform(metrics, lambda m: m['metric']) ||
             CASE WHEN coalesce(getvariable('score_counts'), false)
                  THEN list_transform(metrics, lambda m: m['count_metric']) ELSE []::VARCHAR[] END,
             lambda name: lower(name))))
         <> 1 + len(metrics) * CASE WHEN coalesce(getvariable('score_counts'), false) THEN 2 ELSE 1 END
        THEN error('score_reference: duplicate generated output column names')
    ELSE array_to_string(list_transform(metrics,
    lambda m: 'max(score) FILTER (WHERE metric = ' || chr(39) || replace(m['metric'], chr(39), chr(39) || chr(39)) || chr(39) ||
         ') AS "' || replace(m['metric'], '"', '""') || '"' ||
         CASE WHEN coalesce(getvariable('score_counts'), false) THEN
             ', max(n) FILTER (WHERE metric = ' || chr(39) || replace(m['metric'], chr(39), chr(39) || chr(39)) || chr(39) ||
             ') AS "' || replace(m['count_metric'], '"', '""') || '"' ELSE '' END ||
         CASE WHEN coalesce(getvariable('score_bounds'), false) THEN
             ', max(bound) FILTER (WHERE metric = ' || chr(39) || replace(m['metric'], chr(39), chr(39) || chr(39)) || chr(39) ||
             ') AS "' || replace(m['metric'], '"', '""') || '_BOUND"' ELSE '' END), ', ') END)[1] ||
$score$
FROM results GROUP BY SAMPLE ORDER BY SAMPLE
$score$ END)
