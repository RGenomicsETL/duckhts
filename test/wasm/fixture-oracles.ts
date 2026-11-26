import { readFile } from "node:fs/promises";
import path from "node:path";
import { gunzipSync } from "node:zlib";

// Independent oracles: parse the physical fixture records in Node, including
// duplicates, without DuckHTS.
export const root = path.resolve(__dirname, "../..");

export type BedRow = { chrom: string; start: number; end: number; name: string };
export type GffRow = { seqname: string; feature: string; start: number; end: number };
export type VcfRow = { chrom: string; pos: number; id: string; ref: string; alt: string };

export async function bedFixture(name = "fixture_mixed_regions.bed") {
  const text = await readFile(path.join(root, "test/data", name), "utf8");
  const rows: BedRow[] = text.trim().split("\n").map((line) => {
    const [chrom, start, end, name] = line.split("\t");
    return { chrom, start: Number(start), end: Number(end), name };
  });
  return { text, rows };
}

export async function gffFixture(name = "gff_attrs.gff3") {
  const text = await readFile(path.join(root, "test/data", name), "utf8");
  const rows: GffRow[] = text.trim().split("\n").filter((line) => !line.startsWith("#")).map((line) => {
    const [seqname, , feature, start, end] = line.split("\t");
    return { seqname, feature, start: Number(start), end: Number(end) };
  });
  return { text, rows };
}

// region_union.vcf.gz with its tabix index. overlapping(chrom, start, end)
// returns the records whose 1-based closed span [POS, END or POS+len(REF)-1]
// overlaps the 1-based closed region, as an htslib region query does.
export async function vcfFixture(name = "region_union.vcf.gz") {
  const bytes = await readFile(path.join(root, "test/data", name));
  const index = await readFile(path.join(root, "test/data", `${name}.tbi`));
  const records = gunzipSync(bytes).toString("utf8").trim().split("\n")
    .filter((line) => !line.startsWith("#")).map((line) => {
      const [chrom, pos, id, ref, alt, , , info] = line.split("\t");
      const end = info.split(";").find((field) => field.startsWith("END="));
      const row: VcfRow = { chrom, pos: Number(pos), id, ref, alt };
      return { row, end: end ? Number(end.slice(4)) : Number(pos) + ref.length - 1 };
    });
  const rows = records.map(({ row }) => row);
  const overlapping = (chrom: string, start: number, end: number) =>
    records.filter(({ row, end: last }) => row.chrom === chrom && row.pos <= end && last >= start)
      .map(({ row }) => row);
  return { bytes, index, rows, overlapping };
}

// Number of FASTQ records in a four-line-per-record file.
export async function fastqRecordCount(relativePath: string) {
  const text = await readFile(path.join(root, relativePath), "utf8");
  return text.trim().split("\n").length / 4;
}
