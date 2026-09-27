CREATE TABLE features AS
WITH input AS (
  SELECT row_number() OVER () AS rid, ID, Parent, seqname, source, feature,
         start, "end", score, strand, frame, attributes_pairs
  FROM read_gff('{input}', attributes_pairs := true, attributes := ['ID', 'Parent'])
), numbered AS (
  SELECT *, row_number() OVER (PARTITION BY ID ORDER BY rid) - 1 AS occurrence
  FROM input
)
SELECT rid,
       CASE WHEN ID IS NULL THEN 'row:' || rid
            WHEN occurrence = 0 THEN ID
            ELSE ID || '_' || occurrence END AS id,
       seqname AS seqid, source, feature AS featuretype, start, "end", score,
       strand, frame, Parent AS parents, attributes_pairs AS attributes
FROM numbered
ORDER BY seqid, start, "end";

CREATE TABLE relations AS
SELECT DISTINCT f.id AS child, trim(p) AS parent, 1 AS level
FROM features f, unnest(string_split(f.parents, ',')) AS u(p)
WHERE f.parents IS NOT NULL;

CREATE TABLE closure AS
WITH RECURSIVE c(ancestor, descendant, level) AS (
  SELECT parent, child, level FROM relations
  UNION ALL
  SELECT c.ancestor, r.child, c.level + 1
  FROM c JOIN relations r ON r.parent = c.descendant
)
SELECT DISTINCT ancestor, descendant, level FROM c;

CREATE INDEX closure_ancestor ON closure(ancestor);
CHECKPOINT;
