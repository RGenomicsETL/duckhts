-- Features keep their raw ID first; GFFBase-compatible unique IDs are assigned below.
CREATE TABLE features AS
WITH input AS (
  SELECT row_number() OVER () AS rid, ID, Parent, seqname, source, feature,
         start, "end", score, strand, frame, attributes_pairs
  FROM read_gff('{input}', attributes_pairs := true, attributes := ['ID', 'Parent'])
)
SELECT rid, coalesce(ID, 'row:' || rid) AS id,
       seqname AS seqid, source, feature AS featuretype, start, "end", score,
       strand, frame, Parent AS parents, attributes_pairs AS attributes
FROM input
ORDER BY seqid, start, "end";

-- @resolve-duplicate-ids
-- GFFBase create_unique claims names in file order: a row whose ID is already
-- claimed takes the next free <ID>_<n> from a per-ID counter, so a later physical
-- "e1_1" can itself be renamed ("e1_1_1"). Only rows whose ID shares a stem (the ID
-- without trailing _<digits> groups) with a duplicated ID can interact; the driver
-- replays the rule over those rows and writes the renames back.

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
