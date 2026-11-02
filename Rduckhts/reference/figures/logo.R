# Generate man/figures/logo.svg, the DuckHTS hex sticker.
# Run from r/Rduckhts: Rscript man/figures/logo.R
#
# A duck drawn by rows of beads: runs of beads join into read-like capsules,
# and capsules cluster at two peak positions. Below it, a coverage track built
# from bead stacks and a row of copy-number segments.

set.seed(20)

ink <- "#b3a300"
dark <- "#1b1f2a"
accent <- "#e8492a"

inside_ellipse <- function(x, y, cx, cy, rx, ry) ((x - cx) / rx)^2 + ((y - cy) / ry)^2 <= 1
inside_duck <- function(x, y) {
  body <- inside_ellipse(x, y, 232, 272, 128, 56)
  tail <- x >= 84 & x <= 180 & y <= 280 & y >= 196 + (x - 84) * 0.25 & y >= 300 - (x - 84) * 1.1
  head <- inside_ellipse(x, y, 322, 162, 58, 56)
  neck <- x >= 282 & x <= 352 & y >= 190 & y <= 236
  beak <- x >= 366 & x <= 446 & abs(y - 166) <= (446 - x) * 0.24
  body | tail | head | neck | beak
}
peak_weight <- function(x) exp(-((x - 196) / 26)^2) + exp(-((x - 300) / 20)^2)

bead <- function(x, y, r) sprintf('<circle cx="%.1f" cy="%.1f" r="%.1f"/>', x, y, r)
capsule <- function(x0, x1, y, h) {
  sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="%.1f"/>',
          x0 - h / 2, y - h / 2, x1 - x0 + h, h, h / 2)
}

shapes <- character()
step <- 7
for (y in seq(108, 326, by = 9)) {
  xs <- seq(84, 446, by = step)
  xs <- xs[inside_duck(xs, y)]
  k <- 1L
  while (k <= length(xs)) {
    run <- 1L
    while (k + run <= length(xs) && xs[[k + run]] - xs[[k + run - 1L]] == step) run <- run + 1L
    remaining <- run
    x <- xs[[k]]
    p_read <- 0.18 + 0.7 * peak_weight(x)
    if (runif(1) < p_read && remaining >= 3L) {
      len <- min(remaining, sample(3:7, 1L))
      shapes <- c(shapes, capsule(x, x + (len - 1L) * step, y, 4.5 + 3.5 * peak_weight(x)))
      k <- k + len
    } else {
      if (runif(1) > 0.08) shapes <- c(shapes, bead(x, y, 2.2))
      k <- k + 1L
    }
  }
}
eye <- sprintf('<circle cx="336" cy="148" r="7" fill="%s"/>', dark)

track <- character()
for (x in seq(64, 456, by = 8)) {
  depth <- 2L + round(9 * (exp(-((x - 196) / 30)^2) + 0.8 * exp(-((x - 300) / 22)^2)) + runif(1))
  for (d in seq_len(depth)) track <- c(track, bead(x, 404 - (d - 1L) * 6, 2.2))
}

segments <- c(
  capsule(64, 160, 424, 7),
  capsule(176, 236, 416, 7),
  capsule(252, 330, 424, 7),
  capsule(346, 456, 432, 7)
)
segment_accent <- sub("<rect ", sprintf('<rect fill="%s" ', accent), segments[[2L]])

svg <- c(
  '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 520 600" width="520" height="600" role="img" aria-labelledby="title">',
  '  <title id="title">duckhts</title>',
  '  <defs>',
  '    <linearGradient id="bg" x1="0" y1="0" x2="0" y2="1">',
  '      <stop offset="0" stop-color="#fff100"/>',
  '      <stop offset="1" stop-color="#fff7a8"/>',
  '    </linearGradient>',
  '  </defs>',
  '  <polygon points="260,6 514,153 514,447 260,594 6,447 6,153" fill="url(#bg)"/>',
  sprintf('  <g fill="%s">', ink),
  paste0("    ", c(shapes, track, segments[-2L])),
  '  </g>',
  paste0("  ", segment_accent),
  paste0("  ", eye),
  sprintf('  <text x="260" y="512" text-anchor="middle" fill="%s" font-family="\'DejaVu Sans Mono\', \'Menlo\', \'Consolas\', monospace" font-size="50" font-weight="700">duckhts</text>', dark),
  sprintf('  <polygon points="260,6 514,153 514,447 260,594 6,447 6,153" fill="none" stroke="%s" stroke-width="12" stroke-linejoin="round"/>', dark),
  '</svg>'
)
writeLines(svg, "man/figures/logo.svg")
