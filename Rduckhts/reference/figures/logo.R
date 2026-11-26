# Generate man/figures/logo.svg, the DuckHTS hex sticker.
# Run from r/Rduckhts: Rscript man/figures/logo.R
#
# A duck drawn as a halftone of horizontal strokes, in the style of the DuckDB
# 2.0 artwork: closely spaced rows in which a shade map sets the mark, from
# sparse dots through thin lines to thick read-like capsules. Below it, a
# coverage track in the same marks and a row of copy-number segments.

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

# Fraction of a ring of probes around (x, y) that falls outside the duck: near
# its outline the shade rises, so the silhouette reads as a darker contour.
edge_weight <- function(x, y, radius = 9) {
  angles <- seq(0, 2 * pi, length.out = 13L)[-13L]
  outside <- vapply(angles, function(a) {
    !inside_duck(x + radius * cos(a), y + radius * sin(a))
  }, logical(length(x)))
  rowMeans(matrix(outside, nrow = length(x)))
}

# Shade in [0, 1]: a light body, a dark folded wing, a darker crown and nape,
# the outline, and a little grain.
shade <- function(x, y) {
  wing <- exp(-(((x - 214) / 74)^2 + ((y - 262) / 24)^2))
  wing_fold <- exp(-(((x - 168) / 40)^2 + ((y - 252) / 12)^2))
  crown <- exp(-(((x - 300) / 34)^2 + ((y - 132) / 22)^2))
  belly <- exp(-(((x - 250) / 90)^2 + ((y - 312) / 14)^2))
  value <- 0.36 + 0.45 * wing + 0.25 * wing_fold + 0.45 * crown +
    0.25 * belly + 0.85 * edge_weight(x, y) + stats::rnorm(length(x), sd = 0.06)
  pmin(pmax(value, 0), 1)
}

bead <- function(x, y, r) sprintf('<circle cx="%.1f" cy="%.1f" r="%.2f"/>', x, y, r)
capsule <- function(x0, x1, y, h) {
  sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.2f" rx="%.2f"/>',
          x0 - h / 2, y - h / 2, x1 - x0 + h, h, h / 2)
}

# Mark class of a shade: 0 blank, 1 dot, 2 hairline, 3 line, 4 capsule.
mark_class <- function(s) findInterval(s, c(0.12, 0.3, 0.5, 0.72))
mark_height <- c(0, 0, 1.0, 1.9, 3.4)

# One halftone row: samples every `step` pixels inside `keep`, consecutive
# samples of one class merge into a stroke, broken now and then into dashes.
halftone_row <- function(xs, y, s) {
  cls <- mark_class(s)
  out <- character()
  k <- 1L
  while (k <= length(xs)) {
    c0 <- cls[[k]]
    run <- 1L
    max_run <- sample(4:14, 1L)
    while (k + run <= length(xs) && cls[[k + run]] == c0 &&
           xs[[k + run]] - xs[[k + run - 1L]] < 3 && run < max_run) {
      run <- run + 1L
    }
    x0 <- xs[[k]]
    x1 <- xs[[k + run - 1L]]
    if (c0 == 1L) {
      for (x in seq(x0, x1, by = 3)) out <- c(out, bead(x, y, 0.7))
    } else if (c0 >= 2L) {
      out <- c(out, capsule(x0, max(x1 - 1.5, x0), y, mark_height[[c0 + 1L]]))
    }
    k <- k + run
  }
  out
}

step <- 2
shapes <- character()
for (y in seq(104, 330, by = 3.3)) {
  xs <- seq(80, 450, by = step)
  xs <- xs[inside_duck(xs, y)]
  if (length(xs)) shapes <- c(shapes, halftone_row(xs, y, shade(xs, y)))
}
eye <- sprintf('<circle cx="336" cy="148" r="6.5" fill="%s"/>', dark)

# Coverage track: per column a stack of rows whose height follows depth, the
# top rows lighter, drawn with the same marks.
depth_at <- function(x) {
  2 + 13 * exp(-((x - 196) / 30)^2) + 10 * exp(-((x - 300) / 22)^2)
}
track <- character()
for (y in seq(404, 330, by = -3.3)) {
  xs <- seq(60, 460, by = step)
  height <- depth_at(xs) * 3.3 + stats::rnorm(length(xs), sd = 2)
  keep <- 404 - y <= height
  if (!any(keep)) next
  top <- pmin(1, (height[keep] - (404 - y)) / 14)
  track <- c(track, halftone_row(xs[keep], y, 0.2 + 0.6 * top))
}

segments <- c(
  capsule(64, 160, 424, 6),
  capsule(176, 236, 416, 6),
  capsule(252, 330, 424, 6),
  capsule(346, 456, 432, 6)
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
