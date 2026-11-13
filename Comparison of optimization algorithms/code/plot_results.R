library(ggplot2)
library(patchwork)
library(ragg)

aa <- commandArgs(FALSE)
script <- sub('^--file=','',aa[grepl('^--file=',aa)][1])
root <- dirname(dirname(normalizePath(script,winslash='/')))
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args))args[1] else file.path(root,'figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(file.path(root, 'results/all_replicates_results.csv'), check.names = FALSE)
s <- read.csv(file.path(root, 'results/algorithm_summary.csv'), check.names = FALSE)
algorithms <- c('PSO', 'Tabu', 'OSNB', 'SA', 'GA', 'DE')
stopifnot(nrow(d) == 6000L, setequal(unique(d$Algorithm), algorithms),
          !anyDuplicated(d[c('Replicate', 'Algorithm')]),
          all(table(d$Algorithm) == 1000L),
          all(is.finite(d$BestObjective)), all(is.finite(d$RuntimeSec)),
          all(d$RuntimeSec > 0), length(unique(d$Replicate)) == 1000L)
stopifnot(all(table(d$Replicate) == 6L))
audit <- do.call(rbind, lapply(algorithms, function(a) {
  z <- d[d$Algorithm == a, ]; ref <- s[s$Algorithm == a, ]
  stopifnot(nrow(ref) == 1L, ref$N == nrow(z))
  values <- c(mean(z$BestObjective), sd(z$BestObjective),
              mean(z$RuntimeSec), median(z$RuntimeSec))
  reference <- as.numeric(ref[1, c('MeanBestObjective', 'SD', 'MeanRuntimeSec', 'MedianRuntimeSec')])
  stopifnot(max(abs(values - reference)) < 1e-9)
  data.frame(Algorithm = a, N = nrow(z), MeanBestObjective = values[1],
             SD = values[2], MeanRuntimeSec = values[3], MedianRuntimeSec = values[4],
             MaxSummaryDifference = max(abs(values - reference)))
}))
write.csv(audit, file.path(out, 'figure_data_check.csv'), row.names = FALSE)
d$Algorithm <- factor(d$Algorithm, levels = algorithms)
colors <- c(PSO = '#7222A8', Tabu = '#4264C4', OSNB = '#00B89C',
            SA = '#23C64B', GA = '#82D31A', DE = '#D3DF00')
labels <- c(PSO = 'Particle swarm\noptimization', Tabu = 'Tabu search',
            OSNB = 'Proposed method\n(OSNB)', SA = 'Simulated\nannealing',
            GA = 'Genetic\nalgorithm', DE = 'Differential\nevolution')
base <- theme_classic(base_size = 9, base_family = 'Arial') +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = .55),
        axis.line = element_blank(), axis.ticks = element_line(linewidth = .4),
        axis.ticks.length = grid::unit(1.3, 'mm'),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = 'bold', size = 8.5),
        axis.title.y = element_text(face = 'bold', size = 10, margin = margin(r = 5)),
        plot.title = element_text(face = 'bold', size = 11, hjust = .5, margin = margin(b = 7)),
        plot.margin = margin(7, 9, 4, 6), legend.position = 'none')
make_panel <- function(y, title, ylab) {
  ggplot(d, aes(x = Algorithm, y = .data[[y]], fill = Algorithm)) +
    geom_boxplot(width = .70, linewidth = .45, outlier.size = .50,
                 outlier.alpha = .65, outlier.stroke = 0, coef = 1.5) +
    scale_fill_manual(values = colors, drop = FALSE) +
    scale_x_discrete(labels = labels, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(.04, .06))) +
    labs(x = NULL, y = ylab, title = title) + base
}
p1 <- make_panel('BestObjective', 'A. Total expected progeny value', 'Total expected progeny value')
stopifnot(min(d$RuntimeSec) > .01, max(d$RuntimeSec) < 20)
# Transform coordinates after calculating box statistics on the raw seconds.
# Tick labels stay in seconds rather than displaying logarithm values.
p2 <- make_panel('RuntimeSec', 'B. Runtime', 'Runtime (s)') +
  scale_y_continuous(breaks = c(.01, .02, .05, .1, .2, .5, 1, 2, 5, 10),
                     labels = c('0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1', '2', '5', '10'),
                     expand = expansion(mult = 0)) +
  coord_transform(y = 'log10', ylim = c(.01, 20)) +
  theme(panel.grid.major.y = element_line(color = 'grey93', linewidth = .22, linetype = 'dotted'))
fig <- p1 + p2 + plot_layout(ncol = 2, widths = c(1, 1))
ggsave(file.path(out, 'Figure3.png'), fig, width = 8, height = 3.65,
       units = 'in', dpi = 220, device = ragg::agg_png, bg = 'white')
ggsave(file.path(out, 'Figure3.pdf'), fig, width = 8, height = 3.65,
       units = 'in', device = grDevices::cairo_pdf, bg = 'white')
ggsave(file.path(out, 'Figure3_1000dpi.tiff'), fig, width = 8, height = 3.65,
       units = 'in', dpi = 1000, device = ragg::agg_tiff, compression = 'lzw', bg = 'white')
writeLines(c(
  'Figure 3. Performance of six upper-level optimization algorithms across 1,000 independent simulated populations.',
  '(A) Total expected progeny value returned by each algorithm, calculated as the sum of expected progeny genetic values across the 12 mating pairs, not the per-progeny mean.',
  '(B) Method runtime (RuntimeSec) on a base-10 logarithmic axis. Tick labels show the original times in seconds: 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, and 10. Common input preparation is excluded.',
  'Runtime box statistics are computed on the untransformed observations before the display coordinates are transformed. All observations are visible; there is no axis break or inset.',
  'Each box contains 1,000 results. The central line is the median, the box spans the 25th to 75th percentiles,',
  'and whiskers extend to the most extreme observations within 1.5 interquartile ranges. Observations beyond the whiskers are plotted individually.',
  'All six algorithms were evaluated on the same set of populations. No observations were excluded or winsorized.',
  'OSNB is the proposed method. The reference figure included hill climbing, which was not part of this experiment and is omitted here.',
  'Color and panel structure follow the supplied Figure2.tiff; all plotted values come from all_replicates_results.csv.',
  'The supplied algorithm_summary.csv was independently reconciled against all 6,000 replicate-level records.'
), file.path(out, 'Figure3_caption.txt'))
capture.output(sessionInfo(), file = file.path(out, 'R_sessionInfo.txt'))
print(audit)
