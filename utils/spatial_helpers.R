## Helpers of Fig. 1D and Fig. S1F.

suppressPackageStartupMessages({ library(grid) })

## -2 log10(max(p_a, p_b)), p = upper normal tail of the z-score; 2.6 = both in the top 5%.
coexpression_score <- function(a, b) {
   p_upper <- function(x) 1 - pnorm(as.numeric(scale(x)))
   -log10(pmax(p_upper(a), p_upper(b))^2)
}

## 0 or 180 degrees only.
rotate_section <- function(object, slide, degrees) {
   if (degrees == 0) return(object)
   if (degrees != 180) stop("only 0 and 180 degrees are supported")
   im <- object@images[[slide]]
   for (b in names(im@boundaries)) {
      bd <- im@boundaries[[b]]
      co <- bd@coords
      ctr <- (apply(co, 2, min) + apply(co, 2, max)) / 2
      new <- sweep(-sweep(co, 2, ctr), 2, ctr, FUN = "+")
      dimnames(new) <- dimnames(co)
      bd@coords <- new
      if (.hasSlot(bd, "bbox")) bd@bbox[] <- rbind(range(new[, 1]), range(new[, 2]))
      im@boundaries[[b]] <- bd
   }
   object@images[[slide]] <- im
   object
}

## tissue_positions_list columns 5-6 = pixel row, col.
he_crop <- function(raw_folder, barcodes, degrees = 0, max_px = 1600) {
   tif <- list.files(raw_folder, "\\.tif(\\.gz)?$", full.names = TRUE)[1]
   pos <- list.files(raw_folder, "tissue_positions_list\\.csv(\\.gz)?$", full.names = TRUE)[1]
   if (grepl("\\.gz$", tif)) {
      tmp <- tempfile(fileext = ".tif")
      R.utils::gunzip(tif, destname = tmp, remove = FALSE)
      tif <- tmp
   }
   xy <- read.csv(pos, header = FALSE)
   xy <- xy[xy[[1]] %in% barcodes, ]
   pitch <- 0.5 * median(diff(sort(unique(xy[[6]]))))
   x0 <- max(0, min(xy[[6]]) - pitch); y0 <- max(0, min(xy[[5]]) - pitch)
   w <- max(xy[[6]]) + pitch - x0;     h <- max(xy[[5]]) + pitch - y0
   img <- magick::image_read(tif)
   img <- magick::image_crop(img, magick::geometry_area(w, h, x0, y0))
   if (degrees != 0) img <- magick::image_rotate(img, degrees)
   magick::image_scale(img, paste0(max_px, "x", max_px, ">"))
}

column_header_grob <- function(labels, group_label = "Post-MI", n_ungrouped = 1L, fs = 1) {
   n <- length(labels)
   xc <- (seq_len(n) - 0.5) / n
   fsz <- max(5.5, 8 * fs)
   kids <- gList()
   for (i in seq_len(n))
      kids <- gList(kids, textGrob(labels[i], x = xc[i], y = 0.16, vjust = 0,
                                   gp = gpar(fontsize = fsz, col = "grey15")))
   j <- (n_ungrouped + 1L):n
   x0 <- xc[j[1]] - 0.5 / n + 0.012; x1 <- xc[j[length(j)]] + 0.5 / n - 0.012
   gTree(children = gList(kids,
      segmentsGrob(x0 = x0, x1 = x1, y0 = 0.56, y1 = 0.56, gp = gpar(col = "grey25", lwd = 0.7)),
      textGrob(group_label, x = (x0 + x1) / 2, y = 0.66, vjust = 0, gp = gpar(fontsize = fsz, col = "grey15"))))
}

colourbar_grob <- function(cols, title, bar_width_in = 2.2, title_side = c("top", "bottom"), note = NULL) {
   title_side <- match.arg(title_side)
   bw <- unit(bar_width_in, "inches"); bh <- unit(0.115, "inches"); y0 <- unit(0.56, "npc")
   edge <- function(p) unit(0.5, "npc") + (p - 0.5) * bw
   title_y <- if (title_side == "top") y0 + bh * 0.5 + unit(0.085, "inches")
              else y0 - bh * 0.5 - unit(0.115, "inches")
   kids <- gList(
      rasterGrob(matrix(cols, nrow = 1), x = 0.5, y = y0, width = bw, height = bh, interpolate = TRUE),
      textGrob(title, x = 0.5, y = title_y, vjust = if (title_side == "top") 0.5 else 1,
               gp = gpar(fontsize = 8, col = "grey15")),
      textGrob("Min", x = edge(0) - unit(0.06, "inches"), y = y0, hjust = 1, gp = gpar(fontsize = 8, col = "grey15")),
      textGrob("Max", x = edge(1) + unit(0.06, "inches"), y = y0, hjust = 0, gp = gpar(fontsize = 8, col = "grey15")))
   if (!is.null(note))
      kids <- gList(kids, textGrob(note, x = 0.5, y = y0 - bh * 0.5 - unit(0.115, "inches"), vjust = 1,
                                   gp = gpar(fontsize = 6.5, col = "grey40")))
   gTree(children = kids)
}

## `genes`: data.frame with sym, side ("upstream"/"downstream"), note.
S1C_COL <- c(linc = "#E26235", upstream = "#BBE1D6", downstream = "#E2DBBF")
s1c_label_grob <- function(linc_label, linc_note, genes, greyed = FALSE, letter = NULL, fs = 1) {
   n  <- nrow(genes)
   yr <- 1 - (seq_len(n) - 0.5) / n
   fsz <- function(size, floor = 5) max(floor, size * fs)
   al <- if (greyed) 0.40 else 1
   kids <- gList(
      roundrectGrob(x = 0.160, y = 0.5, width = 0.130, height = 0.86, r = unit(1.4, "mm"),
                    gp = gpar(fill = S1C_COL[["linc"]], col = NA, alpha = al)),
      textGrob(linc_label, x = 0.139, y = 0.5, rot = 90,
               gp = gpar(fontsize = fsz(8.5, 7), fontface = "bold", col = if (greyed) "grey25" else "white")),
      textGrob(linc_note, x = 0.200, y = 0.5, rot = 90,
               gp = gpar(fontsize = fsz(5.8, 5.2), col = if (greyed) "grey30" else "white")))
   if (!is.null(letter))
      kids <- gList(kids, textGrob(letter, x = unit(0.012, "npc"), y = unit(1, "npc") - unit(0.10 * fs, "inches"),
                                   hjust = 0, vjust = 1, gp = gpar(fontsize = 12 * fs, fontface = "bold")))
   hbox <- min(0.34, 0.62 / n)
   for (i in seq_len(n)) {
      up <- identical(genes$side[i], "upstream")
      col <- if (greyed) "grey60" else "grey35"
      flat <- abs(yr[i] - 0.5) < 1e-9
      xs <- if (flat) c(0.227, 0.398) else c(0.227, 0.305, 0.305, 0.398)
      ys <- if (flat) c(0.5, 0.5) else c(0.5, 0.5, yr[i], yr[i])
      kids <- gList(kids,
         linesGrob(x = xs, y = ys, arrow = arrow(length = unit(1.4 * fs, "mm"), type = "closed", angle = 20),
                   gp = gpar(col = col, fill = col, lwd = 0.9, linejoin = "round")),
         textGrob(if (up) "5'" else "3'", x = 0.35, y = unit(yr[i], "npc") + unit(0.4 * fs, "mm"), vjust = 0,
                  gp = gpar(fontsize = max(4.8, 5.2 * fs), fontface = "italic", col = col)),
         roundrectGrob(x = 0.700, y = yr[i], width = 0.560, height = hbox, r = unit(1.2, "mm"),
                       gp = gpar(fill = S1C_COL[[genes$side[i]]], col = NA, alpha = al)),
         textGrob(genes$sym[i], x = 0.700, y = yr[i] + 0.20 * hbox,
                  gp = gpar(fontsize = fsz(8, 6.2), fontface = "bold", col = if (greyed) "grey55" else "grey20")),
         textGrob(genes$note[i], x = 0.700, y = yr[i] - 0.22 * hbox,
                  gp = gpar(fontsize = fsz(6, 5.2), col = if (greyed) "grey65" else "grey45")))
   }
   gTree(children = kids)
}
