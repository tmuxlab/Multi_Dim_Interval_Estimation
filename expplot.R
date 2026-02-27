# Example data
df <- data.frame(
        i        = 1:3,
        moduleID = c(0, 2, 3),
        LOCnorm  = c(0.745, 0.225, 0.320),
        m        = c(3.0, 4.5, 5.9)
)

# point size (optional)
cex_pts <- 1.1 + (df$m - min(df$m)) / (max(df$m) - min(df$m) + 1e-9) * 0.8

# ---- key fix: smaller margins ----
op <- par(mar = c(4.2, 4.0, 2.0, 0.8), xpd = NA)

plot(
        df$moduleID, df$LOCnorm,
        xlim = c(-0.1, 3.1), ylim = c(0, 1),
        xaxs = "i", yaxs = "i",
        xaxt = "n",
        xlab = "moduleID",
        ylab = "LOCnorm",
        pch  = 16,
        cex  = cex_pts
)
axis(1, at = 0:3, labels = 0:3)

text(df$moduleID + 0.07, df$LOCnorm + 0.04, labels = paste0("i=", df$i), cex = 0.9)

mtext("", side = 3, line = 0.2, cex = 0.9)

par(op)
