# Load function (התאימי את הנתיב אם שמרת בקובץ אחר)
source("R/normalize.R")

# --- צעצוע נתונים: counts 5 גנים x 6 דגימות ---
set.seed(1)
counts <- matrix(rpois(5*6, lambda = 50), nrow = 5,
                 dimnames = list(paste0("gene", 1:5),
                                 paste0("S", 1:6)))

# מטא־דאטה לדוגמא (בשביל VST חייבים rownames = sample IDs)
meta <- data.frame(group = rep(c("A","B"), each = 3),
                   row.names = colnames(counts))

cat("Counts dim: ", dim(counts), "\n")

# --- 1) TMM + logCPM ---
x_tmm <- normalize_counts(counts, method = "TMMlogCPM")
cat("\n[TMMlogCPM] dim:", dim(x_tmm), " attr(method) =", attr(x_tmm, "method"), "\n")
print(round(summary(as.numeric(x_tmm)), 3))

# --- 2) VST (DESeq2) ---
x_vst <- normalize_counts(counts, meta = meta, method = "VST")
cat("\n[VST] dim:", dim(x_vst), " attr(method) =", attr(x_vst, "method"), "\n")
print(round(summary(as.numeric(x_vst)), 3))

# הצלחה בסיסית: שני המטריצות קיימות ובגודל תואם
stopifnot(nrow(x_tmm) >= 1, ncol(x_tmm) == ncol(counts))
stopifnot(nrow(x_vst)  >= 1, ncol(x_vst)  == ncol(counts))

cat("\n✅ Smoke test passed.\n")
