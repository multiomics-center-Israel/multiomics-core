# write_pheatmap_png() draws a pheatmap onto its own PNG device, so the report
# no longer depends on which graphics device happens to be current.
# Synthetic data.

test_that("write_pheatmap_png writes a PNG and leaves the current device alone", {
    skip_if_not_installed("pheatmap")
    m <- matrix(c(1, 2, 3, 4, 2, 1, 4, 3, 5, 1, 2, 2), nrow = 3,
                dimnames = list(paste0("g", 1:3), paste0("S", 1:4)))
    # Two devices open, the later one current. Closing the PNG with dev.off()
    # alone hands control to the next device in the list, wrapping round to
    # the first, which is the bug this guards against.
    grDevices::pdf(NULL)
    first <- grDevices::dev.cur()
    grDevices::pdf(NULL)
    second <- grDevices::dev.cur()
    on.exit({ grDevices::dev.off(second); grDevices::dev.off(first) }, add = TRUE)
    ph <- pheatmap::pheatmap(m, silent = TRUE)
    # pheatmap() itself can switch devices while it builds, so take the
    # baseline after it: this test is about write_pheatmap_png() alone.
    grDevices::dev.set(second)
    before <- grDevices::dev.cur()

    f <- write_pheatmap_png(ph, width = 4, height = 3)

    expect_true(file.exists(f))
    expect_identical(readBin(f, "raw", 4), as.raw(c(0x89, 0x50, 0x4e, 0x47)))
    expect_identical(grDevices::dev.cur(), before)
})
