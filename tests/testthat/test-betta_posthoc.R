data("Lee")
lee_phy <- phyloseq::tax_glom(Lee, taxrank = "Phylum")

test_that("betta_posthoc function works as expected", {
  div_obj <- divnet(lee_phy)
  bett_obj <- betta(data = data.frame(sample_data(lee_phy),
                                      chats = sapply(div_obj$shannon, function(x) x$estimate),
                                      ses = sapply(div_obj$shannon, function(x) x$error)),
                    formula = chats ~ char,
                    ses = "ses")
  bett_post <- betta_posthoc(test_obj = bett_obj,
                             metadata = data.frame(sample_data(lee_phy)))
  expect_true(all.equal(bett_obj$table[-1, 1] %>% unname(), 
                        -bett_post$results_df$Estimates[1:4] %>% unname()))
  # we expect estimates to align when comparing baseline level to other levels, not 
  # p-values because of differences in how these are calculated in betta and betta_lincom
})
