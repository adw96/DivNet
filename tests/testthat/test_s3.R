library(DivNet)
context("Test output")

data(Lee)
divnet_phylum_char <- phyloseq::tax_glom(Lee, taxrank="Phylum") %>%
  divnet(X = "char", tuning = "test")

test_that("hypothesis testing works", {
  expect_is(testDiversity(divnet_phylum_char), "matrix")
})

sample_data(Lee)$sample_id <- 1:nrow(sample_data(Lee))

Lee_subset <- phyloseq::tax_glom(Lee, taxrank="Phylum") %>%
  phyloseq::subset_samples(char %in% c("glassy", "altered"))


divnet_phylum_sample <-   divnet(Lee_subset,
                                 X = "sample_id", tuning = "test")

ss_mat = diag(nrow(sample_data(Lee_subset)))
colnames(ss_mat) <- sample_data(Lee_subset)$sample_id


test_that("beta diversity hypothesis testing works", {
  expect_is(testBetaDiversity(dv = divnet_phylum_sample, h0 = "aitchison",
                              groups = sample_data(Lee_subset)$char,
                              sample_specimen_matrix = ss_mat,
                              n_boot = 10), "list")
})

test_that("beta diversity hypothesis testing works for bray-curtis", {
  expect_is(testBetaDiversity(dv = divnet_phylum_sample, h0 = "bray-curtis",
                              groups = sample_data(Lee_subset)$char,
                              sample_specimen_matrix = ss_mat,
                              n_boot = 10), "list")
})

test_that("beta diversity hypothesis testing works with two even groups", {
  Lee_even <- Lee_subset %>%
    phyloseq::subset_samples(sample_id %in% c(3:4, 6:9, 15:16))
  divnet_phylum_even <-   divnet(Lee_even,
                                   X = "sample_id", tuning = "test")
  ss_mat_even = diag(nrow(sample_data(Lee_even)))
  colnames(ss_mat_even) <- sample_data(Lee_even)$sample_id
  expect_is(testBetaDiversity(dv = divnet_phylum_even, h0 = "bray-curtis",
                              groups = sample_data(Lee_even)$char,
                              sample_specimen_matrix = ss_mat_even,
                              n_boot = 10), "list")
  expect_is(testBetaDiversity(dv = divnet_phylum_even, h0 = "euclidean",
                              groups = sample_data(Lee_even)$char,
                              sample_specimen_matrix = ss_mat_even,
                              n_boot = 10), "list")
  expect_is(testBetaDiversity(dv = divnet_phylum_even, h0 = "aitchison",
                              groups = sample_data(Lee_even)$char,
                              sample_specimen_matrix = ss_mat_even,
                              n_boot = 10), "list")
})

test_that("beta diversity hypothesis testing works for euclidean", {
  expect_is(testBetaDiversity(dv = divnet_phylum_sample, h0 = "euclidean",
                              groups = sample_data(Lee_subset)$char,
                              sample_specimen_matrix = ss_mat,
                              n_boot = 10), "list")
})

test_that("Bray-Curtis pseudo-F function works",
          {
            expect_equal(get_bc_test_statistic(bc_mat = divnet_phylum_sample$`bray-curtis`,
                                               groups = sample_data(Lee_subset)$char,
                                               unique_groups = c("altered","glassy"),
                                               n_groups = 2,
                                               n_specimens = nrow(sample_data(Lee_subset))
            ),
            25.88456,
            tolerance = 1e-2)
            
          })

test_that("Euclidean pseudo-F function works",
          {
            expect_equal(get_euc_test_statistic(euc_mat = divnet_phylum_sample$euclidean,
                                                groups = sample_data(Lee_subset)$char,
                                                unique_groups = c("altered","glassy"),
                                                n_groups = 2,
                                                n_specimens = nrow(sample_data(Lee_subset))
            ),
            26.55571,
            tolerance = 1e-2)
            
          })

test_that("Aitchison distance calculation works",
          {
            expect_is(get_aitchison_distance(divnet_phylum_sample$fitted_z),
                      "matrix")
          })


test_that("Aitchison pseudo-F function works",
          {
            expect_equal(get_euc_test_statistic(euc_mat = get_aitchison_distance(divnet_phylum_sample$fitted_z),
                                                groups = sample_data(Lee_subset)$char,
                                                unique_groups = c("altered","glassy"),
                                                n_groups = 2,
                                                n_specimens = nrow(sample_data(Lee_subset))
            ),
            26.09607,
            tolerance = 1e-2)
            
          })

fake_composition_mat <- rbind(rep(1/10,10),
                              exp(-4:5)/sum(exp(-4:5)))

test_that("CLR function works",
          {
            expect_equal(log_ratio(fake_composition_mat),
                         rbind(0,
                               seq(-4.5,4.5,by = 1)))
          })


test_that("plotting works", {
  expect_is(plot(divnet_phylum_char), "ggplot")
})

test_that("printing works", {
  expect_null(print(divnet_phylum_char))
})