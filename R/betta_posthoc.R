#' Helper function to generate compact letter display
#'
#' @param pvals Named vector of p-values for pairwise comparisons
#' @param alpha Significance threshold
#' @param letter_set Character vector of letters to use
#' @return Named vector of compact letter display codes
#' @keywords internal
generate_cld <- function(pvals, alpha = 0.05, letter_set = LETTERS) {
  # Use multcompLetters with p-values
  cld_letters <- multcompView::multcompLetters(
    pvals,
    threshold = alpha,
    Letters = letter_set
  )$Letters

  return(cld_letters)
}

#' Post-hoc pairwise comparisons and compact letter display from diversity test object
#'
#' Performs pairwise post-hoc linear combination tests based on a fitted betta model object
#' (`test_obj`) that contains a formula of the form `estimate ~ ...`. The function
#' extracts the model formula and metadata variables, constructs all pairwise
#' comparisons between levels, runs `betta_lincom()` on each linear contrast,
#' adjusts p-values, computes significance codes, and returns compact letter display (CLD)
#' indicating significantly different groups.
#'
#' @param test_obj An object returned from the initial diversity test function (`betta()`).
#'   This object must contain a formula accessible as `test_obj$function.args$formula`,
#'   where the RHS defines the model terms used in `model.matrix`.
#' @param metadata A data frame containing the variables used in the model formula.
#'   Each variable should correspond to a column used in `test_obj$function.args$formula`.
#' @param p_adjust_method Character string specifying the p-value adjustment method
#'   (default `"BH"`). See `?p.adjust` for available methods.
#' @param alpha Numeric significance level used for compact letter display grouping (default `0.05`).
#' @param letter_set Character vector of letters or numbers used for compact letter display labels
#'   (default `LETTERS`).
#' @param human_sep Character string used as separator for output group labels (default `"."`).
#'
#' @return A list with two components:
#' \describe{
#'   \item{results_df}{A data frame containing pairwise comparison results with columns:
#'     `group_1`, `group_2`, `Estimates`, `Standard Errors`, `Lower CIs`, `Upper CIs`,
#'     `p_value`, `p_adjusted`, and `significance`.}
#'   \item{cld_df}{A data frame containing group-level compact letter display labels.
#'     This data frame contains one row per group with columns corresponding to
#'     the original variables from the formula, a `cld` column for the letter code,
#'     and group-level estimates: `estimate`, `std_error`, `lower_ci`, `upper_ci`.}
#' }
#'
#' @importFrom dplyr select distinct arrange mutate across everything all_of left_join bind_rows
#' @importFrom purrr map_dfr
#' @importFrom tidyr separate
#' @importFrom stringr str_replace_all
#' @importFrom multcompView multcompLetters
#' @importFrom tibble tibble
#'
#' @export
betta_posthoc <- function(
  test_obj,
  metadata,
  p_adjust_method = "BH",
  alpha = 0.05,
  letter_set = LETTERS,
  human_sep = "."
) {
  
  # Input validation
  # Check for required components
  if (
    !is.list(test_obj) ||
      is.null(test_obj$function.args) ||
      is.null(test_obj$function.args$formula)
  ) {
    stop(
      "test_obj must be a fitted betta model with accessible formula in test_obj$function.args$formula"
    )
  }

  # Additional check to ensure it has betta-like structure for betta_lincom
  if (is.null(test_obj$table)) {
    stop(
      "test_obj must be a fitted betta model with a 'table' component (required for betta_lincom)"
    )
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }

  # Internal separator for combining group names
  internal_sep <- "__SEP__"

  # Extract the right-hand side of the formula
  full_formula <- test_obj$function.args$formula
  rhs_formula <- as.formula(paste("~", as.character(full_formula)[3]))
  rhs_vars <- all.vars(rhs_formula)

  # Check that all variables in formula are in metadata
  missing_vars <- setdiff(rhs_vars, names(metadata))
  if (length(missing_vars) > 0) {
    stop(paste(
      "Variables not found in metadata:",
      paste(missing_vars, collapse = ", ")
    ))
  }

  # Create all unique combinations of factor levels
  combos <- metadata %>%
    dplyr::select(dplyr::all_of(rhs_vars)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(dplyr::across(dplyr::everything()))

  # Create group labels
  combos$group <- apply(combos, 1, paste, collapse = internal_sep)

  # Create model matrix for linear combinations
  mm <- model.matrix(rhs_formula, data = combos)

  # Generate all pairwise comparisons
  pairwise_indices <- utils::combn(nrow(combos), 2)

  # Perform pairwise comparisons using betta_lincom
  results_df <- purrr::map_dfr(seq_len(ncol(pairwise_indices)), function(i) {
    idx1 <- pairwise_indices[1, i]
    idx2 <- pairwise_indices[2, i]

    g1 <- combos$group[idx1]
    g2 <- combos$group[idx2]

    # Linear combination: difference between groups
    lin_com <- mm[idx1, ] - mm[idx2, ]

    # Perform linear combination test
    res <- betta_lincom(test_obj, lin_com)

    # Extract results
    res %>%
      dplyr::mutate(p_value = as.numeric(`p-values`)) %>%
      dplyr::select(
        Estimates,
        `Standard Errors`,
        `Lower CIs`,
        `Upper CIs`,
        p_value
      ) %>%
      dplyr::mutate(group_1 = g1, group_2 = g2, .before = Estimates)
  })

  # Add adjusted p-values and significance codes
  results_df <- results_df %>%
    dplyr::mutate(
      p_adjusted = p.adjust(p_value, method = p_adjust_method)#,
      #significance = dplyr::case_when(
      #  p_adjusted <= 0.001 ~ "***",
      #  p_adjusted <= 0.01 ~ "**",
      #  p_adjusted <= 0.05 ~ "*",
      #  p_adjusted <= 0.1 ~ ".",
      #  TRUE ~ " "
      #)
    )

  # Create p-value vector for compact letter display
  # multcompLetters needs a named vector of p-values for all pairwise comparisons
  pvals <- results_df$p_adjusted
  names(pvals) <- paste(results_df$group_1, results_df$group_2, sep = "-")

  # Generate compact letter display using p-value vector approach
  cld_letters <- generate_cld(pvals, alpha, letter_set)

  cld_df <- tibble::tibble(
    group_label = names(cld_letters),
    cld = as.character(cld_letters)
  )

  # Add group-level estimates and confidence intervals to CLD table
  # Calculate predicted values for each group using the model matrix
  group_estimates <- tibble::tibble()

  for (i in seq_len(nrow(combos))) {
    group_name <- combos$group[i]

    # Get linear combination for this group (intercept + effects)
    lin_com_group <- mm[i, ]

    # Calculate group estimate using betta_lincom
    group_result <- betta_lincom(test_obj, lin_com_group)

    group_estimates <- dplyr::bind_rows(
      group_estimates,
      tibble::tibble(
        group_label = group_name,
        estimate = as.numeric(group_result$Estimates),
        std_error = as.numeric(group_result$`Standard Errors`),
        lower_ci = as.numeric(group_result$`Lower CIs`),
        upper_ci = as.numeric(group_result$`Upper CIs`)
      )
    )
  }

  # Join estimates with CLD letters
  cld_df <- cld_df %>%
    dplyr::left_join(group_estimates, by = "group_label")

  # Safely escape internal separator for tidyr::separate()
  escaped_sep <- gsub(
    "([\\^\\$\\.|\\?\\*\\+\\(\\)\\[\\]\\{\\}\\\\])",
    "\\\\\\1",
    internal_sep
  )

  cld_df <- cld_df %>%
    tidyr::separate(
      group_label,
      into = rhs_vars,
      sep = escaped_sep,
      remove = TRUE
    )

  # Replace internal separator with readable one
  results_df <- results_df %>%
    dplyr::mutate(
      group_1 = stringr::str_replace_all(group_1, internal_sep, human_sep),
      group_2 = stringr::str_replace_all(group_2, internal_sep, human_sep)
    )

  cld_df <- cld_df %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(rhs_vars),
        ~ stringr::str_replace_all(., internal_sep, human_sep)
      )
    ) %>%
    # Also clean up the group_label if it still exists
    dplyr::mutate(
      dplyr::across(
        dplyr::any_of("group_label"),
        ~ stringr::str_replace_all(., internal_sep, human_sep)
      )
    )

  # Return results as a proper betta_posthoc object
  new_betta_posthoc(results_df, cld_df)
}


#' Create a betta_posthoc object with proper class
#'
#' @param results_df Pairwise comparison results
#' @param cld_df Compact letter display results
#' @return A betta_posthoc object
#' @keywords internal
new_betta_posthoc <- function(results_df, cld_df) {
  structure(
    list(
      results_df = results_df,
      cld_df = cld_df
    ),
    class = "betta_posthoc"
  )
}
