# source(here::here("scripts/config.R"))

########################
# VARIABLE definitions #
########################
source(here::here("scripts/variable_names.R"))

######################
# FUNCTION defitions #
######################
get_counts <- function(DF, COLNAME) {
  checkmate::assert_data_frame(DF)
  checkmate::assert_string(COLNAME)
  DF <- as.data.frame(table(DF[[COLNAME]])) %>%
    rename(antibody = Var1)

  return(DF)
}



get_med_stats <- function(DF) {
  ab_counts <- get_counts(DF %>% distinct(id_sample, antibody), "antibody")
  med_counts <- get_counts(DF, "drugdesc")

  combined_ab <- full_join(ab_counts, med_counts) %>%
    arrange(antibody)
  return(combined_ab)
}

get_ab_stats <- function(DF) {
  ab_counts <- get_ab_counts(DF)
  bleed_counts <- get_bleed_counts(DF)

  combined_ab <- full_join(ab_counts, bleed_counts)
  combined_ab <- combined_ab %>%
    arrange(antibody)
  return(combined_ab)
}

print_ab_table <- function(DF) {
  DF <- DF %>%
    # as_tabyl() %>%
    flextable::flextable() %>%
    flextable::autofit() %>%
    flextable::theme_booktabs() %>%
    flextable::bold(i = 1, bold = TRUE, part = "header") # %>%
  # vline(border = fp_border_default(), part = "all") %>%
  # hline(border = fp_border_default(), part = "all")
  return(DF)
}

print_ab_table_color <- function(DF) {
  color_ild <- "#56B4E9"
  color_questionable <- "#CC79A7"
  color_ild_related <- "#F0E442"
  color_non_ild <- "#999999"

  DF <- DF %>%
    flextable::flextable() %>%
    flextable::autofit() %>%
    flextable::theme_booktabs() %>%
    flextable::bold(i = 1, bold = TRUE, part = "header") %>%
    flextable::vline(border = flextable::fp_border_default(), part = "all") %>%
    flextable::hline(border = flextable::fp_border_default(), part = "all") %>%
    flextable::bg(~ antibody %in% ASYS_MDA5, bg = color_ild) %>%
    flextable::bg(~ antibody %in% MAA_ab, bg = color_ild_related) %>%
    flextable::bg(~ antibody == "UNDEFINED", bg = color_questionable) %>%
    flextable::bg(~ antibody == "NEGATIVE", bg = color_questionable) %>%
    flextable::bg(~ antibody %in% MYOSITIS_NO_ILD, bg = color_non_ild)
  return(DF)
}

##################################
# background med filter function #
##################################

assert_prior_row_same <- function(DF, REF, COMP_1, i, j) {
  if (nrow(DF) == 0) {
    return(TRUE)
  }
  if (REF[i] == DF$bleed_date[nrow(DF)] && (COMP_1[i] == DF$startdate[nrow(DF)] || COMP_1[j] == DF$startdate[nrow(DF)])) {
    # if(DF$bleed_date[nrow(DF) - 1] == DF$bleed_date[nrow(DF)] && (COMP_1[i] == DF$startdate[nrow(DF)] || COMP_1[j] == DF$startdate[nrow(DF)])){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

insert_bg_days <- function(DF, REF, COMP_1, COMP_2, i, j, INSERT_ROW, MED) {
  checkmate::assert_string(MED)

  on_med <- sprintf("on_%s", MED)
  days_on_med <- sprintf("days_on_%s", MED)
  days_off_med <- sprintf("days_off_%s", MED)

  DF[INSERT_ROW, on_med] <- TRUE

  # CASE: enddate EXISTS
  if (!is.na(COMP_2[j])) {
    # CASE: bleed >= enddate, so calc days_on and days_off
    if (REF[i] >= COMP_2[j]) {
      DF[INSERT_ROW, days_on_med] <- COMP_2[j] - COMP_1[j]
      DF[INSERT_ROW, days_off_med] <- REF[i] - COMP_2[j]
    }
    # CASE: enddate is later than bleed, so calc with startdate
    else if (REF[i] < COMP_2[j]) {
      # CASE: ensure startdate exists and use startdate to calc
      if (!is.na(COMP_1[j]) && REF[i] >= COMP_1[j]) {
        DF[INSERT_ROW, days_on_med] <- REF[i] - COMP_1[j]
      }
    }
  }
  # CASE: enddate does NOT exist
  if (is.na(COMP_2[j])) {
    # TODO: if bleeddate >= startdate
    if (REF[i] >= COMP_1[j]) {
      DF[INSERT_ROW, days_on_med] <- REF[i] - COMP_1[j]
    }
    # TODO: if bleeddate < startdate SKIP
    if (REF[i] >= COMP_1[j]) {}
  }

  return(DF)
}

add_med <- function(DF, DRUGCODE, REF, COMP_1, COMP_2, i, j, INSERT_ROW) {
  checkmate::assert_int(INSERT_ROW)
  DF[INSERT_ROW, "on_bg_med"] <- TRUE

  if (DRUGCODE[j] == as.integer(drugnames["azathioprine"])) {
    DF <- insert_bg_days(DF, REF, COMP_1, COMP_2, i, j, INSERT_ROW, "aza")
  } else if (DRUGCODE[j] == as.integer(drugnames["ivig"])) {
    DF <- insert_bg_days(DF, REF, COMP_1, COMP_2, i, j, INSERT_ROW, "ivig")
  } else if (DRUGCODE[j] == as.integer(drugnames["methotrexate"])) {
    DF <- insert_bg_days(DF, REF, COMP_1, COMP_2, i, j, INSERT_ROW, "mtx")
  } else if (DRUGCODE[j] == as.integer(drugnames["mycophenolate_mofetil"])) {
    DF <- insert_bg_days(DF, REF, COMP_1, COMP_2, i, j, INSERT_ROW, "mmf")
  } else if (DRUGCODE[j] == as.integer(drugnames["rituximab"])) {
    DF <- insert_bg_days(DF, REF, COMP_1, COMP_2, i, j, INSERT_ROW, "rtx")
  } else {
    DF <- insert_bg_days(DF, REF, COMP_1, COMP_2, i, j, INSERT_ROW, "other")
  }


  return(DF)
}

bool_bg_med <- function(DF) {
  df <- DF %>%
    mutate(
      on_bg_med = FALSE,
      on_mtx = "",
      on_mmf = "",
      on_ivig = "",
      on_aza = "",
      on_pred = "",
      on_rtx = "",
      on_other = "",
      days_on_mtx = "",
      days_on_mmf = "",
      days_on_ivig = "",
      days_on_aza = "",
      days_on_pred = "",
      days_on_rtx = "",
      days_on_other = "",
      days_off_mtx = "",
      days_off_mmf = "",
      days_off_ivig = "",
      days_off_aza = "",
      days_off_pred = "",
      days_off_rtx = "",
      days_off_other = "",
    )

  return(df)
}

insert_background_med <- function(DF, UNI_DF, DRUGCODE, REF, COMP_1, COMP_2, i, j) {
  checkmate::assert_character(DRUGCODE[i])
  if (nrow(DF) > 0) {
    if (REF[i] != REF[j] && DF[nrow(DF), "drugcode"] == DRUGCODE[i]) {
      # if DF[last_row] == UNI_DF[checking_row], prevent adding
      # a second copy of UNI_DF[checking_row] to DF
      # even so, UNI_DF[checking_row] will call add_med()
      # to overwrite the row bg_med status
      if (assert_prior_row_same(DF, REF, COMP_1, i, j) == TRUE) {
      }
      # case assert_prior_row_same(DF, REF, COMP_1, i, j) == FALSE
      else {
        DF <- rbindlist(list(DF, UNI_DF[i, ]))
      }
    }
    # else if (REF[i] == DF$bleed_date[nrow(DF)] && (COMP_1[i] == DF$startdate[nrow(DF)] || COMP_1[j] == DF$startdate[nrow(DF)])){
    else if (assert_prior_row_same(DF, REF, COMP_1, i, j) == TRUE) {
    } else {
      DF <- rbindlist(list(DF, UNI_DF[j, ]))
      DF <- rbindlist(list(DF, UNI_DF[i, ]))
    }
  } else {
    DF <- rbindlist(list(DF, UNI_DF[j, ]))
    DF <- rbindlist(list(DF, UNI_DF[i, ]))
  }

  n_temp <- nrow(DF)
  DF <- add_med(DF, DRUGCODE, REF, COMP_1, COMP_2, i, j, n_temp)

  return(DF)
}

not_within_time <- function(DF, UNI_DF, DRUGCODE, REF, COMP_1, COMP_2, i, j, DAYS) {

}

# if two med start dates are within DAY, add to DF
is_within_time <- function(DF, UNI_DF, DF_ID, DRUGCODE, REF, COMP_1, COMP_2, i, j, DAYS) {
  checkmate::assert_date(REF[i])
  checkmate::assert_date(COMP_1[j])
  checkmate::assert_date(COMP_2[j])
  checkmate::assert_int(DAYS)
  if (i > j) {
    # if bleeddate[i] occurs after startdate[j] and there is no enddate[j]
    if (DF_ID[i] == DF_ID[j]) {
      if (REF[i] >= COMP_1[j]) {
        if (is.na(COMP_2[j])) {
          # is bleeddate[i] - startdate[j] <= DAYS?
          if (REF[i] - COMP_1[j] <= DAYS) {
            DF <- insert_background_med(DF, UNI_DF, DRUGCODE, REF, COMP_1, COMP_2, i, j)
          }
        }
        # if bleeddate[i] occurs after startdate[j] and there IS an endate[j]
        else if (!is.na(COMP_2[j])) {
          # check if bleed - start <= DAYS
          if (REF[i] - COMP_1[j] <= DAYS) {
            DF <- insert_background_med(DF, UNI_DF, DRUGCODE, REF, COMP_1, COMP_2, i, j)
          }
          # check if bleed - end <= DAYS
          else if (REF[i] - COMP_2[j] <= DAYS) {
            DF <- insert_background_med(DF, UNI_DF, DRUGCODE, REF, COMP_1, COMP_2, i, j)
          }
          # if bleed - end is out of DAYS range, ensure row[j] is not duplicate
          # of row[i]; if not duplicate, add current row[i]
          else if (assert_prior_row_same(DF, REF, COMP_1, i, j) == FALSE) {
            DF <- rbindlist(list(DF, UNI_DF[i, ]))
          }
          # if none of above, do not add row[i]
          else {}
        }
      }
      # case when prior row's startdate is AFTER bleed date
      # thus, shouldn't call insert_background_med() but
      # should just add the row
      if (REF[i] < COMP_1[j]) {
        if (assert_prior_row_same(DF, REF, COMP_1, i, j) == FALSE) {
          DF <- rbindlist(list(DF, UNI_DF[i, ]))
        }
      }
    }
    # if DF_ID[i] and DF_ID[j] are not the same...
    if (DF_ID[i] != DF_ID[j]) {
      # do NOT add DF_ID[i] if it already exists in DF
      if (assert_prior_row_same(DF, REF, COMP_1, i, j) == TRUE) {}
      # if DF_ID[i] hasn't been added to DF yet, add it but
      # don't call insert_bg_med(), since IDs don't match
      else {
        DF <- rbindlist(list(DF, UNI_DF[i, ]))
      }
    }
  }
  # returns TRUE if DF[i] and DF[j] are same rows and nrow(DF) >=0
  if (i <= j && assert_prior_row_same(DF, REF, COMP_1, i, j) == TRUE && nrow(DF) > 0) {
  }
  # adds first row if nrow(DF) = 0
  else if (nrow(DF) == 0) {
    DF <- rbindlist(list(DF, UNI_DF[i, ]))
  }
  return(DF)
}

create_background <- function(UNI_DF, DF_ID, DRUGCODE, REF, COMP_1, COMP_2, DAYS) {
  checkmate::assert_character(DF_ID)
  checkmate::assert_date(REF)
  checkmate::assert_date(COMP_1)
  checkmate::assert_date(COMP_2)
  checkmate::assert_int(DAYS)

  n_uni <- nrow(UNI_DF)
  temp_df <- UNI_DF[FALSE, ]

  for (i in 1:n_uni) {
    for (j in 1:n_uni) {
      if (i == j) {
        if (i != 1) {
          break()
        }
        if (i == 1) {
          temp_df <- is_within_time(temp_df, UNI_DF, DF_ID, DRUGCODE, REF, COMP_1, COMP_2, i, j, DAYS)
          break()
        }
      }
      if (i > j) {
        # ensures past meds are for the same patient, else moves to next iter j
        if (DF_ID[i] == DF_ID[j]) {
          temp_df <- is_within_time(temp_df, UNI_DF, DF_ID, DRUGCODE, REF, COMP_1, COMP_2, i, j, DAYS)
        }
        if (DF_ID[i] != DF_ID[j]) {
          # below should should that background med state is FALSE
          # and will just end up rbind() the row in question
          # without setting on_bg_med to be TRUE
          temp_df <- is_within_time(temp_df, UNI_DF, DF_ID, DRUGCODE, REF, COMP_1, COMP_2, i, j, DAYS)
        }
      }

      # after above checks are performed, decide how to iter
      if (i > j) {
        next()
      }
      # in theory, it should never be that i < j; here for safety
      else if (i < j) {
        break()
      }
    }
    # counts which iter row i out of total rows (n_uni) are done
    # cat(sprintf("\rIter row: %s/%s", i, n_uni))
  }

  temp_df <- temp_df %>%
    arrange(id_patient_combined, startdate)
  # %>%
  # distinct()

  temp_df <- as.data.frame(temp_df)
  return(temp_df)
}


###################
# DEFINE FUNCTION #
###################
get_earliest <- function(DF) {
  temp_df <- DF[FALSE, ]
  uni_compiled <- as.character(unlist(distinct(DF, id_sample)))
  n_rows <- length(uni_compiled)

  for (i in 1:n_rows) {
    temp <- DF %>% filter(id_sample == uni_compiled[i])
    if (nrow(temp) > 1) {
      temp <- DF %>%
        filter(id_sample == uni_compiled[i]) %>%
        filter(bleed_date == min(bleed_date))
      temp_df <- rbindlist(list(temp_df, temp))
    } else {
      temp_df <- rbindlist(list(temp_df, temp))
    }
  }
  return(temp_df)
}

get_bleed_counts <- function(DF) {
  DF <- as.data.frame(table(distinct(DF, .data[["antibody"]], .data[["bleed_date"]])[["antibody"]]))

  DF <- DF %>%
    rename(
      antibody = Var1,
      number_of_bleeds = Freq
    )
  return(DF)
}

get_ab_counts <- function(DF) {
  RESULT <- get_counts(DF %>% distinct(id_sample, antibody), "antibody")
  return(RESULT)
}

get_med_counts <- function(DF, COLNAME) {
  return(get_counts(DF %>% distinct(id_sample, antibody), COLNAME))
}

ab_dx <- query_db(DB_NORM, "SELECT * FROM antibody_diagnosis")

filter_out_bleeds_on_hand <- function(DF) {
  ban <- import_file(here::here("reports/query_results"), "bleeds_at_nih", "csv") %>%
    mutate(across(everything(), as.character))

  DF <- DF %>%
    distinct(id_sample, antibody, bleed_date) %>%
    mutate(across(everything(), as.character))

  lf_not_in_ban <- anti_join(lf, ban)


  return(lf_not_in_ban)
}

ensure_fu_days_exists <- function(DF, LIST) {
  # LIST is a list of numeric days that MUST exist in the passed dataframe DF
  # ALL days in LIST must be present for it to be counted in the returned DF
  DF <- DF %>%
    filter(followup_day %in% LIST)

  temp_df <- DF[FALSE, ]

  n_list <- length(LIST)
  df_id <- as.character(unlist(DF[["id_sample"]]))
  uni_df_id <- as.character(unlist(distinct(DF, .data[["id_sample"]])))
  n_id <- length(uni_df_id)

  for (i in 1:n_id) {
    temp <- DF[DF[["id_sample"]] == uni_df_id[i], ]
    if (nrow(temp) < n_list) {
      next()
    } else {
      if (length(unique(temp[["startdate"]])) == 1) {
        day_counter <- 0

        for (list_iter in 1:n_list) {
          if (nrow(temp[temp[["followup_day"]] == LIST[list_iter], ]) == 1) {
            day_counter <- day_counter + 1
          } else {
            next()
          }
        }
        if (day_counter == n_list) {
          temp_df <- rbindlist(list(temp_df, temp))
        } else {
          next()
        }
      } else {
        next()
      }
    }
  }
  return(temp_df)
}

get_bleeds_at_nih <- function() {
  ab_dx <- query_db(DB_NORM, "SELECT DISTINCT * FROM antibody_diagnosis")
  emr_ids <- query_db(DB_NORM, "SELECT DISTINCT id_sample, id_patient_combined FROM emr_ids")
  ab_dx <- left_join(ab_dx, emr_ids)
  samples_at_nih <- import_file(here::here("reports/query_results"), "samples_at_nih", "csv") %>%
    mutate(across(everything(), as.character)) %>%
    rename(id_sample = rosen_sera_number)

  for_eric <- import_file(here::here("data/to_do"), "for_eric_july_28", "csv") %>%
    rename(id_sample = "Sample ID", antibody = "Suspected Antibody") %>%
    mutate(antibody = str_to_upper(antibody))
  for_eric <- filter_antibodies(for_eric, ASYS_MDA5) %>%
    rename(used_in_old_SOMA = `Samples used in SOMA`, 
           used_in_old_olink = `Samples used in Olink`
    ) %>%
    mutate(
      across(c(antibody, used_in_old_SOMA, used_in_old_olink), str_to_upper),
      across(everything(), as.character)
    ) %>%
    mutate(
      used_in_old_olink = fct_recode(used_in_old_olink, "TRUE" = "1"),
      used_in_old_SOMA = fct_recode(used_in_old_SOMA, "TRUE" = "MDA5"),
      used_in_old_SOMA = fct_recode(used_in_old_SOMA, "TRUE" = "JO1"),
      used_in_old_SOMA = fct_recode(used_in_old_SOMA, "TRUE" = "EJ"),
      used_in_old_SOMA = fct_recode(used_in_old_SOMA, "TRUE" = "OJ"),
      used_in_old_SOMA = fct_recode(used_in_old_SOMA, "TRUE" = "PL7"),
      used_in_old_SOMA = fct_recode(used_in_old_SOMA, "TRUE" = "PL12"),
    ) %>%
    filter(used_in_old_SOMA == "TRUE" | used_in_old_olink == "TRUE") %>%
    select(id_sample, antibody, used_in_old_SOMA, used_in_old_olink)
  
  # TODO: add bleed dates from Euroimmun to for_eric
  soma_dates <- add_ab_dx(
    add_emr_ids(
      query_db(DB_NORM, "SELECT DISTINCT id_patient_combined, bleed_date as bleed_date FROM antibodies_detailed WHERE source = 'Euroimmun' ")  
    )
  ) %>%
    select(-c(id_patient_combined, id_emr))


  # TODO: join for_eric and samples_at_nih, which accounts for all blood samples we have on hand and thus don't need to send
  for_eric_dates <- left_join(for_eric, soma_dates)

  nih_plus_soma <- full_join(for_eric_dates, samples_at_nih) %>%
    mutate(bleed_date = as.Date(bleed_date))

  return(nih_plus_soma)
}

print_bleeds_at_nih <- function(DF_CANDIDATES) {
  DF_AT_NIH <- query_db(DB_NORM, "SELECT * FROM bleeds_at_nih") %>%
    mutate(across(everything(), as.character))

  for_eric <- import_file(here::here("data/to_do"), "for_eric_july_28", "csv") %>%
    rename(id_sample = "Sample ID", antibody = "Suspected Antibody") %>%
    mutate(antibody = str_to_upper(antibody))
  for_eric <- filter_antibodies(for_eric, ASYS_MDA5)
  for_eric <- for_eric[for_eric[["Samples used in SOMA"]] != "", ] %>%
    filter("Samples used in SOMA" != "") %>%
    select(id_sample, antibody)

  # TODO: add bleed dates from Euroimmun to for_eric
  soma_dates <- query_db(DB_NORM, "SELECT DISTINCT id_patient_combined, bleed_date as bleed_date FROM antibodies_detailed WHERE source = 'Euroimmun' ")
  soma_dates <- left_join(soma_dates, query_db(DB_NORM, "SELECT id_patient_combined, id_sample FROM emr_ids"), relationship = "many-to-many")

  soma_dates <- inner_join(soma_dates, ab_dx, relationship = "many-to-many") %>%
    distinct(id_sample, antibody, bleed_date)

  # TODO: join for_eric and samples_at_nih, which accounts for all blood samples we have on hand and thus don't need to send
  for_eric_dates <- left_join(for_eric, soma_dates) %>% filter(!is.na(bleed_date))

  DF_CANDIDATES <- DF_CANDIDATES %>%
    dplyr::mutate(across(everything(), as.character)) %>%
    distinct(.data[["id_sample"]], .data[["bleed_date"]])

  DF_AT_NIH <- DF_AT_NIH %>%
    dplyr::mutate(across(everything(), as.character))

  n_bleeds_on_hand <- nrow(inner_join(DF_CANDIDATES, DF_AT_NIH))
  n_candidates <- nrow(DF_CANDIDATES)
  n_soma <- nrow(inner_join(DF_CANDIDATES, for_eric_dates))
  n_not_soma <- nrow(inner_join(DF_CANDIDATES, DF_AT_NIH)) - nrow(inner_join(DF_CANDIDATES, for_eric_dates))

  cat(sprintf("We already have n = %s / %s bleeds on hand.\n", n_bleeds_on_hand, n_candidates))

  cat(sprintf(
    "Of that:\n\t- n = %s / %s are from SOMA testing\n\t- n = %s / %s are from other non-SOMA testing at the NIH\n",
    n_soma,
    n_bleeds_on_hand,
    n_not_soma,
    n_bleeds_on_hand
  ))

  cat(sprintf(
    "\nThus, we need to request n = %d / %d samples",
    n_candidates - n_bleeds_on_hand,
    n_candidates
  ))
}

return_soma_runs_at_nih <- function(DF_CANDIDATES) {
  DF_CANDIDATES <- DF_CANDIDATES %>%
    mutate(across(everything(), as.character))
  DF_AT_NIH <- query_db(DB_NORM, "SELECT * FROM bleeds_at_nih") %>%
    mutate(across(everything(), as.character))

  for_eric <- import_file(here::here("data/to_do"), "for_eric_july_28", "csv") %>%
    rename(id_sample = "Sample ID", antibody = "Suspected Antibody") %>%
    mutate(antibody = str_to_upper(antibody))
  for_eric <- filter_antibodies(for_eric, ASYS_MDA5)
  for_eric <- for_eric[for_eric[["Samples used in SOMA"]] != "", ] %>%
    filter("Samples used in SOMA" != "") %>%
    select(id_sample, antibody)

  # TODO: add bleed dates from Euroimmun to for_eric
  soma_dates <- query_db(DB_NORM, "SELECT DISTINCT id_patient_combined, bleed_date as bleed_date FROM antibodies_detailed WHERE source = 'Euroimmun' ")
  soma_dates <- left_join(soma_dates, query_db(DB_NORM, "SELECT id_patient_combined, id_sample FROM emr_ids"), relationship = "many-to-many")

  soma_dates <- inner_join(soma_dates, ab_dx, relationship = "many-to-many") %>%
    distinct(id_sample, antibody, bleed_date)

  # TODO: join for_eric and samples_at_nih, which accounts for all blood samples we have on hand and thus don't need to send
  for_eric_dates <- left_join(for_eric, soma_dates) %>% filter(!is.na(bleed_date))

  RETURN <- inner_join(DF_CANDIDATES, for_eric_dates) %>%
    distinct(id_sample, antibody, bleed_date)

  return(RETURN)
}

print_proposed_by_ab <- function(DF) {
  cat(sprintf("N = %s unique patients\n", nrow(DF %>% distinct(id_sample))))

  n_total_bleeds <- nrow(distinct(DF, id_sample, bleed_date))
  cat(sprintf("N = %d total bleeds\n\n", n_total_bleeds))

  for (iter in 1:nrow(distinct(DF, antibody))) {
    match_this_ab <- distinct(DF, antibody)[iter, ]
    n_ab_matches <- nrow(DF %>% filter(antibody %in% match_this_ab) %>% distinct(id_sample))
    n_total <- nrow(DF %>% distinct(id_sample))

    cat(sprintf(
      "n = %d / %d %s patients\n",
      n_ab_matches,
      n_total,
      match_this_ab
    ))
  }
}

get_ild_bleeds_in_range <- function(AB_LIST, MONTHS_ULN, MONTHS_INTERVAL, MONTHS_BUFFER){
  DAYS_ULN <- (MONTHS_ULN + MONTHS_BUFFER) * 30
  
  ab_dx <- query_db(DB_NORM, "SELECT DISTINCT * FROM antibody_diagnosis") %>%
    filter(chart_confirmed == "TRUE")
  all_bleeds <- query_db(DB_NORM, "SELECT DISTINCT * FROM blood_samples_ebo") %>% distinct(id_sample, bleed_date)
  
  ild_bleeds <- filter_antibodies(
    full_join(all_bleeds, ab_dx, relationship = "many-to-many"),
    AB_LIST) %>%
    distinct(id_sample, bleed_date, antibody)
  
  
  list_day_0 <- get_earliest(ild_bleeds) %>%
    rename(startdate = bleed_date)
    
  compiled <- full_join(list_day_0, ild_bleeds) %>%
    mutate(across(c(startdate, bleed_date), as.Date)) %>%
    mutate(date_diff = bleed_date - startdate) %>%
    filter(date_diff <= DAYS_ULN, id_sample != "NULL") %>%
    arrange(id_sample, startdate, bleed_date)
  
  # debug(shadow_set_of_ranges)
  all_bleeds_with_index <- shadow_set_of_ranges(compiled, MONTHS_ULN, MONTHS_BUFFER, MONTHS_INTERVAL) %>%
    arrange(id_sample, startdate, bleed_date, followup_day)
  
  # Ensures that listed bleeds have at least one follow up within range
  RESULT <- show_multiple_occurences(all_bleeds_with_index, "id_sample") %>%
    filter(!id_sample %in% long_soma_list[["id_sample"]])
  
  return(RESULT)
}

print_nums_table_graph <- function(DF) {
  print_proposed_by_ab(DF)
  
  DF %>%
    tabyl(followup_day, antibody) %>%
    adorn_totals(where = "col") %>%
    adorn_totals(where = "row") %>%
    adorn_percentages() %>%
    adorn_pct_formatting() %>%
    adorn_ns(position = "front") %>%
    flextable::flextable() %>%
    flextable::autofit()
  
  
  bleed_graph <- data.frame(table(DF[["antibody"]], DF[["followup_day"]])) %>%
    rename(antibody = Var1, followup_day = Var2) %>%
    arrange(antibody, followup_day)
  
  bleed_graph  %>%
    # filter(followup_day != 0) %>%
    ggplot(mapping = aes(x = antibody, y = Freq, fill = followup_day)) +
    geom_bar(position = "dodge", stat = "identity")+
    labs(title = "Frequencies, Including Day 0 Bleed")
}