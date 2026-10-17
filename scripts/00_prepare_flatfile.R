# Reduces the full flat file to the records and columns used in this study.

source("config.R", local = TRUE)
source("R/utils.R", local = TRUE)

raw  <- read.csv(RAW_FLATFILE, stringsAsFactors = FALSE)
keep <- select_records(raw, ANALYSIS_TYPES)
slim <- keep[, c(FLATFILE_COLS, OBS_COLS)]

write.csv(slim, DATA_FILE, row.names = FALSE)
message(sprintf("%s: %d records, %d events, %d stations, %d columns.",
                DATA_FILE, nrow(slim), length(unique(slim$EQ_ID)),
                length(unique(slim$Sta)), ncol(slim)))
