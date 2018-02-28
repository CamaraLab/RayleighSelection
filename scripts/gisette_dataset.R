generate_names <- function(dataset)
{
  counts <- c(1, 1)
  names(counts) <- c(4, 9)
  labels <- character(length(dataset))

  for (idx in seq(dataset))
  {
    num <- toString(dataset[[idx]])
    count <- counts[[num]]
    lbl <- paste(c(num, toString(count)), collapse='.')
    labels[idx] <- lbl
    counts[[num]] <- count + 1
  }

  return (labels)
}

gisette_train <- read.table(file='gisette/gisette_train.data')
gisette_train_labels <- read.table(file='gisette/gisette_train.labels')
gisette_train_labels[gisette_train_labels == 1] <- 4
gisette_train_labels[gisette_train_labels == -1] <- 9
gisette_train <- t(gisette_train)
colnames(gisette_train) <- generate_names(unlist(gisette_train_labels))
save(gisette_train, file='gisette_train.RData')

gisette_valid <- read.table(file='gisette/gisette_valid.data')
gisette_valid_labels <- read.table(file='gisette/gisette_valid.labels')
gisette_valid_labels[gisette_valid_labels == 1] <- 4
gisette_valid_labels[gisette_valid_labels == -1] <- 9
gisette_valid <- t(gisette_valid)
colnames(gisette_valid) <- generate_names(unlist(gisette_valid_labels))
save(gisette_valid, file='gisette_valid.RData')
