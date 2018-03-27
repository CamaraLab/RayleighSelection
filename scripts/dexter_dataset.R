dexter_train <- read.table(file='dexter_train.csv', header=TRUE, sep=',')
dexter_valid <- read.table(file='dexter_valid.csv', header=TRUE, sep=',')

save(dexter_train, file='dexter_train.RData')
save(dexter_valid, file='dexter_valid.RData')
