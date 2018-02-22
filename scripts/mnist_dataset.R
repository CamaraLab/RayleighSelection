# Load the MNIST digit recognition dataset into R
# http://yann.lecun.com/exdb/mnist/
# assume you have all 4 files and gunzip'd them
# creates train$n, train$x, train$y  and test$n, test$x, test$y
# e.g. train$x is a 60000 x 784 matrix, each row is one digit (28x28)
# call:  show_digit(train$x[5,])   to see a digit.
# brendan o'connor - gist.github.com/39760 - anyall.org

load_mnist <- function() {
  load_image_file <- function(filename) {
    ret = list()
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    ret$n = readBin(f,'integer',n=1,size=4,endian='big')
    nrow = readBin(f,'integer',n=1,size=4,endian='big')
    ncol = readBin(f,'integer',n=1,size=4,endian='big')
    x = readBin(f,'integer',n=ret$n*nrow*ncol,size=1,signed=F)
    ret$x = matrix(x, ncol=nrow*ncol, byrow=T)
    close(f)
    ret
  }
  load_label_file <- function(filename) {
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    n = readBin(f,'integer',n=1,size=4,endian='big')
    y = readBin(f,'integer',n=n,size=1,signed=F)
    close(f)
    y
  }
  train <<- load_image_file('mnist/train-images-idx3-ubyte')
  test <<- load_image_file('mnist/t10k-images-idx3-ubyte')
  
  train$y <<- load_label_file('mnist/train-labels-idx1-ubyte')
  test$y <<- load_label_file('mnist/t10k-labels-idx1-ubyte')  
}


show_digit <- function(arr784, col=gray(12:1/12), ...) {
  image(matrix(arr784, nrow=28)[,28:1], col=col, ...)
}

generate_names <- function(dataset)
{
  counts <- rep(1, 10)
  names(counts) <- 0:9
  labels <- character(length(dataset))

  for (idx in seq(dataset))
  {
    num <- dataset[[idx]]
    count <- counts[[num+1]]
    lbl <- paste(c(toString(num), toString(count)), collapse='.')
    labels[idx] <- lbl
    counts[[num+1]] <- count + 1
  }
  return (labels)
}

load_mnist()

mnist_train <- data.frame(train$x)
row.names(mnist_train) <- generate_names(train$y)
mnist_train <- t(mnist_train)
save(mnist_train, file = "mnist_train.RData")

mnist_test <- data.frame(test$x)
row.names(mnist_test) <- generate_names(test$y)
mnist_test <- t(mnist_test)
save(mnist_test, file = "mnist_test.RData")
