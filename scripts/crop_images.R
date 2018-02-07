library(magick)

lfw <- data.frame(rep(NA, 100*100))
lfw_attributes <- read.table('lfw_attributes.txt', header = TRUE, sep = '\t', skip = 1)
trim_attributes <- rep(FALSE, dim(lfw_attributes)[1])

people <- list.files('LFW/lfw2', full.names = TRUE)

for(person in people)
{
    images <- list.files(person, full.names = TRUE)

    if (length(images) < 20) next

    i <- 1
    name <- strsplit(person, '/')[[1]][3]
    attrib_name<- gsub("_", " ", name)

    for (image in images)
    {

        img <- image_read(image)
        img <- image_crop(img, "100x100+150+150")
        img <- as.integer(img[[1]])
        dim(img) <- NULL

        attrib <- (lfw_attributes['person'] == attrib_name &
                   lfw_attributes['imagenum'] == i)
        trim_attributes <- trim_attributes | attrib

        if(!any(attrib))
        {
            print(label)
        }
        else
        {
            label <- paste(name, i, sep = '.')
            lfw[label] <- img
        }

        i <- i + 1
    }
}

lfw <- Filter(function(x)!all(is.na(x)), lfw)
save(lfw, file = "lfw.RData")

lfw_attributes <- lfw_attributes[trim_attributes, ]
lfw_attributes$person <- gsub(" ", "_", lfw_attributes$person)
row.names(lfw_attributes) <- paste(lfw_attributes$person, lfw_attributes$imagenum, sep='.')
lfw_attributes$person <- NULL
lfw_attributes$imagenum <- NULL
save(lfw_attributes, file = "lfw_attributes.RData")
