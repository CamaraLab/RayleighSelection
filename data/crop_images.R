library(magick)

pixels = list()

people <- list.files('LFW/lfw2', full.names = TRUE)

save_images = FALSE

for(person in people)
{
    images <- list.files(person, full.names = TRUE)

    if (length(images) < 20) next

    for (image in images)
    {
        img <- image_read(image)
        img <- image_crop(img, "100x100+150+150")

        pixels <- cbind(pixels, as.integer(img[[1]]))

        if (save_images)
        {
            pth <- strsplit(image, '/')[[1]]
            person <- pth[[3]]
            image <- pth[[4]]

            pth <- file.path('cropped_images', person)

            if (!dir.exists(pth)) dir.create(pth, recursive = TRUE)

            pth <- file.path(pth, image)

            image_write(img, path = pth)
        }
    }
}

lfw <- as.matrix(pixels, nrows = 100*100)
save(lfw, file = "lfw.RData")

image_distances <- dist(t(lfw))
image_distances <- as.matrix(image_distances)
save(image_distances, file = "image_distances.RData")
