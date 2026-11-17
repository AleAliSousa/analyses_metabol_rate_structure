# Optional acquisition, separate from the offline analysis runner.
# Existing files are checksum-verified; missing files are downloaded and checked.
library(digest)
m <- read.csv("metadata/nov14_source_downloads.csv",stringsAsFactors=FALSE)
options(timeout=max(600,getOption("timeout")))
cache <- new.env(parent=emptyenv())
for (i in seq_len(nrow(m))) {
  r <- m[i,]
  if (!file.exists(r$local_path)) {
    if (!exists(r$url,envir=cache,inherits=FALSE)) {
      tmp <- tempfile(fileext=if(r$container=="zip") ".zip" else ".download")
      download.file(r$url,tmp,mode="wb",quiet=FALSE)
      assign(r$url,tmp,envir=cache)
    }
    tmp <- get(r$url,envir=cache,inherits=FALSE)
    source <- tmp
    if (r$container!="file") {
      unpack <- tempfile(); dir.create(unpack)
      if (r$container=="zip") unzip(tmp,files=r$archive_member,exdir=unpack)
      else untar(tmp,files=r$archive_member,exdir=unpack)
      source <- file.path(unpack,r$archive_member)
    }
    if (!identical(digest(source,algo="sha256",file=TRUE),r$sha256))
      stop("Downloaded content changed; inspect before using: ",r$local_path)
    dir.create(dirname(r$local_path),recursive=TRUE,showWarnings=FALSE)
    stopifnot(file.copy(source,r$local_path,overwrite=FALSE))
  }
  stopifnot(identical(digest(r$local_path,algo="sha256",file=TRUE),r$sha256))
}
message("All ",nrow(m)," Nov 14 source files verified.")
