# fixes a bug in DECIPHER
# contruibuted by Eric Scott Wright, author of DECIPHER

.subset <- function(x, i) {
	ans <- new(class(x))
	g <- x@ranges@group[i]
	u <- unique(g)
	ans@pool <- x@pool[u]
	ans@ranges@group <- match(g, u)
	ans@ranges@start <- x@ranges@start[i]
	ans@ranges@width <- x@ranges@width[i]
	ans
}

.replace <- function(x, y, index) {
	ans <- new(class(x))
	ans@pool <- c(x@pool,
		y@pool)
	g <- .Call("getPools",
		ans@pool@xp_list,
		PACKAGE="DECIPHER")
	w <- which(!duplicated(g))
	ans@pool <- ans@pool[w]
	m <- match(g, g[w])
	ans@ranges@group <- x@ranges@group
	ans@ranges@group[index] <- y@ranges@group + length(x@pool)
	g <- m[ans@ranges@group]
	u <- unique(g)
	ans@pool <- ans@pool[u]
	ans@ranges@group <- match(g, u)
	ans@ranges@start <- x@ranges@start
	ans@ranges@start[index] <- y@ranges@start
	ans@ranges@width <- x@ranges@width
	ans@ranges@width[index] <- y@ranges@width
	ans
}

.append <- function(x, y) {
	ans <- new(class(x))
	ans@pool <- c(x@pool,
		y@pool)
	ans@ranges@group <- c(x@ranges@group,
		y@ranges@group + length(x@pool))
	ans@ranges@start <- c(x@ranges@start,
		y@ranges@start)
	ans@ranges@width <- c(x@ranges@width,
		y@ranges@width)
	ans
}

.switch <- function(x) {
	if (class(x)=="DNAStringSet") {
		ans <- new("RNAStringSet")
	} else {
		ans <- new("DNAStringSet")
	}
	ans@pool <- x@pool
	ans@ranges@group <- x@ranges@group
	ans@ranges@start <- x@ranges@start
	ans@ranges@width <- x@ranges@width
	ans
}
