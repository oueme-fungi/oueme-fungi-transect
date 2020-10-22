# updated version of AlignSeqs from DECIPHER package,
# contributed by Eric Scott Wright (author of DECIPHER)

AlignSeqs <- function(myXStringSet,
	guideTree=NULL,
	iterations=2,
	refinements=1,
	gapOpening=c(-18, -16),
	gapExtension=c(-2, -1),
	useStructures=TRUE,
	structures=NULL,
	FUN=AdjustAlignment,
	levels=c(0.9, 0.7, 0.7, 0.4, 10, 5, 5, 2),
	alphabet=AA_REDUCED[[1]],
	processors=1,
	verbose=TRUE,
	...) {

	# error checking
	type <- switch(class(myXStringSet),
		`DNAStringSet` = 1L,
		`RNAStringSet` = 2L,
		`AAStringSet` = 3L,
		stop("pattern must be an AAStringSet, DNAStringSet, or RNAStringSet."))
	l <- length(myXStringSet)
	if (l < 2)
		stop("At least two sequences are required in myXStringSet.")
	if (!is.null(guideTree)) {
		if (class(guideTree)=="dendrogram") {
			if (attr(guideTree, "height") > 0.501)
				stop("Total height of guideTree can be at most 0.5.")
			heights <- rapply(guideTree,
				function(x)
					attr(x, "height"))
			if (any(heights > 0.001) || any(heights < -0.001))
				stop("All tips of guideTree must have a height of zero.")
			labels <- rapply(guideTree,
				function(x)
					attr(x, "label"))
			if (any(duplicated(labels)))
				stop("Leaf labels in guideTree must be unique.")
			labels <- match(labels, names(myXStringSet))
			if (any(is.na(labels)))
				stop("Leaf labels in guideTree must match names of myXStringSet.")
			guideTree <- rapply(guideTree,
				function(x) {
					if (is.leaf(x))
						x[1] <- match(attr(x, "label"),
							names(myXStringSet))
					return(x)
				},
				how="replace")
		} else {
			stop("guideTree must be a dendrogram object.")
		}
	}
	if (!is.numeric(iterations))
		stop("iterations must be a numeric.")
	if (iterations != floor(iterations))
		stop("iterations must be a whole number.")
	if (iterations < 0)
		stop("iterations must be at least zero.")
	if (!is.numeric(refinements))
		stop("refinements must be a numeric.")
	if (refinements != floor(refinements))
		stop("refinements must be a whole number.")
	if (refinements < 0)
		stop("refinements must be at least zero.")
	if (!is.logical(useStructures))
		stop("useStructures must be a logical.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') must be removed before alignment.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') must be removed before alignment.")
	FUN <- match.fun(FUN)
	if (!is.numeric(levels))
		stop("levels must be a numeric vector.")
	if (length(levels) != 8)
		stop("levels must be length 8.")
	if (levels[1] < 0 || levels[1] > 1)
		stop("levels[1] must be between zero and one (inclusive).")
	if (levels[2] < 0 || levels[2] > 1)
		stop("levels[2] must be between zero and one (inclusive).")
	if (levels[3] < 0 || levels[3] > 1)
		stop("levels[3] must be between zero and one (inclusive).")
	if (levels[4] < 0 || levels[4] > 1)
		stop("levels[4] must be between zero and one (inclusive).")
	if (levels[5] <= 0 || is.infinite(levels[5]))
		stop("levels[5] must be between zero and Inf (exclusive).")
	if (levels[5] != floor(levels[5]))
		stop("levels[5] must be a whole number.")
	if (levels[6] != floor(levels[6]))
		stop("levels[6] must be a whole number.")
	if (levels[6] < 0)
		stop("levels[6] must be at least zero.")
	if (levels[7] <= 0)
		stop("levels[7] must be at least zero.")
	if (levels[8] <= 0)
		stop("levels[8] must be at least zero.")
	if (any(alphabet==""))
		stop("No elements of alphabet can be empty.")
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors)!=processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- detectCores()
	} else {
		processors <- as.integer(processors)
	}

	args <- list(...)
	n <- names(args)
	m <- character(length(n))
	for (i in seq_along(n)) {
		m[i] <- match.arg(n[i],
			names(formals(AlignProfiles)))
	}

	if (length(gapOpening)==2) {
		gapOpeningMin <- gapOpening[1]
		gapOpeningMax <- gapOpening[2]
	} else if (length(gapOpening)==1) {
		gapOpeningMin <- gapOpeningMax <- gapOpening
	} else {
		stop("gapOpening must be a vector of length 1 or 2.")
	}
	if (length(gapExtension)==2) {
		gapExtensionMin <- gapExtension[1]
		gapExtensionMax <- gapExtension[2]
	} else if (length(gapExtension)==1) {
		gapExtensionMin <- gapExtensionMax <- gapExtension
	} else {
		stop("gapExtension must be a vector of length 1 or 2.")
	}
	if (!is.numeric(gapOpeningMax))
		stop("gapOpening must be a numeric.")
	if (gapOpeningMin > gapOpeningMax)
		stop("gapOpening[1] must be less than or equal to gapOpening[2].")
	if (!is.numeric(gapExtensionMax))
		stop("gapExtension must be a numeric.")
	if (gapExtensionMin > gapExtensionMax)
		stop("gapExtension[1] must be less than or equal to gapExtension[2].")
	gapOpeningSlope <- gapOpeningMax - gapOpeningMin
	gapExtensionSlope <- gapExtensionMax - gapExtensionMin

	# prepare structures and structure matrix
	if (useStructures) {
		if (type==3L) {
			if (is.null(structures)) {
				structures <- PredictHEC(myXStringSet,
					type="probabilities")
			} else {
				if (length(structures) != l)
					stop("structures is not the same length as myXStringSet.")
				if (typeof(structures) != "list")
					stop("structures must be a list.")
			}
			if (refinements > 0) {
				w <- which(m=="structureMatrix")
				if (length(w) > 0) {
					structureMatrix <- args[[w]]
					# assume structures and matrix are ordered the same
					if (!is.double(structureMatrix))
						stop("structureMatrix must contain numerics.")
					if (!is.matrix(structureMatrix))
						stop("structureMatrix must be a matrix.")
					if (dim(structureMatrix)[1] != dim(structureMatrix)[2])
						stop("structureMatrix is not square.")
				} else {
					structureMatrix <- matrix(c(6, 1, -2, 1, 13, 0, -2, 0, 1),
						nrow=3) # order is H, E, C
				}
				if (dim(structureMatrix)[1] != dim(structures[[1]])[1])
					stop("Dimensions of structureMatrix are incompatible with the structures.")
			}
		} else {
			if (!is.null(structures)) {
				if (length(structures) != l)
					stop("structures is not the same length as myXStringSet.")
				if (typeof(structures) != "list")
					stop("structures must be a list.")

				w <- which(m=="structureMatrix")
				if (length(w) > 0) {
					structureMatrix <- args[[w]]
					# assume structures and matrix are ordered the same
					if (!is.double(structureMatrix))
						stop("structureMatrix must contain numerics.")
					if (!is.matrix(structureMatrix))
						stop("structureMatrix must be a matrix.")
					if (dim(structureMatrix)[1] != dim(structureMatrix)[2])
						stop("structureMatrix is not square.")
				} else {
					stop("structureMatrix must be specified when structures are provided.")
				}
				if (dim(structureMatrix)[1] != dim(structures[[1]])[1])
					stop("Dimensions of structureMatrix are incompatible with the structures.")
				replace <- TRUE
			} else if (type==2L) {
				w <- which(m=="structureMatrix")
				if (length(w) > 0) {
					structureMatrix <- args[[w]]
					# assume the matrix is in the correct order
					if (!is.double(structureMatrix))
						stop("structureMatrix must contain numerics.")
					if (!is.matrix(structureMatrix))
						stop("structureMatrix must be a matrix.")
					if (dim(structureMatrix)[1] != dim(structureMatrix)[2])
						stop("structureMatrix is not square.")
					if (dim(structureMatrix)[1] != 3)
						stop("structureMatrix must be 3 x 3 when structures is NULL.")
				} else { # use the default structureMatrix
					structureMatrix <- matrix(c(6, 0, 0, 0, 30, -6, 0, -6, 30),
						nrow=3) # order is ., (, )
				}
				replace <- FALSE
			} else {
				structureMatrix <- numeric()
				useStructures <- FALSE
			}
		}
	} else {
		w <- which(m=="structureMatrix")
		if (length(w) > 0)
			stop("structureMatrix provided when useStructures is FALSE.")
		structureMatrix <- numeric() # needed for colScores
		if (!is.null(structures))
			stop("structures provided when useStructures is FALSE.")
	}

	# prepare substitution matrix
	if (type==3L) {
		r <- strsplit(alphabet, "", fixed=TRUE)
		alphabet <- setNames(rep(0L, 20),
			AA_STANDARD)
		for (i in seq_along(r)) {
			w <- which(!(r[[i]] %in% AA_STANDARD))
			if (length(w) > 0)
				stop("Unrecognized letter(s) found in alphabet:  ",
					paste(r[[i]][w], collapse=", "),
					".")
			w <- which(alphabet[r[[i]]] != 0L)
			if (length(w) > 0)
				stop("Repeated amino acids found in alphabet:  ",
					paste(r[[i]][w], collapse=", "),
					".")
			alphabet[r[[i]]] <- i
		}
		w <- which(alphabet==0L)
		if (length(w) > 0)
			stop("Standard amino acids missing from alphabet:  ",
				paste(names(w), collapse=", "),
				".")
		sizeAA <- max(alphabet)
		if (sizeAA==1)
			stop("More than one grouping of amino acids is required in the alphabet.")
		sizeAA <- as.integer(floor(log(4294967295, sizeAA)))
		alphabet <- alphabet - 1L

		subM <- TRUE
		w <- which(m=="substitutionMatrix")
		AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
		if (length(w)==1L) {
			sM <- args[[w]]
			args[w] <- NULL
			if (is.character(sM)) {
				if (!(sM %in% c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100", "PAM30", "PAM40", "PAM70", "PAM120", "PAM250", "MIQS")))
					stop("Invalid substitutionMatrix.")
			}
			AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
				"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
			if (is.matrix(sM)) {
				if (any(!(AAs %in% dimnames(sM)[[1]])) ||
					any(!(AAs %in% dimnames(sM)[[2]])))
					stop("substitutionMatrix is incomplete.")
			} else {
				sM <- eval(parse(text=data(list=sM, envir=environment(), package=ifelse(substitutionMatrix=="MIQS", "DECIPHER", "Biostrings"))))
			}
			sM <- sM[AAs, AAs]
			sM <- sM + 0 # convert to numeric matrix
		} else {
			# use PFASUM50
			sM <- matrix(c(4.1181,-1.1516,-1.3187,-1.4135,0.4271,-0.5467,-0.6527,0.1777,-1.6582,-1.1243,-1.1843,-1.0235,-0.5685,-1.9515,-0.6072,0.8284,0.0361,-2.5368,-2.1701,0.0661,-11,-1.1516,6.341,0.0543,-0.6628,-3.2085,1.6006,0.5067,-1.961,0.7706,-3.5053,-3.0357,2.938,-1.9894,-3.7846,-1.3455,-0.4194,-0.5594,-2.1629,-1.7957,-2.9403,-11,-1.3187,0.0543,6.4672,2.3024,-2.5179,0.8192,0.5566,0.1585,1.104,-4.1629,-4.0977,0.8743,-2.6216,-3.805,-1.0904,1.1291,0.3253,-3.7763,-1.874,-3.6076,-11,-1.4135,-0.6628,2.3024,6.8156,-4.358,0.6705,2.582,-0.5667,-0.196,-5.475,-5.1661,0.226,-3.9595,-5.3456,-0.5662,0.4273,-0.5218,-4.7691,-3.4644,-4.5477,-11,0.4271,-3.2085,-2.5179,-4.358,13.5349,-3.3641,-4.3086,-2.1614,-1.8945,-0.7546,-0.9453,-3.8239,-0.5923,-0.8182,-3.6019,-0.3927,-0.801,-1.9317,-1.1607,0.0673,-11,-0.5467,1.6006,0.8192,0.6705,-3.3641,5.5795,2.1372,-1.5923,1.0862,-3.3001,-2.7545,1.872,-1.1216,-3.6631,-1.0426,0.1982,-0.0434,-3.061,-1.9214,-2.6993,-11,-0.6527,0.5067,0.5566,2.582,-4.3086,2.1372,5.5684,-1.6462,-0.2488,-4.1849,-4.0275,1.4821,-2.7964,-4.8311,-0.7028,0.0283,-0.312,-4.1969,-2.9489,-3.281,-11,0.1777,-1.961,0.1585,-0.5667,-2.1614,-1.5923,-1.6462,7.6508,-1.8185,-4.7058,-4.4215,-1.5991,-3.2786,-3.9992,-1.4409,0.184,-1.4823,-3.8328,-3.7343,-3.7264,-11,-1.6582,0.7706,1.104,-0.196,-1.8945,1.0862,-0.2488,-1.8185,9.7543,-3.3812,-2.8685,0.1425,-1.8724,-1.2545,-1.5333,-0.4285,-0.8896,-0.9385,1.6476,-2.8729,-11,-1.1243,-3.5053,-4.1629,-5.475,-0.7546,-3.3001,-4.1849,-4.7058,-3.3812,5.1229,2.5319,-3.5454,1.8309,0.9346,-3.4603,-3.0985,-1.2543,-1.5006,-1.117,3.3961,-11,-1.1843,-3.0357,-4.0977,-5.1661,-0.9453,-2.7545,-4.0275,-4.4215,-2.8685,2.5319,4.7049,-3.4581,2.5303,1.687,-3.365,-3.1578,-1.8626,-0.5308,-0.6881,1.4829,-11,-1.0235,2.938,0.8743,0.226,-3.8239,1.872,1.4821,-1.5991,0.1425,-3.5454,-3.4581,5.5476,-2.164,-4.3516,-0.7583,0.0275,-0.1516,-3.5889,-2.4422,-3.0453,-11,-0.5685,-1.9894,-2.6216,-3.9595,-0.5923,-1.1216,-2.7964,-3.2786,-1.8724,1.8309,2.5303,-2.164,7.0856,1.2339,-3.0823,-1.7587,-0.7402,-0.5841,-0.3946,0.9477,-11,-1.9515,-3.7846,-3.805,-5.3456,-0.8182,-3.6631,-4.8311,-3.9992,-1.2545,0.9346,1.687,-4.3516,1.2339,7.4322,-3.6222,-3.0316,-2.2851,2.6305,3.8302,0.1942,-11,-0.6072,-1.3455,-1.0904,-0.5662,-3.6019,-1.0426,-0.7028,-1.4409,-1.5333,-3.4603,-3.365,-0.7583,-3.0823,-3.6222,9.1796,-0.0652,-0.8587,-3.3634,-3.3006,-2.5443,-11,0.8284,-0.4194,1.1291,0.4273,-0.3927,0.1982,0.0283,0.184,-0.4285,-3.0985,-3.1578,0.0275,-1.7587,-3.0316,-0.0652,4.2366,1.8491,-3.1454,-2.1838,-2.1839,-11,0.0361,-0.5594,0.3253,-0.5218,-0.801,-0.0434,-0.312,-1.4823,-0.8896,-1.2543,-1.8626,-0.1516,-0.7402,-2.2851,-0.8587,1.8491,4.8833,-2.8511,-1.8993,-0.2699,-11,-2.5368,-2.1629,-3.7763,-4.7691,-1.9317,-3.061,-4.1969,-3.8328,-0.9385,-1.5006,-0.5308,-3.5889,-0.5841,2.6305,-3.3634,-3.1454,-2.8511,13.6485,3.3017,-1.851,-11,-2.1701,-1.7957,-1.874,-3.4644,-1.1607,-1.9214,-2.9489,-3.7343,1.6476,-1.117,-0.6881,-2.4422,-0.3946,3.8302,-3.3006,-2.1838,-1.8993,3.3017,8.7568,-1.2438,-11,0.0661,-2.9403,-3.6076,-4.5477,0.0673,-2.6993,-3.281,-3.7264,-2.8729,3.3961,1.4829,-3.0453,0.9477,0.1942,-2.5443,-2.1839,-0.2699,-1.851,-1.2438,4.6928,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,14),
				nrow=21,
				ncol=21,
				dimnames=list(AAs, AAs))
		}
	} else {
		bases <- c("A", "C", "G",
			ifelse(type==2L, "U", "T"))
		w <- which(m=="substitutionMatrix")
		if (length(w) > 0) {
			subM <- TRUE
			sM <- args[[w]]
			args[w] <- NULL
			if (is.matrix(sM)) {
				if (any(!(bases %in% dimnames(sM)[[1]])) ||
					any(!(bases %in% dimnames(sM)[[2]])))
					stop("substitutionMatrix is incomplete.")
				sM <- sM[bases, bases]
				sM <- sM + 0 # convert to numeric matrix
				if (type==2L) { # RNAStringSet
					sM2 <- sM
					dimnames(sM2) <- list(DNA_BASES, DNA_BASES)
				}
			} else {
				stop("substitutionMatrix must be NULL or a matrix.")
			}
		} else if (type==2L &&
			!("misMatch" %in% m) &&
			!("perfectMatch" %in% m)) {
			subM <- TRUE
			sM <- matrix(c(14, 4, 7, 4, 4, 15, 4, 7, 7, 4, 15, 4, 4, 7, 4, 14),
				nrow=4,
				ncol=4,
				dimnames=list(bases, bases))
			sM2 <- matrix(c(14, 4, 7, 4, 4, 15, 4, 7, 7, 4, 15, 4, 4, 7, 4, 14),
				nrow=4,
				ncol=4,
				dimnames=list(DNA_BASES, DNA_BASES))
		} else {
			subM <- FALSE
			w <- which(m=="perfectMatch")
			if (length(w) > 0) {
				PM <- args[[w]]
				if (!is.numeric(PM))
					stop("perfectMatch must be a numeric.")
			} else {
				PM <- 5
			}
			w <- which(m=="misMatch")
			if (length(w) > 0) {
				MM <- args[[w]]
				if (!is.numeric(MM))
					stop("misMatch must be a numeric.")
			} else {
				MM <- 0
			}
			sM <- matrix(as.numeric(MM),
				nrow=4,
				ncol=4,
				dimnames=list(bases, bases))
			diag(sM) <- PM
		}
	}

	if (length(args)==0)
		args <- NULL

	if (verbose)
		time.1 <- Sys.time()

	w.x <- width(myXStringSet)
	if (type==3L) {
		wordSize <- ceiling(log(100*quantile(w.x, 0.99),
			.Call("alphabetSizeReducedAA",
				myXStringSet,
				alphabet,
				PACKAGE="DECIPHER")))
		if (wordSize > sizeAA)
			wordSize <- sizeAA
		if (wordSize < 1)
			wordSize <- 1
	} else {
		wordSize <- ceiling(log(100*quantile(w.x, 0.99),
			.Call("alphabetSize",
				myXStringSet,
				PACKAGE="DECIPHER")))
		if (wordSize > 15)
			wordSize <- 15
		if (wordSize < 2)
			wordSize <- 2
	}

	if (is.null(guideTree)) {
		if (verbose) {
			cat("Determining distance matrix based on shared ",
				wordSize,
				"-mers:\n",
				sep="")
			flush.console()
			pBar <- txtProgressBar(max=100, style=ifelse(interactive(), 3, 1))
		} else {
			pBar <- NULL
		}

		if (type==3L) {
			v <- .Call("enumerateSequenceReducedAA",
				myXStringSet,
				wordSize,
				alphabet,
				PACKAGE="DECIPHER")
		} else { # DNAStringSet or RNAStringSet
			v <- .Call("enumerateSequence",
				myXStringSet,
				wordSize,
				PACKAGE="DECIPHER")
		}

		d <- .Call("matchOrder",
			v,
			verbose,
			pBar,
			processors,
			PACKAGE="DECIPHER")
		attr(d, "Size") <- l
		attr(d, "Diag") <- TRUE
		attr(d, "Upper") <- TRUE
		class(d) <- "dist"

		if (verbose) {
			setTxtProgressBar(pBar, 100)
			close(pBar)
			cat("\n")
			time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
		}

		if (verbose) {
			cat("\nClustering into groups by similarity:\n")
			flush.console()
		}
		dimnames(d) <- NULL
		suppressWarnings(guideTree <- IdClusters(d,
			method="single",
			type="dendrogram",
			verbose=verbose,
			processors=processors))
	}

	if (verbose) {
		time.1 <- Sys.time()
		cat("Aligning Sequences:\n")
		flush.console()
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1), max=100)
		before <- steps <- 0L
		nsteps <- l - 1L
	}

	.align <- function(guideTree) { # iteratively align sequences
		# initialize a stack of maximum length (l)
		stack <- vector("list", l)
		visit <- logical(l) # node already visited
		parent <- integer(l) # index of parent node
		index <- integer(l) # index in parent node
		pos <- 1L # current position in the stack
		stack[[pos]] <- guideTree
		while (pos > 0L) { # more nodes to visit
			if (visit[pos]) { # ascending tree
				visit[pos] <- FALSE # reset visit
				dend <- stack[[pos]]

				# align subtrees
				treeLengths <- numeric(length(dend))
				inherit <- attr(dend, "inherit")
				if (!is.null(inherit)) {
					if (inherit) {
						u <- unlist(dend)
						seqs <<- .replace(seqs,
							.Call("removeCommonGaps",
								.subset(seqs_prev, u),
								type,
								processors,
								PACKAGE="DECIPHER"),
							u)

						if (verbose) {
							steps <<- steps + 1L
							percentComplete <- as.integer(100L*steps/nsteps)
							if (percentComplete > before) {
								setTxtProgressBar(pBar, percentComplete)
								before <<- percentComplete
							}
						}
					}

					if (length(dend) > 1) {
						for (i in seq_len(length(dend))) {
							h <- attr(dend[[i]], "treeLength")
							if (!is.null(h))
								treeLengths[i] <- h
						}

						h <- attr(dend, "height")
						for (i in seq_len(length(dend))) {
							m <- unlist(dend[i])
							treeLengths[i] <- treeLengths[i] + h - heights[m[1]]
							weights[m] <<- weights[m] + (h - heights[m])/length(m)
							heights[m] <<- h
						}
					} else if (is.leaf(dend)) {
						heights[unlist(dend)] <- attr(dend, "height")
					} else { # inherit from subtree
						treeLengths[1] <- attr(dend[[1]], "treeLength")
					}
				} else if (length(dend) > 1) {
					for (i in seq_len(length(dend))) {
						h <- attr(dend[[i]], "treeLength")
						if (!is.null(h))
							treeLengths[i] <- h
					}

					h <- attr(dend, "height")
					members <- vector("list", length(dend))
					for (i in seq_len(length(dend))) {
						m <- unlist(dend[i])
						treeLengths[i] <- treeLengths[i] + h - heights[m[1]]
						weights[m] <<- weights[m] + (h - heights[m])/length(m)
						heights[m] <<- h
						members[[i]] <- m
					}

					h <- h*2 # total length of sub-tree
					GO <- h*gapOpeningSlope + gapOpeningMin
					GE <- h*gapExtensionSlope + gapExtensionMin

					for (i in 2:length(dend)) {
						x <- unlist(members[1:(i - 1)])
						y <- members[[i]]
						combo <- c(x, y)

						p.weight <- weights[x]
						w <- which(p.weight <= 0)
						if (length(w) > 0)
							p.weight[w] <- 1
						p.weight <- p.weight/mean(p.weight)
						s.weight <- weights[y]
						w <- which(s.weight <= 0)
						if (length(w) > 0)
							s.weight[w] <- 1
						s.weight <- s.weight/mean(s.weight)

						pattern <- .subset(seqs, x)
						subject <- .subset(seqs, y)

						if (subM) {
							if (useStructures) {
								if (is.null(structures)) {
									if (type==2L && h < LEVEL) {
										# align as DNAStringSet (faster because no pairs matrix)
										temp <- .switch(do.call(AlignProfiles,
											args=c(list(pattern=.switch(pattern),
													subject=.switch(subject),
													p.weight=p.weight,
													s.weight=s.weight,
													substitutionMatrix=sM2,
													processors=processors,
													gapOpening=GO,
													gapExtension=GE),
												args)))
									} else {
										temp <- do.call(AlignProfiles,
											args=c(list(pattern=pattern,
													subject=subject,
													p.weight=p.weight,
													s.weight=s.weight,
													substitutionMatrix=sM,
													processors=processors,
													gapOpening=GO,
													gapExtension=GE),
												args))
									}
								} else {
									if (type==2L && h < LEVEL) {
										# align as DNAStringSet (faster because no pairs matrix)
										temp <- .switch(do.call(AlignProfiles,
											args=c(list(pattern=.switch(pattern),
													subject=.switch(subject),
													p.weight=p.weight,
													s.weight=s.weight,
													p.struct=structures[x],
													s.struct=structures[y],
													substitutionMatrix=sM2,
													processors=processors,
													gapOpening=GO,
													gapExtension=GE),
												args)))
									} else {
										temp <- do.call(AlignProfiles,
											args=c(list(pattern=pattern,
													subject=subject,
													p.weight=p.weight,
													s.weight=s.weight,
													p.struct=structures[x],
													s.struct=structures[y],
													substitutionMatrix=sM,
													processors=processors,
													gapOpening=GO,
													gapExtension=GE),
												args))
									}
								}
							} else {
								if (type==2L && h < LEVEL) {
									# align as DNAStringSet (faster because no pairs matrix)
									temp <- .switch(do.call(AlignProfiles,
										args=c(list(pattern=.switch(pattern),
												subject=.switch(subject),
												p.weight=p.weight,
												s.weight=s.weight,
												substitutionMatrix=sM2,
												processors=processors,
												gapOpening=GO,
												gapExtension=GE),
											args)))
								} else {
									temp <- do.call(AlignProfiles,
										args=c(list(pattern=pattern,
												subject=subject,
												p.weight=p.weight,
												s.weight=s.weight,
												substitutionMatrix=sM,
												processors=processors,
												gapOpening=GO,
												gapExtension=GE),
											args))
								}
							}
						} else {
							if (useStructures) {
								temp <- do.call(AlignProfiles,
									args=c(list(pattern=pattern,
											subject=subject,
											p.weight=p.weight,
											s.weight=s.weight,
											p.struct=structures[x],
											s.struct=structures[y],
											processors=processors,
											gapOpening=GO,
											gapExtension=GE),
										args))
							} else {
								temp <- do.call(AlignProfiles,
									args=c(list(pattern=pattern,
											subject=subject,
											p.weight=p.weight,
											s.weight=s.weight,
											processors=processors,
											gapOpening=GO,
											gapExtension=GE),
										args))
							}
						}

						if (h > LEVEL &&
							length(temp) >= levels[6]) {
							weight <- weights[combo]
							w <- which(weight <= 0)
							if (length(w) > 0)
								weight[w] <- 1
							weight <- weight/mean(weight)

							if (subM) {
								temp <- FUN(temp,
									substitutionMatrix=sM,
									weight=weight,
									processors=processors)
							} else {
								temp <- FUN(temp,
									weight=weight,
									processors=processors)
							}
						}

						seqs <<- .replace(seqs,
							temp,
							combo)

						if (verbose) {
							steps <<- steps + 1L
							percentComplete <- as.integer(100L*steps/nsteps)
							if (percentComplete > before) {
								setTxtProgressBar(pBar, percentComplete)
								before <<- percentComplete
							}
						}
					}
				} else if (is.leaf(dend)) {
					heights[unlist(dend)] <- attr(dend, "height")
				} else { # inherit from subtree
					treeLengths[1] <- attr(dend[[1]], "treeLength")
				}

				attr(stack[[pos]], "treeLength") <- sum(treeLengths)
				# replace self in parent
				if (parent[pos] > 0)
					stack[[parent[pos]]][[index[pos]]] <- stack[[pos]]
				pos <- pos - 1L # pop off of stack
			} else { # descending tree
				visit[pos] <- TRUE
				p <- pos
				for (i in seq_along(stack[[p]])) {
					if (!is.leaf(stack[[p]][[i]])) {
						# push subtree onto stack
						pos <- pos + 1L
						stack[[pos]] <- stack[[p]][[i]]
						parent[[pos]] <- p
						index[[pos]] <- i
					}
				}
			}
		}

		return(attr(stack[[1]], "treeLength"))
	}

	.reorder <- function(dend) {
		# initialize a stack of maximum length (l)
		stack <- vector("list", l)
		visit <- logical(l) # node already visited
		parent <- integer(l) # index of parent node
		index <- integer(l) # index in parent node
		pos <- 1L # current position in the stack
		stack[[pos]] <- dend
		while (pos > 0L) { # more nodes to visit
			if (visit[pos]) { # ascending tree
				visit[pos] <- FALSE # reset visit

				# sort tree by descending width
				members <- lapply(stack[[pos]], unlist)
				o <- order(sapply(members,
						function(x)
							max(w.x[x])),
					lengths(members),
					sapply(members, min),
					decreasing=TRUE)
				stack[[pos]][] <- stack[[pos]][o]

				# replace self in parent
				if (parent[pos] > 0)
					stack[[parent[pos]]][[index[pos]]] <- stack[[pos]]
				pos <- pos - 1L # pop off of stack
			} else { # descending tree
				visit[pos] <- TRUE
				p <- pos
				for (i in seq_along(stack[[p]])) {
					if (!is.leaf(stack[[p]][[i]])) {
						# push subtree onto stack
						pos <- pos + 1L
						stack[[pos]] <- stack[[p]][[i]]
						parent[[pos]] <- p
						index[[pos]] <- i
					}
				}
			}
		}
		return(stack[[1L]])
	}
	guideTree <- .reorder(guideTree)

	ns <- names(myXStringSet)
	seqs <- myXStringSet
	weights <- heights <- numeric(l)
	if (type==3L) {
		LEVEL <- levels[1]
	} else {
		LEVEL <- levels[3]
	}
	treeLength <- .align(guideTree)

	if (refinements==0 && iterations==0) {
		if (verbose) {
			setTxtProgressBar(pBar, 100)
			close(pBar)
			cat("\n")
			time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			flush.console()
			cat("\n")
		}

		names(seqs) <- ns
		return(seqs)
	}

	.RNAStructures <- function(seqs, weights) {
		w <- which(weights <= 0)
		if (length(w) > 0)
			weights[w] <- 1
		weights <- weights/mean(weights)

		if (replace) {
			structureMatrix <- matrix(c(6, 0, 0, 0, 30, -6, 0, -6, 30),
				nrow=3) # order is ., (, )
			replace <- FALSE
		}

		PredictDBN(seqs,
			type="search",
			weight=weights,
			processors=processors,
			verbose=verbose)
	}

	minTreeLength <- levels[7]

	if (iterations > 0) {
		if (type==3L) {
			LEVEL <- levels[2]
		} else {
			LEVEL <- levels[4]
		}

		.order <- function(dend) {
			# initialize a stack of maximum length (l)
			stack <- vector("list", l)
			visit <- logical(l) # node already visited
			pos <- 1L # current position in the stack
			stack[[pos]] <- dend
			while (pos > 0L) { # more nodes to visit
				if (visit[pos]) { # ascending tree
					visit[pos] <- FALSE # reset visit

					j <<- j + 1L
					orders[[j]] <<- unlist(stack[[pos]])

					pos <- pos - 1L # pop off of stack
				} else { # descending tree
					visit[pos] <- TRUE
					p <- pos
					for (i in seq_along(stack[[p]])) {
						if (!is.leaf(stack[[p]][[i]])) {
							# push subtree onto stack
							pos <- pos + 1L
							stack[[pos]] <- stack[[p]][[i]]
						}
					}
				}
			}
		}

		# record the alignment order in the original tree
		j <- 0L
		orders <- vector("list", l - 1L)
		.order(guideTree)

		.compare <- function(dend) {
			# initialize a stack of maximum length (l)
			stack <- vector("list", l)
			visit <- logical(l) # node already visited
			parent <- integer(l) # index of parent node
			index <- integer(l) # index in parent node
			found <- logical(l) # exact subtree found
			pos <- 1L # current position in the stack
			stack[[pos]] <- dend
			while (pos > 0L) { # more nodes to visit
				if (visit[pos]) { # ascending tree
					visit[pos] <- FALSE # reset visit

					# replace self in parent
					if (parent[pos] > 0)
						stack[[parent[pos]]][[index[pos]]] <- stack[[pos]]
					pos <- pos - 1L # pop off of stack
				} else { # descending tree
					if (found[pos]) {
						attr(stack[[pos]], "inherit") <- FALSE
						if (verbose)
							nsteps <<- nsteps - 1L
					} else {
						o <- unlist(stack[[pos]])
						j <<- j + 1L
						orders[[j]] <<- o

						w <- which(ls==length(o))
						for (i in seq_along(w)) {
							if (all(orders_prev[[w[i]]]==o)) {
								found[pos] <- TRUE
								break
							}
						}

						if (found[pos])
							attr(stack[[pos]], "inherit") <- TRUE
					}
					visit[pos] <- TRUE
					p <- pos
					for (i in seq_along(stack[[p]])) {
						if (!is.leaf(stack[[p]][[i]])) {
							# push subtree onto stack
							pos <- pos + 1L
							stack[[pos]] <- stack[[p]][[i]]
							parent[[pos]] <- p
							index[[pos]] <- i
							found[pos] <- found[p]
						}
					}
				}
			}
			return(stack[[1L]])
		}
	}

	for (it in seq_len(iterations)) {
		seqs_prev <- seqs
		ls <- lengths(orders)
		orders_prev <- orders

		if (verbose) {
			setTxtProgressBar(pBar, 100)
			close(pBar)
			cat("\n")
			time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			if (iterations > 1)
				 cat("\nIteration ",
					it,
					" of ",
					iterations,
					":\n",
					sep="")
		}

		if (type==2L &&
			treeLength >= minTreeLength &&
			useStructures) {
			structures <- .RNAStructures(seqs, weights)
		} else if (verbose) {
			cat("\n")
		}
		minTreeLength <- levels[8]

		if (verbose) {
			cat("Determining distance matrix based on alignment:\n")
			flush.console()
		}

		d <- DistanceMatrix(seqs,
			type="dist",
			verbose=verbose,
			processors=processors,
			includeTerminalGaps=TRUE)

		if (verbose) {
			cat("Reclustering into groups by similarity:\n")
			flush.console()
		}

		orgTree <- guideTree
		dimnames(d) <- NULL
		suppressWarnings(guideTree <- IdClusters(d,
			method="UPGMA",
			type="dendrogram",
			collapse=0,
			verbose=verbose,
			processors=processors))

		if (verbose) {
			time.1 <- Sys.time()
			cat("Realigning Sequences:\n")
			flush.console()
			pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1), max=100)
			before <- steps <- 0L
			nsteps <- l - 1L
		}

		j <- 0L
		orders <- vector("list", l - 1L)
		guideTree <- .reorder(guideTree)
		guideTree <- .compare(guideTree)

		seqs <- myXStringSet
		weights <- heights <- numeric(l)
		treeLength <- .align(guideTree)

		if (it < iterations && all(seqs==seqs_prev))
			break
	}

	myXStringSet <- seqs
	w <- which(weights <= 0)
	if (length(w) > 0)
		weights[w] <- 1
	weights <- weights/mean(weights)

	if (verbose) {
		setTxtProgressBar(pBar, 100)
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))

		if (iterations > 0 && it < iterations) {
			cat("\nAlignment converged - skipping remaining",
				ifelse(iterations - it > 1,
					"iterations.\n",
					"iteration.\n"))
		}
	}

	if (refinements > 0) {
		if (type==3L) {
			functionCall <- "colScoresAA"
		} else {
			functionCall <- "colScores"
		}

		GO <- gapOpeningMax/2 # applied at both ends
		colScores <- function(seqs, structs, weights) {
			scores <- .Call(functionCall,
				seqs,
				sM,
				GO,
				gapExtensionMax,
				weights,
				structs,
				structureMatrix)
			return(sum(scores))
		}

		# define trustworthy groups
		if (type==3L) { # AAStringSet
			cutoff <- ifelse(iterations > 0,
				levels[2]/2, # (fraction identical)/2
				levels[1]/2) # (fraction ordered k-mers)/2
		} else { # DNAStringSet or RNAStringSet
			cutoff <- ifelse(iterations > 0,
				levels[4]/2, # (fraction identical)/2
				levels[3]/2) # (fraction ordered k-mers)/2
		}
		guideTree <- cut(guideTree, cutoff)$lower

		# refinement
		n <- length(guideTree)
		if (n > 2) { # more than 2 groups
			if (type==2L &&
				treeLength >= minTreeLength &&
				(iterations==0 || (iterations > 0 && !all(seqs==seqs_prev))) &&
				useStructures) {
				structures <- .RNAStructures(seqs, weights)
			} else if (verbose) {
				cat("\n")
			}

			if (verbose) {
				time.1 <- Sys.time()
				cat("Refining the alignment:\n")
				flush.console()
				pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
			}

			score <- colScores(myXStringSet, structures, weights)
			vec <- seq_along(myXStringSet)

			for (ref in seq_len(refinements)) {
				org_score <- score
				count <- 0L
				for (i in seq_len(n)) {
					x <- unlist(guideTree[[i]])
					y <- vec[-x]
					o <- c(x, y)

					pattern <- .subset(myXStringSet, x)
					pattern <- .Call("removeCommonGaps",
						pattern,
						type,
						processors,
						PACKAGE="DECIPHER")
					subject <- .subset(myXStringSet, y)
					subject <- .Call("removeCommonGaps",
						subject,
						type,
						processors,
						PACKAGE="DECIPHER")

					p.weight <- weights[x]
					p.weight <- p.weight/mean(p.weight)
					s.weight <- weights[y]
					s.weight <- s.weight/mean(s.weight)

					if (subM) {
						if (useStructures) {
							if (is.null(structures)) {
								temp <- do.call(AlignProfiles,
									args=c(list(pattern=pattern,
											subject=subject,
											p.weight=p.weight,
											s.weight=s.weight,
											substitutionMatrix=sM,
											processors=processors,
											gapOpening=gapOpeningMax,
											gapExtension=gapExtensionMax),
										args))
							} else {
								temp <- do.call(AlignProfiles,
									args=c(list(pattern=pattern,
											subject=subject,
											p.weight=p.weight,
											s.weight=s.weight,
											p.struct=structures[x],
											s.struct=structures[y],
											substitutionMatrix=sM,
											processors=processors,
											gapOpening=gapOpeningMax,
											gapExtension=gapExtensionMax),
										args))
							}
						} else {
							temp <- do.call(AlignProfiles,
								args=c(list(pattern=pattern,
										subject=subject,
										p.weight=p.weight,
										s.weight=s.weight,
										substitutionMatrix=sM,
										processors=processors,
										gapOpening=gapOpeningMax,
										gapExtension=gapExtensionMax),
									args))
						}
					} else {
						if (useStructures) {
							temp <- do.call(AlignProfiles,
								args=c(list(pattern=pattern,
										subject=subject,
										p.weight=p.weight,
										s.weight=s.weight,
										p.struct=structures[x],
										s.struct=structures[y],
										processors=processors,
										gapOpening=gapOpeningMax,
										gapExtension=gapExtensionMax),
									args))
						} else {
							temp <- do.call(AlignProfiles,
								args=c(list(pattern=pattern,
										subject=subject,
										p.weight=p.weight,
										s.weight=s.weight,
										processors=processors,
										gapOpening=gapOpeningMax,
										gapExtension=gapExtensionMax),
									args))
						}
					}

					if (useStructures) {
						temp_score <- colScores(temp, structures[o], weights[o])
					} else {
						temp_score <- colScores(temp, NULL, weights[o])
					}

					if (temp_score > score) {
						score <- temp_score
						myXStringSet <- .subset(temp, order(o))

						count <- count + 1L
						if (count %% levels[5] ||
							l < levels[6])
							next # refine every nth change

						if (subM) {
							temp <- FUN(myXStringSet,
								substitutionMatrix=sM,
								weight=weights,
								processors=processors)
						} else {
							temp <- FUN(myXStringSet,
								weight=weights,
								processors=processors)
						}
						temp_score <- colScores(temp, structures, weights)

						if (temp_score > score) {
							score <- temp_score
							myXStringSet <- temp
						}
					}

					if (verbose)
						setTxtProgressBar(pBar,
							(i + n*(ref - 1))/(n*refinements))
				}

				if (org_score==score) # no changes
					break
			}

			if (verbose) {
				setTxtProgressBar(pBar, 1)
				close(pBar)
				cat("\n")
				time.2 <- Sys.time()
				print(round(difftime(time.2,
					time.1,
					units='secs'),
					digits=2))
				cat("\n")

				if (refinements > 0 && ref < refinements)
					cat("Alignment converged - skipping remaining",
						ifelse(refinements - ref > 1,
							"refinements.\n\n",
							"refinement.\n\n"))
			}
		} else if (verbose) {
			cat("\n")
		}
	} else if (verbose) {
		cat("\n")
	}

	if (l >= levels[6]) {
		# apply a final adjustment
		if (subM) {
			myXStringSet <- FUN(myXStringSet,
				substitutionMatrix=sM,
				weight=weights,
				processors=processors)
		} else {
			myXStringSet <- FUN(myXStringSet,
				weight=weights,
				processors=processors)
		}
	}

	names(myXStringSet) <- ns

	return(myXStringSet)
}
