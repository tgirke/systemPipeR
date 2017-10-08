######################################################################################
## Functions to convert R Markdown vignettes to paged markdown for Jekyll Doc Theme ##
######################################################################################
## Author: Thomas Girke
## Last update: Mar 6, 2016

md2Jekyll <- function(mdfile="Rbasics.knit.md", sidebartitle=NULL, sidebarpos, outfilebasename=NULL, outpath="./", sidebar_url_path="./", fenced2highlight=TRUE, image_dir=NULL) {
    ## (1) Import md file 
    md <- readLines(mdfile)

    ## Remove lines/patterns that are not appropriate for the Jekyll Doc Theme
    md <- md[!grepl("\\[Back to Table of Contents\\]", md)]
    md <- gsub("</br>", "", md) # Removes orphan breaks

    ## Convert underline-section tags to comment-section tags (some *.Rmd/*.md files use this format)
    sectionunderlineindex <- grepl("^==={1,}", md)
    subsectionunderlineindex <- grepl("^---{1,}", md)
    md[which(sectionunderlineindex) - 1] <- paste0("# ", md[which(sectionunderlineindex) - 1])
    md[which(subsectionunderlineindex) - 1] <- paste0("## ", md[which(subsectionunderlineindex) - 1])
    md <- md[!(sectionunderlineindex | subsectionunderlineindex)]

    ## (2) Parse specific entries in front matter
    mymaindoctitle <- md[1:20][grepl("^title:", tolower(md[1:20]))]
    mymaindoctitle <- gsub("(^.*? )|_|\\*|\"", "", mymaindoctitle)
    mymaindoctitle <- gsub(":", " -", mymaindoctitle)
    myauthor <- md[1:20][grepl("^author:", tolower(md[1:20]))]
    myauthor <- gsub("(^.*? )|_|\\*|\"", "", myauthor)
    mydate <- md[1:20][grepl("^date:", tolower(md[1:20]))]
    mydate <- gsub("(^.*? )|_|\\*|\"", "", mydate)
    mybibliography <- gsub("(^.* )|( .*)", "", md[1:20][grepl("bibliography:", md[1:20])])
    mypath <- gsub("(^.*\\/).*", "\\1", mdfile) # Note, this expects *.bib in the same directory as *.md
    if(grepl("\\/$", mypath)) mybibliography <- paste0(mypath, mybibliography)
    
    ## (3) Reformat citations and generate bibliography
    if(file.exists(mybibliography)) {
        biblist <- renderBib(x=md, bibtex=mybibliography)
        md <- biblist[[1]] # Returns text with reformatted citations
        bibliography_section <- biblist[[2]] # Returns references for bibliography section
    } else {
        md <- md
        bibliography_section <- NULL
    }

    ## (4) Optionally, convert fenced backtick code tags to 'highlight' tags
    ## Get code chunk ranges to protect them from being split if they contain #s at beginning of lines
    chunkpos <- grep("^```", md)
    ma <- matrix(chunkpos, length(chunkpos)/2, 2, byrow=TRUE)
    codechunk_ranges <- unlist(sapply(seq_along(ma[,1]), function(x) ma[x,1]:ma[x,2]))
    
    if(fenced2highlight==TRUE) {
        if((length(chunkpos) %% 2) != 0) stop("Odd number of chunk tags.")
        if(length(chunkpos) != 0) {
            ## Process code chunk positions
            codestartpos <- grep("^```\\w{1,}", md)
            languagetag <- gsub("^`{3,}", "", md[codestartpos])
            languagetag <- paste0("{% highlight ", languagetag, " %}")
            md[codestartpos] <- languagetag
            codeendpos <- chunkpos[which(chunkpos %in% codestartpos) + 1]
            md[codeendpos] <- "{% endhighlight %}"

            ## Process output chunk positions
            outputpos <- chunkpos[!chunkpos %in% c(codestartpos, codeendpos)]
            if(length(outputpos) != 0) {
                outputstartpos <- matrix(outputpos, (length(outputpos)/2), 2, byrow=TRUE)[,1]
                outputendpos <- matrix(outputpos, (length(outputpos)/2), 2, byrow=TRUE)[,2]
                md[outputstartpos] <- "{% highlight txt %}" # note, just highlight without language tag gives error in Jekyll
                md[outputendpos] <- "{% endhighlight %}"
            }
        }
    }
   
    ## (5) Collect all images in one directory and adjust image file paths in md file accordingly
    ## Generate new image_dir
    if(is.null(image_dir)) {
        image_dir <- gsub(".knit.md$", "_images", mdfile)
        image_dir <- gsub("^.*/", "", image_dir)
        ## If image directory *_images does not exist then use *_files instead
        if(!file.exists(paste0(mypath, image_dir))) {
            image_dir <- gsub(".knit.md$", "_files", mdfile)
            image_dir <- gsub("^.*/", "", image_dir)
        }
    } else {
        image_dir <- image_dir
    }
    image_dir <- paste0(outpath, "/", image_dir)
    image_dir <- gsub("/{1,}", "/", image_dir)
    dir.create(image_dir, showWarnings=FALSE)
    ## Get html image file paths
    htmlimgindex <- grepl("src=\".*\"", md)
    htmlimgindex <- htmlimgindex & !grepl("iframe {1,} src=", md) # To exclude, iframe tags
    if(sum(htmlimgindex) > 0) {
        htmlimgtag <- md[htmlimgindex]
        htmlimgpath <- gsub("^.*?src=\"(.*?)\".*", "\\1", htmlimgtag)
        ## Correct html image paths in md accordingly
        image_dir2 <- gsub("^.*/", "", image_dir) # Note, image path is relative in html source
        newhtmlimgpath <- paste0(image_dir2, "/", gsub("^.*/", "", htmlimgpath))
        newhtmlimgpath <- paste0(gsub("(^.*?src=\").*", "\\1", htmlimgtag), 
                                 "./pages/mydoc/", newhtmlimgpath, 
                                 gsub("^.*?src=\".*?(\".*)", "\\1", htmlimgtag))
        md[htmlimgindex] <- newhtmlimgpath
        ## Copy image files into image_dir
        file.copy(from=htmlimgpath[file.exists(htmlimgpath)], to=image_dir, overwrite=TRUE)
    } 
    
    ## Get md image file paths
    mdimgindex <- grepl("\\!\\[.*{0,}\\]\\(.*\\)", md)
    if(sum(mdimgindex) > 0) {
        mdimgtag <- md[mdimgindex]
        mdimgpath <- gsub("<!.*>", "", mdimgtag)  
        mdimgpath <- gsub("\\!\\[.*{0,}\\]\\((.*)\\)\\\\{0,}", "\\1", mdimgpath)
        ## Correct image paths in md accordingly
        image_dir2 <- gsub("^.*/", "", image_dir) # Note, image path is relative in html source
        newmdimgpath <- paste0(image_dir2, "/", gsub("^.*/", "", mdimgpath))
        newmdimgpath <- paste0(gsub("(\\!\\[.*{0,}\\]\\().*", "\\1", mdimgtag),
                               "./pages/mydoc/", newmdimgpath,
                               ")")
        md[mdimgindex] <- newmdimgpath
        ## Copy image files into image_dir
        file.copy(from=mdimgpath[file.exists(mdimgpath)], to=image_dir, overwrite=TRUE)
    }
    cat(paste("Saved image files to directory:", image_dir, "\n")) 

    ## (6) Split on main sections (^# )
    ## Identify split positions
    splitpos <- grep("^# {1,}", md)
    splitpos <- splitpos[!splitpos %in% codechunk_ranges]
    if(length(splitpos) != 0) {
        splitdist <- c(splitpos, length(md)+1) - c(1, splitpos) 
        mdlist <- split(md, factor(rep(c(1, splitpos), splitdist)))
    
        ## Get alternative format info and remove front matter of R Markdown 
        altformatpath <- paste0("_alternativeFormatlinks.md")
        if(file.exists(altformatpath)) {
            altformats <- readLines(altformatpath)
        } else {
            altformats <- ""
        }
        mdlist <- mdlist[-1] # Removes R Markdown front matter
    
        ## Write sections stored in list components to separate files
        titles <- sapply(seq_along(mdlist), function(x) mdlist[[x]][1])
        mytitles <- gsub("# {1,}", "", titles)
        mytitles <- paste0(1:length(mytitles), ". ", mytitles)
    } else {
        stop("mdfile is expected to contain at least one section tag.")
    }
    ## Add references if section exists
    reflistpos <- which(sapply(seq_along(mdlist), function(x) mdlist[[x]][1]) %in% "# References")
    mdlist[[reflistpos]] <-  mdlist[[reflistpos]][nchar(mdlist[[reflistpos]])>0]
    if(length(reflistpos)==1 & length(mdlist[[reflistpos]])==1) {
        mdlist[[reflistpos]] <- c("# References", " ", paste0(seq_along(bibliography_section), ". ", bibliography_section))
    }
    
    ## (7) Add Jekyll Doc front matter to each list component
    for(i in seq_along(mdlist)) { 
        frontmatter <- c(starttag="---", 
                         title=paste0("title: ", mytitles[i]), 
                         last_updated=paste0("last_updated: ", date()), 
                         sidebar="sidebar: mydoc_sidebar",
                         permalink="permalink: ",
                         endtag="---")
        mdlist[[i]] <- c(as.character(frontmatter), mdlist[[i]][-1])
    }
    
    ## Special handling of first page
    mdlist[[1]][2] <- paste0("title: ", mymaindoctitle, " <br> <br> ", mytitles[1]) # Uses in front matter main title of source document
    mdlist[[1]] <- c(mdlist[[1]][1:6], myauthor, "", mydate, "", altformats, "", mdlist[[1]][7:length(mdlist[[1]])])    

    ## (8) Write sections to files named after input files with numbers appended
    filenumbers <- sprintf(paste0("%0", as.character(nchar(length(mdlist))), "d"), seq_along(mdlist))
    ## If no outfilebasename is provided the use default where they are named after input file with ".knit.md" stripped off
    if(is.null(outfilebasename)) {
        outfilebasename <- gsub(".knit.md$", "", mdfile)
    } else {
        outfilebasename <- outfilebasename    
    }
    filenames <- paste0(gsub("/$", "", outpath), "/mydoc_", outfilebasename, "_", filenumbers, ".md")
    filenames <- gsub("/{1,}", "/", filenames)

    ## Add permalink info to front matter 
    permalink <- paste0(gsub("(^.*/)|(md$)", "", filenames), "html")
    for(i in seq_along(mdlist)) mdlist[[i]][5] <- paste0(mdlist[[i]][5], permalink[i])
   
    ## Next/previous page image links
    nextpageurl <- paste0('<a href=\"', permalink, '\"><img src=\"images/right_arrow.png\" alt="Next page."></a>')
    nextpageurl <- nextpageurl[c(2:length(nextpageurl), 1)]
    nextpageurl <- paste0(nextpageurl, "</center>")
    previouspageurl <- paste0('<a href=\"', permalink, '\"><img src=\"images/left_arrow.png\" alt="Previous page."></a>')
    previouspageurl <- previouspageurl[c(1,1:(length(previouspageurl)-1))]
    previouspageurl <- paste0("<br><br><center>", previouspageurl, "Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page")
    for(i in seq_along(mdlist)) mdlist[[i]] <- c(mdlist[[i]], previouspageurl[i], nextpageurl[i])

    ## Write subpages to files
    for(i in seq_along(mdlist)) {
        writeLines(mdlist[[i]], filenames[i])
        cat(paste("Created file:", filenames[i]), "\n") 
    }
    
    ## (9) Copy R Markdown *.html and *.pdf to pages/mydoc. This makes it easier to render files on GitHub
    rmdhtmlfile <- paste0(outpath, "/", gsub(".knit.md", ".html", mdfile))
    file.copy(from=basename(rmdhtmlfile), to=rmdhtmlfile, overwrite = TRUE)
    cat(paste("Copied R Markdown HTML file to:", rmdhtmlfile), "\n")
    rmdpdffile <- paste0(outpath, "/", gsub(".knit.md", ".pdf", mdfile))
    file.copy(from=basename(rmdpdffile), to=rmdpdffile, overwrite = TRUE)
    cat(paste("Copied R Markdown PDF file to:", rmdpdffile), "\n")

    ## (10) Register new files in sidebar (_data/sidebars/mydoc_sidebar.yml)
    sb <- readLines("../../_data/sidebars/mydoc_sidebar.yml") 
    splitFct <- function(sb, pattern) {    
        splitpos <- grep(pattern, sb)
        if(length(splitpos) != 0) {
            splitdist <- c(splitpos, length(sb)+1) - c(1, splitpos) 
            sblist <- split(sb, factor(rep(c(1, splitpos), splitdist)))
            mynames <- gsub("^.*title: {1,}", "", sb[splitpos])
            mynames <- gsub(" {1,}", "_", mynames)
            names(sblist) <- c("header", mynames)
            return(sblist)
        } else {
            stop("mydoc_sidebar.yml is expected to contain at least one component.")
        }
    }
    ## Split on first level
    sblist <- splitFct(sb=sb, pattern="^ {2,2}- title: \\w{1,}")
    ## Split on second level resulting in nested list
    header <- sblist[1]
    sblist <- sblist[-1]
    sblist <- sapply(names(sblist), function(x) splitFct(sb= sblist[[x]], pattern="^ {4,4}- title: \\w{1,}"), simplify=FALSE)
    if(is.null(sidebartitle)) {
        sidebartitle <- gsub("(^.*?) {1,}.*", "\\1", mymaindoctitle)
    } else {
        sidebartitle <- sidebartitle
    }
    sidebartitle <- gsub("(^ {1,})|( {1,}$)", "", sidebartitle)
    sblist <- sblist[!names(sblist) %in% sidebartitle] # Removes existing section entry
    sblist <- c(header, sblist)        
    ## Construct new sidebar entries
    myurls <- paste0("", basename(filenames))
    myurls <- gsub(".md$", ".html", myurls)
    sectionheader <- c(paste0("  - title: ", sidebartitle),
                       "    output: web, pdf",
                       "    folderitems:",
                       "")
    subsections <- c("    - title: ",
                     "      url: /",
                     "      output: web",
                     "")
    subsectionlist <- lapply(seq_along(mytitles), function(x) subsections)
    for(i in seq_along(subsectionlist)) {
        subsectionlist[[i]][1] <- paste0(subsectionlist[[i]][1], mytitles[i])
        subsectionlist[[i]][2] <- paste0(subsectionlist[[i]][2], myurls[i])
    }
    names(subsectionlist) <- mytitles
    sectionlist <- list(c(list(header=sectionheader), subsectionlist))
    names(sectionlist) <- sidebartitle
    sblist <- c(sblist[1:(sidebarpos)], sectionlist, sblist[(sidebarpos+1):length(sblist)])
    sidebarfile <- paste0(sidebar_url_path, "/", "mydoc_sidebar.yml")
    sidebarfile <- gsub("/{1,}", "/", sidebarfile)
    writeLines(unlist(sblist), sidebarfile)
    cat(paste("Created file", sidebarfile), "\n")
}

## Usage:
# setwd("~/Dropbox/Websites/manuals/vignettes/Rbasics")
# source("../md2jekyll.R")
# md2Jekyll(mdfile="bioassayR.knit.md", sidebartitle=NULL, sidebarpos=12, outfilebasename=NULL, outpath="../../mydoc", sidebar_url_path="../../_data/mydoc/", fenced2highlight=TRUE, image_dir=NULL)

##############################################################
## Functions to reformat citatons and generate bibliography ##
##############################################################
## Import BibTeX file
## Import into list
bibtexImp <- function(file) {
	bibtexv <- readLines(file)
	bibtexv <- gsub("^(\\s|\\t){1,}", "", bibtexv) # Removes spaces/tabs at beginning of lines
	bibtexv <- bibtexv[!grepl("^(\\s|\\t){0,}%", bibtexv)] # removes all comment lines
    bibtexv <- bibtexv[!grepl("^(\\s|\\t){0,}$", bibtexv)] # removes all empty fields
    ## Concatenate entries wrapping over several lines to single one
    splitAt2 <- function(x, pos) {
        out <- list()
        pos2 <- c(1, pos, length(x)+1)
        for (i in seq_along(pos2[-1])) {
            out[[i]] <- x[pos2[i]:(pos2[i+1]-1)]
        }
        return(out)
    }
    index <- grep("(^@)|(^.* {0,}=)|(^})", bibtexv)
    indexlist <- splitAt2(seq_along(bibtexv), index[-1])
    bibtexv <- sapply(seq_along(indexlist), function(x) paste(bibtexv[indexlist[[x]]], collapse=" "))
	## Split into list 
    start <- regexpr("^@", bibtexv, perl=T)
	start <- which(start==1)
	end <- regexpr("^}", bibtexv, perl=T)
	end <- which(end==1)
	if(length(start)!=length(end)) cat("Error in BibTex source: number of '^@' differes from number of '^}' \n") 
	index <- data.frame(start=start, end=end)
	myfactor <- rep(1:length(index[,1]), index[,2]-index[,1]+1)
	bibtexList <- split(bibtexv, myfactor) # generates BibTex list
	names(bibtexList) <- unlist(lapply(names(bibtexList), function(x) gsub(".*\\{|,$", "", bibtexList[[x]][1])))
	bibtexList
}

## Convert to character or list formats
formatBibtex <- function(bibTexSubDB=bibTexSubDB, printurl=TRUE, format="list") {
	bibTexSubDB <- lapply(bibTexSubDB, function(x) gsub("(= {0,})\\{{1,}|\\}{1,} {0,}(,)$", "\\1\"\\2", x)) # replaces all start and end braces by quotation marks
	authors <- lapply(names(bibTexSubDB), function(x) gsub(".* \\\"(.*)\".*", "\\1", bibTexSubDB[[x]][grep("^author", bibTexSubDB[[x]])]))
	names(authors) <- 1:length(authors)
	authors <- lapply(names(authors), function(x) gsub(" {1,}and {1,}", ", ", authors[x]))
	authors <- lapply(authors, function(x) x[1]) # necessary to turn empty list components into vector with "NA"
	year <- lapply(names(bibTexSubDB), function(x) gsub(".*?(\\d{2,}).*", "\\1", bibTexSubDB[[x]][grep("^year", bibTexSubDB[[x]])]))
	year <- lapply(year, function(x) x[1]) # necessary to turn empty list components into vector with "NA"
	title <- lapply(names(bibTexSubDB), function(x) gsub(".*\"(.*)\".*", "\\1", bibTexSubDB[[x]][grep("^title", bibTexSubDB[[x]])]))
	names(title) <- 1:length(title)
	title <- lapply(names(title), function(x) gsub("\\{|\\}|,$", "", title[x]))
	title <- lapply(title, function(x) x[1]) # necessary to turn empty list components into vector with "NA"
	journal <- lapply(names(bibTexSubDB), function(x) gsub(".*\"(.*)\".*", "\\1", bibTexSubDB[[x]][grep("^journal", bibTexSubDB[[x]])]))
	journal <- lapply(journal, function(x) x[1]) # necessary to turn empty list components into vector with "NA"
	volume <- lapply(names(bibTexSubDB), function(x) gsub(".*?(\\d{1,}).*", "\\1", bibTexSubDB[[x]][grep("^volume", bibTexSubDB[[x]])]))
	volume <- lapply(volume, function(x) x[1]) # necessary to turn empty list components into vector with "NA"
	pages <- lapply(names(bibTexSubDB), function(x) gsub(".*\"(.*)\".*", "\\1", bibTexSubDB[[x]][grep("^pages", bibTexSubDB[[x]])]))
	pages <- lapply(pages, function(x) x[1]) # necessary to turn empty list components into vector with "NA"
	myurl <- lapply(names(bibTexSubDB), function(x) gsub(".*\"(.*)\".*", "\\1", bibTexSubDB[[x]][grep("^url", bibTexSubDB[[x]])]))
	myurl <- lapply(myurl, function(x) x[1]) # necessary to turn empty list components into vector with "NA"
    myurl <- lapply(seq_along(myurl), function(x) paste0("[URL](", myurl[[x]], ")"))
    if(format=="character") {
        if(printurl==TRUE) {
	    	refv <- paste(unlist(authors), " (", unlist(year), ") ", unlist(title), ". ", unlist(journal), ", ", unlist(volume), ": ", unlist(pages), "; ", unlist(myurl), sep="")	
	    } else {
	    	refv <- paste(unlist(authors), " (", unlist(year), ") ", unlist(title), ". ", unlist(journal), ", ", unlist(volume), ": ", unlist(pages), sep="")	
	    }
    names(refv) <- names(bibTexSubDB)
    return(refv)
    }
    if(format=="list") {
        if(printurl==TRUE) {
	        m <- cbind(authors, year, title, journal, volume, pages, myurl)
	        m <- sapply(seq_along(m[,1]), function(x) m[x,], simplify=FALSE)
        } else {
	        m <- cbind(authors, year, title, journal, volume, pages)
	        m <- sapply(seq_along(m[,1]), function(x) m[x,], simplify=FALSE)
        }
    names(m) <- names(bibTexSubDB)
    return(m)
    }
}
## Usage:
## Import BibTeX file into list
# bibTexDB <- bibtexImp(file="bibtex.bib")
# bibTexDB <- bibtexImp(file="04_longevityTools_eQTL/bibtex.bib")
## Reformat BibTeX list
# formatBibtex(bibTexSubDB=bibTexDB, printurl=TRUE, format="character")
# formatBibtex(bibTexSubDB=bibTexDB, printurl=TRUE, format="list")

## Render citations in text and return bibliography section
renderBib <- function(x, bibtex="bibtex.bib") {
    ## Import md file
    ## Import BibTeX
    bibtexDB <- bibtexImp(bibtex)
    ## Collapse text to one line
    xc <- paste(x, collapse=" ")
    ## Locate citations in text
    pos <- gregexpr("(\\[@.*?\\])|(\\[-@.*?\\])", xc)
    mpos <- as.numeric(pos[[1]]) 
    mlength <- attributes(pos[[1]])$match.length
    ## Extract citations
    citation <- substring(xc, mpos, mpos + (mlength - 1))
    citation <- gsub("(\\[)|(\\])|(^ {1,})|( {1,}$)", "", citation)    
    citation <- gsub("(, )|(; )", " ", citation)
    citation <- gsub(" {0,}, {0,}", "_", citation) # Support for split on comma 
    citation <- gsub(" {1,}", "___", citation) # Support for split on space
    citation <- strsplit(citation, "___")
    ## Creat bibliography section for citations with bibtexDB
    bibliography <- formatBibtex(bibTexSubDB=bibtexDB, printurl=TRUE, format="character")
    ids <- unique(gsub("(@)|(-@)", "", unlist(citation)))    
    bibliographysub <- bibliography[ids]
    names(bibliographysub)[is.na(names(bibliographysub))] <- ids[!ids %in% names(bibliographysub)] # Maintains names of orphan citations 
    bibliographysub <- bibliographysub[order(bibliographysub)]
    ## Render citations for insertion into source text
    bibtexlist <- formatBibtex(bibTexSubDB=bibtexDB, printurl=FALSE, format="list")
    genCitations <- function(x, bibtexlist) {
        .genCitations <- function(x, bibtexlist) {
            x2 <- x
            if(grepl("^-@", x)) {
                suppressauthor <- TRUE 
            } else {
                suppressauthor <- FALSE
            }    
            x <- gsub("(@)|(-@)", "", x)
            if(!is.null(bibtexlist[[x]])) {
                myauthors <- unlist(strsplit(bibtexlist[[x]]$authors, ", "))
                if(length(myauthors) > 2) citationstring <- paste(myauthors[1], "et al.,", bibtexlist[[x]]$year)
                if(length(myauthors) <= 2) citationstring <- paste(myauthors[1], ",", bibtexlist[[x]]$year)
                if(suppressauthor==TRUE) citationstring <- bibtexlist[[x]]$year
            } else {
                citationstring <- x2
            }
            return(citationstring)
        }
        citationstring <- sapply(x, function(x) .genCitations(x, bibtexlist))
        return(paste(citationstring, collapse="; "))
    }
    citation <- unique(unlist(citation))
    citationstrings <- sapply(citation, genCitations, bibtexlist)
    citationstrings <- citationstrings[order(names(citationstrings))] # Important to process suppress author citations first
    ## Insert citations into source text. Note, the following simply uses a pattern-based 
    ## replacement to insert the citations. Using a position-based replacement would
    ## be more flexible with respect to formatting. This could be implemented later
    ## The index required for this is available in 'mpos' and 'mlength'.
    if(length(citationstrings)!=0) {
        for(i in seq_along(citationstrings)) {
            x <- gsub(paste0("\\[ {0,}", names(citationstrings[i]), " {0,}\\]"), paste0("(", citationstrings[i], ")"), x)
            x <- gsub(paste0("\\[ {0,}", names(citationstrings[i])), paste0("(", citationstrings[i]), x)
            x <- gsub(paste0(names(citationstrings[i]), " {0,}\\]"), paste0(citationstrings[i], ")"), x)
            x <- gsub(names(citationstrings[i]), citationstrings[i], x)
        }
    } else {
        x <- x
    }
    ## Return results as list
    mylist <- list(md=x, bibliography=bibliographysub)
    return(mylist)
}

## Usage:
# x <- c("Some text", "with [@Peters2015-fc, @12234] citation [@Sood2015-pb]", "and anoter [-@Sood2015-pb]", "some complication [@Sood2015-pb, @Peters2015-fc]")
# x <- readLines("05_longevityTools_eDRUG/longevityTools_eDRUG.knit.md") 
# renderBib(x=x, bibtex="05_longevityTools_eDRUG/bibtex.bib")


## Run from command-line with arguments
myargs <- commandArgs()
md2Jekyll(mdfile=myargs[6], sidebartitle=NULL, sidebarpos=as.numeric(myargs[7]), outfilebasename=NULL, outpath="../../pages/mydoc", sidebar_url_path="../../_data/sidebars/", fenced2highlight=FALSE, image_dir=NULL)
# $ Rscript ../md2jekyll.R bioassayR.knit.md 8


