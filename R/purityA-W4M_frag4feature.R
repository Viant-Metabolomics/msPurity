#' @title Using a purityA object, link MS/MS datas to XCMS features, respecting the grouping step
#'
#' @description
#'
#' **General**
#'
#' This function is a deviate function from frag4feature. It makes the same thing with dividing each file into its group. So it 
#' assigns fragmentation spectra (MS/MS) stored within a purityA class object to grouped features within an XCMS xset object.
#'
#' XCMS calculates individual chromatographic peaks for each mzML file (saved in xset@@peaks), these are then grouped together
#' (using xcms.group). Ideally the mzML files that contain the MS/MS spectra also contain sufficient MS1 scans for XCMS to detect
#' MS1 chromatographic features. If this is the case, to determine if a MS2 spectra is to be linked to an XCMS grouped feature,
#' the associated acquisition time of the MS/MS event has to be within the retention time window defined for the individual peaks
#' associated for each file. The precursor m/z value also has to be within the user ppm tolerance to XCMS feature.
#'
#' See below for representation of the linking (the \*------\* represent a many-to-many relationship) e.g. 1 or more MS/MS events can be
#' linked to 1 or more individual feature and an individual XCMS feature can be linked to 1 or more grouped XCMS features
#'
#' * \[grouped XCMS feature - across files\] \*------\*  \[individual XCMS feature - per file\] \*------\*  \[MS/MS spectra\]
#'
#' Alternatively, if the "useGroup" argument is set to TRUE, the full width of the grouped peak (determined as the minimum rtmin
#' and maximum rtmax of the all associated individual peaks) will be used. This option should be used if the mzML file with
#' MS/MS has very limited MS1 data and so individual chromatographic peaks might not be detected with the mzML files containing the
#' MS/MS data. However, it should be noted this may lead to potential inaccurate linking.
#'
#' * \[grouped XCMS peaks\] \*------\* \[MS/MS spectra\]
#'
#'
#' **Example LC-MS/MS processing workflow**
#'
#' The purityA object can be used for further processing including linking the fragmentation spectra to XCMS features, averaging fragmentation, database creation and spectral matching (from the created database). See below for an example workflow
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.xcmsSet -> xcms.merge -> xcms.group -> xcms.retcor -> xcms.group -> (xset)
#'  * Fragmentation processing
#'    + (xset, pa, sampleMetadata) -> **W4M_frag4feature** -> filterFragSpectra -> averageAllFragSpectra -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#' **Additional notes**
#'
#' * If using only a single file, then grouping still needs to be performed within XCMS before frag4feature can be used.
#' * Fragmentation spectra below a certain precursor ion purity can be be removed (see plim argument).
#' * A SQLite database can be created directly here but the functionality has been deprecated and the createDatabase function should now be used
#' * Can experience some problems when using XCMS version < 3 and obiwarp retention time correction.
#'
#'
#' @aliases W4M_frag4feature
#'
#' @param pa object; purityA object
#' @param xset object; xcmsSet object derived from the same files as those used to create the purityA object
#' @param sampleMetadata character; table which contains filename link to a class/group (colnames should contains "" and "")
#' @param ppm numeric; ppm tolerance between precursor mz and XCMS feature mz
#' @param plim numeric; minimum purity of precursor to be included
#' @param intense boolean; If TRUE the most intense precursor will be used. If FALSE the precursor closest to the center of the isolation window will be used
#' @param useGroup boolean; Ignore individual peaks and just find matching fragmentation spectra within the (full) rtmin rtmax of each grouped feature
#' @param convert2RawRT boolean; If retention time correction has been used in XCMS set this to TRUE
#' @param create_db boolean; (Deprecated, to be removed - use createDatabase function) SQLite database will be created of the results
#' @param db_name character; (Deprecated, to be removed - use createDatabase function) If create_db is TRUE, a custom database name can be used, default is a time stamp
#' @param out_dir character; (Deprecated, to be removed - use createDatabase function) Path where database will be created
#' @param grp_peaklist dataframe; (Deprecated, to be removed - use createDatabase function) Can use any peak dataframe to add to databse. Still needs to be derived from the xset object though
#' @param use_group boolean; (Deprecated, to be removed - replaced with useGroup argument for style consistency)
#' @return Returns one or more purityA object (pa), depends of the number of classes, with the following slots populated:
#'
#' * pa@@grped_df: A dataframe of the grouped XCMS features linked to the associated fragmentation spectra precursor details is recorded here
#' * pa@@grped_ms2: A list of fragmentation spectra associated with each grouped XCMS feature is recorded here
#' * pa@@f4f_link_type: The linking method is recorded here (e.g. individual peaks or grouped - "useGroup=TRUE")
#'
#'
#' @examples
#'
#' sampleMetadata <- read.csv("my_sampleMetadata_file.csv", sep="\t", header=TRUE)
#' msmsPths <- list.files(system.file("extdata", "lcms",
#'                        "mzML", package="msPurityData"), full.names = TRUE,
#'                        pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths)
#' pa <- W4M_frag4feature(pa, xset,sampleMetadata)
#' @md
#' @export

#To add when we will build it in a package
#setMethod(f="W4M_frag4feature", signature="purityA",
#          definition= function(pa, xset, sampleMetadata=NULL, ppm=5, plim=NA, intense=TRUE, convert2RawRT=TRUE, useGroup=FALSE, create_db=FALSE,
#                                out_dir='.', db_name=NA, grp_peaklist=NA, use_group=NA){})
W4M_frag4feature <- function(pa, xset, sampleMetadata=NULL, ppm=5, plim=NA, intense=TRUE, convert2RawRT=TRUE, useGroup=FALSE, create_db=FALSE,
                                out_dir='.', db_name=NA, grp_peaklist=NA, use_group=NA){
  
    #Packages to load (or to add in NAMESPACE ?)
    pkgs <- c("R.utils","plyr")
    loadAndDisplayPackages(pkgs)

    #Need it to run the latest msPurity
    #sourceDirectory("/home/jsaintvanne/W4M/msPurity/R")

    #Works on xset and pa object before process them
    xset <- getxcmsSetObject(xset)
    xset <- egal_pa_xset_filenames(xset, pa)
    #Checking and building a good sampleMetadata containing only MSMS files
    if(is.null(sampleMetadata)){
	   sampleMetadataMSMS <- buildSamplemetadataFromXCMS(xset,pa)
    }else{
        #Verify if we have a file
        if(class(sampleMetadata) == "character"){
            sampleMetadataMSMS <- buildSamplemetadataFromFile(sampleMetadata,xset,pa)
        }else{
            sampleMetadataMSMS <- buildSamplemetadataFromTable(sampleMetadata,xset,pa)
        }
    }

    #Verify if sampleMetadataMSMS is null !
    if(nrow(sampleMetadataMSMS) == 0){
	   error_message <- "No file has msLevel = 2 in your sampleMetadata...\n"
	   cat(error_message)
	   stop(error_message)
    }
    #Order it by class names
    sampleMetadataMSMS <- sampleMetadataMSMS[order(sampleMetadataMSMS[,"class"]),]
    cat("\nsampleMetadata :\n\n")
    print(sampleMetadataMSMS)

    #Run MSMS file by MSMS file to find matches features
    outputdata <- lapply(unique(sampleMetadataMSMS$class),function(x){
    
        sameClass <- sampleMetadataMSMS[which(sampleMetadataMSMS$class == x),]
        fileOfSameClass <- as.character(sameClass[,"MSMS"])

        # 1 - Save class
        class <- as.character(sameClass[1,"class"])

        # 2 - Modify the XcmsSet object to have only peaks and groups from the good class
        xsettempo <- modifyXsetObject(xset,class)

        # 3 - Modify the pa object to have only the file we are working on!
        patempo <- modifyPaObject(pa,fileOfSameClass,class)

        # 4 - Find if retention time have been processed for this file (maybe it can be add somewhere in package ?)
        if(NA %in% match(xsettempo@rt$raw,xsettempo@rt$corrected)){
            convert2RawRT= TRUE
        }

        # 6 - Run f4f function from msPurity
        cat(paste("\n--------------- Run frag4feature in msPurity package ---------------\n"))
        cat(paste("Searching in",nrow(patempo@puritydf),"MS/MS datas...\n"))
        patempo <- frag4feature(pa=patempo, 
                       xset=xsettempo, 
                       ppm=ppm, 
                       plim=plim,
                       intense=intense,
                       convert2RawRT=convert2RawRT,
                       useGroup=useGroup,
                       create_db=create_db,
                       out_dir=out_dir,
                       db_name='alldata.sqlite',
                       grp_peaklist=grp_peaklist,
                       use_group=useGroup)
        cat("\nFrag4feature finish !\n\n")

        patempo@grped_df$filename <- sapply(patempo@grped_df$fileid, function(x) names(patempo@fileList)[as.integer(x)])

    return(patempo)
    })

    names(outputdata) <- as.character(unique(sampleMetadataMSMS$class))
    return(outputdata)
}


###################
#### FUNCTIONS ####
###################

#@author G. Le Corguille
# This function will
# - load the packages
# - display the sessionInfo
loadAndDisplayPackages <- function(pkgs) {
    for(pkg in pkgs) suppressPackageStartupMessages( stopifnot( library(pkg, quietly=TRUE, logical.return=TRUE, character.only=TRUE)))

    sessioninfo = sessionInfo()
    cat(sessioninfo$R.version$version.string,"\n")
    cat("Main packages:\n")
    for (pkg in names(sessioninfo$otherPkgs)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
    cat("Other loaded packages:\n")
    for (pkg in names(sessioninfo$loadedOnly)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
}

# This function retrieve a xset like object
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getxcmsSetObject <- function(xobject) {
    # XCMS 1.x
    if (class(xobject) == "xcmsSet")
        return (xobject)
    # XCMS 3.x
    if (class(xobject) == "XCMSnExp") {
        # Get the legacy xcmsSet object
        suppressWarnings(xset <- as(xobject, 'xcmsSet'))
        if (!is.null(xset@phenoData$sample_group))
            sampclass(xset) <- xset@phenoData$sample_group
        else
            sampclass(xset) <- "."
        return (xset)
    }
}

#This function is to have a correction between xset filenames and pa filenames
## xset filenames are absolute path of each file
## pas filenames are exactly what the user enter as input during purityA function
egal_pa_xset_filenames <- function(xset, pa){
    if(all(xset@filepaths != pa@fileList)){
        if(all(names(pa@fileList) == basename(xset@filepaths))){
            for(x in 1:length(xset@filepaths)) {
                xset@filepaths[x] <- pa@fileList[[x]]
            }
            return(xset)
        }
    }
    return(xset)
}

# 3 functions to have a sampleMetadata looks like this 
#(with no extension for files and just basename) :
# MSMS 	   class
#file1    class1
#file2    class2
#file3    class1
#First from xcmsSet object with all MSMS files possible
buildSamplemetadataFromXCMS <- function(xset, pa){
    if(!(is.null(xset))){
        MSMS <- NULL
  		class <- NULL
  		filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                       "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
        filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
  		for(i in 1:length(xset@filepaths)) {
            
            #files <- names(pa@fileList[match(files,gsub(filepattern,"",pa@fileList))])

            #Select only MSMS files
    		if(2 %in% unique(readMSData(xset@filepaths[i],mode="onDisk")@featureData@data$msLevel)) {
            	sampname <- gsub(filepattern, "",basename(xset@filepaths[i])) 
        		MSMS <- c(MSMS, sampname)
                #Verify the colname containing groups/classes
                if("sample_group" %in% colnames(xset@phenoData)){
                    class <- c(class,make.names(as.character(xset@phenoData[i,"sample_group"])))
                }else{
                    if("sampleGroups" %in% colnames(xset@phenoData)){
                        class <- c(class,make.names(as.character(xset@phenoData[i,"sampleGroups"])))
                    }else{
                        if("class" %in% colnames(xset@phenoData)){
                            class <- c(class,make.names(as.character(xset@phenoData[i,"class"])))
                        }else{
                            print("error to find groups/class")
                        }
                    }
                }
    		}
  		}
        class <- make.names(class) #Correction of spaces in class names
  		MSMSclassfile <- data.frame(MSMS,class)
  		return(MSMSclassfile)
	}else{
		error_message <- "No xset available in buildSamplemetadataFromXCMS\n"
		cat(error_message)
		stop(error_message)
		return(NULL)
    }	
}
#Second from file
buildSamplemetadataFromFile <- function(sampleMetadata,xset,pa){
    finaleMSMSclasses <- NULL
	if(!(is.null(sampleMetadata))){
  		#Read the CSV file
  		MSMSclasses <- read.csv(file=sampleMetadata, sep="\t", header=TRUE)
  		finaleMSMSclasses <- buildSamplemetadataFromTable(MSMSclasses, xset@filepaths,pa)
        return(finaleMSMSclasses)
	}else{
  		error_message <- "No sampleMetadata file enter\n"
  		cat(error_message)
		stop(error_message)
		return(NULL)
	}
}
#Third from table
buildSamplemetadataFromTable <- function(sampleMetadata, files,pa){
    finaleMSMSclasses <- NULL
    #remove .cdf, .mzXML filepattern
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                         "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
    #Verify the colnames
    if("MSMS" %in% colnames(sampleMetadata)){
        if("class" %in% colnames(sampleMetadata)){
            MSMS <- sampleMetadata[,"MSMS"]
            class <- sampleMetadata[,"class"]
        }else if("sample_group" %in% colnames(sampleMetadata)){
            MSMS <- sampleMetadata[,"MSMS"]
            class <- sampleMetadata[,"sample_group"]
        }else{
            error_message <- "Can't find any class or sample_group column from sampleMetadata ! \n"
            cat(error_message)
            stop(error_message)
            return(NULL)
        }
    }else if("sample_name" %in% colnames(sampleMetadata)){
        if("class" %in% colnames(sampleMetadata)){
            MSMS <- sampleMetadata[,"sample_name"]
            class <- sampleMetadata[,"class"]
        }else if("sample_group" %in% colnames(sampleMetadata)){
            MSMS <- sampleMetadata[,"sample_name"]
            class <- sampleMetadata[,"sample_group"]
        }else{
            error_message <- "Can't find any class or sample_group column from sampleMetadata ! \n"
            cat(error_message)
            stop(error_message)
            return(NULL)
        }
    }else if(NA %in% as.numeric(rownames(sampleMetadata))){ #Find rownames as non numerics (= filenames??)
        if("class" %in% colnames(sampleMetadata)){
            MSMS <- gsub(filepattern,"",basename(rownames(sampleMetadata)))
            class <- sampleMetadata[,"class"]
        }else if("sample_group" %in% colnames(sampleMetadata)){
            MSMS <- gsub(filepattern,"",basename(rownames(sampleMetadata)))
            class <- sampleMetadata[,"sample_group"]
        }else{
            error_message <- "Can't find any class or sample_group column from sampleMetadata ! \n"
            cat(error_message)
            stop(error_message)
            return(NULL)
        }
    }else{
        error_message <- "Can't find any file name ! Please add a column named \"sample_name\" with your filenames \n"
        cat(error_message)
        stop(error_message)
        return(NULL)
    }
            
    #Rebuild sampleMetadata
    class <- make.names(class) #Correction of spaces in class names
    #MSMS <- make.names(MSMS)
    sampleMetadata <- data.frame(MSMS,"class" = class)
    
    #Keep only MSMS files
    files <- names(pa@fileList[match(files,gsub(filepattern,"",pa@fileList))])
    for(i in 1:length(files)){
        sampname<-gsub(filepattern, "",basename(files[i]))  
        if(2 %in% unique(readMSData(files[i],mode="onDisk")@featureData@data$msLevel)){
            finaleMSMSclasses <- rbind(finaleMSMSclasses,sampleMetadata[which(sampleMetadata[,"MSMS"] == sampname),])
        }
    }

    return(finaleMSMSclasses)
}

#This function create a tempo xcmsSet object with only peaks from class we are looking at
modifyXsetObject <- function(xset, class){
    cat(paste("\n--------------- Modify xcmsSet object for \"",class,"\" class file(s) ---------------\n"))
    cat("STEP 1 : find good groups from our class\n")
    #STEP 1 : select groups in xset@groups for the good class and save ids of groups deleted
    #Select groups which contain peaks in the same class as file class
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
    #Save line where we have one peak in one file of my class
    allsavedgroups <- xset@groups[which(xset@groups[,class] > 0),]
    #Save line number where 0 peaks is in my files of my class
    alldeletedgroups <- which(xset@groups[,class] == 0)
    cat("\tDeleting",length(alldeletedgroups),"group(s) from groups tabledata and stock",nrow(allsavedgroups),"group(s) from class",class,"from xset\n")

    cat("STEP 2 : add grpid for each peak and delete when they are not in the good group\n")
    #STEP 2 : Delete peaks from the wrong groups (groups to delete)   
    #Add new col for grpid number for each peak
    grpidcol <- rep(NA, nrow(xset@peaks))
    xset@peaks <- cbind(xset@peaks,grpidcol)
    #Complete grpid for each peak
    for(x in 1:length(xset@groupidx)){
        for(y in 1:length(xset@groupidx[[x]])){
            xset@peaks[xset@groupidx[[x]][y],"grpidcol"] <- x
        }
    }
    #Delete peaks where groupidx have to be delete
    peaktodelete <- NULL
    peaktodelete <- which(xset@peaks[,"grpidcol"] %in% alldeletedgroups)
    if(length(peaktodelete) > 0){
        xset@peaks <- xset@peaks[-peaktodelete,]
    }
    cat("\tWe now have",nrow(xset@peaks),"peaks and deleting",length(peaktodelete),"peaks which were not in groups where we have peak from class \"",class,"\" !\n")
    
    cat("STEP 3 : Rebuild groupidx and groups\n")
    #STEP 3 : Rebuild groupidx with new row for groups and peaks (usefull for groupval after...)
    newgrpidx <- xset@groupidx
    newgrpidx <- lapply(newgrpidx, function(a) a <- NULL)
    for(p in 1:nrow(xset@peaks)){
        if(!(is.na(xset@peaks[p,"grpidcol"]))){
        	newgrpidx[[xset@peaks[p,"grpidcol"]]] <- c(newgrpidx[[xset@peaks[p,"grpidcol"]]],p)
        }else{
        	next
        }   
    }
    if(length(alldeletedgroups) > 0){
        newgrpidx <- newgrpidx[-alldeletedgroups]
    }
    xset@groupidx <- newgrpidx
    cat("\tRebuild groupidx with",length(xset@groupidx),"groupidx ! \n")
    xset@groups <- allsavedgroups
    cat("\tRebuild groups with",nrow(xset@groups),"groups ! \n")

    #Delete line grpidcol after order them 
    xset@peaks <- subset(xset@peaks, select = -(which(colnames(xset@peaks)=="grpidcol")))

    # cat("STEP 4 : Keep only MS files and MS/MS files from our class\n")
    # keep_file <- NULL
    # for(m in 1:length(xset@filepaths)){
    #     if(!(2 %in% unique(readMSData(xset@filepaths[m],mode="onDisk")@featureData@data$msLevel))) {
    #         keep_file <- c(keep_file, xset@filepaths[m])
    #         print(keep_file)
    #     }else{
    #         if(gsub(filepattern,"",basename(xset@filepaths[m])) 
    #             %in% sampleMetadataMSMS[which(sampleMetadataMSMS[,"class"] == class),"MSMS"]){
    #             keep_file <- c(keep_file,xset@filepaths[m])
    #             print(keep_file)
    #         }
    #     }
    # }
    # xset@filepaths <- xset@filepaths[match(keep_file,xset@filepaths)]

    return(xset)
}

#This function keep MSMS informations only from our files we are working with
modifyPaObject <- function(pa, files,class){
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")

	cat(paste("\n--------------- Modify pa object for \"",class,"\" class file(s) ---------------\n"))
	#cat("STEP 1 : modify pa@fileList\n")
	# keep_file <- NULL
 #    for(m in 1:length(pa@fileList)){
 #        if(!(2 %in% unique(readMSData(pa@fileList[m],mode="onDisk")@featureData@data$msLevel))) {
 #            keep_file <- c(keep_file, pa@fileList[m])
 #            print(keep_file)
 #        }else{
 #            if(gsub(filepattern,"",basename(pa@fileList[m])) 
 #                %in% sampleMetadataMSMS[which(sampleMetadataMSMS[,"class"] == class),"MSMS"]){
 #                keep_file <- c(keep_file,pa@fileList[m])
 #                print(keep_file)
 #            }
 #        }
 #    }
 #    pa@fileList <- pa@fileList[match(keep_file,pa@fileList)]

    #cat("STEP 1 : modify pa@puritydf\n")
    cat("Select datas for file(s) : ")
    cat(files,sep=", ")
    cat("\n")
	fileIndex <- which(gsub(filepattern,"",names(pa@fileList)) %in% files)
	pa@puritydf <- pa@puritydf[which(pa@puritydf[,"fileid"] %in% fileIndex),]
	return(pa)
}