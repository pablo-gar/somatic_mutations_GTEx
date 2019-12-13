library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
if(!file.exists("./refcds_hg19_fastAccess.Rdata"))
    stop("RefCDS Rdata file not found, are you sure you are doing 'source(x, chdir = T)'")
load("./refcds_hg19_fastAccess.Rdata") # Contains RefCDS and gr_genes object(s)
#mutations <- read.table("../mutationCount/map/Liver/n6_0.0_0.5/SRR1100991.txt", sep = "\t", stringsAsFactors = F)

#' General function to query the RefCDS database
#' @param x - matching values for reference, only queries from this values will be returned
#' @param reference - the field matching x
#' @param query - the field for which value will be returned
#' @examples
#' queryRefCDS(reference = "gene_id", query = "CDS_length")
queryRefCDS <- function(x = NULL, referenceField, queryField, RefCDS = .RefCDS) {
    
    referenceField <- as.character(referenceField)
    queryField <- as.character(queryField)
    
    referenceInd <- setNames(1:length(RefCDS), sapply(RefCDS,function(x) x[[referenceField]]))
    queryValues <- setNames(sapply(RefCDS, function(x) x[[queryField]]), names(referenceInd))
    
    if(is.null(x)) {
        results <- queryValues
        names(results) <- names(referenceInd)
    } else {
        x <- as.character(x)
        results <- setNames(rep(NA, length(x)), x)
        referenceInd <- referenceInd[names(referenceInd) %in% x]
        if (length(referenceInd) == 0) {
            warning("No elements from x where found in the field ", referenceField)
            return (results)
        }
        
        queryValues <- queryValues[referenceInd]
        
        referenceSubs <- x %in% names(referenceInd)
        results[referenceSubs] <- queryValues[x[referenceSubs]]
    }
    
    return(results)
}
                        
    

#' Overlaps mutation table with genes
#' IMPORTANT, FOR GENE IN THE CRICK STRAND it returns the complementary of each mutation and its reference
#' @param mutations - data.frame, with positional columns: chromosome, position, ref base, alt base, context, coverage, n_supporting
#' @param RefCDS - object modified from dndscv containing info for all CDS in hg19 (https://github.com/im3sanger/dndscv)
#' @param gr_genes - GRanges object modified from dndscv containing the genomic ranges for all CDS in hg19
#' @return mutation table with two extra columns, gene_name and gene_id (ensembl)
overlapMutationsGenes <- function(mutations, RefCDS = .RefCDS, gr_genes = .gr_genes) {
    
    nt = c("A","C","G","T")
    
    colnames(mutations)[1:7] = c("chr","pos","ref","mut", "context", "coverage", "n_support")
    mutations$chr <- gsub("chr", "", mutations$chr)
    
    ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
    gr_genes_ind = ind[gr_genes$names]
    
    # Mapping mutations to genes
    gr_muts = GRanges(mutations$chr, IRanges(mutations$pos,mutations$pos))
    ol = as.matrix(findOverlaps(gr_muts, gr_genes, type="any", select="all"))
    mutations = mutations[ol[,1],] # Duplicating subs if they hit more than one gene
    mutations$geneind = gr_genes_ind[ol[,2]]
    mutations$gene_id = sapply(RefCDS,function(x) x$gene_id)[mutations$geneind]
    mutations$gene_name = sapply(RefCDS,function(x) x$gene_name)[mutations$geneind]
    mutations$cds_length = sapply(RefCDS,function(x) x$CDS_length)[mutations$geneind]
    
    # Annotating strand and strand specific mutations
    mutations$strand = sapply(RefCDS,function(x) x$strand)[mutations$geneind]
    mutations$ref_cod = mutations$ref
    mutations$mut_cod = mutations$mut
    isminus = (mutations$strand==-1)
    compnt = setNames(rev(nt), nt)
    mutations$ref_cod[isminus] = compnt[mutations$ref[isminus]]
    mutations$mut_cod[isminus] = compnt[mutations$mut[isminus]]
    
    return (mutations)
}
    
    

#' TODO
#' @param mutations - data.frame, with positional columns: chromosome, position, ref base, alt base, context, coverage, n_supporting
#' @param RefCDS - object modified from dndscv containing info for all CDS in hg19 (https://github.com/im3sanger/dndscv)
#' @param gr_genes - GRanges object modified from dndscv containing the genomic ranges for all CDS in hg19
#' @return TODO
annotateMutations <- function(mutations, RefCDS = .RefCDS, gr_genes = .gr_genes) {

    mutations <- overlapMutationsGenes(mutations, RefCDS, gr_genes)
    
    # Combinations of nucleotide, trinucleotides and substitutions
    nt = c("A","C","G","T")
    trinucs = paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
    trinucinds = setNames(1:64, trinucs)
    
    trinucsubs = NULL
    for (j in 1:length(trinucs)) {
        trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j],1,1), setdiff(nt,substr(trinucs[j],2,2)), substr(trinucs[j],3,3), sep=""), sep=">"))
    }
    trinucsubsind = setNames(1:192, trinucsubs)
    
    # Additional annotation of substitutions
    snv = (mutations$ref %in% nt & mutations$mut %in% nt)
    indels = mutations[!snv,]
    mutations = mutations[snv,]
    
    for (j in 1:length(RefCDS)) {
        RefCDS[[j]]$N = array(0, dim=c(192,4)) # Initialising the N matrices
    }
    
    # Subfunction: obtaining the codon positions of a coding mutation given the exon intervals
    
    chr2cds = function(pos,cds_int,strand) {
        if (strand==1) {
            return(which(pos==unlist(apply(cds_int, 1, function(x) x[1]:x[2]))))
        } else if (strand==-1) {
            return(which(pos==rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2])))))
        }
    }
    
    # Annotating the functional impact of each substitution and populating the N matrices
    
    ref3_cod = mut3_cod = wrong_ref = aachange = ntchange = impact = array(NA, nrow(mutations))
    
    for (j in 1:nrow(mutations)) {
    
        geneind = mutations$geneind[j]
        pos = mutations$pos[j]
        if (any(pos == RefCDS[[geneind]]$intervals_splice)) { # Essential splice-site substitution
        
            impact[j] = "Essential_Splice"; impind = 4
            pos_ind = (pos==RefCDS[[geneind]]$intervals_splice)
            cdsnt = RefCDS[[geneind]]$seq_splice[pos_ind]
            # CHANGED FROM ORIGINAL
            ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], mutations$ref_cod[j], RefCDS[[geneind]]$seq_splice1down[pos_ind])
            mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_splice1down[pos_ind])
            aachange[j] = ntchange[j] = "-"
        
        } else { # Coding substitution
        
            pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, RefCDS[[geneind]]$strand)
            cdsnt = RefCDS[[geneind]]$seq_cds[pos_ind]
            # CHANGED FROM ORIGINAL
            ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], mutations$ref_cod[j], RefCDS[[geneind]]$seq_cds1down[pos_ind])
            mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_cds1down[pos_ind])
            codon_pos = c(ceiling(pos_ind/3)*3-2, ceiling(pos_ind/3)*3-1, ceiling(pos_ind/3)*3)
            old_codon = as.character(as.vector(RefCDS[[geneind]]$seq_cds[codon_pos]))
            pos_in_codon = pos_ind-(ceiling(pos_ind/3)-1)*3
            #ADDED ON TOP OF ORIGINAL
            old_codon[pos_in_codon] = as.character(mutations$ref_cod[j])
            new_codon = old_codon; 
            new_codon[pos_in_codon] = as.character(mutations$mut_cod[j])
            old_aa = seqinr::translate(old_codon)
            new_aa = seqinr::translate(new_codon)
            aachange[j] = sprintf('%s%0.0f%s',old_aa,ceiling(pos_ind/3),new_aa)
            ntchange[j] = sprintf('%s%0.0f%s',mutations$ref_cod[j],pos_ind,mutations$mut_cod[j])
        
            # Annotating the impact of the mutation
            if (new_aa == old_aa){ 
                impact[j] = "Synonymous"; impind = 1
            } else if (new_aa == "*"){
                impact[j] = "Nonsense"; impind = 3
            } else if (old_aa != "*"){
                impact[j] = "Missense"; impind = 2
            } else if (old_aa=="*") {
                impact[j] = "Stop_loss"; impind = NA
            }
        }
        
        if (!is.na(impind)) { # Correct base annotation in the input mutation file
            trisub = trinucsubsind[ paste(ref3_cod[j], mut3_cod[j], sep=">") ]
            RefCDS[[geneind]]$N[trisub,impind] = RefCDS[[geneind]]$N[trisub,impind] + 1 # Adding the mutation to the N matrices
        }
      
        if (round(j/1e4)==(j/1e4)) { message(sprintf('    %0.3g %%...', round(j/nrow(mutations),2)*100)) }
    }
    
    mutations$ref3_cod = ref3_cod
    mutations$mut3_cod = mut3_cod
    mutations$aachange = aachange
    mutations$ntchange = ntchange
    mutations$impact = impact
    mutations$pid = sapply(RefCDS,function(x) x$protein_id)[mutations$geneind]

    return(mutations)
}

# Get last exon coordinates
get_last_exon_per_genes <- function() {
    result <- lapply(.RefCDS, FUN = function(x) {
                      if (x$strand == -1) { 
                          interval <- x$intervals_cds[1,] 
                      } else { 
                          interval <- x$intervals_cds[nrow(x$intervals_cds),]
                      }
                      
                      data.frame(gene_name = x$gene_name, strand = x$strand, chr = x$chr, last_exon_start = interval[1], last_exon_end = interval[2])
                    })

    result  <- do.call(rbind, result)
    return(result)
}
