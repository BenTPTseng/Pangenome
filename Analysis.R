#==============================================================================#
#-----------------------------------[MAIN]-------------------------------------#
#                    CpG (CCGG) Motif Stability Analysis                       #
#                     in a GBZ-derived Pangenome Graph                         #
#==============================================================================#

#-----------------------------------[GOAL]-------------------------------------#
# This script quantifies the stability of the CpG motif (CCGG) across
# vg-derived graph paths extracted from a GBZ chunk (chr6).
#
# A reference genome is used to identify CpG loci and construct a
# flanking-sequence context index (left + right flanks).
#
# These reference-defined CpG sites are then projected onto graph path
# sequences to assess whether the motif is preserved or mutated.
#
# NOTE:
# This is a reference-indexed, context-based motif conservation analysis
# performed on vg graph paths (from a GBZ chunk), NOT a coordinate-based
# alignment between reference and haplotypes.
#------------------------------------------------------------------------------#

#-----------------------------------[STEPS]------------------------------------#
# 1. Load reference genome (FASTA)
#
# 2. Identify all CpG (CCGG) occurrences in the reference genome
#    and extract flanking sequences (±15 bp)
#
# 3. Build a reference index:
#    flank_key = (left_flank + right_flank) → CpG locus ID(s)
#
# 4. Load vg graph paths from a GBZ chunk (chr6)
#    (paths represent graph-derived sequence traversals)
#
# 5. Extract full nucleotide sequences for each vg path using vg
#
# 6. Perform sliding-window scan over each path sequence:
#    - extract left flank, center motif, right flank
#    - construct flank_key for each window
#
# 7. Match flank_key against reference index to identify CpG loci
#
# 8. For matched loci:
#    - count observation (m)
#    - check whether center sequence equals "CCGG"
#    - count stability (s)
#
# 9. Aggregate results across all paths:
#    - total CpG loci detected
#    - stable vs unstable occurrences
#    - locus-specific stability (s/m)
#    - loci with no observations in graph paths
#------------------------------------------------------------------------------#
#==============================================================================#

#===============================================================================
# Load required packages
#===============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("data.table", "parallel", "stringdist", "stringi")

for (p in pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

#===============================================================================
# Initial settings
#===============================================================================
# Define the working directory:
wd_dir <- "" #FIXME
setwd(dir = wd_dir)

# Path to the vg executable on the computer
vg_path <- "" #FIXME

#===============================================================================
# Initial setup
#===============================================================================
ref_file   <- "GRCh38_chr6.fa"
gbz_file   <- "chunk_0_chr6.gbz"
path_file  <- "chr6_paths.txt"

flank      <- 15
motif      <- "CCGG"
motif_len  <- nchar(motif)
width      <- motif_len + 2 * flank

# cores      <- max(1, (parallel::detectCores() - 1) )
cores      <- 5
cat("Using cores:", cores, "\n")

block_size <- 5000000  # sliding window batch size to prevents RAM blowups on huge paths

#===============================================================================
# Custom function(s)
#===============================================================================
# FASTA LOADER #
read_fasta <- function(x) {
  cat("[1/6] Reading FASTA...\n")
  lines <- readLines(x)
  seq <- paste(lines[!grepl(pattern = "^>", x = lines)], collapse = "")
  cat(" FASTA loaded. Length:", nchar(seq), "\n")
  return(seq)
}

# VG PATH EXTRACTOR #
get_path_seq <- function(x) {
  out <- tryCatch(system2(vg_path,
                          args   = c("paths", "-x", gbz_file, "-F", "-Q", x),
                          stdout = TRUE,
                          stderr = FALSE),
                  error = function(e) NULL)

  if (is.null(out) || length(out) == 0) { return(NULL) }

  paste(out, collapse = "")
}

# STREAMING WINDOW PROCESSOR # (DO NOT RUN!)
# # V3
# process_window_chunk <- function(seq, start_idx, end_idx) {
# 
#   starts <- start_idx:end_idx
# 
#   wins <- substring(text  = seq, first = starts, last  = starts + width - 1)
# 
#   lefts   <- substr(x = wins, start = 1, stop = flank)
#   rights  <- substr(x = wins, start = flank + motif_len + 1, stop = width)
#   centers <- substr(x = wins, start = flank + 1, stop = flank + motif_len)
# 
#   flank_key <- paste0(lefts, rights)
# 
#   hit_idx <- match(flank_key, flank_names)
#   valid <- which(!is.na(hit_idx))
# 
#   if (!length(valid)) { return(NULL) }
# 
#   list(keys    = flank_names[hit_idx[valid]],
#        centers = centers[valid])
# }

#===============================================================================
# Main code
#===============================================================================

cat("Starting main code...\n")

#-------------------------------------------------------------------------------
# Step 1 — Read reference fasta
#-------------------------------------------------------------------------------
ref_seq <- read_fasta(x = ref_file)

#-------------------------------------------------------------------------------
# Step 2 — Read paths
#-------------------------------------------------------------------------------
cat("[2/6] Loading paths...\n")
paths <- readLines(path_file)
cat(" TOTAL PATHS: ", length(paths), "\n")

#-------------------------------------------------------------------------------
# Step 3 — Run motif scan
#-------------------------------------------------------------------------------
cat("[3/6] Scanning reference for motif...\n")

t0 <- Sys.time()
ref_pos <- stringi::stri_locate_all_fixed(pattern = motif,
                                          str     = ref_seq,
                                          overlap = FALSE)[[1]][,1]

cat(" Reference motif hits: ", length(ref_pos), "\n")
cat(" Scan time: ", round(difftime(time1 = Sys.time(),
                                   time2 = t0,
                                   units = "secs" ), 2),
    "sec\n")

ref_dt <- data.table::data.table(locus = seq_along(ref_pos),
                                 start = ref_pos)

ref_dt <- ref_dt[start > flank & (start + width - 1) <= nchar(ref_seq)]
# ref_dt <- ref_dt[ref_dt$start > flank & (ref_dt$start + width - 1) <= nchar(ref_seq)]

ref_dt[, kmer := substring(text  = ref_seq,
                           first = start - flank,
                           last  = start + flank + motif_len - 1)]

ref_dt[, left_flank  := substring(text  = kmer,
                                  first = 1,
                                  last  = flank)]

ref_dt[, right_flank := substring(text  = kmer,
                                  first = flank + motif_len + 1,
                                  last  = width)]

ref_dt[, flank_key := paste0(left_flank, right_flank)]

ref_dt[, locus_id := .I]

cat(" Reference loci:", nrow(ref_dt), "\n")

#-------------------------------------------------------------------------------
# Step 4 — Look up reference structures
#-------------------------------------------------------------------------------
cat("[4/6] Preparing reference structures...\n")

ref_map     <- split(x = ref_dt$locus_id, f = ref_dt$flank_key)

flank_index <- stats::setNames(object = seq_along(names(ref_map)), 
                               nm     = names(ref_map))

# ACCUMULATORS
m <- ref_dt |> nrow() |> integer()
s <- ref_dt |> nrow() |> integer()

#-------------------------------------------------------------------------------
# Step 5 — Main loop
#-------------------------------------------------------------------------------
cat("[5/6] Processing paths...\n")

for (p in seq_along(paths)) {
  
  cat("Path:", p, "/", length(paths), "\n")
  utils::flush.console()
  
  seq <- get_path_seq(paths[p])
  
  if (is.null(seq)) { next }
  if (nchar(seq) < width) { next }
  
  nwin       <- nchar(seq) - width + 1
  starts_all <- seq_len(nwin)
  
  #=========================================================
  # CHUNKED SLIDING WINDOW
  #=========================================================
  for (i in seq(1, nwin, by = block_size)) {
    
    j <- min(i + block_size - 1, nwin)
    
    starts <- starts_all[i:j]
    
    #-----------------------------------------------------
    # DEFINE WINDOW HERE
    #-----------------------------------------------------
    lefts <- substring(text  = seq, 
                       first = starts,
                       last  = starts + flank - 1)
    
    centers <- substring(text  = seq,
                         first = starts + flank,
                         last  = starts + flank + motif_len - 1)
    
    rights <- substring(text  = seq,
                        first = starts + flank + motif_len,
                        last  = starts + width - 1)
    
    flank_key <- paste0(lefts, rights)
    
    #-----------------------------------------------------
    # FAST HASH LOOKUP
    #-----------------------------------------------------
    hit_idx <- flank_index[flank_key]
    valid   <- which(!is.na(hit_idx))
    
    if (!length(valid)) { next }
    
    keys        <- names(ref_map)[hit_idx[valid]]
    centers_hit <- centers[valid]
    loci_list   <- ref_map[keys]
    lens        <- lengths(loci_list)
    
    if (all(lens == 0)) { next }
    
    all_loci <- unlist(loci_list, use.names = FALSE)
    
    # Expand centers to loci mapping
    center_rep <- rep(centers_hit, lens)
    
    #-----------------------------------------------------
    # ACCUMULATION
    #-----------------------------------------------------
    m[all_loci] <- m[all_loci] + 1L
    s[all_loci] <- s[all_loci] + (center_rep == motif)
  }
  
  # if (p %% 10 == 0) {
  #   cat("Completed paths:", p, "\n")
  #   flush.console()
  # }
}
#-------------------------------------------------------------------------------
# Step 6 — Final output
#-------------------------------------------------------------------------------
cat("[6/6] Returning final output...\n")

#===============================================================================
# FINAL OUTPUT
#===============================================================================

results <- data.table::data.table(locus = ref_dt$locus,
                                  kmer  = ref_dt$kmer,
                                  m     = m,
                                  s     = s)

results[, stability := data.table::fifelse(test = m > 0, 
                                           yes  = s / m, 
                                           no   = NA_real_)]

#===============================================================================
# SUMMARY
#===============================================================================
cat("\n================ FINAL SUMMARY ================\n")
cat("TOTAL LOCI:", nrow(results), "\n")
cat("MEAN STABILITY:", mean(results$stability, na.rm = TRUE), "\n")
cat("UNMAPPED LOCI:", sum(results$m == 0), "\n")

cat("\nTOP LOW STABILITY:\n")
print(results[order(stability)][1:10])
cat("\nEND OF SCRIPT\n")

cat("                                                         
                                                                                   
                     @@@                                                           
              *%    @=--@                                                          
             @--@:  %@=%@                                                          
             @%-=@.   @%@+                                                         
               .@=@.   @=@    .                                                    
                .@=@   @#@     ..                                                  
              .  @+%%  +@@.                                                        
             ..   @=@  @#*@                                                        
                  @=@==*==@@  .                                                    
               .@@=-------==@:                                                     
          .   @-@@-==::-----=@                                                     
             @%-*==*:@@----=-@                                                     
             @-----@-@=-=%=-=@.                                                    
            .@--------=-@----@.                                                    
             @=-=-----=@=----@=                                                    
              @@%==+@@@=%----#%                                                    
             .     @@==@=-==@@@ .                                                  
     .  @@. .    @@@-=-------==@@@+                                                
      @##=@=*- +@=%--=---------=---=@@                                             
     @@==@=*@  @-------=-------==---===-=@                                         
      @===%. .@=-@=----=-----@@@@@@@@@@=+@                                         
        @==@@@=-==-----------@@      @===@                   .                     
         @=--@-=@=-=----------@.     #%-@                       .                  
          @==-===-==----=-----#%    .@*-@.                                         
           @@=@@==---------=--=@    @=-=-@                                         
             #@===-=------=----@   @=@=@@                                          
              @=---------=--=-=%%    @@                    ..                      
             #%---=----==----===@@@@@@@@@@@@@.                                    .
             @-=--==--===--=-====@=-=----=-=-==@@         .                        
       ....  @===-=--=-=-----==-=-==----=-=--==-=@@       .                        
            .@=-------=======-=------=-----=-=====-@@.@@@@.                        
           .@---=---=---=-=----=--=----=--==-=-=-=-==+@@====@%                     
            @===-==----====-=---=----=---=--==-======--====--=@.                   
          . %*===----===-=--------=-==--=---==----=-========-==%@        .@@@      
             @--=-=-=-=-=------=---------====--==--=========-====@       .@==+@    
             .@-----=----=---=---=--=----=-===-=--=-==-===========%@      #+==@    
              .@=-------------=-----=*@@@@@==-=-====@@@@@@@@========+@@@@@=+==@    
                #@---=------------#@.  @@@@@@@@@@@@@@.    @@@@%==============@.    
                  .@@===-=====%@@.                          :@@@@=========%@.      
                                                                ..=++=....         
              
")  
  
