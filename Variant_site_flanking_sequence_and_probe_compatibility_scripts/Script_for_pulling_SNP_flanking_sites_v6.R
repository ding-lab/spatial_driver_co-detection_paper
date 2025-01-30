suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
input_variants <- read.table("/diskmnt/Projects/Users/austins2/tools/xenium_snvs/version_4/input_variant_targets.tsv", header = T, sep = '\t')
ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 100)
#listFilterOptions(mart = ensembl, filter = 'ensembl_transcript_id')
row_names <- row.names(input_variants)
if (!("Ensembl_transcript_no_version" %in% colnames(input_variants))) {
    if ("ensembl_gene_ID" %in% colnames(input_variants)) {
        for (row_name in row_names) {
            ensembl_gene_ID <- input_variants[row_name,"ensembl_gene_ID"]
            ensembl_transcript_table = getBM(attributes=c('ensembl_gene_id',
                                                          'ensembl_transcript_id',
                                                          'ensembl_transcript_id_version',
                                                          'transcript_is_canonical'),
                                             filters=c('ensembl_gene_id'),
                                             values=list(ENSG=ensembl_gene_ID), 
                                             mart=ensembl, 
                                             checkFilters=FALSE)
            transcript_canon <- subset(ensembl_transcript_table, transcript_is_canonical == 1)
            input_variants[row_name,"ensembl_transcript_id"] <- ensembl_transcript_table[(ensembl_transcript_table$transcript_is_canonical == 1), "ensembl_transcript_id"][is.na(ensembl_transcript_table[(ensembl_transcript_table$transcript_is_canonical == 1), "ensembl_transcript_id"]) == FALSE]
            input_variants[row_name,"ensembl_transcript_id_version"] <- ensembl_transcript_table[(ensembl_transcript_table$transcript_is_canonical == 1), "ensembl_transcript_id_version"][is.na(ensembl_transcript_table[(ensembl_transcript_table$transcript_is_canonical == 1), "ensembl_transcript_id_version"]) == FALSE]
        }
    } else if ("hugo_symbol" %in% colnames(input_variants)) {
        for (row_name in row_names) {
            hgnc_symbol <- input_variants[row_name,"hugo_symbol"]
            ensembl_transcript_table = getBM(attributes=c('ensembl_gene_id',
                                                          'ensembl_transcript_id',
                                                          'ensembl_transcript_id_version',
                                                          'transcript_is_canonical'),
                                             filters=c('hgnc_symbol'),
                                             values=list(HGNC=hgnc_symbol), 
                                             mart=ensembl, 
                                             checkFilters=FALSE)
            transcript_canon <- subset(ensembl_transcript_table, transcript_is_canonical == 1)
            input_variants[row_name,"Ensembl_transcript_no_version"] <- ensembl_transcript_table[(ensembl_transcript_table$transcript_is_canonical == 1), "ensembl_transcript_id"][is.na(ensembl_transcript_table[(ensembl_transcript_table$transcript_is_canonical == 1), "ensembl_transcript_id"]) == FALSE]
            input_variants[row_name,"Ensembl_transcript_id_version"] <- ensembl_transcript_table[(ensembl_transcript_table$transcript_is_canonical == 1), "ensembl_transcript_id_version"][is.na(ensembl_transcript_table[(ensembl_transcript_table$transcript_is_canonical == 1), "ensembl_transcript_id_version"]) == FALSE]
        }
    }
}
#ensembl_transcripts = c('ENST00000256078')
#ensembl_transcripts = c('ENST00000263967')
ensembl_transcripts <- unique(input_variants$Ensembl_transcript_no_version)
ensembl_exons_table = getBM(attributes=c('ensembl_gene_id',
                                         'ensembl_transcript_id',
                                         'ensembl_transcript_id_version',
                                         'ensembl_exon_id',
                                         'exon_chrom_start',
                                         'exon_chrom_end',
                                         'strand',
                                         'rank'),
                            filters=c('ensembl_transcript_id'),
                            values=list(ENST=ensembl_transcripts), 
                            mart=ensembl, 
                            checkFilters=FALSE)
ensembl_exons_table = as.data.frame(ensembl_exons_table)
# > ensembl_exons_table # PIK3CA-201
#    ensembl_gene_id ensembl_transcript_id ensembl_transcript_id_version                                                                                                                                                                                               
# 1  ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 2  ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 3  ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 4  ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 5  ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 6  ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 7  ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 8  ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 9  ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 10 ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 11 ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 12 ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 13 ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 14 ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 15 ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 16 ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 17 ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 18 ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 19 ENSG00000121879       ENST00000263967             ENST00000263967.4                                                                                                                                                                                               
# 20 ENSG00000121879       ENST00000263967             ENST00000263967.4
# 21 ENSG00000121879       ENST00000263967             ENST00000263967.4
#    ensembl_exon_id exon_chrom_start exon_chrom_end strand rank
# 1  ENSE00001493081        179148357      179148603      1    1
# 2  ENSE00001139995        179198750      179199177      1    2
# 3  ENSE00000997375        179199690      179199899      1    3
# 4  ENSE00001077693        179201290      179201540      1    4
# 5  ENSE00001077692        179203544      179203789      1    5
# 6  ENSE00001077694        179204503      179204588      1    6
# 7  ENSE00001077691        179209595      179209700      1    7
# 8  ENSE00001128470        179210186      179210338      1    8
# 9  ENSE00001128465        179210431      179210565      1    9
# 10 ENSE00001077674        179218210      179218334      1   10
# 11 ENSE00000826291        179219196      179219277      1   11
# 12 ENSE00000826292        179219571      179219735      1   12
# 13 ENSE00003485038        179219949      179220052      1   13
# 14 ENSE00003568097        179220986      179221157      1   14
# 15 ENSE00003489671        179224081      179224187      1   15
# 16 ENSE00003485539        179224700      179224821      1   16
# 17 ENSE00000826297        179225962      179226040      1   17
# 18 ENSE00000826298        179229272      179229442      1   18
# 19 ENSE00000826299        179230004      179230121      1   19
# 20 ENSE00000826300        179230225      179230376      1   20
# 21 ENSE00001139987        179234094      179240093      1   21
#   ensembl_gene_id ensembl_transcript_id ensembl_transcript_id_version # KRAS-201
# 1 ENSG00000133703       ENST00000256078            ENST00000256078.10
# 2 ENSG00000133703       ENST00000256078            ENST00000256078.10
# 3 ENSG00000133703       ENST00000256078            ENST00000256078.10
# 4 ENSG00000133703       ENST00000256078            ENST00000256078.10
# 5 ENSG00000133703       ENST00000256078            ENST00000256078.10
# 6 ENSG00000133703       ENST00000256078            ENST00000256078.10
#   ensembl_exon_id exon_chrom_start exon_chrom_end strand rank
# 1 ENSE00001644818         25225614       25225773     -1    4
# 2 ENSE00002477035         25205246       25209911     -1    6
# 3 ENSE00003903543         25250751       25250929     -1    1
# 4 ENSE00000936617         25245274       25245395     -1    2
# 5 ENSE00001719809         25227234       25227412     -1    3
# 6 ENSE00001189807         25215437       25215560     -1    5
# These provides a table with just 1 row per transcript sequence and all of the info is in comma separated lists for each sequence
# ensembl_exons_for_variant=getBM(attributes=c('ensembl_gene_id',
#                                              'ensembl_transcript_id',
#                                              'ensembl_transcript_id_version',
#                                              'ensembl_exon_id',
#                                              'exon_chrom_start',
#                                              'exon_chrom_end',
#                                              'strand','transcript_flank','rank'),
#                                 filters=c('ensembl_transcript_id','upstream_flank'),
#                                 values=list(ENST='ENST00000256078',Upstream=100), 
#                                 mart=ensembl, 
#                                 checkFilters=FALSE)
# ensembl_exons_for_variant=getBM(attributes=c('ensembl_gene_id',
#                                              'ensembl_transcript_id',
#                                              'ensembl_transcript_id_version',
#                                              'ensembl_exon_id',
#                                              'exon_chrom_start',
#                                              'exon_chrom_end',
#                                              'strand','transcript_flank'),
#                                 filters=c('ensembl_transcript_id','downstream_flank'),
#                                 values=list(ENST='ENST00000256078',Downstream=100), 
#                                 mart=ensembl, 
#                                 checkFilters=FALSE)
hgnc_symbol_table = getBM(attributes=c('ensembl_transcript_id',
                                       'hgnc_symbol'),
                          filters=c('ensembl_transcript_id'),
                          values=list(ENST=ensembl_transcripts), 
                          mart=ensembl, 
                          checkFilters=FALSE)
#   ensembl_transcript_id hgnc_symbol
# 1       ENST00000263967      PIK3CA
# ensembl_exons_for_variant=getBM(attributes=c('ensembl_gene_id',
#                                              'ensembl_transcript_id',
#                                              'ensembl_transcript_id_version','transcript_mane_select'),
#                                 filters=c('ensembl_transcript_id'),
#                                 values=list(ENST='ENST00000311936'), 
#                                 mart=ensembl, 
#                                 checkFilters=FALSE)
row_names <- row.names(input_variants)
input_variants$Variant_front_flank_Start <- NA
input_variants$Variant_front_flank_Stop <- NA
input_variants$Variant_rear_flank_Start <- NA
input_variants$Variant_rear_flank_Stop <- NA
input_variants$Variant_front_flank_sequence <- NA
input_variants$Variant_rear_flank_sequence <- NA
input_variants$Variant_complete_flank_sequence <- NA
input_variants$rear_flank_spans_introns <- FALSE
input_variants$front_flank_spans_introns <- FALSE
input_variants$Front_flank_full_coordinates <- NA
input_variants$Front_flank_coordinates_based_sequence <- NA
input_variants$Rear_flank_full_coordinates <- NA
input_variants$Rear_flank_coordinates_based_sequence <- NA
input_variants$Variant_is_5bp_or_closer_to_a_splice_junction <- FALSE
flank_length = 21
ensembl_exons_table$exon_length = abs(ensembl_exons_table$exon_chrom_start - ensembl_exons_table$exon_chrom_end)+1
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringi))
for (row_name in row_names) {
    var_ensembl_transcript_id = input_variants[row_name,"Ensembl_transcript_no_version"]
    print(var_ensembl_transcript_id)
    var_start = input_variants[row_name,"Start"]
    var_stop = input_variants[row_name,"Stop"]
    var_ref = input_variants[row_name,"Ref_allele"]
    var_alt = input_variants[row_name,"Alt_allele"]
    var_chr = input_variants[row_name,"Chromosome"]
    transcript_exons_table = ensembl_exons_table[(ensembl_exons_table$ensembl_transcript_id == var_ensembl_transcript_id),]
    transcript_exons_table = transcript_exons_table[order(transcript_exons_table$rank, decreasing=F), ]
    rownames(transcript_exons_table) <- transcript_exons_table$rank
    transcript_hgnc_table = hgnc_symbol_table[(hgnc_symbol_table$ensembl_transcript_id == var_ensembl_transcript_id),]
    print(transcript_hgnc_table)
    strand = unique(transcript_exons_table$strand)
    # trasnscript_txdb <- txdb<-makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
    #                     dataset="hsapiens_gene_ensembl",
    #                     transcript_ids=var_ensembl_transcript_id)
    # exons_txdb <- exons(txdb)
    # check if the transcript is on the forwards strand (1 is forwards, -1 is reverse)
    if (strand > 0) {
        #forward strand is 1
        strand_char = "+"
        #transcript_exons_start_table = transcript_exons_table[((transcript_exons_table$exon_chrom_start < var_start) & (transcript_exons_table$exon_chrom_stop > var_start)),]
        #transcript_exons_start_initial = transcript_exons_start_table$exon_chrom_start
        #transcript_exons_stop_initial = transcript_exons_start_table$exon_chrom_stop
        #exon_start_distance = abs(transcript_exons_start_initial - var_start) #for insertions or deletions this logic will need to be changed
        #exon_stop_distance = abs(transcript_exons_stop_initial - var_stop) #for insertions or deletions this logic will need to be changed
        # order the exons by exon_chrom_start site (create a column for exon order (ordered_exon number)) # <- this information is already present in the rank column
        # for each exon calculate the length from the exon start to the exon stop and save to new value in column #done above and saved to ensembl_exons_table$exon_length
        # for each exon use the exon start to the exon stop to pull the sequence from Granges. 
        exon_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                 paste(rep('chr', times = length(transcript_exons_table$ensembl_exon_id)), 
                                       rep(var_chr,times=length(transcript_exons_table$ensembl_exon_id)),sep=''),
                                 start=transcript_exons_table$exon_chrom_start,
                                 end=transcript_exons_table$exon_chrom_end,
                                 strand=rep(strand_char, times = length(transcript_exons_table$ensembl_exon_id)))
        if (length(transcript_exons_table$ensembl_exon_id) == 1) {
            print("single exon transcript")
            transcript_exons_table$sequence <- as.character(exon_sequences)
            print(transcript_exons_table$sequence)
        } else {
            exon_sequences_df <- as.data.frame(exon_sequences)
            transcript_exons_table$sequence <- exon_sequences_df$x
        }
        #print(transcript_exons_table)
        print(var_start)
        print(var_stop)
        # concat all of the Granges sequences together.
        transcript_sequence <- paste(transcript_exons_table$sequence, collapse='')
        # figure out which exon contains the snv and save to variable (make it so a single line can handle forwards or reverse that way you can just copy it below)
        exon_with_snv <- transcript_exons_table$ensembl_exon_id[(((transcript_exons_table$exon_chrom_start <= var_start) & (transcript_exons_table$exon_chrom_end >= var_stop)) | ((transcript_exons_table$exon_chrom_start >= var_start) & (transcript_exons_table$exon_chrom_end <= var_stop)))]
        print(exon_with_snv)
        # save the ordered exon number for the exon that the snv to new variable
        exon_with_snv_rank <- transcript_exons_table$rank[ transcript_exons_table$ensembl_exon_id == exon_with_snv ]
        print(exon_with_snv_rank)
        # find the distance from the var_start to the start of the exon that the SNV is present in  # do a test to make sure this is correct 
        snv_to_exon_start_distance = abs(var_start - transcript_exons_table$exon_chrom_start[transcript_exons_table$ensembl_exon_id == exon_with_snv])
        print(snv_to_exon_start_distance)
        # find the distance from the var_stop to the end of the exon that the SNV is present in  # do a test to make sure this is correct 
        snv_to_exon_stop_distance = abs(var_stop - transcript_exons_table$exon_chrom_end[transcript_exons_table$ensembl_exon_id == exon_with_snv])
        print(snv_to_exon_stop_distance)
        # find the distance for the start of the transcript to the snv
        snv_to_transcript_start_distance = snv_to_exon_start_distance + sum(transcript_exons_table$exon_length[transcript_exons_table$rank < exon_with_snv_rank])
        # find the distance from the end of the transcript to the snv
        snv_to_transcript_stop_distance = snv_to_exon_stop_distance + sum(transcript_exons_table$exon_length[transcript_exons_table$rank > exon_with_snv_rank])
        # calculate length of leftover flank in Upstream (will be positive if the SNV is closer to the start than flank length)
        upstream_flank_distance = flank_length - snv_to_transcript_start_distance
        # calculate length of leftover flank in Downstream (will be positive if the SNV is closer to the stop than flank length)
        downstream_flank_distance = flank_length - snv_to_transcript_stop_distance
        # check if the snv is the upstream flank distance is greater than 0
        print(upstream_flank_distance)
        print(flank_length)
        print(snv_to_transcript_start_distance)
        if (upstream_flank_distance > 0) {
            # check if the snv is the upstream flank distance is less than 0
            if (downstream_flank_distance > 0){
                print("This transcript is very short. the distance from the transcript start to the SNV start is less than the flank length (",flank_length," bp). Also, the distance from the transcript stop to the SNV stop is also less than the flank length (",flank_length," bp). As a result, some of the flanking sequence includes the 5' sequence upstream of the transcript and the 3' sequence downstream of the transcript. These upstream and downstream sequences are not necessarily a part of the UTRs. This means they may or may not be transcribed and should be manually reviewed. They are denoted by a '|' in the output sequence.",sep='')
                upstream_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                            paste('chr', var_chr, sep=''),
                                            start=transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == 1] - upstream_flank_distance, # extract only the needed length of the leftover flank from the Upstream sequenced
                                            end=transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == 1] - 1, #subtract 1 to offset since bounds are inclusive and exon start nucleotide is already included above.
                                            strand=strand_char)
                upstream_sequence <- as.character(upstream_sequence) # change to character for pasting.
                downstream_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                              paste('chr', var_chr, sep=''),
                                              start=transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == max(transcript_exons_table$rank)] + 1 , # add 1 to offset since bounds are inclusive and exon stop nucleotide is already included above.
                                              end=transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == max(transcript_exons_table$rank)] + downstream_flank_distance, # extract only the needed length of the leftover flank from the downstream sequence (since bounds are inclusive don't add 1 here)
                                              strand=strand_char)
                downstream_sequence <- as.character(downstream_sequence) # change to character for pasting.
                front_flank_remainder_in_transcript = flank_length - upstream_flank_distance
                rear_flank_remainder_in_transcript = flank_length - downstream_flank_distance
                front_transcript_flank <- stringr::str_sub(transcript_sequence, 1, snv_to_transcript_start_distance)
                rear_transcript_flank <- stringr::str_sub(transcript_sequence, snv_to_transcript_start_distance+1+nchar(var_ref), snv_to_transcript_start_distance+1+nchar(var_ref)+rear_flank_remainder_in_transcript-1)
                flank_sequence <- paste(upstream_sequence,'|',front_transcript_flank,'[',var_ref,'/',var_alt,']',rear_transcript_flank,'|',downstream_sequence,sep='')
                # record the flank sequence on the snv_table
                input_variants[row_name,"Variant_front_flank_sequence"] <- paste(upstream_sequence,'|',front_transcript_flank,sep='')
                input_variants[row_name,"Variant_rear_flank_sequence"] <- paste(rear_transcript_flank,'|',downstream_sequence,sep='')
                input_variants[row_name,"Variant_complete_flank_sequence"] <- flank_sequence
            } else {
                # warn the user that the distance from snv to the start of exon 1 is less than the flank length. and the flank will include some of the 5' sequence upstream of exon 1
                print(paste("Warning: SNV distance to the start of the first exon is less than flank length (",flank_length," bp). As a result, some of the flanking sequence will include the 5' upstream sequence. This will be denoted by a '|' in the output sequence. The downstream sequence may not be part of the the 5' UTR.", sep=''))
                # pull the transcript_flank in a new biomart query for that transcript.
                upstream_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                            paste('chr', var_chr, sep=''),
                                            start=transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == 1] - upstream_flank_distance, # extract only the needed length of the leftover flank from the Upstream sequenced
                                            end=transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == 1] - 1, #subtract 1 to offset since bounds are inclusive and exon start nucleotide is already included above.
                                            strand=strand_char)
                upstream_sequence <- as.character(upstream_sequence)
                # calculate length of flank in exon_1 # this is just the snv_to_exon_start_distance needed
                # slice the length of front flank from Granges sequence.
                front_flank_remainder_in_transcript = flank_length - upstream_flank_distance
                front_transcript_flank <- stringr::str_sub(transcript_sequence, 1, snv_to_transcript_start_distance)
                # slice the length of rear flank from Granges sequence.
                rear_transcript_flank <- stringr::str_sub(transcript_sequence, snv_to_transcript_start_distance+1+nchar(var_ref), snv_to_transcript_start_distance+1+nchar(var_ref)+flank_length-1) # both start and stop are inclusive so the +1+nchar(var_ref) is the offest so the SNV site isn't included and the -1 is the offset so we don't include an extra base on the rear flank
                # paste Upstream flank, front_exon flank, "[",ref,"/",alt,"]", rear flank and save to variant table
                flank_sequence <- paste(upstream_sequence,'|',front_transcript_flank,'[',var_ref,'/',var_alt,']',rear_transcript_flank,sep='')
                # record the flank sequence on the snv_table
                input_variants[row_name,"Variant_front_flank_sequence"] <- paste(upstream_sequence,'|',front_transcript_flank,sep='')
                input_variants[row_name,"Variant_rear_flank_sequence"] <- rear_transcript_flank
                input_variants[row_name,"Variant_complete_flank_sequence"] <- flank_sequence
            }
            # check if just the downstream flank distance is greater than 0
        } else if (downstream_flank_distance > 0) {
            print(paste("Warning: SNV distance to the end of the last exon is less than flank length (",flank_length," bp). As a result, some of the flanking sequence will include the 3' downstream sequence. This will be denoted by a '|' in the output sequence. The downstream sequence may not be part of the the 5' UTR.", sep=''))
            downstream_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                          paste('chr', var_chr, sep=''),
                                          start=transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == max(transcript_exons_table$rank)] + 1 , # add 1 to offset since bounds are inclusive and exon stop nucleotide is already included above.
                                          end=transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == max(transcript_exons_table$rank)] + downstream_flank_distance, # extract only the needed length of the leftover flank from the downstream sequence (since bounds are inclusive don't add 1 here)
                                          strand=strand_char)
            downstream_sequence <- as.character(downstream_sequence) # change to character for pasting.
            front_transcript_flank <- stringr::str_sub(transcript_sequence, snv_to_transcript_start_distance+1-flank_length, snv_to_transcript_start_distance)
            rear_flank_remainder_in_transcript = flank_length - downstream_flank_distance
            rear_transcript_flank <- stringr::str_sub(transcript_sequence, snv_to_transcript_start_distance+1+nchar(var_ref), snv_to_transcript_start_distance+1+nchar(var_ref)+rear_flank_remainder_in_transcript-1)
            flank_sequence <- paste(front_transcript_flank,'[',var_ref,'/',var_alt,']',rear_transcript_flank,'|',downstream_sequence,sep='')
            # record the flank sequence on the snv_table
            input_variants[row_name,"Variant_front_flank_sequence"] <- front_transcript_flank
            input_variants[row_name,"Variant_rear_flank_sequence"] <- paste(rear_transcript_flank,'|',downstream_sequence,sep='')
            input_variants[row_name,"Variant_complete_flank_sequence"] <- flank_sequence
        } else {
            print("SNV is located further from the start and end of the transcript than the flank length")
            front_transcript_flank <- stringr::str_sub(transcript_sequence, snv_to_transcript_start_distance+1-flank_length, snv_to_transcript_start_distance)
            rear_transcript_flank <- stringr::str_sub(transcript_sequence, snv_to_transcript_start_distance+1+nchar(var_ref), snv_to_transcript_start_distance+1+nchar(var_ref)+flank_length-1) # both start and stop are inclusive so the +1+nchar(var_ref) is the offest so the SNV site isn't included and the -1 is the offset so we don't include an extra base on the rear flank
            flank_sequence <- paste(front_transcript_flank,'[',var_ref,'/',var_alt,']',rear_transcript_flank,sep='')
            # record the flank sequence on the snv_table
            input_variants[row_name,"Variant_front_flank_sequence"] <- front_transcript_flank
            input_variants[row_name,"Variant_rear_flank_sequence"] <- rear_transcript_flank
            input_variants[row_name,"Variant_complete_flank_sequence"] <- flank_sequence
        }
    } else {
        # reverse strand is -1
        #forward strand is 1
        strand_char = "-"
        library(stringi)
        var_ref <- stri_reverse(chartr("ATCG","TAGC",var_ref))
        var_alt <- stri_reverse(chartr("ATCG","TAGC",var_alt))
        # order the exons by exon_chrom_start site (create a column for exon order (ordered_exon number)) # <- this information is already present in the rank column
        # for each exon calculate the length from the exon start to the exon stop and save to new value in column #done above and saved to ensembl_exons_table$exon_length
        # for each exon use the exon start to the exon stop to pull the sequence from Granges. 
        exon_sequences <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                 paste(rep('chr', times = length(transcript_exons_table$ensembl_exon_id)), 
                                       rep(var_chr,times=length(transcript_exons_table$ensembl_exon_id)),sep=''),
                                 start=transcript_exons_table$exon_chrom_start,
                                 end=transcript_exons_table$exon_chrom_end,
                                 strand=rep(strand_char, times = length(transcript_exons_table$ensembl_exon_id)))
        if (length(transcript_exons_table$ensembl_exon_id) == 1) {
            print("single exon transcript")
            transcript_exons_table$sequence <- as.character(exon_sequences)
            print(transcript_exons_table$sequence)
        } else {
            exon_sequences_df <- as.data.frame(exon_sequences)
            transcript_exons_table$sequence <- exon_sequences_df$x
        }
        # concat all of the Granges sequences together.
        transcript_sequence <- paste(transcript_exons_table$sequence, collapse='')
        # figure out which exon contains the snv and save to variable (make it so a single line can handle forwards or reverse that way you can just copy it below)
        exon_with_snv <- transcript_exons_table$ensembl_exon_id[(((transcript_exons_table$exon_chrom_start <= var_start) & (transcript_exons_table$exon_chrom_end >= var_stop)) | ((transcript_exons_table$exon_chrom_start >= var_start) & (transcript_exons_table$exon_chrom_end <= var_stop)))]
        # save the ordered exon number for the exon that the snv to new variable
        exon_with_snv_rank <- transcript_exons_table$rank[ transcript_exons_table$ensembl_exon_id == exon_with_snv ]
        # find the distance from the var_start to the start of the exon that the SNV is present in  # do a test to make sure this is correct 
        snv_to_exon_start_distance = abs(var_start - transcript_exons_table$exon_chrom_end[transcript_exons_table$ensembl_exon_id == exon_with_snv]) # flip start and stop in the reverse strand because the coordinate system is reversed
        # find the distance from the var_stop to the end of the exon that the SNV is present in  # do a test to make sure this is correct 
        snv_to_exon_stop_distance = abs(var_stop - transcript_exons_table$exon_chrom_start[transcript_exons_table$ensembl_exon_id == exon_with_snv]) # flip start and stop in the reverse strand because the coordinate system is reversed
        # find the distance for the start of the transcript to the snv
        snv_to_transcript_start_distance = snv_to_exon_start_distance + sum(transcript_exons_table$exon_length[transcript_exons_table$rank < exon_with_snv_rank])
        # find the distance from the end of the transcript to the snv
        snv_to_transcript_stop_distance = snv_to_exon_stop_distance + sum(transcript_exons_table$exon_length[transcript_exons_table$rank > exon_with_snv_rank])
        # calculate length of leftover flank in Upstream (will be positive if the SNV is closer to the start than flank length)

        print(exon_with_snv)
        #print(315)
        print(snv_to_exon_start_distance)
        print(snv_to_transcript_start_distance)
        upstream_flank_distance = flank_length - snv_to_transcript_start_distance
        print(upstream_flank_distance)
        # calculate length of leftover flank in Downstream (will be positive if the SNV is closer to the stop than flank length)
        downstream_flank_distance = flank_length - snv_to_transcript_stop_distance
        if (upstream_flank_distance > 0) {
            # check if the snv is the upstream flank distance is less than 0
            if (downstream_flank_distance > 0){
                print("This transcript is very short. the distance from the transcript start to the SNV start is less than the flank length (",flank_length," bp). Also, the distance from the transcript stop to the SNV stop is less than the flank length (",flank_length," bp). As a result, some of the flanking sequence includes the 5' sequence upstream of the transcript and the 3' sequence downstream of the transcript. These upstream and downstream sequences are not necessarily a part of the UTRs. This means they may or may not be transcribed and should be manually reviewed. They are denoted by a '|' in the output sequence.",sep='')
                upstream_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                            paste('chr', var_chr, sep=''),
                                            start=transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == 1] + 1, # extract only the needed length of the leftover flank from the Upstream sequenced
                                            end=transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == 1] + upstream_flank_distance, #subtract 1 to offset since bounds are inclusive and exon start nucleotide is already included above.
                                            strand=strand_char)
                upstream_sequence <- as.character(upstream_sequence) # change to character for pasting.
                downstream_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                              paste('chr', var_chr, sep=''),
                                              start=transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == max(transcript_exons_table$rank)] - downstream_flank_distance , # add 1 to offset since bounds are inclusive and exon stop nucleotide is already included above.
                                              end=transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == max(transcript_exons_table$rank)] - 1 , # extract only the needed length of the leftover flank from the downstream sequence (since bounds are inclusive don't add 1 here)
                                              strand=strand_char)
                downstream_sequence <- as.character(downstream_sequence) # change to character for pasting.
                front_flank_remainder_in_transcript = flank_length - upstream_flank_distance
                rear_flank_remainder_in_transcript = flank_length - downstream_flank_distance
                front_transcript_flank <- stringr::str_sub(transcript_sequence, 1, snv_to_transcript_start_distance)
                rear_transcript_flank <- stringr::str_sub(transcript_sequence, snv_to_transcript_start_distance+1+nchar(var_ref), snv_to_transcript_start_distance+1+nchar(var_ref)+rear_flank_remainder_in_transcript-1)
                flank_sequence <- paste(upstream_sequence,'|',front_transcript_flank,'[',var_ref,'/',var_alt,']',rear_transcript_flank,'|',downstream_sequence,sep='')
                # record the flank sequence on the snv_table
                input_variants[row_name,"Variant_front_flank_sequence"] <- paste(upstream_sequence,'|',front_exon_flank,sep='')
                input_variants[row_name,"Variant_rear_flank_sequence"] <- paste(rear_exon_flank,'|',downstream_sequence,sep='')
                input_variants[row_name,"Variant_complete_flank_sequence"] <- flank_sequence
            } else {
                # warn the user that the distance from snv to the start of exon 1 is less than the flank length. and the flank will include some of the 5' sequence upstream of exon 1
                print(paste("Warning: SNV distance to the start of the first exon is less than flank length (",flank_length," bp). As a result, some of the flanking sequence will include the 5' upstream sequence. This will be denoted by a '|' in the output sequence. The downstream sequence may not be part of the the 5' UTR.", sep=''))
                # pull the transcript_flank in a new biomart query for that transcript.
                upstream_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                            paste('chr', var_chr, sep=''),
                                            start=transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == 1] + 1, # extract only the needed length of the leftover flank from the Upstream sequenced
                                            end=transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == 1] + upstream_flank_distance, #subtract 1 to offset since bounds are inclusive and exon start nucleotide is already included above.
                                            strand=strand_char)
                upstream_sequence <- as.character(upstream_sequence)
                # calculate length of flank in exon_1 # this is just the snv_to_exon_start_distance needed
                # slice the length of front flank from Granges sequence.
                front_flank_remainder_in_transcript = flank_length - upstream_flank_distance
                front_transcript_flank <- stringr::str_sub(transcript_sequence, 1, snv_to_transcript_start_distance)
                # slice the length of rear flank from Granges sequence.
                rear_transcript_flank <- stringr::str_sub(transcript_sequence, snv_to_transcript_start_distance+1+nchar(var_ref), snv_to_transcript_start_distance+1+nchar(var_ref)+flank_length-1) # both start and stop are inclusive so the +1+nchar(var_ref) is the offest so the SNV site isn't included and the -1 is the offset so we don't include an extra base on the rear flank
                # paste Upstream flank, front_exon flank, "[",ref,"/",alt,"]", rear flank and save to variant table
                flank_sequence <- paste(upstream_sequence,'|',front_transcript_flank,'[',var_ref,'/',var_alt,']',rear_transcript_flank,sep='')
                # record the flank sequence on the snv_table
                input_variants[row_name,"Variant_front_flank_sequence"] <- paste(upstream_sequence,'|',front_transcript_flank,sep='')
                input_variants[row_name,"Variant_rear_flank_sequence"] <- rear_transcript_flank
                input_variants[row_name,"Variant_complete_flank_sequence"] <- flank_sequence
            }
            # check if just the downstream flank distance is greater than 0
        } else if (downstream_flank_distance > 0) {
            print(paste("Warning: SNV distance to the end of the last exon is less than flank length (",flank_length," bp). As a result, some of the flanking sequence will include the 3' downstream sequence. This will be denoted by a '|' in the output sequence. The downstream sequence may not be part of the the 5' UTR.", sep=''))
            downstream_sequence <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                          paste('chr', var_chr, sep=''),
                                          start=transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == max(transcript_exons_table$rank)] - downstream_flank_distance , # add 1 to offset since bounds are inclusive and exon stop nucleotide is already included above.
                                          end=transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == max(transcript_exons_table$rank)] - 1 , # extract only the needed length of the leftover flank from the downstream sequence (since bounds are inclusive don't add 1 here)
                                          strand=strand_char)
            downstream_sequence <- as.character(downstream_sequence) # change to character for pasting.
            front_transcript_flank <- stringr::str_sub(transcript_sequence, snv_to_transcript_start_distance+1-flank_length, snv_to_transcript_start_distance)
            rear_flank_remainder_in_transcript = flank_length - downstream_flank_distance
            rear_transcript_flank <- stringr::str_sub(transcript_sequence, snv_to_transcript_start_distance+1+nchar(var_ref), snv_to_transcript_start_distance+1+nchar(var_ref)+rear_flank_remainder_in_transcript-1)
            flank_sequence <- paste(front_transcript_flank,'[',var_ref,'/',var_alt,']',rear_transcript_flank,'|',downstream_sequence,sep='')
            # record the flank sequence on the snv_table
            input_variants[row_name,"Variant_front_flank_sequence"] <- front_transcript_flank
            input_variants[row_name,"Variant_rear_flank_sequence"] <- paste(rear_transcript_flank,'|',downstream_sequence,sep='')
            input_variants[row_name,"Variant_complete_flank_sequence"] <- flank_sequence
        } else {
            print("SNV is located further from the start and end of the transcript than the flank length")
            front_transcript_flank <- stringr::str_sub(transcript_sequence, snv_to_transcript_start_distance+1-flank_length, snv_to_transcript_start_distance)
            rear_transcript_flank <- stringr::str_sub(transcript_sequence, snv_to_transcript_start_distance+1+nchar(var_ref), snv_to_transcript_start_distance+1+nchar(var_ref)+flank_length-1) # both start and stop are inclusive so the +1+nchar(var_ref) is the offest so the SNV site isn't included and the -1 is the offset so we don't include an extra base on the rear flank
            flank_sequence <- paste(front_transcript_flank,'[',var_ref,'/',var_alt,']',rear_transcript_flank,sep='')
            # record the flank sequence on the snv_table
            input_variants[row_name,"Variant_front_flank_sequence"] <- front_transcript_flank
            input_variants[row_name,"Variant_rear_flank_sequence"] <- rear_transcript_flank
            input_variants[row_name,"Variant_complete_flank_sequence"] <- flank_sequence
        }
    }
    # do front flank separately from rear flank to keep this from being too insanely long.
    # this is just the front flank so we use distance from the SNV start to the Exon start
    front_exons <- c(exon_with_snv_rank)
    front_exons <- c(front_exons, rev(transcript_exons_table$rank[transcript_exons_table$rank < exon_with_snv_rank]))
    front_coordinates_string <- ""
    count = flank_length
    front_exon_flank_sequence <- ""
    print("front")
    if (strand > 0) { # strand is positive
        print("forwards")
        for (exon_rank in front_exons) {
            exon_start <- transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == exon_rank]
            if (exon_rank == exon_with_snv_rank) {
                count_remainder = count - abs(var_start - exon_start)
                if (abs(var_start - exon_start) <= 5) {
                    input_variants[row_name,"Variant_is_5bp_or_closer_to_a_splice_junction"] <- TRUE
                    print(paste0("Warning: Variant is 5 b or closer to a splice junction. Manually review the variant at chromosome ",var_chr," start ",var_start," stop ",var_stop," reference allele ",var_ref,". This variant is flagged as TRUE in the output column Variant_is_5bp_or_closer_to_a_splice_junction. This variant and flanking site NEEDS to be manually reviewed."))
                }
                if (count_remainder < 1) {
                    print(1)
                    front_flank_start = var_start - count
                    front_flank_stop = var_start - 1
                    front_coordinates_string = paste0(as.character(front_flank_start),"-",as.character(front_flank_stop))
                    front_exon_flank_sequence <- paste0(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=front_flank_start, end=front_flank_stop, strand=strand_char)),front_exon_flank_sequence)
                    count = count_remainder
                    break
                } else {
                    print(2)
                    input_variants[row_name,"front_flank_spans_introns"] <- TRUE
                    front_flank_stop = var_start - 1
                    front_coordinates_string = paste0(as.character(exon_start),"-",as.character(front_flank_stop))
                    front_exon_flank_sequence <- paste0(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=exon_start, end=front_flank_stop, strand=strand_char)),front_exon_flank_sequence)
                    count = count_remainder
                }
            } else {
                exon_stop <- transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == exon_rank]
                exon_len <- transcript_exons_table$exon_length[transcript_exons_table$rank == exon_rank]
                count_remainder = count - exon_len
                if (count_remainder < 1) {
                    print(3)
                    front_flank_start = exon_stop - count + 1 # This may or may not need a +/- 1 offset. I have one here right now but it might need to be changed to the other.
                    front_coordinates_string = paste0(as.character(front_flank_start),"-",as.character(exon_stop),";",front_coordinates_string)
                    front_exon_flank_sequence <- paste0(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=front_flank_start, end=exon_stop, strand=strand_char)),front_exon_flank_sequence)
                    count = count_remainder
                    break
                } else {
                    print(4)
                    front_coordinates_string = paste0(as.character(exon_start),"-",as.character(exon_stop),";",front_coordinates_string)
                    front_exon_flank_sequence <- paste0(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=exon_start, end=exon_stop, strand=strand_char)),front_exon_flank_sequence)
                    count = count_remainder
                }
            }
        }
        # get the upstream sequence coordinates if the flank extends past the first Exon. This only happens when count is still greater than 0 after reaching the end of the vector of exon ranks
        if (count > 0) {
            print(5)
            upstream_start = transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == 1] - count
            upstream_stop = transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == 1] - 1
            front_coordinates_string = paste0(as.character(upstream_start),"-",as.character(upstream_stop),";",front_coordinates_string)
            front_exon_flank_sequence <- paste0(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=upstream_start, end=upstream_stop, strand=strand_char)),"|",front_exon_flank_sequence)
        }
    } else { # strand is negative
        print("reverse")
        for (exon_rank in front_exons) {
            exon_start <- transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == exon_rank]
            if (exon_rank == exon_with_snv_rank) {
                count_remainder = count - abs(exon_start - var_start)
                if (abs(exon_start - var_start) <= 5) {
                    input_variants[row_name,"Variant_is_5bp_or_closer_to_a_splice_junction"] <- TRUE
                    print(paste0("Warning: Variant is 5 b or closer to a splice junction. Manually review the variant at chromosome ",var_chr," start ",var_start," stop ",var_stop," reference allele ",var_ref,". This variant is flagged as TRUE in the output column Variant_is_5bp_or_closer_to_a_splice_junction. This variant and flanking site NEEDS to be manually reviewed."))
                }
                if (count_remainder < 1) {
                    print(1)
                    front_flank_start = var_start + count
                    front_flank_stop = var_start + 1
                    front_coordinates_string = paste0(as.character(front_flank_start),"-",as.character(front_flank_stop)) # the way this is currently will result in the bigger number showing first which may be confusing to some people.
                    front_exon_flank_sequence <- paste0(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=front_flank_stop, end=front_flank_start, strand=strand_char)),front_exon_flank_sequence) # check the order here
                    count = count_remainder
                    break
                } else {
                    print(2)
                    input_variants[row_name,"front_flank_spans_introns"] <- TRUE
                    front_flank_stop = var_start + 1
                    front_coordinates_string = paste0(as.character(exon_start),"-",as.character(front_flank_stop)) # the way this is currently will result in the bigger number showing first which may be confusing to some people.
                    front_exon_flank_sequence <- paste0(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=front_flank_stop, end=exon_start, strand=strand_char)), front_exon_flank_sequence)
                    count = count_remainder
                }
            } else {
                exon_stop <- transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == exon_rank]
                exon_len <- transcript_exons_table$exon_length[transcript_exons_table$rank == exon_rank]
                count_remainder = count - exon_len
                #print(exon_len)
                if (count_remainder < 1) {
                    print(3)
                    print(count)
                    front_flank_start = exon_stop + count - 1 # This might need a +/- 1 offset
                    print(front_flank_start)
                    print(exon_stop)
                    front_coordinates_string = paste0(as.character(front_flank_start),"-",as.character(exon_stop),";",front_coordinates_string) # the way this is currently will result in the bigger number showing first which may be confusing to some people.
                    front_exon_flank_sequence <- paste0(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=exon_stop, end=front_flank_start, strand=strand_char)), front_exon_flank_sequence)
                    count = count_remainder
                    break
                } else {
                    print(4)
                    front_coordinates_string = paste0(as.character(exon_start),"-",as.character(exon_stop),";",front_coordinates_string) # the way this is currently will result in the bigger number showing first which may be confusing to some people.
                    front_exon_flank_sequence <- paste0(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=exon_stop, end=exon_start, strand=strand_char)), front_exon_flank_sequence)
                    count = count_remainder
                }
            }
        }
        if (count > 0) {
            print(5)
            upstream_start = transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == 1] + count
            upstream_stop = transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == 1] + 1
            front_coordinates_string = paste0(as.character(upstream_start),"-",as.character(upstream_stop),";",front_coordinates_string) # the way this is currently will result in the bigger number showing first which may be confusing to some people.
            front_exon_flank_sequence <- paste0(as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=upstream_start, end=upstream_stop, strand=strand_char),"|",front_exon_flank_sequence))
        }
    }
    rear_exons <- c(exon_with_snv_rank)
    rear_exons <- c(rear_exons, transcript_exons_table$rank[transcript_exons_table$rank > exon_with_snv_rank])
    rear_coordinates_string <- ""
    count = flank_length
    rear_exon_flank_sequence <- ""
    print("rear")
    if (strand > 0) { # forwards strand
        print("forwards")
        for (exon_rank in rear_exons) {
            exon_stop <- transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == exon_rank]
            if (exon_rank == exon_with_snv_rank) {
                count_remainder = count - abs(var_stop - exon_stop)
                if (abs(var_stop - exon_stop) <= 5) {
                    input_variants[row_name,"Variant_is_5bp_or_closer_to_a_splice_junction"] <- TRUE
                    print(paste0("Warning: Variant is 5 b or closer to a splice junction. Manually review the variant at chromosome ",var_chr," start ",var_start," stop ",var_stop," reference allele ",var_ref,". This variant is flagged as TRUE in the output column Variant_is_5bp_or_closer_to_a_splice_junction. This variant and flanking site needs to be manually reviewed."))
                }
                if (count_remainder < 1) {
                    print(1)
                    rear_flank_start = var_stop + 1
                    rear_flank_stop = var_stop + count
                    rear_coordinates_string = paste0(as.character(rear_flank_start), "-", as.character(rear_flank_stop))
                    rear_exon_flank_sequence <- paste0(rear_exon_flank_sequence, as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=rear_flank_start, end=rear_flank_stop, strand=strand_char)))
                    count = count_remainder
                    break
                } else {
                    print(2)
                    input_variants[row_name,"rear_flank_spans_introns"] <- TRUE
                    rear_flank_start = var_stop + 1
                    rear_coordinates_string = paste0(as.character(rear_flank_start), "-", as.character(exon_stop))
                    rear_exon_flank_sequence <- paste0(rear_exon_flank_sequence, as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=rear_flank_start, end=exon_stop, strand=strand_char)))
                    count = count_remainder
                }
            } else {
                exon_start <- transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == exon_rank]
                exon_len <- transcript_exons_table$exon_length[transcript_exons_table$rank == exon_rank]
                count_remainder = count - exon_len
                if (count_remainder < 1) {
                    print(3)
                    rear_flank_stop = exon_start + count - 1  # This might need a +/- 1 offset
                    rear_coordinates_string = paste0(rear_coordinates_string,";",as.character(exon_start),"-",as.character(rear_flank_stop))
                    rear_exon_flank_sequence <- paste0(rear_exon_flank_sequence, as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=exon_start, end=rear_flank_stop, strand=strand_char)))
                    count = count_remainder
                    break
                } else {
                    print(4)
                    print(exon_start)
                    print(exon_stop)
                    rear_coordinates_string = paste0(rear_coordinates_string,";",as.character(exon_start),"-",as.character(exon_stop))
                    rear_exon_flank_sequence <- paste0(rear_exon_flank_sequence, as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=exon_start, end=exon_stop, strand=strand_char)))
                    count = count_remainder
                }
            }
        }
        if (count > 0) {
            print(5)
            downstream_start = transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == max(transcript_exons_table$rank)] + 1
            downstream_stop = transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == max(transcript_exons_table$rank)] + count
            rear_coordinates_string = paste0(rear_coordinates_string,';',as.character(downstream_start),"-",as.character(downstream_stop)) # the way this is currently will result in the bigger number showing first which may be confusing to some people.
            rear_exon_flank_sequence <- paste0(rear_exon_flank_sequence,'|',as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=downstream_start, end=downstream_stop, strand=strand_char)))
        }
    } else { # strand is negative
        print("reverse")
        for (exon_rank in rear_exons) {
            exon_stop <- transcript_exons_table$exon_chrom_start[transcript_exons_table$rank == exon_rank]
            if (exon_rank == exon_with_snv_rank) {
              count_remainder = count - abs(var_stop - exon_stop)
                print(paste0("count_remainder is ",count_remainder))
                if (abs(var_stop - exon_stop) <= 5) {
                    input_variants[row_name,"Variant_is_5bp_or_closer_to_a_splice_junction"] <- TRUE
                    print(paste0("Warning: Variant is 5 b or closer to a splice junction. Manually review the variant at chromosome ",var_chr," start ",var_start," stop ",var_stop," reference allele ",var_ref,". This variant is flagged as TRUE in the output column Variant_is_5bp_or_closer_to_a_splice_junction. This variant and flanking site needs to be manually reviewed."))
                }
                if (count_remainder < 1) {
                    print(1)
                    rear_flank_start = var_stop - 1
                    rear_flank_stop = var_stop - count 
                    rear_coordinates_string = paste0(as.character(rear_flank_start),"-",as.character(rear_flank_stop)) # the way this is currently will result in the bigger number showing first which may be confusing to some people.
                    rear_exon_flank_sequence <- paste0(rear_exon_flank_sequence, as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=rear_flank_stop, end=rear_flank_start, strand=strand_char))) # check the order here
                    count = count_remainder
                    break
                } else {
                    print(2)
                    input_variants[row_name,"rear_flank_spans_introns"] <- TRUE
                    rear_flank_start = var_stop - 1
                    rear_coordinates_string = paste0(as.character(rear_flank_start),"-",as.character(exon_stop)) # the way this is currently will result in the bigger number showing first which may be confusing to some people.
                    rear_exon_flank_sequence <- paste0(rear_exon_flank_sequence, as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=exon_stop, end=rear_flank_start, strand=strand_char)))
                    count = count_remainder
                }
            } else {
                exon_start <- transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == exon_rank]
                exon_len <- transcript_exons_table$exon_length[transcript_exons_table$rank == exon_rank]
                count_remainder = count - exon_len
                if (count_remainder < 1) {
                    print(3)
                    rear_flank_stop = exon_start - count + 1 # This might need a +/- 1 offset
                    rear_coordinates_string = paste0(rear_coordinates_string,";",as.character(exon_start),"-",as.character(rear_flank_stop))
                    rear_exon_flank_sequence <- paste0(rear_exon_flank_sequence, as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=rear_flank_stop, end=exon_start, strand=strand_char)))
                    count = count_remainder
                    break
                } else {
                    print(4)
                    rear_coordinates_string = paste0(rear_coordinates_string,";",as.character(exon_start),"-",as.character(exon_stop))
                    rear_exon_flank_sequence <- paste0(rear_exon_flank_sequence, as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=exon_stop, end=exon_start, strand=strand_char)))
                    count = count_remainder
                }
            }
        }
        if (count > 0) {
            print(5)
            downstream_start = transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == max(transcript_exons_table$rank)] - 1
            downstream_stop = transcript_exons_table$exon_chrom_end[transcript_exons_table$rank == max(transcript_exons_table$rank)] - count
            rear_coordinates_string = paste0(rear_coordinates_string,';',as.character(downstream_start),"-",as.character(downstream_stop)) # the way this is currently will result in the bigger number showing first which may be confusing to some people.
            rear_exon_flank_sequence <- paste0(rear_exon_flank_sequence,'|',as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, paste0('chr',var_chr), start=downstream_start, end=downstream_stop, strand=strand_char)))
        }
    }
    input_variants[row_name,"Variant_front_flank_Start"] <- front_flank_start
    input_variants[row_name,"Variant_front_flank_Stop"] <- front_flank_stop
    input_variants[row_name,"Variant_rear_flank_Start"] <- rear_flank_start
    input_variants[row_name,"Variant_rear_flank_Stop"] <- rear_flank_stop
    input_variants[row_name,"Front_flank_full_coordinates"] <- front_coordinates_string
    input_variants[row_name,"Front_flank_coordinates_based_sequence"] <- front_exon_flank_sequence
    input_variants[row_name,"Rear_flank_full_coordinates"] <- rear_coordinates_string
    input_variants[row_name,"Rear_flank_coordinates_based_sequence"] <- rear_exon_flank_sequence
}
# now check for preferred, neutral, and non-preferred junctions.
# Preferred is AT, TA, GA, AG
# Neutral is TT, CT, CA, TC, AC, CC, TG, AA
# Non-preferred is CG, GT, GG, GC
# checking if there is an insertion/deletion event or if everything is a SNV.


# this block only covers when the variants are SNVs. If a variant is an indel it may not be accurate and should be manually reviewed.
# front/Ref_rear ligation means the junction is the last base of the front flank and the first base of the reference allele.
input_variants$ref_junction_v1[input_variants$Strand == "+"] <- paste0(substr(gsub("[|]","",input_variants$Variant_front_flank_sequence[input_variants$Strand == "+"]),nchar(gsub("[|]","",input_variants$Variant_front_flank_sequence[input_variants$Strand == "+"])),nchar(gsub("[|]","",input_variants$Variant_front_flank_sequence[input_variants$Strand == "+"]))), substr(input_variants$Ref_allele[input_variants$Strand == "+"], 1,1))
# front_Ref/rear ligation means the junction is the last base of reference allele and the first base of the rear flank.
input_variants$ref_junction_v2[input_variants$Strand == "+"] <- paste0(substr(input_variants$Ref_allele[input_variants$Strand == "+"], nchar(input_variants$Ref_allele[input_variants$Strand == "+"]),nchar(input_variants$Ref_allele[input_variants$Strand == "+"])),substr(gsub("[|]","",input_variants$Variant_rear_flank_sequence[input_variants$Strand == "+"]),1,1))
# front/Alt_rear ligation means the junction is the last base of the front flank and the first base of the alternate allele.
input_variants$alt_junction_v1[input_variants$Strand == "+"] <- paste0(substr(gsub("[|]","",input_variants$Variant_front_flank_sequence[input_variants$Strand == "+"]),nchar(gsub("[|]","",input_variants$Variant_front_flank_sequence[input_variants$Strand == "+"])),nchar(gsub("[|]","",input_variants$Variant_front_flank_sequence[input_variants$Strand == "+"]))), substr(input_variants$Alt_allele[input_variants$Strand == "+"], 1,1))
# front_Alt/rear ligation means the junction is the last base of the alternate allele and the first base of the rear flank.
input_variants$alt_junction_v2[input_variants$Strand == "+"] <- paste0(substr(input_variants$Alt_allele[input_variants$Strand == "+"], nchar(input_variants$Alt_allele[input_variants$Strand == "+"]),nchar(input_variants$Alt_allele[input_variants$Strand == "+"])),substr(gsub("[|]","",input_variants$Variant_rear_flank_sequence[input_variants$Strand == "+"]),1,1))
# front/Ref_rear ligation means the junction is the last base of the front flank and the first base of the reference allele.
input_variants$ref_junction_v1[input_variants$Strand == "-"] <- paste0(substr(gsub("[|]","",input_variants$Variant_front_flank_sequence[input_variants$Strand == "-"]),nchar(gsub("[|]","",input_variants$Variant_front_flank_sequence[input_variants$Strand == "-"])),nchar(gsub("[|]","",input_variants$Variant_front_flank_sequence[input_variants$Strand == "-"]))), stri_reverse(chartr("ATCG","TAGC",substr(input_variants$Ref_allele[input_variants$Strand == "-"], 1,1))))
# front_Ref/rear ligation means the junction is the last base of reference allele and the first base of the rear flank.
input_variants$ref_junction_v2[input_variants$Strand == "-"] <- paste0(stri_reverse(chartr("ATCG","TAGC",substr(input_variants$Ref_allele[input_variants$Strand == "-"], nchar(input_variants$Ref_allele[input_variants$Strand == "-"]),nchar(input_variants$Ref_allele[input_variants$Strand == "-"])))),substr(gsub("[|]","",input_variants$Variant_rear_flank_sequence[input_variants$Strand == "-"]),1,1))
# front/Alt_rear ligation means the junction is the last base of the front flank and the first base of the alternate allele.
input_variants$alt_junction_v1[input_variants$Strand == "-"] <- paste0(substr(gsub("[|]","",input_variants$Variant_front_flank_sequence[input_variants$Strand == "-"]),nchar(gsub("[|]","",input_variants$Variant_front_flank_sequence[input_variants$Strand == "-"])),nchar(gsub("[|]","",input_variants$Variant_front_flank_sequence[input_variants$Strand == "-"]))), stri_reverse(chartr("ATCG","TAGC",substr(input_variants$Alt_allele[input_variants$Strand == "-"], 1,1))))
# front_Alt/rear ligation means the junction is the last base of the alternate allele and the first base of the rear flank.
input_variants$alt_junction_v2[input_variants$Strand == "-"] <- paste0(stri_reverse(chartr("ATCG","TAGC",substr(input_variants$Alt_allele[input_variants$Strand == "-"], nchar(input_variants$Alt_allele[input_variants$Strand == "-"]),nchar(input_variants$Alt_allele[input_variants$Strand == "-"])))),substr(gsub("[|]","",input_variants$Variant_rear_flank_sequence[input_variants$Strand == "-"]),1,1))

input_variants$Ref_Ligation_junction_quality <- "Manual_review_needed"
input_variants$Alt_Ligation_junction_quality <- "Manual_review_needed"
preferred_junctions = c("AT","TA","GA","AG")
neutral_junctions = c("TT","CT","CA","TC","AC","CC","TG","AA")
nonpreferred_junctions = c("CG","GT","GG","GC")
input_variants$Ref_Ligation_junction_quality[(input_variants$ref_junction_v1 %in% nonpreferred_junctions) | (input_variants$ref_junction_v2 %in% nonpreferred_junctions)] <- "Non-preferred"
input_variants$Ref_Ligation_junction_quality[(input_variants$ref_junction_v1 %in% neutral_junctions) | (input_variants$ref_junction_v2 %in% neutral_junctions)] <- "Neutral"
input_variants$Ref_Ligation_junction_quality[(input_variants$ref_junction_v1 %in% preferred_junctions) | (input_variants$ref_junction_v2 %in% preferred_junctions)] <- "Preferred"
input_variants$Alt_Ligation_junction_quality[(input_variants$alt_junction_v1 %in% nonpreferred_junctions) | (input_variants$alt_junction_v2 %in% nonpreferred_junctions)] <- "Non-preferred"
input_variants$Alt_Ligation_junction_quality[(input_variants$alt_junction_v1 %in% neutral_junctions) | (input_variants$alt_junction_v2 %in% neutral_junctions)] <- "Neutral"
input_variants$Alt_Ligation_junction_quality[(input_variants$alt_junction_v1 %in% preferred_junctions) | (input_variants$alt_junction_v2 %in% preferred_junctions)] <- "Preferred"
input_variants$Ref_Ligation_junction_quality[(input_variants$Ref_allele == "-")] <- "Manual_review_needed"
input_variants$Alt_Ligation_junction_quality[(input_variants$Alt_allele == "-")] <- "Manual_review_needed"

write.table(input_variants, file = "/diskmnt/Projects/Users/austins2/tools/xenium_snvs/version_4/input_variant_targets_flanks.tsv",sep='\t', quote = F)
# The indel flank sites in the output are wrong (at least the rear flank for them is anyways). 
# For now there are few enough and I don't have time to figure out what is going on with them right now so I'm just going to manually review and generate them using the Ensembl transcript sequence and the UCSC genome browser.
# If I have time I'll come back and fix this. 
