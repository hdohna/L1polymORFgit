# The script below explore a catalog of full-length L1


# Read repeat table and subset to get only L1HS rows with fragment size below 
# MaxFragLength
RepeatTable <- read.delim("D:/L1polymORF/Data/repeatsHg38")
RepeatTable <- RepeatTable[nchar(as.character(RepeatTable$genoName)) <= 5, ]
RepeatTable <- RepeatTable[RepeatTable$repName == "L1HS",]

# Save table as csv file
write.csv(RepeatTable, "D:/L1polymORF/Data/repeatsHg38_L1HS.csv")