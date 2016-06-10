# The following script creates per chromosome one shell scipt to run isPcr and
# submits all the scripts to the pbs batch system

# Create all chromosome paths
Chroms <- paste("chr", c(1:22, "X", "Y"), sep = "")
ChromPaths <- paste("/home/hzudohna/L1polymORF/Data/", Chroms,
                    ".fas", sep = "")

# Create all output paths
OutputPaths <- paste("/share/diskarray3/hzudohna/L1Primers/L1catalog.isPcr.output.", 
                     Chroms, ".txt", sep = "")

# Loop through paths, create a pbs script that calls isPcr and submit the
# script
for (i in 1:length(Chroms)){
  CMDLines <- c('/share/diskarray1/isPcr/isPcr \\',
                paste(ChromPaths[i], '\\'),
                '/share/diskarray3/hzudohna/L1Primers/L1catalog.isPcr.input.txt \\',
                OutputPaths[i])
  ScriptPath <- paste("/home/hzudohna/qsubScript_isPcr", Chroms[i], sep = "_")
  pbsName    <- paste("isPcr", Chroms[i], sep = "_")
  CreateAndCallqsubScript(file = ScriptPath, qsubCommandLines = CMDLines,
                          scriptName = pbsName)
}